# ==============================================================================
# 1. CORE UTILITY FUNCTIONS
# ==============================================================================

#' @keywords internal
#' @noRd
slice_to_matrix <- function(arr, k) {
  m <- arr[,,k, drop = FALSE]
  dim(m) <- dim(arr)[1:2]
  dimnames(m) <- dimnames(arr)[1:2]
  return(m)
}

#' Iterated Grid Search for Function Optimization
iterated_grid_search <- function(min_param, max_param, func, step_length = 100,
                                 max_iterations = 5, precision_threshold = 1e-3,
                                 exp_transform = FALSE, ...) {
  for (iteration in 1:max_iterations) {
    search_grid <- seq(min_param, max_param, length.out = step_length)
    if (exp_transform) search_grid <- exp(search_grid)
    
    func_values <- sapply(search_grid, func, ...)
    min_index <- which.min(func_values)
    
    best_param <- search_grid[min_index]
    best_value <- func_values[min_index]
    
    # Calculate the new range safely
    new_min_param <- search_grid[max(1, min_index - 1)]
    new_max_param <- search_grid[min(length(search_grid), min_index + 1)]
    
    if (abs(new_max_param - new_min_param) < precision_threshold) break
    
    if (exp_transform) {
      min_param <- log(new_min_param)
      max_param <- log(new_max_param)
    } else {
      min_param <- new_min_param
      max_param <- new_max_param
    }
  }
  
  return(list(best_param = best_param, best_value = best_value))
}

# Helper to combine mnl and exo priors
#' @keywords internal
#' @noRd
combine_utilities <- function(utilitites_mnl, utilitites_exo, weights) {
  intersecting_cols <- intersect(colnames(utilitites_mnl), colnames(utilitites_exo))
  if (length(intersecting_cols) > 0) {
    temp_utilities_comb <- (1 - weights[, intersecting_cols, drop = FALSE]) * utilitites_mnl[, intersecting_cols, drop = FALSE] +
      weights[, intersecting_cols, drop = FALSE] * utilitites_exo[, intersecting_cols, drop = FALSE]
    colnames(temp_utilities_comb) <- intersecting_cols
    return(temp_utilities_comb)
  }
  return(NULL)
}

# ==============================================================================
# 2. OPTIMIZATION CONFIGURATION & ADAPTIVE WRAPPERS
# ==============================================================================

#' @keywords internal
#' @noRd
make_nloptr_opts <- function(options) {
  list(
    algorithm = options$algorithm,
    objective = get(options$obj_fun, mode = "function", inherits = TRUE),
    gradient  = if (!is.null(options$grad_obj_fun)) get(options$grad_obj_fun, mode = "function", inherits = TRUE) else NULL,
    ineqc     = if (!is.null(options$ineqc_fun) && options$algorithm == "NLOPT_LN_COBYLA") get(options$ineqc_fun, mode = "function", inherits = TRUE) else NULL,
    xtol_rel = options$xtol_rel,
    xtol_abs = options$xtol_abs,
    maxeval  = options$max_eval,
    stopval  = options$max_diff
  )
}

#' @keywords internal
#' @noRd
adaptive_nloptr <- function(initial_solution, x0, lb, ub, opts_base, priors_sub, 
                            targets_sub, restrictions_sub, options, nloptr_extra_args) {
  
  res <- list(
    objective = options$max_diff * 1000,
    solution = initial_solution,
    lower_bounds = lb,
    upper_bounds = ub,
    options = opts_base
  )
  runs <- list()
  
  current_try <- 1
  max_tries <- max(1, options$max_retry)
  previous_objective <- Inf # For early bailout tracking
  
  # Combine lists once before the loop to skip do.call() list concatenation overhead
  full_args <- c(
    list(
      eval_f = opts_base$objective, eval_grad_f = opts_base$gradient,
      eval_g_ineq = opts_base$ineqc, opts = res$options, 
      mu0 = priors_sub, targets = targets_sub, 
      restrictions = restrictions_sub, cutoff = options$cutoff, 
      lambda = options$pen_lambda, log_scale = options$optimize_log_scale
    ),
    nloptr_extra_args
  )
  
  while ((current_try <= max_tries) && (res$objective > options$max_diff)) {
    # Dynamically scale evaluation limits
    res$options$maxeval <- max(options$min_eval, ceiling(options$max_eval * (current_try / max_tries)))
    
    # Expand bounds based on scale
    if (options$optimize_log_scale) {
      res$lower_bounds <- lb - (current_try - 1) * 5
      res$upper_bounds <- ub + (current_try - 1) * 5
      res$solution <- pmax(pmin(res$solution, res$upper_bounds - 1e-6), res$lower_bounds + 1e-6)
    } else {
      res$lower_bounds <- lb / (5^(current_try - 1))
      res$upper_bounds <- ub * (5^(current_try - 1))
      res$solution <- pmax(pmin(res$solution, res$upper_bounds * (1 - 1e-6)), res$lower_bounds * (1 + 1e-6))
    }
    
    names(res$lower_bounds) <- names(res$upper_bounds) <- names(x0)
    
    # JITTER: Inject 1% random noise to kick it out of local minima
    if (current_try > 1) {
      noise <- runif(length(res$solution), min = 0.99, max = 1.01)
      if (options$optimize_log_scale) {
        res$solution <- res$solution + log(noise) 
      } else {
        res$solution <- res$solution * noise
      }
      res$solution <- pmax(pmin(res$solution, res$upper_bounds), res$lower_bounds)
    }
    
    # Inject the dynamic variables directly into the pre-packed list
    full_args$x0 <- res$solution
    full_args$lb <- res$lower_bounds
    full_args$ub <- res$upper_bounds
    full_args$opts$maxeval <- res$options$maxeval
    
    nloptr_out <- do.call(nloptr::nloptr, full_args)
    
    nloptr_out$message <- paste0(sprintf("Adaptive try %d\n", current_try), nloptr_out$message)
    runs[[paste0("try_", current_try)]] <- nloptr_out
    
    res$solution  <- nloptr_out$solution
    res$objective <- if (!is.null(nloptr_out$objective)) nloptr_out$objective else res$objective
    res$message   <- nloptr_out$message
    names(res$solution) <- names(x0)
    
    # DEAD END BAILOUT: Break if improved by less than 0.01%
    if (current_try > 1 && (previous_objective - res$objective) < (previous_objective * 0.0001)) {
      if (!is.null(options$err.txt) && options$err.txt != "") {
        cat(sprintf("[%s] NLOPTR trapped in dead end at Try %d. Bailing out early.\n", options$err.txt, current_try))
      }
      break
    }
    
    previous_objective <- res$objective
    current_try <- current_try + 1
  }
  
  return(list(final = res, runs = runs))
}

# ==============================================================================
# 3. STATE MANAGEMENT & REPORTING
# ==============================================================================

#' @keywords internal
#' @noRd
calculate_target_deviation <- function(pred_totals, target_totals, options) {
  # Call the centralized function from optim_aux_func.R
  calculate_objective_error(pred_totals, target_totals, options$obj_fun, options$optimize_log_scale)
}

#' @keywords internal
#' @noRd
update_residual_state <- function(targets, current_mu = NULL, areas = NULL, 
                                  restrictions = NULL, pixel_avail, type = "lu") {
  
  # Initialize a zero-matrix for tracking allocations at the pixel level
  allocated_tracker_matrix <- matrix(0, 
                                     nrow = length(pixel_avail), 
                                     ncol = length(targets),
                                     dimnames = list(pixel_avail, names(targets)))
  
  res <- list(
    targets = targets,
    areas = if (!is.null(areas)) areas[pixel_avail] else NULL,
    already_allocated = allocated_tracker_matrix  
  )
  
  # 1. Update based on what is ALREADY allocated (if a state is provided)
  if (!is.null(current_mu)) {
    # Robustly handle matrix or vector inputs
    if (is.matrix(current_mu) || is.data.frame(current_mu)) {
      res$already_allocated <- current_mu[pixel_avail, , drop = FALSE]
    } else {
      # Handle vector input by converting to single-column matrix if it matches pixel count
      # This addresses "incorrect number of dimensions" errors when single-crop states are passed
      current_mu_vec <- as.vector(current_mu)
      res$already_allocated[, 1] <- current_mu_vec[pixel_avail]
    }
    allocated_totals <- colSums(res$already_allocated)
    
    # ---- TARGET OVERSHOOT PROTECTION ----
    overshoot_ratios <- targets / pmax(allocated_totals, 1e-12)
    bad_targets <- overshoot_ratios < (1 - 1e-10) 
    
    if (any(bad_targets)) {
      for (tname in names(targets)[bad_targets]) {
        res$already_allocated[, tname] <- res$already_allocated[, tname] * overshoot_ratios[tname]
      }
      allocated_totals <- colSums(res$already_allocated)
    }
    
    res$targets <- pmax(targets - allocated_totals, 0) # Prevent negative targets
    
    if (type == "lu" && !is.null(areas)) {
      used_land <- rowSums(res$already_allocated)
      res$areas <- pmax(areas[pixel_avail] - used_land, 0) # Prevent negative land
    }
  }
  
  return(res)
}

#' @keywords internal
#' @noRd
summarize_allocations <- function(out_mu, active_targets) {
  allocated_totals <- colSums(out_mu)
  
  summary_df <- data.frame(
    Target_Name = names(active_targets),
    Target_Amount = as.numeric(active_targets),
    Allocated_Amount = as.numeric(allocated_totals),
    Difference = as.numeric(allocated_totals - active_targets)
  )
  
  summary_df$Percent_Error <- ifelse(
    summary_df$Target_Amount == 0, 
    0, 
    round((summary_df$Difference / summary_df$Target_Amount) * 100, 4)
  )
  
  return(list(
    breakdown = summary_df,
    grand_total = sum(allocated_totals)
  ))
}

#' @keywords internal
#' @noRd
normalize_overshoot <- function(out_mu, areas, pixel_avail) {
  actual_pixel_sums <- rowSums(out_mu[pixel_avail, , drop = FALSE])
  overshoot_ratios <- actual_pixel_sums / pmax(areas[pixel_avail], 1e-12)
  bad_pixels <- overshoot_ratios > (1 + 1e-10)
  
  if (any(bad_pixels)) {
    out_mu[pixel_avail[bad_pixels], ] <- out_mu[pixel_avail[bad_pixels], ] / overshoot_ratios[bad_pixels]
  }
  return(out_mu)
}

# ==============================================================================
# 4. SOLVER PHASES
# ==============================================================================

#' @keywords internal
#' @noRd
generalized_ipf_fit <- function(x0 = NULL, mu0, targets, alloc_func, alloc_args,
                                 max_iter = 1000, tolerance = 1e-4, damping = 0.5, 
                                 verbose = FALSE, log_scale = TRUE, log_mu0 = NULL,
                                 options = list()) {
  n_targets <- length(targets)
  current_iter <- 0
  error_change <- Inf
  
  if (log_scale) {
    x <- if (!is.null(x0) && length(x0) == n_targets) x0 else rep(0, n_targets)
    log_targ <- log(pmax(targets, 1e-12))
    if (is.null(log_mu0)) log_mu0 <- log(pmax(mu0, 1e-12))
    alloc_call <- as.call(c(list(alloc_func, x = quote(x), log_scale = TRUE, log_mu0 = log_mu0), alloc_args))
  } else {
    x <- if (!is.null(x0) && length(x0) == n_targets) pmax(x0, 1e-12) else rep(1, n_targets)
    alloc_call <- as.call(c(list(alloc_func, x = quote(x), mu0 = mu0, log_scale = FALSE), alloc_args))
  }
  
  while (current_iter < max_iter && !is.na(error_change) && error_change > tolerance) {
    if (log_scale) {
      mu_current <- eval(alloc_call)
      log_pred <- matrixStats::colLogSumExps(mu_current)
      
      if (isTRUE(verbose) && (current_iter %% 50 == 0)) {
        pred_linear <- exp(log_pred)
        cat(sprintf("   IPF Iter %d: Max Dev %.4f, Total Pred %.2f, Total Targ %.2f\n", 
                    current_iter, error_change, sum(pred_linear), sum(targets)))
      }
      
      delta_logx <- damping * (log_targ - log_pred)
      delta_logx[is.na(delta_logx)] <- 0
      delta_logx <- pmin(pmax(delta_logx, -50), 50)
      
      x <- pmin(pmax(x + delta_logx, -400), 400)
      error_change <- max(abs(delta_logx))
    } else {
      mu_current <- eval(alloc_call)
      pred <- colSums(mu_current)
      
      if (isTRUE(verbose) && (current_iter %% 50 == 0)) {
        cat(sprintf("   IPF Iter %d: Max Dev %.4e, Total Pred %.2f, Total Targ %.2f\n", 
                    current_iter, error_change, sum(pred), sum(targets)))
      }
      
      ratio <- targets / pmax(pred, 1e-12)
      ratio[is.na(ratio)] <- 1
      ratio <- pmin(pmax(ratio, 1e-12), 1e12)
      
      x_new <- x * (ratio ^ damping)
      error_change <- max(abs(log(x_new / pmax(x, 1e-12))))
      x <- pmax(x_new, 1e-12)
    }
    current_iter <- current_iter + 1
  }
  
  converged <- (error_change <= tolerance)
  
  if (log_scale) {
    mu_final <- eval(alloc_call)
    final_pred <- exp(matrixStats::colLogSumExps(mu_final))
  } else {
    mu_final <- eval(alloc_call)
    final_pred <- colSums(mu_final)
  }
  
  final_dev <- calculate_target_deviation(final_pred, targets, options)
  converged <- !any(final_dev > tolerance)
  
  if (isTRUE(verbose)) {
    cat(sprintf("   IPF Final: Max Dev %.4e, Total Pred %.2f, Total Targ %.2f, Converged: %s\n", 
                max(final_dev), sum(final_pred), sum(targets), converged))
  }
  
  return(list(solution = x, iterations = current_iter, converged = converged, final_deviation = final_dev))
}

#' @keywords internal
#' @noRd
run_ipf_phase <- function(master_out_mu, x0, active_targets, priors, pixel_avail, 
                          alloc_func, alloc_args, log_priors, options) {
  
  res <- generalized_ipf_fit(
    x0 = x0, mu0 = priors[pixel_avail, , drop = FALSE], targets = active_targets,
    alloc_func = alloc_func, alloc_args = alloc_args, tolerance = options$max_diff,
    damping = if(options$optimize_log_scale) 0.7 else 1.0, 
    max_iter = if(isTRUE(options$use_ipf) && isFALSE(options$use_multi_nloptr) && isFALSE(options$use_single_nloptr)) options$max_eval else options$min_eval,
    log_scale = options$optimize_log_scale, log_mu0 = log_priors, verbose = options$verbose,
    options = options
  )
  
  phase_mu <- do.call(alloc_func, c(list(x = res$solution, log_scale = options$optimize_log_scale, 
                                         log_mu0 = log_priors, mu0 = priors[pixel_avail, , drop = FALSE]), alloc_args))
  if (options$optimize_log_scale) phase_mu <- exp(phase_mu)
  
  deviation_vec <- calculate_target_deviation(colSums(phase_mu), active_targets, options)
  
  list(
    phase_mu = phase_mu,
    out_mu = master_out_mu + phase_mu,
    solution = res$solution,
    converged = res$converged, 
    runs = list(
      converged = res$converged,
      max_dev = max(deviation_vec),
      iters = res$iterations
    )
  )
}

#' @keywords internal
#' @noRd
run_multi_nloptr_phase <- function(master_out_mu, curr_sol, x0, active_targets, priors, restrictions, pixel_avail, 
                             alloc_func, alloc_args, log_priors, nloptr_extra_args, options) {
  
  if (options$optimize_log_scale) {
    lb <- curr_sol - 20
    ub <- curr_sol + 20
  } else {
    lb <- pmin(curr_sol, options$MIN_EXP) / 20
    ub <- pmax(curr_sol, x0) * 20
  }
  names(lb) <- names(ub) <- names(x0)
  
  adaptive_res <- adaptive_nloptr(
    initial_solution = curr_sol, x0 = x0, lb = lb, ub = ub, opts_base = make_nloptr_opts(options),
    priors_sub = priors[pixel_avail, , drop = FALSE], targets_sub = active_targets,
    restrictions_sub = restrictions[pixel_avail, , drop = FALSE], options = options, nloptr_extra_args = nloptr_extra_args
  )
  
  new_sol <- adaptive_res$final$solution
  
  phase_mu <- do.call(alloc_func, c(list(x = new_sol, log_scale = options$optimize_log_scale, 
                                         log_mu0 = log_priors, mu0 = priors[pixel_avail, , drop = FALSE]), alloc_args))
  if (options$optimize_log_scale) phase_mu <- exp(phase_mu)
  
  deviation_vec <- calculate_target_deviation(colSums(phase_mu), active_targets, options)
  
  list(
    phase_mu = phase_mu,
    out_mu = master_out_mu + phase_mu,
    solution = new_sol,
    converged = !any(deviation_vec > options$max_diff),
    runs = list(
      converged = !any(deviation_vec > options$max_diff),
      max_dev = max(deviation_vec),
      adaptive_summary = adaptive_res$final
    )
  )
}

#' @keywords internal
#' @noRd
run_single_nloptr_phase <- function(master_out_mu, curr_sol, x0, active_targets, active_areas, priors, restrictions, 
                            pixel_avail, type, alloc_func, alloc_args, nloptr_extra_args, log_priors, options, res_main_final) {
  
  phase_mu <- matrix(0, nrow = length(pixel_avail), ncol = length(active_targets), dimnames = list(pixel_avail, names(active_targets)))
  solver_runs <- list()
  boost_opts_base <- make_nloptr_opts(options)
  
  failing_targets <- names(active_targets[active_targets > options$max_diff])

  if (length(failing_targets) > 1) {
    deviation_vec <- calculate_target_deviation(rep(0, length(active_targets)), active_targets, options)
    # 1. Safely extract the exact deviation values for ONLY the failing targets
    devs <- deviation_vec[failing_targets]
    
    if (!is.null(restrictions)) {
      # 2. Count restricted pixels
      restr_counts <- colSums(restrictions[pixel_avail, failing_targets, drop = FALSE], na.rm = TRUE)
      
      # 3. Sort Primary: Most restricted. Secondary: Highest deviation.
      failing_targets <- failing_targets[order(restr_counts, devs, decreasing = TRUE)]
    } else {
      # If no restrictions, sort purely by the highest deviation
      failing_targets <- failing_targets[order(devs, decreasing = TRUE)]
    }
  }  
  
  for (tname in failing_targets) {
    potential_mask <- if (!is.null(restrictions)) (1 - restrictions[pixel_avail, tname]) else rep(1, length(pixel_avail))
    available_and_allowed <- active_areas * potential_mask
    max_possible_area <- sum(available_and_allowed)
    
    if (max_possible_area > 1e-9) {
      curr_lb <- res_main_final$lower_bounds[[tname]]
      curr_ub <- res_main_final$upper_bounds[[tname]]
      x0_single <- pmax(pmin(x0[[tname]], curr_ub - 1e-6), curr_lb + 1e-6) 
      
      boost_call_args <- list(
        x0 = x0_single, lb = curr_lb, ub = curr_ub, opts = res_main_final$options, 
        targets = active_targets[[tname]], mu0 = priors[pixel_avail, tname, drop = FALSE],
        eval_f = boost_opts_base$objective, eval_grad_f = boost_opts_base$gradient,
        restrictions = if (!is.null(restrictions)) restrictions[pixel_avail, tname, drop = FALSE] else NULL,
        cutoff = options$cutoff, lambda = options$pen_lambda, log_scale = options$optimize_log_scale
      )
      
      current_boost_extra_args <- nloptr_extra_args
      if (type == "lu") current_boost_extra_args$areas <- as.matrix(available_and_allowed)
      
      res_single <- do.call(nloptr::nloptr, c(boost_call_args, current_boost_extra_args))
      solver_runs[[paste0("boost_", tname)]] <- res_single
      
      temp_alloc_args <- alloc_args
      if (type == "lu") temp_alloc_args$areas <- active_areas
      if (!is.null(restrictions)) temp_alloc_args$restrictions <- restrictions[pixel_avail, tname, drop = FALSE]
      
      temp_mu <- do.call(alloc_func, c(list(
        x = res_single$solution, log_scale = options$optimize_log_scale, 
        log_mu0 = if(options$optimize_log_scale) log_priors[, tname, drop=FALSE] else NULL, 
        mu0 = if(!options$optimize_log_scale) priors[pixel_avail, tname, drop = FALSE] else NULL
      ), temp_alloc_args))
      
      if (options$optimize_log_scale) temp_mu <- exp(temp_mu)
      phase_mu[pixel_avail, tname] <- as.vector(temp_mu)
      curr_sol[[tname]] <- res_single$solution 
    }
  }
  
  new_out_mu <- master_out_mu + phase_mu
  new_deviation_vec <- calculate_target_deviation(colSums(new_out_mu), colSums(master_out_mu) + active_targets, options)
  
  list(
    phase_mu = phase_mu,
    out_mu = new_out_mu,
    solution = curr_sol,
    converged = !any(new_deviation_vec > options$max_diff),
    runs = list(
      converged = !any(new_deviation_vec > options$max_diff),
      max_dev = max(new_deviation_vec),
      failing_targets = failing_targets
    )
  )
}

#' @keywords internal
#' @noRd
run_mopup_phase <- function(master_out_mu, curr_sol, active_targets, active_areas, restrictions, pixel_avail, type, options) {
  
  phase_mu <- matrix(0, nrow = length(pixel_avail), ncol = length(active_targets), dimnames = list(pixel_avail, names(active_targets)))
  solver_runs <- list()
  failing_targets <- names(active_targets[active_targets > options$max_diff])
  
  if (isTRUE(options$verbose)) {
    cat(sprintf("   [Mop-up] Phase called for %s downscaling. Failing targets count: %d\n", 
                type, length(failing_targets)))
  }
  
  if (length(failing_targets) > 1) {
    # Extract the remaining amounts for the failing targets
    amounts <- active_targets[failing_targets]
    
    if (!is.null(restrictions)) {
      # Count restricted pixels
      restr_counts <- colSums(restrictions[pixel_avail, failing_targets, drop = FALSE], na.rm = TRUE)
      
      # order() evaluates sequentially: Primary = restr_counts, Secondary = amounts
      failing_targets <- failing_targets[order(restr_counts, amounts, decreasing = TRUE)]
    } else {
      # If no restrictions exist, fallback to sorting purely by largest amount
      failing_targets <- failing_targets[order(amounts, decreasing = TRUE)]
    }
  }
  
  for (tname in failing_targets) {
    remainder <- active_targets[[tname]]
    potential_mask <- if (!is.null(restrictions)) (1 - restrictions[pixel_avail, tname]) else rep(1, length(pixel_avail))
    
    if (type == "pop") {
      # Use current allocation weights (fallback to potential_mask if empty)
      weights <- master_out_mu[pixel_avail, tname] * potential_mask
      if (sum(weights) < 1e-12) {
        weights <- potential_mask
      }
      total_weight <- sum(weights)
      if (total_weight > 1e-12) {
        adjustment <- weights * remainder / total_weight
        phase_mu[, tname] <- adjustment
        allocated_amount <- sum(adjustment)
        missing_amount <- remainder - allocated_amount
        solver_runs[[tname]] <- sprintf("Mopped up pop count: %.4f", remainder)
        if (isTRUE(options$verbose)) {
          cat(sprintf("     └─ %s: remainder=%.4f | allocated=%.4f | missing=%.4f\n", 
                      tname, remainder, allocated_amount, missing_amount))
        }
      } else {
        if (isTRUE(options$verbose)) {
          cat(sprintf("     └─ %s: remainder=%.4f | skip (no allowed pixels)\n", tname, remainder))
        }
      }
    } else {
      # For Land Use (lu)
      available_and_allowed <- pmax(active_areas * potential_mask, 0)
      total_land_left <- sum(available_and_allowed)
      if (total_land_left > 1e-12) {
        adjustment <- (available_and_allowed) * min(remainder, total_land_left)/total_land_left
        phase_mu[, tname] <- adjustment
        active_areas <- pmax(active_areas - adjustment, 0)
        allocated_amount <- sum(adjustment)
        missing_amount <- remainder - allocated_amount
        solver_runs[[tname]] <- sprintf("Mopped up area: %.4f", remainder)
        if (isTRUE(options$verbose)) {
          cat(sprintf("     └─ %s: remainder=%.4f | allocated=%.4f | missing=%.4f\n", 
                      tname, remainder, allocated_amount, missing_amount))
        }
      } else {
        if (isTRUE(options$verbose)) {
          cat(sprintf("     └─ %s: remainder=%.4f | skip (no available land area)\n", tname, remainder))
        }
      }
    }
  }
  
  new_out_mu <- master_out_mu + phase_mu
  new_deviation_vec <- calculate_target_deviation(colSums(new_out_mu), colSums(master_out_mu) + active_targets, options)
  
  list(
    phase_mu = phase_mu,
    out_mu = new_out_mu,
    solution = curr_sol,
    converged = !any(new_deviation_vec > options$max_diff),
    runs = list(
      converged = !any(new_deviation_vec > options$max_diff),
      max_dev = max(new_deviation_vec),
      failing_targets = failing_targets
    )
  )
}

# ==============================================================================
# 5. MAIN ORCHESTRATOR
# ==============================================================================

#' @keywords internal
#' @noRd
optim_allocation <- function(targets, priors, restrictions = NULL, pixel_avail, 
                             type = c("lu", "pop"), areas = NULL, grid = NULL, 
                             current_lu_from = NULL, options = list()) {
  type <- match.arg(type)
  pixel_avail <- as.character(pixel_avail)
  
  # Map population control option names if present
  if (!is.null(options$use_nloptr) && is.null(options$use_multi_nloptr)) {
    options$use_multi_nloptr <- options$use_nloptr
  }
  if (!is.null(options$use_sequential_nloptr) && is.null(options$use_single_nloptr)) {
    options$use_single_nloptr <- options$use_sequential_nloptr
  }
  
  solver_runs <- list()
  # ---- 1. SACRED BASELINES ----
  # These are never altered. They are the reference points for everything.
  original_targets <- targets 
  original_areas <- areas
  
  # ---- 2. PHASE 0: PRE-ALLOCATION ----
  current_state <- update_residual_state(
    original_targets, current_mu = NULL, areas = original_areas, 
    restrictions = restrictions, pixel_avail = pixel_avail, type = type
  )
  
  master_out_mu <- current_state$already_allocated # The absolute cumulative matrix
  active_targets <- current_state$targets          # The residual targets for Phase 1/2
  active_areas <- current_state$areas              # The residual land for Phase 1/2
  
  if (all(active_targets == 0)) {
    return(list(solution = rep(0, length(targets)), out_mu = master_out_mu, 
                solver_runs = list(main = "Early return"), summary = summarize_allocations(master_out_mu, original_targets)))
  }
  
  # ---- 3. SETUP SOLVERS ----
  if (type == "lu") {
    alloc_func <- allocate_mnl
    alloc_args <- list(areas = active_areas, # Solvers strictly see residual land
                       restrictions = restrictions[pixel_avail, , drop = FALSE], 
                       cutoff = options$cutoff, area_cutoff = options$area_cutoff)
    nloptr_extra_args <- list(areas = active_areas, area_cutoff = options$area_cutoff)
  } else {
    alloc_func <- allocate_poisson
    
    # Subset count bounds if they are spatial vectors (length > 1)
    c_lb <- if (length(options$count_lb_cutoff) > 1) options$count_lb_cutoff[pixel_avail] else options$count_lb_cutoff
    c_ub <- if (length(options$count_ub_cutoff) > 1) options$count_ub_cutoff[pixel_avail] else options$count_ub_cutoff
    
    alloc_args <- list(restrictions = restrictions[pixel_avail, , drop = FALSE], 
                       cutoff = options$cutoff, count_lb = c_lb, count_ub = c_ub)
    nloptr_extra_args <- list(count_lb = c_lb, count_ub = c_ub)
    if (!options$optimize_log_scale) alloc_args$grid <- grid[pixel_avail] 
  }
  
  # Adaptive tempering loop state
  tempering_current <- options$tempering_start
  tempering_step <- options$tempering_step
  max_tempering <- options$max_tempering
  
  converged <- FALSE
  phase_1_2_mu <- matrix(0, nrow = length(pixel_avail), ncol = length(active_targets), dimnames = list(pixel_avail, names(active_targets)))
  curr_sol <- rep(0, length(active_targets))
  
  # ---- 4 & 5. INITIALIZE & SIMULTANEOUS SOLVERS (ADAPTIVE TEMPERING LOOP) ----
  while (!converged && tempering_current <= max_tempering) {
    # Inject current tempering into args
    alloc_args$tempering <- tempering_current
    nloptr_extra_args$tempering <- tempering_current
    
    n_targs <- length(active_targets)
    if (options$optimize_log_scale) {
      log_priors <- log(pmax(priors[pixel_avail, , drop = FALSE], 1e-12))
      log_mu_init <- do.call(alloc_func, c(list(x = rep(0, n_targs), log_scale = TRUE, log_mu0 = log_priors), alloc_args))
      log_denom <- matrixStats::colLogSumExps(log_mu_init)
      x0 <- log(pmax(active_targets, 1e-12)) - log_denom
    } else {
      log_priors <- NULL
      mu_init <- do.call(alloc_func, c(list(x = rep(1, n_targs), log_scale = FALSE, mu0 = priors[pixel_avail, , drop = FALSE]), alloc_args))
      denom <- colSums(mu_init)
      denom[denom == 0] <- .Machine$double.eps
      x0 <- active_targets / denom
    }
    names(x0) <- names(active_targets)
    curr_sol <- x0
    
    if (isTRUE(options$use_ipf)) {
      res_ipf <- run_ipf_phase(master_out_mu, x0, active_targets, priors, pixel_avail, alloc_func, alloc_args, log_priors, options)
      curr_sol <- res_ipf$solution
      phase_1_2_mu <- res_ipf$phase_mu
      solver_runs[["ipf"]] <- res_ipf$runs
      converged <- res_ipf$converged
    }
    
    if (isTRUE(options$use_multi_nloptr) && !converged) {
      res_nloptr <- run_multi_nloptr_phase(master_out_mu, curr_sol, x0, active_targets, priors, restrictions, pixel_avail, alloc_func, alloc_args, log_priors, nloptr_extra_args, options)
      curr_sol <- res_nloptr$solution
      phase_1_2_mu <- res_nloptr$phase_mu
      solver_runs[["nloptr"]] <- res_nloptr$runs
      converged <- res_nloptr$converged
    } 
    
    if (!converged && tempering_current < max_tempering) {
      if (isTRUE(options$verbose)) {
        cat(sprintf("   ! Solver stuck at tempering %.2f. Re-annealing with %.2f\n", tempering_current, tempering_current + tempering_step))
      }
      tempering_current <- tempering_current + tempering_step
    } else {
      break
    }
  } 
  
  # ---- 6. CRITICAL MID-POINT STATE UPDATE ----
  # Safely combine Phase 0 + Phase 1/2
  master_out_mu <- master_out_mu + phase_1_2_mu
  
  if (isTRUE(options$use_intermediate_target_update)) {
    
    current_state <- update_residual_state(
      original_targets, current_mu = master_out_mu, areas = original_areas, 
      restrictions = restrictions, pixel_avail = pixel_avail, type = type
    )
    
    # Update the variables for the next phases
    master_out_mu <- current_state$already_allocated # Safely shaved down
    active_targets <- current_state$targets          # Deducted targets
    active_areas <- current_state$areas              # Deducted land
  }
  
  # ---- 7. PHASE 3: SEQUENTIAL SHAVING ----
  if (isTRUE(options$use_single_nloptr) && !converged && any(active_targets > options$max_diff)) {
    res_boost <- run_single_nloptr_phase(
      master_out_mu = master_out_mu, curr_sol = curr_sol, x0 = x0, 
      active_targets = active_targets, active_areas = active_areas, # Pass clean residuals!
      priors = priors, restrictions = restrictions, pixel_avail = pixel_avail, type = type,
      alloc_func = alloc_func, alloc_args = alloc_args, nloptr_extra_args = nloptr_extra_args, 
      log_priors = log_priors, options = options, res_main_final = solver_runs[["nloptr"]]$adaptive_summary
    )
    master_out_mu <- res_boost$out_mu
    curr_sol <- res_boost$solution
    solver_runs[["boost"]] <- res_boost$runs
    converged <- res_boost$converged
  }
  
  # ---- 8. PHASE 4: GLOBAL MOP-UP ----
  if (isTRUE(options$use_mopup) && !converged) {
    # Final state check before brute force
    final_state <- update_residual_state(
      original_targets, current_mu = master_out_mu, areas = original_areas, 
      restrictions = restrictions, pixel_avail = pixel_avail, type = type
    )
    
    if (any(final_state$targets > 1e-9)) {
      res_mopup <- run_mopup_phase(
        master_out_mu = final_state$already_allocated, curr_sol = curr_sol,
        active_targets = final_state$targets, active_areas = final_state$areas, 
        restrictions = restrictions, pixel_avail = pixel_avail, type = type, options = options
      )
      master_out_mu <- res_mopup$out_mu
      solver_runs[["mopup"]] <- res_mopup$runs
    }
  }
  
  # ---- 9. FINAL NORMALIZATION & SUMMARY ----
  if (type == "lu") {
    # master_out_mu <- normalize_overshoot(master_out_mu, original_areas, pixel_avail)
  }
  
  solver_runs[["final_tempering"]] <- tempering_current
  
  return(list(
    solution = curr_sol, 
    out_mu = master_out_mu, 
    solver_runs = solver_runs,
    summary = summarize_allocations(master_out_mu, original_targets) # FIXED: Use original targets
  ))
}

# ==============================================================================
# 6. WRAPPER EXPORTS
# ==============================================================================

#' @keywords internal
#' @noRd
optim_lu_allocation <- function(...) optim_allocation(..., type = "lu")
#' @keywords internal
#' @noRd
optim_pop_allocation <- function(...) optim_allocation(..., type = "pop")
