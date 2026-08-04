# ===============================================================
# optim_pop_allocation.R
# ===============================================================
optim_pop_allocation <-     function(targets = NULL,
                                     priors = NULL,
                                     grid = NULL,
                                     restrictions = NULL,
                                     pixel_avail = NULL,
                                     current_pop = NULL,
                                     current_lu_to = NULL,
                                     options = list()) {
  pixel_avail <- as.character(pixel_avail)
  # ---- 0. Helper: build base nloptr options ----
  make_nloptr_opts <- function(options) {
    list(
      algorithm = options$algorithm,
      objective = get(options$obj_fun, mode = "function", inherits = TRUE),
      gradient  = if (!is.null(options$grad_obj_fun))
        get(options$grad_obj_fun, mode = "function", inherits = TRUE) else NULL,
      ineqc     = if (!is.null(options$ineqc_fun) && options$algorithm == "NLOPT_LN_COBYLA")
        get(options$ineqc_fun, mode = "function", inherits = TRUE) else NULL,
      xtol_rel = options$xtol_rel,
      xtol_abs = options$xtol_abs,
      maxeval  = options$max_eval,
      stopval  = options$max_diff
    )
  }
  
  # ---- 1. Helper: adaptive nloptr loop ----
  adaptive_nloptr <- function(initial_solution, x0, lb, ub, opts_base,
                              priors_sub, grid_sub, targets_sub,
                              restrictions_sub) {
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
    
    while ((current_try <= max_tries) && (res$objective > options$max_diff)) {
      res$options$maxeval <- max(options$min_eval,
                                 ceiling(options$max_eval * (current_try / max_tries)))
      res$lower_bounds <- lb / (10^(current_try * 10))
      res$upper_bounds <- ub * 100 * current_try
      
      names(res$lower_bounds) <- names(x0)
      names(res$upper_bounds) <- names(x0)
      
      nloptr_out <- nloptr::nloptr(
        x0 = res$solution,
        eval_f = opts_base$objective,
        eval_grad_f = opts_base$gradient,
        eval_g_ineq = opts_base$ineqc,
        lb = res$lower_bounds,
        ub = res$upper_bounds,
        opts = res$options,
        mu0 = priors_sub,
        targets = targets_sub,
        restrictions = restrictions_sub,
        cutoff = options$cutoff,
        count_lb = options$count_lb_cutoff,
        count_ub = options$count_ub_cutoff,
        lambda = options$pen_lambda,
        log_scale = options$optimize_log_scale,
        tempering = 1.0
      )
      
      nloptr_out$message <- paste0(sprintf("Adaptive try %d\n", current_try), nloptr_out$message)
      runs[[paste0("try_", current_try)]] <- nloptr_out
      
      res$solution  <- nloptr_out$solution
      res$objective <- if (!is.null(nloptr_out$objective)) nloptr_out$objective else res$objective
      res$message   <- nloptr_out$message
      
      current_try <- current_try + 1
      names(res$solution) <- names(x0)
    }
    
    list(final = res, runs = runs)
  }
  
  # ---- 2. Helper: single-target fallback (boost) ----
  boost_target <- function(tname, error_target_val, res_grid, res_main_bounds, opts_base) {
    x0_single <- x0[[tname]]
    lb_single <- res_main_bounds$lower_bounds[[tname]]
    ub_single <- res_main_bounds$upper_bounds[[tname]]
    
    nloptr_out <- nloptr::nloptr(
      x0 = x0_single,
      eval_f = opts_base$objective,
      eval_grad_f = opts_base$gradient,
      eval_g_ineq = opts_base$ineqc,
      lb = lb_single,
      ub = ub_single,
      opts = res_main_bounds$options,
      mu0 = priors[pixel_avail, tname, drop = FALSE],
      targets = error_target_val,
      restrictions = if (!is.null(restrictions))
        restrictions[pixel_avail, tname, drop = FALSE] else NULL,
      cutoff = options$cutoff,
      count_lb = options$count_lb_cutoff,
      count_ub = options$count_ub_cutoff,
      lambda = options$pen_lambda,
      log_scale = options$optimize_log_scale,
      tempering = 1.0
    )
    
    nloptr_out$message <- paste0("Boosted target: ", tname, "\n", nloptr_out$message)
    nloptr_out$note <- "Individual logit boosted"
    nloptr_out
  }
  
  # ---- 3. Early exit: all-zero targets ----
  solver_runs <- list()
  
  if (all(targets == 0)) {
    out_mu <- matrix(0, nrow = length(pixel_avail), ncol = length(targets),
                     dimnames = list(pixel_avail, names(targets)))
    solver_runs[[curr_lu]] <- res_main$message = "All targets zero: early return"
    
    return(list(
      solution = rep(0, length(targets)),
      out_mu = out_mu,
      solver_runs = solver_runs,
    ))
  }
  
  # ---- 4. Initial guess x0 ----
  denom <- colSums(allocate_poisson(
    1,
    priors[pixel_avail, , drop = FALSE],
    restrictions[pixel_avail, , drop = FALSE],
    options$cutoff,
    count_lb = options$count_lb_cutoff,
    count_ub = options$count_ub_cutoff
  ))
  denom[denom == 0] <- .Machine$double.eps
  x0 <- targets / denom
  names(x0) <- names(targets)
  
  # ---- 5. Quick IPF solver ----
  res_main <- count_ipf_fit(
    x0 = x0,
    mu0 = priors[pixel_avail, , drop = FALSE],
    targets = targets,
    restrictions = restrictions[pixel_avail, , drop = FALSE],
    cutoff = options$cutoff,
    count_lb = options$count_lb_cutoff,
    count_ub = options$count_ub_cutoff,
    tolerance = options$max_diff,
    max_iter = if(options$only_ipf_optim) options$max_eval else options$min_eval,
    log_scale = options$optimize_log_scale,
    verbose = FALSE
  )
  
  if ((res_main$converged && !any(res_main$final_deviation > options$max_diff)) || options$only_ipf_optim) {
    res_main$message <- "IPF optimization."
    solver_runs[["main"]] <- res_main
    out_mu <- allocate_poisson(
      res_main$solution,
      priors[pixel_avail, , drop = FALSE],
      restrictions[pixel_avail, , drop = FALSE],
      options$cutoff,
      count_lb = options$count_lb_cutoff,
      count_ub = options$count_ub_cutoff
    )
    
    return(list(
      solution = res_main$solution,
      out_mu = out_mu,
      solver_runs = solver_runs
    ))
  }
  
  # ---- 6. Fallback to adaptive nloptr ----
  opts_base <- make_nloptr_opts(options)
  lb <- pmin(res_main$solution,options$MIN_EXP) / 10
  ub <- pmax(res_main$solution,x0) * 10
  names(lb) <- names(ub) <- names(x0)
  
  adaptive_res <- adaptive_nloptr(
    initial_solution = res_main$solution,
    x0 = x0,
    lb = lb,
    ub = ub,
    opts_base = opts_base,
    priors_sub = priors[pixel_avail, , drop = FALSE],
    grid_sub = grid[pixel_avail],
    targets_sub = targets,
    restrictions_sub = restrictions[pixel_avail, , drop = FALSE]
  )
  
  res_main_final <- adaptive_res$final
  solver_runs[["main"]] <- c(list(adaptive_summary = res_main_final), adaptive_res$runs)
  
  out_mu <- allocate_poisson(res_main_final$solution,
                             priors[pixel_avail, , drop = FALSE],
                             restrictions[pixel_avail, , drop = FALSE],
                             options$cutoff,
                             count_lb = options$count_lb_cutoff,
                             count_ub = options$count_ub_cutoff
                              )
  
  # ---- 7. Optional per-target boost ----
  deviation_sq <- (colSums(out_mu) - targets)^2
  error_ind <- deviation_sq > options$max_diff
  
  if (any(error_ind)) {
    boost_opts_base <- make_nloptr_opts(options)
    res_grid <- grid[pixel_avail] - rowSums(out_mu[pixel_avail, !error_ind, drop = FALSE])
    res_grid <- as.matrix(res_grid)
    res_main_bounds <- res_main_final
    
    for (tname in names(targets[error_ind])) {
      res_single <- boost_target(
        tname = tname,
        error_target_val = targets[[tname]],
        res_grid = res_grid,
        res_main_bounds = res_main_bounds,
        opts_base = boost_opts_base
      )
      
      solver_runs[[paste0("boost_", tname)]] <- res_single
      
      curr_out_mu <- allocate_poisson(res_single$solution,
                                      priors[pixel_avail, tname, drop = FALSE],
                                      if (!is.null(restrictions)) restrictions[pixel_avail, tname, drop = FALSE] else NULL,
                                      options$cutoff,
                                      count_lb = options$count_lb_cutoff,
                                      count_ub = options$count_ub_cutoff,
                                      log_scale = options$optimize_log_scale
                                      )
      
      out_mu[pixel_avail, tname] <- curr_out_mu[pixel_avail, tname]
    }
  }
  
  # ---- 8. Return unified structure ----
  res_main_final$message <- paste("Finished (with nloptr fallback and optional boosts)",
                                  res_main_final$message)
  
  list(
    solution = res_main_final$solution,
    out_mu = out_mu,
    solver_runs = solver_runs,
    message = res_main_final$message
  )
}


count_ipf_fit <- function(x0 = NULL, mu0, targets, restrictions = NULL,
                          cutoff = 0, count_lb = 0, count_ub = Inf,
                          max_iter = 1000, tolerance = 1e-4,
                          damping = 0.5, verbose = FALSE,
                          log_scale = TRUE) {
  
  # ---- 1. Initialization ----
  n_targets <- length(targets)
  
  if (!is.null(x0)) {
    if (length(x0) != n_targets) {
      warning("Initial guess x0 length does not match targets length. Defaulting to x=1.")
      x <- rep(1, n_targets)
    } else {
      x <- pmax(x0, 1e-12)
    }
  } else {
    x <- rep(1, n_targets)
  }
  
  
  current_iter <- 0
  error_change <- Inf
  
  if (verbose) {
    message("Starting IPF balancing ", ifelse(log_scale, "(log-space)...", "(linear-space)..."))
  }
  
  if(log_scale) log_targ <- log(pmax(targets, 1e-12))
  
  # ---- 2. Iteration Loop ----
  while (current_iter < max_iter && error_change > tolerance) {
    if (log_scale) {
      mu_current <- allocate_poisson(x, mu0, restrictions, cutoff, count_lb, count_ub, log_scale=TRUE)
      # max_c <- apply(mu_current,2,max) 
      # log_pred <- max_c + log(colSums(exp(sweep(mu_current,2,max_c,"-"))))
      log_pred <- matrixStats::colLogSumExps(mu_current)
      # Log-space update (with damping)
      delta_logx <- damping * (log_targ - log_pred)
      x_new <- x * exp(delta_logx)
      # Convergence metric: max relative log-change
      error_change <- max(abs(delta_logx))
      x <- x_new
      
    } else {
      mu_current <- allocate_poisson(x, mu0, grid, restrictions, cutoff, count_lb, count_ub, log_scale=FALSE)
      
      # Linear-space update
      pred <- colSums(mu_current)
      ratio <- targets / pmax(pred, 1e-12)  # avoid division by zero
      x_new <- x * (ratio ^ damping)
      
      # Convergence metric: max relative change
      error_change <- max(abs(log(x_new / pmax(x, 1e-12))))
      x <- pmax(x_new, 1e-12)  # avoid negative/zero values
    }
    
    current_iter <- current_iter + 1
    
    if (verbose && current_iter %% 50 == 0) {
      message(sprintf("Iter %d: Max change = %.6f", current_iter, error_change))
    }
  }
  
  # ---- 3. Finalization ----
  converged <- (error_change <= tolerance)
  
  if (!converged && verbose) {
    warning(sprintf("IPF did not converge within %d iterations (final error=%.6f)", 
                    max_iter, error_change))
  } else if (verbose) {
    message(sprintf("Converged in %d iterations (max change: %.6f).", 
                    current_iter, error_change))
  }
  
  # Final squared differences from targets
  final_sqr_diff <- if(log_scale) (exp(log_pred)-targets)^2 else (pred-targets)^2
  
  return(list(
    solution = x,
    iterations = current_iter,
    converged = converged,
    final_deviation = final_sqr_diff
  ))
}
