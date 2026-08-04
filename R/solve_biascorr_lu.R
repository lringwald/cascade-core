# ===============================================================
# solve_biascorr_lu.R
# ===============================================================
library(data.table)
library(nloptr)
# Main solver: fully data.table optimized
solve_biascorr_lu     <- function(targets,
                               areas,
                               xmat,
                               betas, 
                               priors = NULL, 
                               restrictions = NULL, 
                               timestep = NULL,
                               type = NULL,
                               options = downscale_control()
) {
  # Convert inputs to data.table
  targets      <- if(!is.data.table(targets)) as.data.table(targets) else targets
  areas        <- if(!is.data.table(areas)) as.data.table(areas) else areas
  xmat         <- if (!is.null(xmat) && !is.data.table(xmat)) as.data.table(xmat) else xmat
  betas        <- if (!is.null(betas) && !is.data.table(betas)) as.data.table(betas) else betas
  priors       <- if (!is.null(priors) && !is.data.table(priors)) as.data.table(priors) else priors
  restrictions <- if (!is.null(restrictions) && !is.data.table(restrictions)) as.data.table(restrictions) else restrictions
  out_solver <- list()
  full_out_res <- list()
  # Precompute arrays
  ks_avail <- unique(betas$ks)
  
  target_lu_from_levels <- as.character(unique(targets[J(timestep), from_class]))
  
  pb_inner <- cli_progress_bar(
    total = length(target_lu_from_levels),
    format = "   ├─ {pb_spin} {curr_lu} {pb_bar} {pb_percent} [{pb_current}/{pb_total}]",
    format_done = "   ├─ {pb_spin} {timestep} optimizations complete.",
    clear = FALSE
  )
  for(curr_lu in target_lu_from_levels) {
    if (is.na(curr_lu)) next
    err_txt <- paste0(curr_lu, " ", options$err.txt)
    # if(curr_lu %in% c("Cropland.HI.Whea") && as.character(timestep) %in% c("2050")) 
    
    out_solver_runs <- list()
    
    # Extract targets
    curr_targets_vec <- targets[J(timestep, curr_lu), !c("times", "from_class")][, setNames(value, to_class)]
    curr_lu_to <- names(curr_targets_vec)
    
    # Extract betas
    if (!is.null(betas) && nrow(betas) > 0) {
      b_lu_from <- intersect(c(curr_lu, PLCHOLD_LU), betas$from_class)
      curr_betas_dt <- betas[from_class %in% b_lu_from]
      
      if (nrow(curr_betas_dt) > 0) {
        if (PLCHOLD_LU %in% curr_betas_dt$to_class) {
          any_to <- curr_betas_dt[to_class == PLCHOLD_LU]
          other_to <- curr_betas_dt[to_class != PLCHOLD_LU]
          if (nrow(any_to) > 0) {
            expanded_to <- any_to[rep(1:.N, each = length(curr_lu_to))]
            expanded_to[, to_class := rep(curr_lu_to, times = nrow(any_to))]
            curr_betas_dt <- rbind(other_to, expanded_to)
          }
        }
        
        curr_betas_dt <- curr_betas_dt[to_class %in% curr_lu_to]
        curr_betas_dt <- curr_betas_dt[, .(value = max(value, na.rm = TRUE)), by = .(to_class, ks)]
      }
      
      if (nrow(curr_betas_dt) > 0) {
        # Pivot to matrix: rows = ks, cols = to_class
        tmp_betas <- dcast(curr_betas_dt, ks ~ to_class, value.var = "value", fill = 0)
        curr_betas <- as.matrix(tmp_betas[, -1, with=FALSE])
        rownames(curr_betas) <- tmp_betas$ks
        ks_avail <- rownames(curr_betas)
      } else {
        curr_betas <- NULL
        ks_avail <- NULL
      }
      curr_betas_dt <- NULL
    } else {
      curr_betas <- NULL
      ks_avail <- NULL
    }
    
    
    # Extract areas
    curr_areas_vec <- areas[.(curr_lu), on = "from_class", nomatch = 0L][, setNames(value, ns)]
    curr_areas_vec <- curr_areas_vec[curr_areas_vec>0]
    
    ns_avail <- sort(unique(names(curr_areas_vec)))
    if (length(ns_avail) == 0) {
      cli_progress_update(id = pb_inner)
      next
    }
    
    curr_areas_vec <- curr_areas_vec[ns_avail]
    
    # Extract xmat
    if (!is.null(ks_avail) && !is.null(xmat)) {
      # Subset xmat to ns_avail and ks_avail
      curr_xmat_dt <- xmat[CJ(ns=ns_avail, ks=ks_avail), on = .(ns, ks), nomatch = 0L]
      # Convert to matrix: rows = ns, columns = ks
      curr_xmat <- dcast(curr_xmat_dt, ns ~ ks, value.var = "value")
      curr_xmat_dt <- NULL
      
      # Ensure consistent ordering
      setorder(curr_xmat, ns)
      
      # Convert to numeric matrix
      curr_xmat <- as.matrix(curr_xmat, rownames="ns")
      
      # Ensure proper row order and keep rownames
      curr_xmat <- curr_xmat[match(ns_avail, rownames(curr_xmat)), ks_avail, drop = FALSE]
      curr_xmat[is.na(curr_xmat)] <- 0 # Handle missing combinations as zero
      rownames(curr_xmat) <- ns_avail
    } else {
      curr_xmat <- NULL
    }
    
    
    # Priors / Restrictions selection as matrix (Optimized)
    if(!is.null(priors)) {
      # Use key-based filtering (assuming ns, times, from_class are keys)
      p_times <- intersect(c(timestep, PLCHOLD_T), priors$times)
      p_lu_from <- intersect(c(curr_lu, PLCHOLD_LU), priors$from_class)
      
      curr_priors_dt <- priors[times %in% p_times & ns %in% ns_avail & from_class %in% p_lu_from]
      
      if (nrow(curr_priors_dt) > 0) {
        if (PLCHOLD_T %in% curr_priors_dt$times) {
           curr_priors_dt[times == PLCHOLD_T, times := timestep]
        }
        if (PLCHOLD_LU %in% curr_priors_dt$to_class) {
           any_to <- curr_priors_dt[to_class == PLCHOLD_LU]
           other_to <- curr_priors_dt[to_class != PLCHOLD_LU]
           if (nrow(any_to) > 0) {
             expanded_to <- any_to[rep(1:.N, each = length(curr_lu_to))]
             expanded_to[, to_class := rep(curr_lu_to, times = nrow(any_to))]
             curr_priors_dt <- rbind(other_to, expanded_to)
           }
        }
        
        curr_priors_dt <- curr_priors_dt[to_class %in% curr_lu_to]
        curr_priors_dt <- curr_priors_dt[, .(value = max(value, na.rm=TRUE), weight = max(weight, na.rm=TRUE)), by = .(ns, to_class)]
        
        # Create array with a single dcast per measure if possible, or pivot once
        curr_priors_arr <- array(0, dim = c(length(ns_avail), length(curr_lu_to), 2),
                                 dimnames = list(ns = ns_avail, to_class = curr_lu_to, measure = c("value","weight")))
        
        # Melt and dcast once to fill if multiple measures, or just dcast twice (still faster than nested joins)
        for (m in c("value","weight")) {
          tmp <- dcast(curr_priors_dt, ns ~ to_class, value.var = m, fill = 0)
          setorder(tmp, ns)
          mat_tmp <- as.matrix(tmp, rownames = "ns")
          # Handle potential missing columns in tmp if some lu.to aren't in priors
          common_lu_to <- intersect(colnames(mat_tmp), curr_lu_to)
          common_ns_p  <- intersect(rownames(mat_tmp), ns_avail)
          curr_priors_arr[common_ns_p, common_lu_to, m] <- mat_tmp[common_ns_p, common_lu_to, drop = FALSE]
        }
      } else {
        curr_priors_arr <- NULL
      }
      curr_priors_dt <- NULL  
    } else {      
      curr_priors_arr <- NULL
    }
    
    # Restrictions selection as matrix (Optimized)
    if (!is.null(restrictions)) {
      r_times <- intersect(c(timestep, PLCHOLD_T), restrictions$times)
      r_lu_from <- intersect(c(curr_lu, PLCHOLD_LU), restrictions$from_class)
      
      curr_restr_mat_dt <- restrictions[times %in% r_times & ns %in% ns_avail & from_class %in% r_lu_from]
      
      if (nrow(curr_restr_mat_dt) > 0) {
        if (PLCHOLD_T %in% curr_restr_mat_dt$times) {
           curr_restr_mat_dt[times == PLCHOLD_T, times := timestep]
        }
        if (PLCHOLD_LU %in% curr_restr_mat_dt$to_class) {
           any_to <- curr_restr_mat_dt[to_class == PLCHOLD_LU]
           other_to <- curr_restr_mat_dt[to_class != PLCHOLD_LU]
           if (nrow(any_to) > 0) {
             expanded_to <- any_to[rep(1:.N, each = length(curr_lu_to))]
             expanded_to[, to_class := rep(curr_lu_to, times = nrow(any_to))]
             curr_restr_mat_dt <- rbind(other_to, expanded_to)
           }
        }
        
        curr_restr_mat_dt <- curr_restr_mat_dt[to_class %in% curr_lu_to]
        curr_restr_mat_dt <- curr_restr_mat_dt[, .(value = max(value, na.rm=TRUE)), by = .(ns, to_class)]
        
        tmp_r <- dcast(curr_restr_mat_dt, ns ~ to_class, value.var = "value", fill = 0)
        setorder(tmp_r, ns)
        mat_r <- as.matrix(tmp_r, rownames = "ns")
        common_r_lu_to <- intersect(colnames(mat_r), curr_lu_to)
        
        curr_restr_mat <- matrix(0, length(ns_avail), length(curr_lu_to), dimnames = list(ns_avail, curr_lu_to))
        common_ns <- intersect(rownames(mat_r), ns_avail)
        curr_restr_mat[common_ns, common_r_lu_to] <- mat_r[common_ns, common_r_lu_to, drop = FALSE]
      } else {
        curr_restr_mat <- NULL
      }
      curr_restr_mat_dt <- NULL
    } else {
      curr_restr_mat <- NULL
    }
    
    n <- length(curr_areas_vec)
    p <- length(curr_targets_vec)
    
    # Compute econometric priors
    curr_priors <- matrix(0, n, p)
    rownames(curr_priors) <- ns_avail
    colnames(curr_priors) <- curr_lu_to
    
    mnl_priors_cols <- colnames(curr_priors)[colnames(curr_priors) %in% colnames(curr_betas)]
    exo_priors_cols <- if (!is.null(curr_priors_arr)) colnames(curr_priors)[colnames(curr_priors) %in% dimnames(curr_priors_arr)[[2]]] else character(0)
    
    p_mnl <- length(mnl_priors_cols)
    p_exo <- length(exo_priors_cols)
    
    overlap_priors_cols <- intersect(mnl_priors_cols,exo_priors_cols)
    
    
    if(p_mnl > 0 && !is.null(curr_betas)) {
      curr_priors_mu_hat <- matrix(0, n, p_mnl)
      rownames(curr_priors_mu_hat) <- ns_avail
      colnames(curr_priors_mu_hat) <- mnl_priors_cols
      
      # Subset curr_betas to matched columns if necessary
      curr_betas_sub <- curr_betas[, colnames(curr_priors_mu_hat), drop = FALSE]
      curr_priors_mu_hat[] <- curr_xmat[ns_avail, , drop = FALSE] %*% curr_betas_sub
      
      curr_priors_mu_hat[] <- pmin(curr_priors_mu_hat, options$MAX_EXP)
      curr_priors_mu_hat[] <- pmax(curr_priors_mu_hat, log(options$MIN_EXP))
      
      # Exponentiate MNL utilities to get competitive mass (e^U)
      scale_mnl_utility_exp <- if (length(curr_priors_mu_hat) > 0) mean(exp(curr_priors_mu_hat), na.rm = TRUE) else 1
    } else{
      # Emtpy matrix avoids NULL subsetting errors downstream
      curr_priors_mu_hat <- matrix(0, n, 0)
      rownames(curr_priors_mu_hat) <- ns_avail
      colnames(curr_priors_mu_hat) <- character(0)
      
      # Default competitive mass to 1 (neutral prior) if no MNL transitions exist
      scale_mnl_utility_exp <- 1 
    }
    
    if(p_exo>0) {
      curr_priors_mu_exo <- matrix(0, n, p_exo)
      rownames(curr_priors_mu_exo) <- ns_avail
      colnames(curr_priors_mu_exo) <- exo_priors_cols
      curr_priors_mu_exo[] <- log(slice_to_matrix(curr_priors_arr, 1)[ns_avail, exo_priors_cols, drop = FALSE]+1e-12)
      curr_priors_mu_exo[] <- pmin(curr_priors_mu_exo, options$MAX_EXP)
      curr_priors_mu_exo[] <- pmax(curr_priors_mu_exo, log(options$MIN_EXP))
      #Calculate Exo Probability Scale (Average P_Exo per transition)
      scale_exo_prob_vector <- if (length(curr_priors_mu_exo) > 0) colMeans(exp(curr_priors_mu_exo), na.rm = TRUE) else 1
      # Handle potential division by zero
      scale_exo_prob_vector[scale_exo_prob_vector == 0] <- 1e-12 
    } else{
      curr_priors_mu_exo <- NULL
      scale_exo_prob_vector <- 1
    }
    # Determine Lambda (Vector)
    # Lambda_j = Scale_MNL_eU / Scale_Exo_j
    # Calculate Log-Bridge Term (ln(lambda))
    log_lambda_int_factor <- log(scale_mnl_utility_exp / scale_exo_prob_vector)
    
    
    # --- 4. Integration and Final Normalization ---
    
    # 4a. Apply Log-Bridge Term to Exogenous Utilities
    # U_Exo_Final = U_Exo_raw + log(lambda)
    # This scales the Exo utility to match the competitive scale of the MNL utility.
    if (p_exo > 0) {
      curr_priors_mu_exo <- sweep(
        curr_priors_mu_exo, # U_Exo_raw = log(P_Exo)
        2, 
        log_lambda_int_factor, 
        "+"
      )
    }
    
    
    
    if(length(overlap_priors_cols)>0){
      
      curr_priors[,overlap_priors_cols] <-
        combine_utilities(utilitites_mnl=curr_priors_mu_hat[ns_avail, overlap_priors_cols , drop = FALSE],
                          utilitites_exo=curr_priors_mu_exo[ns_avail, overlap_priors_cols, drop = FALSE],
                          weights=slice_to_matrix(curr_priors_arr, 2)[ns_avail, overlap_priors_cols, drop = FALSE])
      
    }
    
    curr_priors[ns_avail, setdiff(colnames(curr_priors_mu_hat),overlap_priors_cols)] <- curr_priors_mu_hat[ns_avail ,setdiff(colnames(curr_priors_mu_hat),overlap_priors_cols)]
    curr_priors[ns_avail, setdiff(colnames(curr_priors_mu_exo),overlap_priors_cols)] <- curr_priors_mu_exo[ns_avail ,setdiff(colnames(curr_priors_mu_exo),overlap_priors_cols)]
    
    # Exponentiation
    curr_priors <- exp(curr_priors) # e^U_Total
    # Final Renormalization
    # Denominator = 1 (for baseline e^U_Base) + sum(all e^U_Total)
    # Calculate final transition probabilities: P_j = e^U_j / D_Total
    curr_priors[is.na(curr_priors)] <- 0
    curr_priors <- sweep(curr_priors,
                         1, 
                         rowSums(cbind(rep(1,nrow(curr_priors)),
                                       curr_priors)), 
                         "/")    
    curr_priors[is.na(curr_priors)] <- 0
    
    # --- 5. Feasibility Checks (Diagnostics) ---
    total_targets <- sum(curr_targets_vec)
    total_area    <- sum(curr_areas_vec)
    rebalance_active <- isTRUE(options$rebalance_targets) && (total_targets > total_area) && (total_targets > 0)
    
    # Global insufficient area flag
    insufficient_area <- (total_targets > total_area + 1e-9)
    
    # Individual restriction feasibility check
    # Check for each target if the unblocked area is enough
    restriction_blocking <- logical(length(curr_targets_vec))
    names(restriction_blocking) <- names(curr_targets_vec)
    
    if (!is.null(curr_restr_mat)) {
      for (tname in names(curr_targets_vec)) {
        allowed_area <- sum(curr_areas_vec * (1 - curr_restr_mat[, tname]))
        if (allowed_area < curr_targets_vec[tname] - 1e-9) {
          restriction_blocking[tname] <- TRUE
        }
      }
    }
    
    if (rebalance_active) {
      # Scale targets down proportionally to match total_area
      curr_targets_vec <- curr_targets_vec * (total_area / total_targets) * (1 - 1e-9)
      cli_alert_warning("Targets (sum = {round(total_targets,2)}) exceed total area ({round(total_area,2)}) for '{curr_lu}' at timestep {timestep}. Rebalancing proportionally.")
    }

    # --- 6. Optimization ---
    if (isTRUE(options$verbose)) cat(sprintf("   Optimizing %s...\n", curr_lu))
    t_start <- if (isTRUE(options$track_performance)) Sys.time() else NULL
    res <- optim_lu_allocation(targets = curr_targets_vec,
                               areas = curr_areas_vec,
                               priors = curr_priors,
                               restrictions = curr_restr_mat,
                               pixel_avail = ns_avail,
                               current_lu_from = curr_lu,
                               options = options)
    t_end <- if (isTRUE(options$track_performance)) Sys.time() else NULL
    duration_s <- if (!is.null(t_start)) as.numeric(difftime(t_end, t_start, units = "secs")) else NA_real_

    # --- 7. Post-process & Diagnostics ---
    # Construct diagnostic row for each target
    diag_dt <- data.table(
      times = timestep,
      from_class = curr_lu,
      to_class = names(curr_targets_vec),
      initial_area_insufficient = insufficient_area,
      rebalancing_applied = rebalance_active,
      restriction_blocking = restriction_blocking[names(curr_targets_vec)],
      ipf_status = if (!is.null(res$solver_runs$ipf)) res$solver_runs$ipf$converged else NA,
      nloptr_status = if (!is.null(res$solver_runs$nloptr)) res$solver_runs$nloptr$converged else NA,
      mopup_status = if (!is.null(res$solver_runs$mopup)) res$solver_runs$mopup$converged else NA,
      final_deviation = calculate_target_deviation(colSums(res$out_mu), curr_targets_vec, options),
      duration_s = duration_s
    )
    
    # Competition detection: If Mopup phase was explicitly called and failed (isFALSE), 
    # and it wasn't a simple budget or individual restriction issue, it's a competition conflict.
    diag_dt[, restriction_competition := (fcoalesce(mopup_status, TRUE) == FALSE) & 
              (restriction_blocking == FALSE) & 
              (initial_area_insufficient == FALSE)]
    
    res_dt <- data.table(ns = ns_avail, lu_resid_area = curr_areas_vec[ns_avail] - rowSums(res$out_mu[ns_avail, , drop = FALSE]))
    
    setnames(res_dt,"lu_resid_area",curr_lu)
    res_dt <- cbind(res_dt, as.data.table(res$out_mu[ns_avail, , drop = FALSE]))
    res_dt <- melt(res_dt, id.vars = "ns", variable.name = "to_class", value.name = "value")
    res_dt[, from_class := curr_lu]
    res_dt[, rebalancing_applied := rebalance_active]
    setcolorder(res_dt, c("ns","from_class", "to_class", "value", "rebalancing_applied"))
    setkey(res_dt, ns, from_class, to_class)
    
    full_out_res[[curr_lu]] <- res_dt
    out_solver[[curr_lu]] <- list(runs = res$solver_runs, diag = diag_dt)
    cli_progress_update(id = pb_inner)
  }
  cli_progress_done(id = pb_inner)
  
  # Combine results and diagnostics
  final_res <- if (length(full_out_res) > 0) rbindlist(full_out_res) else data.table()
  
  # Filter out NULL diags before rbindlist
  diag_list <- lapply(out_solver, `[[`, "diag")
  diag_list <- diag_list[!sapply(diag_list, is.null)]
  final_diag <- if (length(diag_list) > 0) rbindlist(diag_list) else data.table()
  
  return(list(out_res = final_res, out_diag = final_diag, solver_runs = out_solver))
}
