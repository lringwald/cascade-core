# ===============================================================
# solve_biascorr_pop.R
# ===============================================================
library(data.table)

#' Bias Correction Solver for Population and Livestock
#'
#' @param targets data.table of targets
#' @param grid data.table of spatial grid
#' @param xmat data.table of explanatory variables
#' @param betas data.table of coefficients
#' @param priors data.table of priors
#' @param restrictions data.table of restrictions
#' @param timestep current year/time
#' @param options list of control parameters
#'
#' @return list(out_res, solver_runs, out_diag)
solve_biascorr_pop <- function(targets,
                               grid = NULL,
                               xmat = NULL,
                               betas = NULL,
                               priors = NULL,
                               restrictions = NULL,
                               timestep = NULL,
                               options = downscale_control()) {
  # --- POPULATION / LIVESTOCK LOGIC ---
    
    # 1. Preprocessing Targets
    # Expected: data.table [times, pop_type, value]
    # We need a named vector for the current timestep
    if (!is.data.table(targets)) targets <- as.data.table(targets)
    curr_targets_vec <- targets[times == timestep, setNames(value, pop_type)]
    
    if (length(curr_targets_vec) == 0) {
      return(list(out_res = data.table(), solver_runs = list(), out_diag = data.table()))
    }
    
    # 2. Pixel Availability
    # ns_avail determines the spatial scope of the optimization
    ns_avail <- if (!is.null(grid)) {
      unique(grid$ns)
    } else if (!is.null(xmat)) {
      unique(xmat$ns)
    } else if (!is.null(priors)) {
      unique(priors$ns)
    } else {
      stop("solve_biascorr_unified (pop): Either 'grid', 'xmat', or 'priors' must be provided to determine pixel availability.")
    }
    
    curr_pop_types <- names(curr_targets_vec)
    n <- length(ns_avail)
    p <- length(curr_pop_types)
    
    # ---------------------------------------------------------------
    # 3. Econometric Priors from xmat/betas (mirroring LU solve_biascorr)
    # ---------------------------------------------------------------
    curr_priors_mu_hat <- NULL
    
    if (!is.null(xmat) && !is.null(betas) && is.data.table(xmat) && is.data.table(betas) && nrow(betas) > 0) {
      # 3a. Extract betas for current pop_types
      # betas expected: data.table [pop_type, ks, value, ...]
      curr_betas_dt <- betas[pop_type %in% curr_pop_types]
      
      if (nrow(curr_betas_dt) > 0) {
        # Aggregate to unique (pop_type, ks) pairs
        curr_betas_dt <- curr_betas_dt[, .(value = max(value, na.rm = TRUE)), by = .(pop_type, ks)]
        
        # Pivot to matrix: rows = ks, cols = pop_type
        tmp_betas <- dcast(curr_betas_dt, ks ~ pop_type, value.var = "value", fill = 0)
        curr_betas <- as.matrix(tmp_betas[, -1, with = FALSE])
        rownames(curr_betas) <- tmp_betas$ks
        ks_avail <- rownames(curr_betas)
        mnl_pop_cols <- colnames(curr_betas)
        
        # 3b. Extract xmat for current ns_avail and ks_avail
        # xmat expected: data.table [ns, ks, value, ...]
        curr_xmat_dt <- xmat[ns %in% ns_avail & ks %in% ks_avail]
        
        if (nrow(curr_xmat_dt) > 0) {
          # Pivot to matrix: rows = ns, cols = ks
          curr_xmat <- dcast(curr_xmat_dt, ns ~ ks, value.var = "value", fill = 0, fun.aggregate = sum)
          setorder(curr_xmat, ns)
          curr_xmat <- as.matrix(curr_xmat, rownames = "ns")
          
          # Align rows to ns_avail
          curr_xmat_aligned <- matrix(0, n, length(ks_avail), dimnames = list(ns_avail, ks_avail))
          common_ns <- intersect(rownames(curr_xmat), ns_avail)
          common_ks <- intersect(colnames(curr_xmat), ks_avail)
          curr_xmat_aligned[common_ns, common_ks] <- curr_xmat[common_ns, common_ks]
          
          # 3c. Compute mu_hat = X %*% beta (in log-space)
          curr_betas_sub <- curr_betas[ks_avail, mnl_pop_cols, drop = FALSE]
          curr_priors_mu_hat <- curr_xmat_aligned %*% curr_betas_sub
          
          # Apply MAX_EXP / MIN_EXP caps to prevent overflow
          if (!is.null(options$MAX_EXP)) {
            curr_priors_mu_hat[] <- pmin(curr_priors_mu_hat, options$MAX_EXP)
          } else {
            curr_priors_mu_hat[] <- pmin(curr_priors_mu_hat, 100)
          }
          if (!is.null(options$MIN_EXP)) {
            curr_priors_mu_hat[] <- pmax(curr_priors_mu_hat, log(options$MIN_EXP))
          } else {
            curr_priors_mu_hat[] <- pmax(curr_priors_mu_hat, -100)
          }
        }
      }
    }
    
    # ---------------------------------------------------------------
    # 4. Exogenous Priors from priors argument
    # ---------------------------------------------------------------
    curr_priors_exo <- NULL
    
    if (!is.null(priors)) {
      if (!is.data.table(priors)) priors <- as.data.table(priors)
      curr_priors_dt <- priors[times == timestep]
      
      if (nrow(curr_priors_dt) > 0) {
        # Pivot to matrix: rows = ns, cols = pop_type
        curr_priors_mat <- dcast(curr_priors_dt, ns ~ pop_type, value.var = "value", fill = 0, fun.aggregate = sum)
        curr_priors_mat <- as.matrix(curr_priors_mat, rownames = "ns")
        
        curr_priors_exo <- matrix(0, n, p, dimnames = list(ns_avail, curr_pop_types))
        common_ns <- intersect(rownames(curr_priors_mat), ns_avail)
        common_cols <- intersect(colnames(curr_priors_mat), curr_pop_types)
        curr_priors_exo[common_ns, common_cols] <- curr_priors_mat[common_ns, common_cols]
      }
    }
    
    # ---------------------------------------------------------------
    # 5. Combine Econometric + Exogenous Priors
    # ---------------------------------------------------------------
    curr_priors <- matrix(0, n, p, dimnames = list(ns_avail, curr_pop_types))
    
    has_econ <- !is.null(curr_priors_mu_hat) && ncol(curr_priors_mu_hat) > 0
    has_exo <- !is.null(curr_priors_exo)
    
    if (has_econ && has_exo) {
      # Both econometric and exogenous priors exist
      econ_cols <- intersect(colnames(curr_priors_mu_hat), curr_pop_types)
      exo_only_cols <- setdiff(curr_pop_types, econ_cols)
      
      # For pop_types with econometric priors: exponentiate mu_hat
      if (length(econ_cols) > 0) {
        econ_priors <- exp(curr_priors_mu_hat[, econ_cols, drop = FALSE])
        
        # For overlapping pop_types: weighted blend (50/50 in log-space)
        # Scale exo priors to match econometric scale
        exo_subset <- curr_priors_exo[, econ_cols, drop = FALSE]
        exo_subset[exo_subset == 0] <- 1e-12
        
        # Mean-center the econometric log-odds spatially to prevent arbitrary intercepts from overwhelming log_exo
        log_econ <- curr_priors_mu_hat[, econ_cols, drop = FALSE]
        log_econ_centered <- scale(log_econ, center = TRUE, scale = FALSE)
        
        # Soften extreme spatial variances (limit standard deviation to 2 to prevent winner-takes-all)
        log_econ_sd <- apply(log_econ_centered, 2, sd)
        for (j in seq_along(log_econ_sd)) {
          if (!is.na(log_econ_sd[j]) && log_econ_sd[j] > 2) {
             log_econ_centered[, j] <- log_econ_centered[, j] * (2 / log_econ_sd[j])
          }
        }
        
        # Combine: Exo provides absolute baseline scale (weight 1.0), Econ provides spatial modulation (weight 0.5)
        log_exo <- log(pmax(exo_subset, 1e-12))
        curr_priors[, econ_cols] <- exp(log_exo + 0.5 * log_econ_centered)
      }
      
      # For pop_types with only exogenous priors: use directly
      if (length(exo_only_cols) > 0) {
        curr_priors[, exo_only_cols] <- curr_priors_exo[, exo_only_cols, drop = FALSE]
      }
      
    } else if (has_econ) {
      # Only econometric priors
      econ_cols <- intersect(colnames(curr_priors_mu_hat), curr_pop_types)
      if (length(econ_cols) > 0) {
        log_econ <- curr_priors_mu_hat[, econ_cols, drop = FALSE]
        log_econ_centered <- scale(log_econ, center = TRUE, scale = FALSE)
        log_econ_sd <- apply(log_econ_centered, 2, sd)
        for (j in seq_along(log_econ_sd)) {
          if (!is.na(log_econ_sd[j]) && log_econ_sd[j] > 2) {
             log_econ_centered[, j] <- log_econ_centered[, j] * (2 / log_econ_sd[j])
          }
        }
        curr_priors[, econ_cols] <- exp(log_econ_centered)
      }
      
    } else if (has_exo) {
      # Only exogenous priors
      curr_priors <- curr_priors_exo
      
    } else {
      # No priors at all — uniform fallback
      curr_priors[] <- 1
    }
    
    # Ensure no absolute zeros
    curr_priors[curr_priors == 0] <- 1e-12
    curr_priors[is.na(curr_priors)] <- 1e-12
    
    # ---------------------------------------------------------------
    # 6. Restrictions Preprocessing (unchanged)
    # ---------------------------------------------------------------
    if (!is.null(restrictions)) {
      if (!is.data.table(restrictions)) restrictions <- as.data.table(restrictions)
      curr_restr_dt <- restrictions[times == timestep]
      
      curr_restr_mat <- dcast(curr_restr_dt, ns ~ pop_type, value.var = "value", fill = 0, fun.aggregate = max)
      curr_restr_mat <- as.matrix(curr_restr_mat, rownames = "ns")
      
      curr_restr <- matrix(0, n, p, dimnames = list(ns_avail, curr_pop_types))
      common_ns <- intersect(rownames(curr_restr_mat), ns_avail)
      common_cols <- intersect(colnames(curr_restr_mat), curr_pop_types)
      curr_restr[common_ns, common_cols] <- curr_restr_mat[common_ns, common_cols]
    } else {
      curr_restr <- NULL
    }
    
    # ---------------------------------------------------------------
    # 7. Call the Core Allocation Engine
    # ---------------------------------------------------------------
    res <- optim_pop_allocation(
      targets = curr_targets_vec,
      priors = curr_priors,
      restrictions = curr_restr,
      pixel_avail = ns_avail,
      options = options
    )
    
    # 8. Format and Return Output
    out_res <- as.data.table(res$out_mu, keep.rownames = "ns")
    out_res <- melt(out_res, id.vars = "ns", variable.name = "pop_type", value.name = "value")
    
    return(list(
      out_res = out_res,
      solver_runs = res$solver_runs,
      out_diag = res$out_diag
    ))
}
