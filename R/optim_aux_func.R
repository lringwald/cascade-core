# ===============================================================
# optim_aux_func.R
# ===============================================================
#' Multinomial logit mean calculation with restrictions and cutoff
#'
#' @param x Bias correction parameters
#' @param mu_alloc Means of each class
#' @param areas Vector of areas
#' @param restrictions Optional restrictions. Binary, if 1, the log-odds are set to zero
#' @param cutoff Optional, defaults to zero. If set, all log-odds below cutoff are set to zero
#' @param area_cutoff Boolean, defaults to FALSE. If TRUE, cutoff triggers based on area instead
#'
#' @return Returns the vector of log-odds times the areas
#'
#' @keywords internal
allocate_mnl <- function(x, mu0 = NULL, areas, restrictions = NULL, cutoff = 0, area_cutoff = 0, log_scale = FALSE, log_mu0 = NULL, tempering = 1.0) {
  areas <- as.vector(areas)
  
  if(log_scale){
    # Work in log scale (x is passed as log(x))
    if (is.null(log_mu0)) log_mu0 <- log(pmax(mu0, 1e-12))
    
    # log(P_ij) = log(mu0_ij) + x_j - log(sum(exp(log_mu0_i. + x.)))
    # With tempering: log(P_ij) = (log(mu0_ij)/T) + (x_j/T)
    log_mu_alloc <- sweep(log_mu0 / tempering, 2, x / tempering, "+")
    
    # Add implicit baseline (U_base = 0)
    matrix_log_sums <- matrixStats::rowLogSumExps(log_mu_alloc)
    log_row_sums <- matrixStats::rowLogSumExps(cbind(0, matrix_log_sums))
    
    # Normalize and scale by area: log(P_ij) = log_mu_alloc - log_row_sums + log(areas)
    mu_alloc <- sweep(log_mu_alloc, 1, log_row_sums - log(areas), "-")
    
    # Apply restrictions as absolute area capacity caps
    if (!is.null(restrictions)) {
      mu_alloc <- mu_alloc + log(1 - restrictions)
    }
    
    # Apply cutoffs
    if (cutoff > 0) {
      log_cutoff <- log(cutoff)
      mu_alloc[mu_alloc < log_cutoff] <- -Inf
    }
    
  } else {
    # Linear scale
    # P_ij = mu0_ij^(1/T) * x_j^(1/T)
    mu_alloc <- sweep(mu0^(1/tempering), 2, x^(1/tempering), "*")
    
    # Row normalize: P_ij = (mu0_ij * x_j) / (1 + sum(mu0_i. * x.))
    row_sums <- rowSums(mu_alloc) + 1
    mu_alloc <- sweep(mu_alloc, 1, areas / row_sums, "*")
    
    # Apply restrictions
    if (!is.null(restrictions)) {
      mu_alloc <- mu_alloc * (1 - restrictions)
    }
    
    # Apply cutoffs
    if (cutoff > 0) mu_alloc[mu_alloc < cutoff] <- 0
  }
  
  return(mu_alloc)
}

#' Poisson-type mean allocation with restrictions and cutoff
#'
#' @param x Bias correction parameters (length = number of classes)
#' @param mu0 Baseline mean intensities (matrix: rows = spatial units, cols = classes)
#' @param areas Vector of total area (length = nrow(mu0))
#' @param restrictions Optional binary matrix (same dim as mu0): 1 = restricted
#' @param cutoff Optional numeric cutoff; values below are zeroed
#' @param count_lb Lower bounds for allocated counts
#' @param count_ub Upper bounds for allocated counts
#' @param log_scale Logical, whether to apply calculations in log-space
#'
#' @return Matrix of allocated values (same shape as mu0)
#' @keywords internal
allocate_poisson <- function(x, mu0 = NULL,
                             restrictions = NULL,
                             cutoff = 0,
                             count_lb = 0,
                             count_ub = Inf,
                             log_scale = FALSE,
                             log_mu0 = NULL,
                             tempering = 1.0) {
  
  if (log_scale) {
    # Work in log-space for numerical stability (x is passed as log(x))
    if (is.null(log_mu0)) log_mu0 <- log(pmax(mu0, 1e-12))
    
    # λ_ij = μ_ij * x_j   (in log-space: logλ_ij = logμ_ij + logx_j)
    # With tempering: logλ_ij = (logμ_ij/T) + (logx_j/T)
    log_lambda <- sweep(log_mu0 / tempering, 2, x / tempering, "+")
    
    # Optional restriction mask
    if (!is.null(restrictions)) {
      log_lambda[restrictions == 1] <- -Inf
    }
    
    # Apply cutoffs
    if (cutoff > 0) {
      log_cut <- log(cutoff)
      log_lambda[log_lambda < log_cut] <- -Inf
    }
    
    # Apply count cutoffs
    if (count_lb > 0) {
      log_lb <- log(count_lb)
      log_lambda[log_lambda < log_lb] <- -Inf
    }
    if (any(count_ub < Inf)) {
      # The upper bound is a common bound across all classes in the pixel
      # Sum across classes (in log space)
      log_row_sums <- matrixStats::rowLogSumExps(log_lambda)
      over_cap <- log_row_sums > log(count_ub)
      if (any(over_cap)) {
        sub_count_ub <- if (length(count_ub) > 1) count_ub[over_cap] else count_ub
        # Scale down proportionally: new_log_lambda = log_lambda - log_row_sums + log(count_ub)
        log_lambda[over_cap, ] <- sweep(log_lambda[over_cap, , drop = FALSE], 1, log_row_sums[over_cap] - log(sub_count_ub), "-")
      }
    }
    
    mu_alloc <- log_lambda
    
  } else {
    # Linear version
    # With tempering: λ_ij = μ_ij^(1/T) * x_j^(1/T)
    lambda <- sweep(mu0^(1/tempering), 2, x^(1/tempering), "*")
    
    if (!is.null(restrictions)) {
      lambda <- lambda * (1 - restrictions)
    }
    
    if (cutoff > 0) lambda[lambda < cutoff] <- 0
    if (count_lb > 0) lambda[lambda < count_lb] <- 0
    if (any(count_ub < Inf)) {
      # The upper bound is a common bound across all classes in the pixel
      row_sums <- rowSums(lambda)
      over_cap <- row_sums > count_ub
      if (any(over_cap)) {
        sub_count_ub <- if (length(count_ub) > 1) count_ub[over_cap] else count_ub
        # Scale down proportionally
        lambda[over_cap, ] <- sweep(lambda[over_cap, , drop = FALSE], 1, row_sums[over_cap], "/")
        lambda[over_cap, ] <- sweep(lambda[over_cap, , drop = FALSE], 1, sub_count_ub, "*")
      }
    }
    
    mu_alloc <- lambda
  }
  
  return(mu_alloc)
}

#' Squared differences optimization function with multinomial logit
#'
#' @param x Bias correction parameters
#' @param mu0 Priors for each class (pixel-level priors)
#' @param areas Vector of areas
#' @param targets Targets to match
#' @param restrictions Optional restrictions. Binary, if 1, the log-odds are set to zero
#' @param cutoff Optional, defaults to zero. If set, all log-odds below cutoff are set to zero
#'
# ==============================================================================
# ERROR CALCULATION & OBJECTIVE HELPERS
# ==============================================================================

#' @keywords internal
#' @noRd
calculate_objective_error <- function(pred_totals, target_totals, obj_fun_name, optimize_log_scale = FALSE) {
  # Log-space comparison for specific solvers
  if (isTRUE(optimize_log_scale) && obj_fun_name %in% c("sqr_diff_mnl", "sqr_diff_poisson")) {
    return((log(pmax(pred_totals, 1e-12)) - log(pmax(target_totals, 1e-12)))^2)
  }
  
  # Check for relative error handling
  is_rel <- obj_fun_name %in% c("sqr_diff_rel_mnl", "sqr_kl_mnl", "sqr_kl_penalty")
  
  diff <- pred_totals - target_totals
  if (is_rel) {
    diff <- diff / pmax(target_totals, 1e-12)
  }
  return(diff^2)
}


#' Squared differences optimization function with multinomial logit
#'
#' @param x Bias correction parameters
#' @param mu0 Priors for each class (pixel-level priors)
#' @param areas Vector of areas
#' @param targets Target totals to match (vector, length = number of classes)
#' @param restrictions Optional restrictions (same shape as mu0): 1 = restricted
#' @param cutoff Optional cutoff. All log-odds below cutoff are set to zero
#'
#' @return Squared differences of optimized values and targets
#'
#' @keywords internal
sqr_diff_mnl <- function(x, mu0, areas, targets, restrictions = NULL,
                         cutoff = 0, area_cutoff = 0, log_scale = FALSE,
                         lambda = 1e6, tempering = 1.0) {
  x <- x[1:length(targets)]
  
  # Compute mnl allocation from bias parameters x
  mu_alloc <- allocate_mnl(x, mu0, areas, restrictions, cutoff, area_cutoff, log_scale = log_scale, tempering = tempering)
  
  if (log_scale) {
    pred_totals <- exp(matrixStats::colLogSumExps(mu_alloc))
  } else {
    pred_totals <- colSums(mu_alloc)
  }
  
  # --- Base objective: match target sums ---
  obj <- sum(calculate_objective_error(pred_totals, targets, "sqr_diff_mnl", log_scale))
  
  # --- Constraint penalty: per-pixel capacity ---
  if (!is.null(lambda)) {
    row_sum <- if (log_scale) rowSums(exp(mu_alloc)) else rowSums(mu_alloc)
    violation <- pmax(0, row_sum - areas)
    penalty <- lambda * sum(violation^2)
  } else {
    penalty <- 0
  }
  # Combine both
  obj_total <- obj + penalty
  
  return(obj_total)
}

#' Squared differences optimization function with Poisson-type allocation
#'
#' @param x Bias correction parameters
#' @param mu0 Priors for each class (pixel-level priors)
#' @param targets Target totals to match (vector, length = number of classes)
#' @param restrictions Optional restrictions (same shape as mu0): 1 = restricted
#' @param cutoff Optional numeric cutoff; values below are zeroed
#' @param count_lb Lower bounds for allocated counts
#' @param count_ub Upper bounds for allocated counts
#' @param log_scale Logical, perform calculations in log-space for stability
#' @param lambda Optional penalty multiplier for total mismatch (default 1e6)
#'
#' @return Squared differences between predicted and target totals
#' @keywords internal
sqr_diff_poisson <- function(x, mu0, targets,
                             restrictions = NULL,
                             cutoff = 0,
                             count_lb = 0,
                             count_ub = Inf,
                             log_scale = FALSE,
                             lambda = 1e6,
                             tempering = 1.0) {
  
  # Ensure x has same length as number of classes
  x <- x[1:length(targets)]
  
  # Compute Poisson allocation
  mu_alloc <- allocate_poisson(
    x = x, mu0 = mu0, restrictions = restrictions, cutoff = cutoff,
    count_lb = count_lb, count_ub = count_ub, log_scale = log_scale, tempering = tempering
  )
  
  if (log_scale) {
    pred_totals <- exp(matrixStats::colLogSumExps(mu_alloc))
  } else {
    pred_totals <- colSums(mu_alloc)
  }
  
  # --- Base objective: match target sums ---
  obj <- sum(calculate_objective_error(pred_totals, targets, "sqr_diff_poisson", log_scale))
  
  # Optional penalty: encourage negative/inf values or instability
  if (!is.null(lambda)) {
    bad_mask <- if (log_scale) !is.finite(mu_alloc) else (!is.finite(mu_alloc) | mu_alloc < 0)
    penalty <- lambda * sum(bad_mask)
  } else {
    penalty <- 0
  }
  
  # Combine both
  obj_total <- obj + penalty
  
  return(obj_total)
}

#' Relative squared differences optimization function with multinomial logit
#'
#' @param x Bias correction parameters (length = number of crops)
#' @param mu Means of each class (pixel-level priors)
#' @param areas Vector of areas
#' @param targets Country-level targets to match
#' @param restrictions Optional restrictions (pixels x crops), 1 = allowed, 0 = forbidden
#' @param cutoff Optional cutoff. All log-odds below cutoff are set to zero
#'
#' @return Relative squared differences of optimized values and targets
#' @keywords internal
sqr_diff_rel_mnl = function(x, mu0, areas, targets, restrictions = NULL, cutoff = 0,  area_cutoff = 0, tempering = 1.0 ) {
  
  # Extract bias parameters corresponding to crops
  x <- x[1:length(targets)]
  
  # Compute pixel-level allocation using allocate_mnl
  mu_alloc <- allocate_mnl(x, mu0, areas, restrictions, cutoff, area_cutoff, tempering = tempering)
  
  # Base objective: match target sums
  sum(calculate_objective_error(colSums(mu_alloc), targets, "sqr_diff_rel_mnl", FALSE))
}

#' KL divergence optimization function with multinomial logit
#'
#' @param x Bias correction parameters (length = number of crops)
#' @param mu Means of each class (pixel-level priors)
#' @param areas Vector of areas
#' @param targets Country-level targets to match
#' @param restrictions Optional restrictions (pixels x crops), 1 = allowed, 0 = forbidden
#' @param cutoff Optional cutoff. All log-odds below cutoff are set to zero
#'
#' @return KL divergence between pixel-level allocation and priors
#' @keywords internal
sqr_kl_mnl = function(x, mu0, areas, targets, restrictions = NULL, cutoff = 0,  area_cutoff = 0, tempering = 1.0) {
  
  # Extract bias parameters corresponding to crops
  x <- x[1:length(targets)]
  
  # Compute pixel-level allocation using allocate_mnl
  mu_alloc <- allocate_mnl(x, mu0, areas, restrictions, cutoff, area_cutoff, tempering = tempering)
  
  # Ensure no zeros to avoid log(0)
  mu_alloc_safe <- pmax(mu_alloc, 1e-12)
  mu_prior_safe <- pmax(mu, 1e-12)
  
  # Compute KL divergence at the pixel level
  kl_matrix <- mu_alloc_safe * log(mu_alloc_safe / mu_prior_safe)
  
  # Sum KL divergence across all pixels and crops
  kl_total <- sum(kl_matrix)
  
  # Total objective: KL + weighted penalty
  kl_total + sum(calculate_objective_error(colSums(mu_alloc), targets, "sqr_kl_mnl", FALSE))
}

#' KL-divergence penalized allocation objective
#'
#' @param x Vector of bias parameters per transition (length = n_transitions)
#' @param mu Pixel-level priors (pixels x transitions)
#' @param areas Pixel-level capacities (vector)
#' @param targets Country-level targets per transition
#' @param restrictions Pixel-level restrictions (pixels x transitions)
#' @param cutoff Minimum allocation per pixel (default 0)
#' @param weight_penalty Weight for deviation from targets (default 1)
#'
#' @return Scalar objective: KL divergence + weighted target deviation
sqr_kl_penalty <- function(x, mu, areas, targets, restrictions = NULL, cutoff = 0,  area_cutoff = 0 , weight_penalty = 1) {
  
  n_transitions <- length(targets)
  n_pixels <- nrow(mu)
  
  # Bias parameters per transition
  bias <- x[1:n_transitions]
  
  # Compute pixel-level allocation via multinomial logit-style scaling
  # Exponentiate bias to get positive allocation factors
  alloc_factors <- exp(matrix(bias, nrow = n_pixels, ncol = n_transitions, byrow = TRUE))
  
  # Apply priors and restrictions
  if (!is.null(restrictions)) alloc_factors <- alloc_factors * restrictions
  mu_alloc <- mu * alloc_factors
  
  # Clip to pixel capacities
  row_sums <- rowSums(mu_alloc)
  over_cap <- row_sums > areas
  if (any(over_cap)) {
    mu_alloc[over_cap, ] <- sweep(mu_alloc[over_cap, , drop = FALSE], 1, row_sums[over_cap], "/")
    mu_alloc[over_cap, ] <- sweep(mu_alloc[over_cap, , drop = FALSE], 1, areas[over_cap], "*")
  }
  
  # Apply cutoff
  if (cutoff > 0) mu_alloc[mu_alloc < cutoff] <- 0
  
  # KL divergence between allocation and priors
  mu_safe <- pmax(mu_alloc, 1e-12)
  priors_safe <- pmax(mu, 1e-12)
  kl_term <- sum(mu_safe * log(mu_safe / priors_safe))
  
  # Weighted penalty for deviation from country-level targets
  penalty <- weight_penalty * sum(calculate_objective_error(colSums(mu_alloc), targets, "sqr_kl_penalty", FALSE))
  
  # Total objective
  kl_term + penalty
}

#' Gradient of KL divergence optimization function with multinomial logit
#'
#' @param x Bias correction parameters (length = number of crops)
#' @param mu Means of each class (pixel-level priors)
#' @param areas Vector of areas
#' @param targets Country-level targets to match
#' @param restrictions Optional restrictions (pixels x crops), 1 = allowed, 0 = forbidden
#' @param cutoff Optional cutoff. All log-odds below cutoff are set to zero
#'
#' @return Gradient vector of the KL objective with respect to x
#' @keywords internal
grad_kl_mnl <- function(x, mu, areas, targets, restrictions = NULL, cutoff = 0,  area_cutoff = 0 ) {
  
  n_crops <- length(targets)
  x <- x[1:n_crops]
  
  # Compute pixel-level allocation
  mu_alloc <- allocate_mnl(x, mu, areas, restrictions, cutoff, area_cutoff)
  
  # Safe for log
  mu_alloc_safe <- pmax(mu_alloc, 1e-12)
  mu_prior_safe <- pmax(mu, 1e-12)
  
  # Compute pixel-level gradient of KL wrt mu_alloc
  # d/dmu (mu * log(mu / mu_prior)) = 1 + log(mu / mu_prior)
  grad_mu <- 1 + log(mu_alloc_safe / mu_prior_safe)
  
  # Soft penalty for target deviation
  target_error <- (colSums(mu_alloc) - targets) / pmax(targets, 1e-12)
  
  # Backpropagate target penalty to pixels
  grad_target <- sweep(matrix(target_error, nrow = nrow(mu_alloc), ncol = n_crops, byrow = TRUE),
                       2, rep(1, nrow(mu_alloc)), "*")
  
  # Combine KL gradient + penalty
  grad_total_mu <- grad_mu + grad_target
  
  # Convert pixel-level gradient to crop-level bias gradient
  # Assuming allocate_mnl is differentiable and linear in x
  # We'll approximate using numerical sum over pixel weights per crop
  grad_x <- numeric(n_crops)
  for (j in 1:n_crops) {
    grad_x[j] <- sum(grad_total_mu[, j])
  }
  
  grad_x
}


#' Gradient of KL divergence optimization function with multinomial logit
#'
#' @param x Vector of bias parameters per transition (length = n_transitions)
#' @param mu Pixel-level priors (pixels x transitions)
#' @param areas Pixel-level capacities (vector)
#' @param targets Country-level targets per transition
#' @param restrictions Pixel-level restrictions (pixels x transitions)
#' @param cutoff Minimum allocation per pixel (default 0)
#' @param weight_penalty Weight for deviation from targets (default 1)
#'
#' @return Gradient vector of the KL objective with respect to x
#' @keywords internal

#' @keywords internal
#' @noRd
grad_sqr_kl_penalty <- function(x, mu, areas, targets, restrictions = NULL, cutoff = 0,  area_cutoff = 0 , weight_penalty = 1) {
  
  n_transitions <- length(targets)
  n_pixels <- nrow(mu)
  bias <- x[1:n_transitions]
  
  # Pixel-level allocation factors
  alloc_factors <- exp(matrix(bias, nrow = n_pixels, ncol = n_transitions, byrow = TRUE))
  if (!is.null(restrictions)) alloc_factors <- alloc_factors * restrictions
  mu_alloc <- mu * alloc_factors
  
  # Clip to pixel capacities
  row_sums <- rowSums(mu_alloc)
  over_cap <- row_sums > areas
  if (any(over_cap)) {
    mu_alloc[over_cap, ] <- sweep(mu_alloc[over_cap, , drop = FALSE], 1, row_sums[over_cap], "/")
    mu_alloc[over_cap, ] <- sweep(mu_alloc[over_cap, , drop = FALSE], 1, areas[over_cap], "*")
  }
  
  if (cutoff > 0) mu_alloc[mu_alloc < cutoff] <- 0
  
  # Gradient w.r.t. bias parameters
  grad <- numeric(n_transitions)
  mu_safe <- pmax(mu_alloc, 1e-12)
  priors_safe <- pmax(mu, 1e-12)
  
  # KL gradient: d/dx = (1 + log(mu_alloc / priors))
  grad_kl <- colSums(mu_safe * (1 + log(mu_safe / priors_safe)))
  
  # Penalty gradient: d/dx approx via pixel-level aggregation
  target_error <- (colSums(mu_alloc) - targets) / pmax(targets, 1e-12)
  grad_penalty <- 2 * weight_penalty * target_error * colSums(mu_alloc) / pmax(targets, 1e-12)
  
  grad <- grad_kl + grad_penalty
  grad
}

#' Provides analytical gradients for the squared differences optimization function with multinomial logit
#'
#' @inheritParams sqr_diff_mnl
#'
#' @return Vector of length \eqn{J} containing the gradients for each value of \code{x}.
#'
#' Does not work for \code{cutoff} \eqn{> 0} or with restrictions. Generates a warning if these are supplied.
#'
#' @keywords internal
grad_sqr_diff_mnl = function(x,mu,areas,targets,restrictions = NULL,cutoff = 0,  area_cutoff = 0 ) {
  if (!is.null(restrictions) || cutoff > 0) {
    stop("Gradient optimization not yet implemented for restrictions or cutoff > 0. \n
            Consider switching to a gradient free optimiser (see the nloptr documentation for details).")
  }
  n = nrow(mu); J = length(x)
  mu = cbind(1,mu * matrix(exp(x) ,n,J,byrow = TRUE))
  theta = mu / c(rowSums(mu)); theta = theta[,-1,drop=FALSE]
  diff_targets = targets - colSums(areas * theta)
  
  grad_x = rep(0,J)
  for (jj in 1:J) {
    own_d = sum( areas * theta[,jj] * (1 - theta[,jj])) * diff_targets[jj]
    other_d = 0
    for (jj2 in c(1:J)[-jj]) {
      other_d = other_d + sum(areas * theta[,jj] * theta[,jj2]) * diff_targets[jj2]
    }
    grad_x[jj] = -2*own_d + 2*other_d
  }
  return(grad_x)
}


# Inequality constraint: row sums <= available area
#' @keywords internal
#' @noRd
pixel_capacity_constraint <- function(x, mu0, areas, targets, restrictions = NULL, cutoff = 0, area_cutoff = 0, lambda=0, log_scale=FALSE) {
  # Compute mu given x
  x = x[1:length(targets)]
  mu_alloc <- allocate_mnl(x, mu0, areas=areas, restrictions = restrictions, cutoff = cutoff, area_cutoff = area_cutoff, log_scale = log_scale)
  
  if (log_scale) {
    mu_alloc <- exp(mu_alloc)
  }
  
  # Row sums minus available area (must be <= 0)
  c_ineq <- as.numeric(rowSums(mu_alloc) - areas)
  return(c(c_ineq))
}
