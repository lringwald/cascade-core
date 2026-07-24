#' Fienberg rebalancing for two vectors with a given starting matrix
#'
#' @param start_mat An \eqn{n} by \eqn{n} matrix with non-negative elements
#' and row sums equal to 1.
#' @param target_from An \eqn{n} by 1 vector with values >=0.
#' @param target_to An \eqn{n} by 1 vector with values >=0. The sum of
#' \code{target_from} must be equal to the sum of \code{target_to}.
#'
#' @return The rebalanced \eqn{n} by \eqn{n} matrix \code{start_mat}.
#' @export
fienberg <- function(start_mat, target_from, target_to,
                     max.iter = 1000,
                     tol.optim = 1e-8,
                     tol.gain = 1e-10,
                     min.iter = 10,
                     min_flow = NULL) {
  
  if (!is.matrix(start_mat)) stop("start_mat must be a matrix.")
  if (any(start_mat < 0)) stop("start_mat must have non-negative entries.")
  if (abs(sum(target_from) - sum(target_to)) > 1e-6)
    stop("Sum of target_from and target_to must be equal.")
  
  S <- start_mat
  n <- nrow(S)
  k <- ncol(S)
  
  curr.optimVal <- Inf
  curr.gain <- Inf
  iter <- 0
  
  while ((iter < min.iter) || (curr.gain > tol.gain && curr.optimVal > tol.optim && iter < max.iter)) {
    weights <- target_to / colSums(S * target_from)
    # Handle possible division by zero or NaN if targets are 0 but matrix has non-zero
    weights[is.nan(weights) | is.infinite(weights)] <- 0
    S <- sweep(S, 2, weights, "*")
    
    rs <- rowSums(S)
    rs[rs == 0] <- 1 # Prevent division by zero
    S <- sweep(S, 1, rs, "/")
    
    new_optimVal <- sqrt(mean((target_to - colSums(S * target_from))^2))
    curr.gain <- abs(curr.optimVal - new_optimVal)
    curr.optimVal <- new_optimVal
    
    iter <- iter + 1
  }
  
  stop.reason <- if (curr.gain <= tol.gain) {
    paste0("Gain in objective value smaller than tolerance (", tol.gain, ").")
  } else if (curr.optimVal <= tol.optim) {
    paste0("Objective value smaller than tolerance (", tol.optim, ").")
  } else {
    paste0("Maximum iterations reached (", max.iter, ").")
  }
  
  # Optional thresholding of small flows
  if (!is.null(min_flow)) {
    S[S < min_flow] <- 0
    rs <- rowSums(S)
    rs[rs == 0] <- 1
    S <- sweep(S, 1, rs, "/")
  }
  
  return(list(
    adjusted_matrix = S,
    optimVal = curr.optimVal,
    gain = curr.gain,
    iter = iter,
    stop.reason = stop.reason
  ))
}

#' Convert Net Land Use to Gross Land Use Change Transitions
#'
#' @description
#' Harmonizes differing total land use area vectors (`from_targets` and `to_targets`) into an 
#' \eqn{n \times n} gross transition flow matrix using a cross-entropy (Fienberg) optimization.
#'
#' @param from_targets A `data.table` with columns `lu` and `value` (numeric >= 0) representing the initial areas.
#' @param to_targets A `data.table` with columns `lu` and `value` (numeric >= 0) representing the target areas.
#' @param prior Optional `data.table` with columns `from_class`, `to_class` and `value` (numeric >= 0) acting as a prior preference flow matrix.
#' @param keep_areas Character string indicating which marginal totals to preserve: `"from"`, `"to"`, or `"both"` (default).
#'
#' @return A `data.table` containing the solved transition flows with columns `from_class`, `to_class`, and `value`.
#' @export
solve_gross_transitions <- function(from_targets, to_targets, prior=NULL, keep_areas = "both") {
  require(data.table)
  min_cutoff = 1.0e-8
  prior_off_diag = 0.1

  # Error check of from_targets and to_targets
  if (!all(c("lu", "value") %in% names(from_targets))) stop("Missing columns in from_targets")
  if (!all(c("lu", "value") %in% names(to_targets))) stop("Missing columns in to_targets")
  
  if (!all(from_targets$value >= 0)) stop("All from_targets values must be larger or equal to zero.")
  if (!all(to_targets$value >= 0)) stop("All to_targets values must be larger or equal to zero.")
  
  if (!is.data.table(from_targets)) from_targets <- as.data.table(from_targets)
  if (!is.data.table(to_targets)) to_targets <- as.data.table(to_targets)
  
  # Ensure we don't modify the originals by reference unnecessarily if they are passed by user
  from_targets <- copy(from_targets[, .(lu, value)])
  to_targets <- copy(to_targets[, .(lu, value)])
  
  lu_classes <- unique(c(as.character(from_targets$lu), as.character(to_targets$lu)))
  
  sum_from <- sum(from_targets$value)
  sum_to <- sum(to_targets$value)
  
  if (keep_areas == "from") {
    if (sum_to > 0) to_targets[, value := value / sum_to * sum_from]
    sum_to <- sum_from
  } else if (keep_areas == "to") {
    if (sum_from > 0) from_targets[, value := value / sum_from * sum_to]
    sum_from <- sum_to
  } else if (keep_areas == "both") {
    if (abs(sum_from - sum_to) > 1e-6) {
      if (!"NODATA" %in% lu_classes) lu_classes <- c(lu_classes, "NODATA")
      if (sum_from > sum_to) {
        if ("NODATA" %in% to_targets$lu) {
          to_targets[lu == "NODATA", value := value + (sum_from - sum_to)]
        } else {
          to_targets <- rbindlist(list(to_targets, data.table(lu = "NODATA", value = sum_from - sum_to)), use.names = TRUE)
        }
      } else {
        if ("NODATA" %in% from_targets$lu) {
          from_targets[lu == "NODATA", value := value + (sum_to - sum_from)]
        } else {
          from_targets <- rbindlist(list(from_targets, data.table(lu = "NODATA", value = sum_to - sum_from)), use.names = TRUE)
        }
      }
    }
  } else {
    stop("keep_areas has to be ('from','to','both').")
  }

  n <- length(lu_classes)
  
  # Standardize classes using data.table merge
  base_classes <- data.table(lu = lu_classes)
  
  from_std <- from_targets[base_classes, on = "lu"]
  from_std[is.na(value), value := 0]
  
  to_std <- to_targets[base_classes, on = "lu"]
  to_std[is.na(value), value := 0]
  
  # Initialize Matrix A
  matA <- diag(n) + prior_off_diag
  matA <- sweep(matA, 1, rowSums(matA), "/")
  colnames(matA) <- rownames(matA) <- lu_classes

  # Apply prior matrix using optimized prefix matching
  if (!is.null(prior) && nrow(prior) > 0) {
    if (!is.data.table(prior)) setDT(prior)
    
    rn <- rownames(matA)
    cn <- colnames(matA)
    
    for (jjj in seq_len(nrow(prior))) {
      p_from <- prior$from_class[jjj]
      p_to <- prior$to_class[jjj]
      p_val <- prior$value[jjj]
      
      n_from <- nchar(p_from)
      row.pointer <- which(rn == p_from)
      if (length(row.pointer) == 0) {
        row.pointer <- which(startsWith(rn, p_from) & (nchar(rn) == n_from | substr(rn, n_from + 1, n_from + 1) == "."))
      }
      
      n_to <- nchar(p_to)
      col.pointer <- which(cn == p_to)
      if (length(col.pointer) == 0) {
        col.pointer <- which(startsWith(cn, p_to) & (nchar(cn) == n_to | substr(cn, n_to + 1, n_to + 1) == "."))
      }
      
      if (length(row.pointer) > 0 && length(col.pointer) > 0) {
        matA[row.pointer, col.pointer] <- p_val
      }
    }
    
    # Rescale the prior matrix so row sums equal 1 after prior modifications
    rs <- rowSums(matA)
    rs[rs == 0] <- 1 # Prevent division by zero
    matA <- sweep(matA, 1, rs, "/")
  }

  from_vec <- from_std$value
  to_vec <- to_std$value

  # Filter zero margins
  nz_from <- which(from_vec > 0)
  nz_to <- which(to_vec > 0)

  if (length(nz_from) < n || length(nz_to) < n) {
    matA <- matA[nz_from, nz_to, drop = FALSE]
    from_vec <- from_vec[nz_from]
    to_vec <- to_vec[nz_to]
  }

  if (length(from_vec) == 0 || length(to_vec) == 0) {
    return(data.table(from_class = character(0), to_class = character(0), value = numeric(0)))
  }

  # Run fienberg rebalancing
  res <- fienberg(start_mat = matA, target_from = from_vec, target_to = to_vec)

  # Extract results
  out_mat <- sweep(res$adjusted_matrix, 1, from_vec, "*")
  
  dt_out <- as.data.table(out_mat)
  dt_out[, from_class := lu_classes[nz_from]]
  
  output <- melt(dt_out, id.vars = "from_class", variable.name = "to_class", value.name = "value")
  output <- output[value > min_cutoff]
  
  # Ensure character strings match standard downscalr
  output[, from_class := as.character(from_class)]
  output[, to_class := as.character(to_class)]
  
  return(output)
}
