# ===============================================================
# areas_update.R
# ===============================================================
#' Areas summing up function over land-use changes
#'
#' @param res Result from downscale
#' @param areas Data.table of current areas
#' @param priors Priors for updating
#' @param xmat X-matrix data.table
#' @param xmat_proj Projected x matrix
#'
#' @return A data.table of areas
#' @export
areas_sum_to <- function(res, areas, priors, xmat, xmat_proj) {
  out <- res[
    , .(value = sum(value)), by = .(ns, to_class)
  ][
    , .(ns, from_class = to_class, value)
  ]
  
  tol <- 1e-10   # tolerance for numerical noise
  min_val <- min(out$value)
  # Check for values below tolerance
  if (min_val < -tol) {

      message("Warning: some projected area values are below tolerance of ", -tol,
              ". Smallest value found: ", round(min_val, 6), ". Clipping to 0.")
    
  }
  
  # Clip all negative values to zero
  out$value <- pmax(out$value, 0)
  
  return(out[])
}

#' Areas identity function over land-uses
#'
#' @param res Result from downscale
#' @param areas Data.table of current areas
#' @param priors Priors for updating
#' @param xmat X-matrix data.table
#' @param xmat_proj Projected x matrix
#'
#' @return A data.table of areas
#' @export
areas_identity <- function(res, areas, priors, xmat, xmat_proj) {
  
  # Correct for small numerical mistakes
  if (min(res$value) < 0 && min(res$value) > -1e-10) {
    out[value < 0, value := 0]
  }
  
  return(res[])
}
