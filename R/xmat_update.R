# ===============================================================
# xmat_update.R
# ===============================================================
#' X-matrix summing up function over land-use changes
#'
#' @param res Result from downscale
#' @param curr_areas Data.table of current areas
#' @param priors Priors for updating
#' @param xmat X-matrix data.table
#' @param xmat_proj Projected x matrix
#'
#' @return A data.table of dynamically updated xmat columns
#' @import data.table
#' @export
xmat_sum_to <- function(res, curr_areas, priors, xmat, xmat_proj) {
  out <- res$out_res[
    , .(value = sum(value)), by = .(ns, to_class)
  ][
    , .(ns, ks = to_class, value)
  ]
  
  out[value < 0 & value > -1e-10, value := 0]
  return(out[])
}

#' X-matrix calculating shares function over land-use changes
#'
#' @param res Result from downscale
#' @param curr_areas Data.table of current areas
#' @param priors Priors for updating
#' @param xmat X-matrix data.table
#' @param xmat_proj Projected x matrix
#'
#' @return A data.table of dynamically updated xmat columns
#' @import data.table
#' @export
xmat_share_to <- function(res, curr_areas, priors, xmat, xmat_proj) {
  out <- res$out_res[
    , .(value = sum(value)), by = .(ns, to_class)
  ][
    , .(ns, ks = to_class, value)
  ]
  
  out[, value := value / sum(value), by = .(ns)]
  
  out[value < 0 & value > -1e-10, value := 0]
  return(out[])
}

#' X-matrix identity function
#'
#' This function returns the current projections from res
#'
#' @param res Result from downscale
#' @param curr_areas Data.table of current areas
#' @param priors Priors for updating
#' @param xmat X-matrix data.table
#' @param xmat_proj Projected x matrix
#'
#' @return A data.table of dynamically updated xmat columns
#' @import data.table
#' @export
xmat_identity <- function(res, curr_areas, priors, xmat, xmat_proj) {
  out <- copy(res$out_res)
  
  # Correct for small numerical mistakes
  out[value < 0 & value > -1e-10, value := 0]
  
  return(out[])
}
