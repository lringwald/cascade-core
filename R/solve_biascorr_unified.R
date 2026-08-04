# ===============================================================
# solve_biascorr_unified.R
# ===============================================================
library(data.table)

#' Unified Bias Correction Solver
#' 
#' Dispatches to LU or Population/Livestock specific solvers using the same 
#' data.table optimized engine.
#' 
#' @param type "lu" or "pop"
#' @param targets data.table of targets
#' @param areas data.table of areas (LU only)
#' @param grid data.table of spatial grid (Pop only)
#' @param xmat data.table of explanatory variables
#' @param betas data.table of coefficients
#' @param priors data.table of priors
#' @param restrictions data.table of restrictions
#' @param timestep current year/time
#' @param options list of control parameters
#' 
#' @return list(out_res, solver_runs, out_diag)
#' @export
solve_biascorr_unified <- function(type = c("lu", "pop"),
                                   targets,
                                   areas = NULL,
                                   grid = NULL,
                                   xmat = NULL,
                                   betas = NULL,
                                   priors = NULL,
                                   restrictions = NULL,
                                   timestep = NULL,
                                   options = downscale_control()) {
  type <- match.arg(type)
  
  if (type == "lu") {
    # Dispatch to the Land-Use specific orchestrator
    return(solve_biascorr_lu(targets, areas, xmat, betas, priors, restrictions, timestep, type, options))
    
  } else if (type == "pop") {
    # Dispatch to Population specific orchestrator
    return(solve_biascorr_pop(targets, grid, xmat, betas, priors, restrictions, timestep, options))
  }
}