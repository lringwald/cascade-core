# ===============================================================
# downscale_control.R
# ===============================================================
#' Set Options for Downscaling Solver
#'
#' This function configures the downscalR solver with customizable options. It allows the selection of the solver function and adjustment of solver behavior and stopping conditions, tailored for specific downscaling needs.
#'
#' @param solve_fun The downscaling function to use. Valid options are "\code{solve_biascorr}" for bias correction and "\code{solve_notarget}" for non-targeted downscaling. Defaults to "\code{solve_biascorr}".
#' @param obj_fun The objective function to use. 
#' @param algorithm Specifies the solver algorithm. The default algorithm is "\code{NLOPT_LN_SBPLX}". For more information on available algorithms, refer to the \code{\link[nloptr]{nloptr}} package.
#' @param xtol_rel The relative tolerance level of the solver, controlling solution precision in relative terms. Default is \code{1.0e-20}, where a lower value indicates higher precision.
#' @param xtol_abs The absolute tolerance level of the solver, specifying solution precision in absolute terms. Default is \code{1.0e-20}, with a lower value indicating higher precision.
#' @param maxeval The maximum number of evaluations the solver will perform, limiting computational effort. Default is \code{1600}.
#' @param MAX_EXP Sets the numerical cutoff to prevent extremely high values in the Multinomial Logit (MNL) function. Default is \code{log(.Machine$double.xmax)}.
#' @param cutoff An optional cutoff in the share of a pixel to avoid MNL values close to zero. Default is \code{0}.
#' @param area_cutoff An additonal boolean to switch \code{cutoff} to be applied to predicted areas instead MNL values. Default is \code{FALSE}.
#' @param only_ipf_optim If \code{TRUE}, the solver skips the NLOPTR and sequential boost phases, returning the IPF-only result. Default is \code{FALSE}.
#' @param use_ipf Logical, whether to run the Generalized IPF phase. Default \code{TRUE}.
#' @param use_nloptr Logical, whether to run the accelerated NLOPTR phase if IPF doesn't converge. Default \code{TRUE}.
#' @param use_sequential_nloptr Logical, whether to run the sequential per-target boost phase for failing targets. Default \code{TRUE}.
#' @param use_mopup Logical, whether to run the final residual area distribution (mopping) layer. Default \code{TRUE}.
#' @param err.txt An optional error message for debugging purposes. Empty by default.
#'
#' @return Returns a list with the specified options for the downscaling solver.
#'
#' @details
#' The function returns a list of parameters to be used by the solver. It is essential to provide compatible parameters to ensure accurate downscaling results. The function will stop with an error message if \code{solve_fun} is not correctly specified.
#'
#' @export downscale_control
#'
#' @examples
#' # Default settings
#' opts_default = downscale_control()
#'
#' # Custom settings
#' opts_custom = downscale_control(algorithm = "some_other_algorithm",
#'                                 xtol_rel = 1e-5, xtol_abs = 1e-5,
#'                                 maxeval = 2000)
downscale_control = function(mode = c("exact", "fast"),
                             solve_fun = "solve_biascorr", 
                             algorithm = "NLOPT_LN_BOBYQA",# "NLOPT_LN_COBYLA",
                             obj_fun="sqr_diff_mnl", 
                             grad_obj_fun=NULL,
                             ineqc_fun = NULL,
                             pen_lambda = NULL,
                             optimize_log_scale = TRUE,
                             xtol_rel = 1.0e-5,
                             xtol_abs = 1.0e-5,
                             max_eval = 500,
                             min_eval = 100,
                             use_ipf = TRUE,
                             use_multi_nloptr = TRUE,
                             use_intermediate_target_update = TRUE,
                             use_single_nloptr = TRUE,
                             use_mopup = TRUE,
                             max_retry = 10,
                             MIN_EXP = 1e-12,
                             MAX_EXP = .Machine$double.xmax,
                             cutoff = 0, 
                             area_cutoff = 0,
                             max_diff = 1.0e-8,
                             ref_class_adjust_threshold = 1e-12,
                             rebalance_targets = FALSE,
                             track_performance = FALSE,
                             err.txt = "",
                             verbose = FALSE,
                             tempering_start = 1.0,
                             tempering_step = 0.25,
                             max_tempering = 2.0) {
  
  mode <- match.arg(mode)
  
  # Override settings for fast mode
  if (mode == "fast") {
    xtol_rel <- 1e-4
    xtol_abs <- 1e-4
    max_eval <- 500
    min_eval <- 100
    use_multi_nloptr <- FALSE
    use_single_nloptr <- FALSE
    use_mopup <- TRUE # Rely heavily on mopup for speed
    max_tempering <- 1.0 # Disable tempering
  }
  
  if (!solve_fun %in% c("solve_biascorr","solve_notarget")) {stop("solve_fun not correctly specified.")}
  return(list(solve_fun = solve_fun,
              algorithm = algorithm,
              obj_fun=obj_fun, 
              grad_obj_fun=grad_obj_fun, 
              ineqc_fun=ineqc_fun,
              pen_lambda = pen_lambda,
              optimize_log_scale = optimize_log_scale,
              xtol_rel = xtol_rel,
              xtol_abs = xtol_abs,
              max_eval = max_eval,
              min_eval = min_eval,
              max_retry = max_retry,
              use_ipf = use_ipf,
              use_multi_nloptr = use_multi_nloptr,
              use_intermediate_target_update = use_intermediate_target_update,
              use_single_nloptr = use_single_nloptr,
              use_mopup = use_mopup,
              MIN_EXP = MIN_EXP, 
              MAX_EXP = MAX_EXP,
              cutoff = cutoff,
              area_cutoff=area_cutoff,
              max_diff = max_diff,
              ref_class_adjust_threshold=ref_class_adjust_threshold,
              rebalance_targets=rebalance_targets,
              track_performance=track_performance,
              err.txt = err.txt,
              verbose = verbose,
              tempering_start = tempering_start,
              tempering_step = tempering_step,
              max_tempering = max_tempering
  ))
}


#' Set Options for Population Downscaling Solver
#'
#' This function configures solver settings for population downscaling tasks.
#' It provides defaults optimized for Poisson-type allocation and allows
#' overriding of solver behavior, convergence tolerance, and numerical cutoffs.
#'
#' @param solve_fun The solver function to use. Valid options include
#'   \code{"solve_biascorr"} or \code{"solve_notarget"}. Defaults to \code{"solve_biascorr"}.
#' @param obj_fun The objective function to use. Defaults to \code{"sqr_diff_poisson"}.
#' @param algorithm The optimization algorithm used by \pkg{nloptr}. Defaults to
#'   \code{"NLOPT_LN_SBPLX"} (Subplex algorithm).
#' @param optimize_log_scale Logical, whether to operate in log-space for stability.
#'   Defaults to \code{TRUE}.
#' @param xtol_rel Relative convergence tolerance (default \code{1.0e-20}).
#' @param xtol_abs Absolute convergence tolerance (default \code{1.0e-20}).
#' @param max_eval Maximum number of solver evaluations (default \code{2000}).
#' @param cutoff Optional lower cutoff applied to intensities (default \code{0}).
#' @param count_lb_cutoff Lower bound on allocated counts (default \code{0}).
#' @param count_ub_cutoff Upper bound on allocated counts (default \code{Inf}).
#' @param pen_lambda Optional penalty multiplier for violation constraints.
#' @param rebalance_targets Logical, whether to rebalance targets before solving (default \code{TRUE}).
#' @param MAX_EXP Numerical ceiling for exponential terms to prevent overflow.
#'   Defaults to \code{log(.Machine$double.xmax)}.
#' @param only_ipf_optim If \code{TRUE}, only the IPF initialization is run.
#' @param use_ipf Logical, whether to run the initial iteration phase.
#' @param use_nloptr Logical, whether to use NLOPTR for simultaneous optimization.
#' @param use_sequential_nloptr Logical, whether to use sequential single-target optimization for remaining errors.
#' @param use_mopup Logical, whether to apply the final mop-up layer to failing targets.
#' @param err.txt Optional string for caller identification (used in debug messages).
#'
#' @return A list of control parameters for use in population downscaling.
#' @export downscale_control_pop
#'
#' @examples
#' # Default Poisson-style downscaling settings
#' opts_pop = downscale_control_pop()
#'
#' # Custom tolerance and solver setup
#' opts_pop_custom = downscale_control_pop(
#'   obj_fun = "sqr_diff_poisson",
#'   algorithm = "NLOPT_LN_COBYLA",
#'   xtol_rel = 1e-8,
#'   max_eval = 3000
#' )
downscale_control_pop <- function(mode = c("exact", "fast"),
                                  solve_fun = "solve_biascorr",
                                  obj_fun = "sqr_diff_poisson",
                                  algorithm = "NLOPT_LN_SBPLX",
                                  grad_obj_fun=NULL,
                                  ineqc_fun = NULL,
                                  pen_lambda = NULL,
                                  optimize_log_scale = TRUE,
                                  xtol_rel = 1.0e-20,
                                  xtol_abs = 1.0e-20,
                                  max_eval = 3000,
                                  min_eval = 1000,
                                  only_ipf_optim = FALSE,
                                  use_ipf = TRUE,
                                  use_nloptr = TRUE,
                                  use_sequential_nloptr = TRUE,
                                  use_mopup = TRUE,
                                  max_retry = 10,
                                  MIN_EXP = 1e-12,
                                  MAX_EXP = .Machine$double.xmax,
                                  cutoff = 0, 
                                  count_lb_cutoff = 0,
                                  count_ub_cutoff = Inf,
                                  rebalance_targets = TRUE,
                                  track_performance = FALSE,
                                  max_diff = 1.0e-8,
                                  err.txt = "",
                                  tempering_start = 1.0,
                                  tempering_step = 0.25,
                                  max_tempering = 2.0) {
  
  mode <- match.arg(mode)
  
  if (mode == "fast") {
    xtol_rel <- 1e-4
    xtol_abs <- 1e-4
    max_eval <- 500
    min_eval <- 100
    use_nloptr <- FALSE
    use_sequential_nloptr <- FALSE
    use_mopup <- TRUE
    max_tempering <- 1.0 # Disable tempering
  }
  
  if (!solve_fun %in% c("solve_biascorr", "solve_notarget")) {
    stop("solve_fun not correctly specified.")
  }
  
  return(list(
    solve_fun = solve_fun,
    algorithm = algorithm,
    obj_fun=obj_fun, 
    grad_obj_fun=grad_obj_fun, 
    ineqc_fun=ineqc_fun,
    pen_lambda = pen_lambda,
    optimize_log_scale = optimize_log_scale,
    xtol_rel = xtol_rel,
    xtol_abs = xtol_abs,
    max_eval = max_eval,
    min_eval = min_eval,
    max_retry = max_retry,
    only_ipf_optim = only_ipf_optim,
    use_ipf = use_ipf,
    use_nloptr = use_nloptr,
    use_sequential_nloptr = use_sequential_nloptr,
    use_mopup = use_mopup,
    MIN_EXP = MIN_EXP, 
    MAX_EXP = MAX_EXP,
    cutoff = cutoff,
    count_lb_cutoff = count_lb_cutoff,
    count_ub_cutoff = count_ub_cutoff,
    rebalance_targets = rebalance_targets,
    track_performance = track_performance,
    max_diff = max_diff,
    err.txt = err.txt,
    tempering_start = tempering_start,
    tempering_step = tempering_step,
    max_tempering = max_tempering
  ))
}
