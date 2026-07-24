# ===============================================================
# aux_func.R
# ===============================================================
PLCHOLD_REGION = "__NA_REGION__"
PLCHOLD_LU = "__NA_LU__"
PLCHOLD_POPT = "__NA_POPT__"
PLCHOLD_K = "__NA_K__"
PLCHOLD_T = "__NA_TIME__"

# ------------------------------------------------------------------------------
# INTERNAL HELPERS
# ------------------------------------------------------------------------------

#' Check Data.Table Integrity
#' @keywords internal
check_dt_integrity <- function(dt, nm, required_cols, err_txt = "", allow_na = FALSE) {
  if (is.null(dt)) return()
  if (!is.data.table(dt)) stop(paste0(err_txt, nm, " must be a data.table"))
  if (nrow(dt) == 0) stop(paste0(err_txt, "No observations in ", nm))
  if (!allow_na && anyNA(dt)) stop(paste0(err_txt, nm, " contains NA values"))
  
  missing_cols <- setdiff(required_cols, names(dt))
  if (length(missing_cols) > 0) {
    stop(paste0(err_txt, "Missing columns in ", nm, ": ", paste(missing_cols, collapse = ", ")))
  }
}

#' Extract Global Factor Levels
#'
#' Scans the raw inputs to build a master dictionary of all possible factor 
#' levels (land use classes, spatial units, variables, and chronological times).
#' Ensures consistent factor integer mapping across all data.tables.
#'
#' @keywords internal
extract_global_levels <- function(targets, start_areas, priors = NULL, restrictions = NULL, xmat, xmat_proj = NULL, times = NULL) {
  
  # 1. Chronological Time Extraction and Sorting
  if (is.null(times)) {
    raw_times <- if ("times" %in% names(targets)) unique(as.character(targets$times)) else unique(as.character(targets$year))
    # Try numeric sort first
    times <- if (!anyNA(suppressWarnings(as.numeric(raw_times)))) {
      as.character(sort(as.numeric(raw_times)))
    } else {
      sort(raw_times) 
    }
  } else {
    times <- as.character(times)
  }
  
  # 2. Master Land-Use / Population Dictionary
  master_type_levels <- unique(c(
    as.character(start_areas$from_class),
    as.character(targets$from_class),
    as.character(targets$to_class),
    as.character(start_areas$pop_type),
    as.character(targets$pop_type)
  ))
  if (!is.null(priors)) {
    master_type_levels <- unique(c(master_type_levels, as.character(priors$from_class), as.character(priors$to_class), as.character(priors$pop_type)))
  }
  if (!is.null(restrictions)) {
    master_type_levels <- unique(c(master_type_levels, as.character(restrictions$from_class), as.character(restrictions$to_class), as.character(restrictions$pop_type)))
  }
  master_type_levels <- unique(c(master_type_levels, PLCHOLD_LU, PLCHOLD_POPT))
  master_type_levels <- master_type_levels[!is.na(master_type_levels) & master_type_levels != ""]

  # 3. Variable (ks) Dictionary
  ks_levels <- unique(as.character(xmat$ks))
  if (!is.null(xmat_proj)) ks_levels <- unique(c(ks_levels, as.character(xmat_proj$ks)))
  ks_levels <- unique(c(ks_levels, PLCHOLD_K))
  ks_levels <- ks_levels[!is.na(ks_levels) & ks_levels != ""]
  
  # 4. Pixel (ns) Dictionary
  ns_levels <- unique(c(as.character(start_areas$ns), as.character(xmat$ns)))
  if (!is.null(xmat_proj)) ns_levels <- unique(c(ns_levels, as.character(xmat_proj$ns)))
  if (!is.null(priors)) ns_levels <- unique(c(ns_levels, as.character(priors$ns)))
  if (!is.null(restrictions)) ns_levels <- unique(c(ns_levels, as.character(restrictions$ns)))
  ns_levels <- ns_levels[!is.na(ns_levels)]
  
  # Return the Named List
  list(
    from_class = master_type_levels,
    to_class   = master_type_levels,
    pop_type   = master_type_levels,
    ns         = ns_levels,
    ks         = ks_levels,
    times      = unique(c(times, PLCHOLD_T))
  )
}

#' Standardize Factors using Global Levels
#'
#' Applies predefined factor levels to a data.table by reference.
#' Automatically assigns `ordered = TRUE` to the times column.
#'
#' @keywords internal
standardize_factors <- function(dt) {
  if (is.null(dt)) return(NULL)
  if (!is.data.table(dt)) setDT(dt)
  
  # Fetch the dictionary from the system options
  g_levels <- getOption("downscalr.global_levels")
  
  if (is.null(g_levels)) {
    stop("Global factor levels have not been set. Run extract_global_levels first.")
  }
  
  for (col in names(dt)) {
    if (col %in% names(g_levels)) {
      is_ord <- (col == "times")
      set(dt, j = col, value = factor(as.character(dt[[col]]), levels = g_levels[[col]], ordered = is_ord))
    }
  }
  return(dt)
}

# ------------------------------------------------------------------------------
# COMPLETION FUNCTIONS
# ------------------------------------------------------------------------------

#' Complete Input Targets
#' @keywords internal
#' @import data.table
complete_targets <- function(targets) {
  targets <- as.data.table(targets)
  
  if (!"from_class" %in% names(targets)) targets[, from_class := PLCHOLD_LU]
  else if (any(targets$from_class == PLCHOLD_LU)) stop(sprintf("from_class %s is reserved.", PLCHOLD_LU))
  
  if (!"times" %in% names(targets)) targets[, times := PLCHOLD_T]
  else if (any(targets$times == PLCHOLD_T)) stop(sprintf("times %s is reserved.", PLCHOLD_T))
  
  targets <- CJ(from_class = unique(targets$from_class),
                to_class   = unique(targets$to_class),
                times   = unique(targets$times))[targets, on = .(from_class, to_class, times)]
  targets[is.na(value), value := 0]
  
  targets <- standardize_factors(targets)
  setkey(targets, times, from_class, to_class)
  return(targets)
}

#' Complete Input Targets for Population Downscaling
#' @keywords internal
#' @import data.table
complete_targets_pop <- function(targets) {
  targets <- as.data.table(targets)
  
  if (!"pop_type" %in% names(targets)) targets[, pop_type := PLCHOLD_POPT]
  else if (any(targets$pop_type == PLCHOLD_POPT)) stop(sprintf("pop_type %s is reserved.", PLCHOLD_POPT))
  
  if (!"times" %in% names(targets)) targets[, times := PLCHOLD_T]
  else if (any(targets$times == PLCHOLD_T)) stop(sprintf("times %s is reserved.", PLCHOLD_T))
  
  targets <- CJ(pop_type = unique(targets$pop_type),
                times    = unique(targets$times))[targets, on = .(pop_type, times)]
  targets[is.na(value), value := 0]
  targets <- targets[value != 0]
  
  targets <- standardize_factors(targets)
  setkey(targets, times, pop_type)
  return(targets)
}

#' Complete Input Areas
#' @keywords internal
#' @import data.table
complete_areas <- function(areas) {
  areas <- as.data.table(areas)
  
  if (!"from_class" %in% names(areas)) areas[, from_class := PLCHOLD_LU]
  else if (any(areas$from_class == PLCHOLD_LU)) stop(sprintf("from_class %s is reserved.", PLCHOLD_LU))
  
  areas <- areas[, .(ns, from_class, value)]
  ns_orig_class <- class(areas$ns)[1]
  
  # Extract all possible classes from the global factor levels dictionary
  g_levels <- getOption("downscalr.global_levels")
  if (is.null(g_levels) || is.null(g_levels$from_class)) {
    all_classes <- unique(areas$from_class)
  } else {
    all_classes <- setdiff(unique(c(g_levels$from_class, g_levels$to_class)), c(PLCHOLD_LU, PLCHOLD_POPT))
  }
  
  areas <- areas[CJ(ns = unique(ns), from_class = all_classes), on = .(ns, from_class)]
  areas[is.na(value), value := 0]
  
  factor_cols <- if (ns_orig_class %in% c("numeric", "integer")) "from_class" else c("ns", "from_class")
  areas[, (factor_cols) := lapply(.SD, as.factor), .SDcols = factor_cols]
  
  setkey(areas, ns, from_class)
  return(areas)
}

#' Complete Start Grid for Population Downscaling
#' @keywords internal
#' @import data.table
complete_start_grid_pop <- function(grid, targets=NULL) {
  grid <- as.data.table(grid)
  if (!is.null(targets)) targets <- as.data.table(targets)
  
  if (!"pop_type" %in% names(grid)) {
    if (!is.null(targets) && "pop_type" %in% names(targets)) {
      grid <- CJ(ns = unique(grid$ns), pop_type = unique(targets$pop_type))[grid, on = .(ns, pop_type)]
    } else {
      grid[, pop_type := PLCHOLD_POPT]
    }
  } else if (any(grid$pop_type == PLCHOLD_POPT)) {
    stop(sprintf("pop_type %s is reserved.", PLCHOLD_POPT))
  }
  
  if (!"value" %in% names(grid)) grid[, value := as.integer(NA)]
  
  grid <- CJ(ns = unique(grid$ns), pop_type = unique(grid$pop_type))[grid, on = .(ns, pop_type)]
  grid[is.na(value), value := as.integer(NA)]
  
  grid <- standardize_factors(grid)
  setkey(grid, ns, pop_type)
  return(grid)
}

#' Complete Input Xmat
#' @keywords internal
#' @import data.table
complete_xmat <- function(xmat, areas) {
  xmat <- as.data.table(xmat)[, .(ns, ks, value)]
  
  ns_type <- class(areas$ns)[1]
  if (ns_type == "factor") {
    xmat[, ns := as.factor(as.character(ns))]
  } else if (ns_type %in% c("numeric", "integer")) {
    xmat[, ns := as.numeric(as.character(ns))]
  }
  
  xmat <- xmat[CJ(ns = unique(areas$ns), ks = unique(xmat$ks)), on = .(ns, ks)]
  xmat[is.na(value), value := 0]
  
  xmat[, ks := as.factor(ks)]
  if (ns_type == "factor") xmat[, ns := as.factor(ns)]
  
  setkey(xmat, ns, ks)
  return(xmat)
}

#' Complete Input Betas
#' @keywords internal
#' @import data.table
complete_betas <- function(betas) {
  betas <- as.data.table(betas)
  
  if (!"from_class" %in% names(betas)) betas[, from_class := PLCHOLD_LU]
  else if (any(betas$from_class == PLCHOLD_LU)) stop(sprintf("from_class %s is reserved.", PLCHOLD_LU))
  
  betas <- standardize_factors(betas)
  setkey(betas, from_class, to_class, ks)
  return(na.omit(betas))
}

#' Complete Input Betas for Population Downscaling
#' @keywords internal
#' @import data.table
complete_betas_pop <- function(betas) {
  betas <- as.data.table(betas)
  if (!"pop_type" %in% names(betas)) betas[, pop_type := PLCHOLD_POPT]
  
  betas <- standardize_factors(betas)
  setkey(betas, pop_type, ks)
  return(betas)
}

#' Complete Input Priors
#' @keywords internal
#' @import data.table
complete_priors <- function(priors, xmat, targets) {
  priors <- as.data.table(priors)
  
  if (!"times" %in% names(priors)) {
    priors[, dummy_join := 1L]
    times_dt <- data.table(times = unique(targets$times), dummy_join = 1L)
    priors <- priors[times_dt, on = "dummy_join", allow.cartesian = TRUE][, dummy_join := NULL]
  }
  
  if (!"from_class" %in% names(priors)) priors[, from_class := PLCHOLD_LU]
  if (!"to_class" %in% names(priors)) priors[, to_class := PLCHOLD_LU]
  if (!"weight" %in% names(priors)) priors[, weight := 1]
  
  priors <- standardize_factors(priors)
  ns_avail <- unique(xmat$ns)
  priors <- priors[ns %in% ns_avail]
  
  if (all(unique(priors$from_class) == PLCHOLD_LU)) {
    valid_combos <- CJ(
      times = unique(targets$times),
      from_class = PLCHOLD_LU,
      to_class = if (all(unique(priors$to_class) == PLCHOLD_LU)) PLCHOLD_LU else union(priors$to_class, targets$to_class),
      sorted = FALSE
    )
  } else if (all(unique(priors$to_class) == PLCHOLD_LU)) {
    valid_combos <- CJ(
      times = unique(targets$times),
      from_class = union(priors$from_class, targets$from_class),
      to_class = PLCHOLD_LU,
      sorted = FALSE
    )
  } else {
    combos_targets <- unique(targets[, .(times, from_class, to_class)])
    combos_priors <- unique(priors[, .(times, from_class, to_class)])
    valid_combos <- unique(rbindlist(list(combos_targets, combos_priors)))
  }
  
  needed_combos <- CJ(ns = ns_avail, dummy = 1L)
  valid_combos[, dummy := 1L]
  needed_combos <- needed_combos[valid_combos, on = "dummy", allow.cartesian = TRUE][, dummy := NULL]
  
  missing_combos <- needed_combos[!priors, on = .(ns, from_class, to_class, times)]
  if (nrow(missing_combos) > 0) {
    missing_combos[, `:=`(value = 0, weight = 0)]
    priors <- rbindlist(list(priors, missing_combos), use.names = TRUE, fill = TRUE)
  }
  
  setkey(priors, times, ns, from_class, to_class)
  return(priors)
}

#' Complete Input Priors for Population Downscaling
#' @keywords internal
#' @import data.table
complete_priors_pop <- function(priors, xmat, targets) {
  priors <- as.data.table(priors)
  xmat <- as.data.table(xmat)
  targets <- as.data.table(targets)
  
  if (!"times" %in% names(priors)) {
    priors[, dummy := 1L]
    tmp_times <- data.table(times = unique(targets$times), dummy = 1L)
    priors <- merge(priors, tmp_times, by = "dummy", allow.cartesian = TRUE)
    priors[, dummy := NULL]
  }
  
  if (!"pop_type" %in% names(priors)) {
    priors[, pop_type := PLCHOLD_POPT]
  } else {
    if (any(priors$pop_type == PLCHOLD_POPT)) stop(paste0("The (priors) pop_type '", PLCHOLD_POPT, "' is reserved."))
    if (any(is.na(priors$pop_type))) stop("Priors column pop_type contains NA.")
  }
  
  if (!"weight" %in% names(priors)) priors[, weight := 1]
  if (any(priors$weight < 0 | priors$weight > 1)) stop("All prior weights must be between 0 and 1.")
  
  priors <- standardize_factors(priors)
  ns_avail <- unique(xmat$ns)
  priors <- priors[ns %in% ns_avail]
  
  needed_combos <- CJ(
    times = unique(targets$times),
    ns = ns_avail,
    pop_type = if (all(unique(priors$pop_type) == PLCHOLD_POPT)) PLCHOLD_POPT else union(priors$pop_type, targets$pop_type),
    sorted = FALSE
  )
  
  missing_combos <- needed_combos[!priors, on = .(ns, pop_type, times)]
  if (nrow(missing_combos) > 0) {
    missing_combos[, `:=`(value = 0, weight = 0)]
    priors <- rbindlist(list(priors, missing_combos), use.names = TRUE, fill = TRUE)
  }
  
  setkey(priors, times, ns, pop_type)
  return(priors)
}

#' Complete Input Restrictions
#' @keywords internal
#' @import data.table
complete_restrictions <- function(restrictions, xmat, targets) {
  restrictions <- as.data.table(restrictions)
  
  if (!"times" %in% names(restrictions)) {
    restrictions[, dummy_join := 1L]
    times_dt <- data.table(times = unique(targets$times), dummy_join = 1L)
    restrictions <- restrictions[times_dt, on = "dummy_join", allow.cartesian = TRUE][, dummy_join := NULL]
  }
  
  if (!"from_class" %in% names(restrictions)) restrictions[, from_class := PLCHOLD_LU]
  if (!"to_class" %in% names(restrictions)) restrictions[, to_class := PLCHOLD_LU]
  
  restrictions <- standardize_factors(restrictions)
  ns_avail <- unique(xmat$ns)
  restrictions <- restrictions[ns %in% ns_avail]
  
  if (all(unique(restrictions$from_class) == PLCHOLD_LU)) {
    valid_combos <- CJ(
      times = unique(targets$times),
      from_class = PLCHOLD_LU,
      to_class = if (all(unique(restrictions$to_class) == PLCHOLD_LU)) PLCHOLD_LU else union(restrictions$to_class, targets$to_class),
      sorted = FALSE
    )
  } else if (all(unique(restrictions$to_class) == PLCHOLD_LU)) {
    valid_combos <- CJ(
      times = unique(targets$times),
      from_class = union(restrictions$from_class, targets$from_class),
      to_class = PLCHOLD_LU,
      sorted = FALSE
    )
  } else {
    combos_targets <- unique(targets[, .(times, from_class, to_class)])
    combos_restrictions <- unique(restrictions[, .(times, from_class, to_class)])
    valid_combos <- unique(rbindlist(list(combos_targets, combos_restrictions)))
  }
  
  needed_combos <- CJ(ns = ns_avail, dummy = 1L)
  valid_combos[, dummy := 1L]
  needed_combos <- needed_combos[valid_combos, on = "dummy", allow.cartesian = TRUE][, dummy := NULL]
  
  missing_combos <- needed_combos[!restrictions, on = .(ns, from_class, to_class, times)]
  if (nrow(missing_combos) > 0) {
    missing_combos[, value := 0]
    restrictions <- rbindlist(list(restrictions, missing_combos), use.names = TRUE, fill = TRUE)
  }
  
  setkey(restrictions, times, ns, from_class, to_class)
  return(restrictions)
}

#' Complete Input Restrictions for Population Downscaling
#' @keywords internal
#' @import data.table
complete_restrictions_pop <- function(restrictions, xmat, targets) {
  restrictions <- as.data.table(restrictions)
  xmat <- as.data.table(xmat)
  targets <- as.data.table(targets)
  
  if (!"times" %in% names(restrictions)) {
    restrictions[, dummy := 1L]
    tmp_times <- data.table(times = unique(targets$times), dummy = 1L)
    restrictions <- merge(restrictions, tmp_times, by = "dummy", allow.cartesian = TRUE)
    restrictions[, dummy := NULL]
  }
  
  if (!"pop_type" %in% names(restrictions)) {
    restrictions[, pop_type := PLCHOLD_POPT]
  } else {
    if (any(restrictions$pop_type == PLCHOLD_POPT)) stop(paste0("The (restrictions) pop_type '", PLCHOLD_POPT, "' is reserved."))
    if (any(is.na(restrictions$pop_type))) stop("Restrictions column pop_type contains NA.")
  }
  
  restrictions <- standardize_factors(restrictions)
  ns_avail <- unique(xmat$ns)
  restrictions <- restrictions[ns %in% ns_avail]
  
  needed_combos <- CJ(
    times = unique(targets$times),
    ns = ns_avail,
    pop_type = if (all(unique(restrictions$pop_type) == PLCHOLD_POPT)) PLCHOLD_POPT else union(restrictions$pop_type, targets$pop_type),
    sorted = FALSE
  )
  
  missing_combos <- needed_combos[!restrictions, on = .(ns, pop_type, times)]
  if (nrow(missing_combos) > 0) {
    missing_combos[, value := 0]
    restrictions <- rbindlist(list(restrictions, missing_combos), use.names = TRUE, fill = TRUE)
  }
  
  setkey(restrictions, times, ns, pop_type)
  return(restrictions)
}

#' Complete Input Xmat Column Types
#' @keywords internal
#' @import data.table
complete_xmat_coltypes <- function(xmat_coltypes, xmat) {
  if (is.null(xmat_coltypes)) {
    return(data.table(ks = unique(xmat$ks), value = "static"))
  }
  
  if (is.vector(xmat_coltypes) && !is.list(xmat_coltypes)) {
    warning("Deprecated: xmat_coltypes should be a data.table with columns ks and value.")
    if (is.null(names(xmat_coltypes))) stop("No names in xmat_coltypes vector")
    
    dt <- data.table(ks = unique(xmat$ks), value = "static")
    for (k in names(xmat_coltypes)) dt[ks == k, value := xmat_coltypes[[k]]]
    return(dt)
  }
  
  return(as.data.table(xmat_coltypes))
}

#' Complete Input Xmat Projections
#' @keywords internal
#' @import data.table
complete_xmat_proj <- function(xmat_proj) {
  xmat_proj <- as.data.table(xmat_proj)
  
  if (!"times" %in% names(xmat_proj)) {
    xmat_proj[, times := PLCHOLD_T]
  } else if (any(xmat_proj$times == PLCHOLD_T)) {
    stop(sprintf("The times %s is reserved.", PLCHOLD_T))
  }
  
  xmat_proj <- standardize_factors(xmat_proj)
  setkey(xmat_proj, times, ns, ks)
  return(xmat_proj)
}

#' Check Area Consistency with Land-Use Targets
#' @keywords internal
#' @import data.table
target_area_check <- function(targets, areas, restrictions = NULL) {
  targets <- as.data.table(targets)
  areas <- as.data.table(areas)
  
  report_mismatch <- function(dt, timestep) {
    neg_classes <- dt[value < -1 * .Machine$double.eps]
    setnames(neg_classes, "from_class", "class")
    if (nrow(neg_classes) > 0) {
      cat(sprintf("\nMismatch in timestep %s:\n", timestep))
      print(knitr::kable(neg_classes))
    }
    neg_classes
  }
  
  curr_levels <- areas[, .(value = sum(value)), by = .(from_class)]
  timesteps <- sort(unique(targets$times))
  mismatches <- list()
  
  targets_from <- targets[, .(value.from_class = sum(value)), by = .(times, from_class)]
  
  curr_levels_net <- copy(curr_levels)
  setnames(curr_levels_net, "value", "value_net")
  
  targets_to <- targets[value != 0]
  targets_to <- merge(targets_to, curr_levels_net, by = "from_class", all.x = TRUE)
  targets_to[, value := value + ifelse(value_net < 0, value_net / .N, 0), by = .(times, from_class)]
  targets_to_sum <- targets_to[, .(value.to_class = sum(value)), by = .(times, to_class)]
  
  for (t in timesteps) {
    lu_from_targets_t <- targets_from[times == t, .(from_class, value.from_class)]
    lu_to_targets_t <- targets_to_sum[times == t, .(to_class, value.to_class)]
    
    curr_levels <- merge(curr_levels, lu_from_targets_t, by = "from_class", all.x = TRUE)
    curr_levels[is.na(value.from_class), value.from_class := 0]
    curr_levels[, value := value - value.from_class]
    curr_levels[, value.from_class := NULL]
    
    mismatches[[as.character(t)]] <- report_mismatch(curr_levels, t)
    
    curr_levels[value < 0, value := 0] 
    curr_levels <- merge(curr_levels, lu_to_targets_t, by.x = "from_class", by.y = "to_class", all.x = TRUE)
    curr_levels[is.na(value.to_class), value.to_class := 0]
    curr_levels[, value := value + value.to_class]
    curr_levels[, value.to_class := NULL]
  }
}

# ------------------------------------------------------------------------------
# VALIDATION FUNCTIONS
# ------------------------------------------------------------------------------

#' Validate Inputs for Land-Use Allocation Model
#' @keywords internal
#' @import data.table
err_check_inputs <- function(targets, areas, xmat, betas,
                             areas.update.fun, xmat_coltypes,
                             xmat_proj = NULL, xmat_dyn_fun = NULL,
                             priors = NULL, restrictions = NULL,
                             err.txt = "") {
  
  # Base Conversions & Integrity Checks
  targets <- as.data.table(targets)
  areas <- as.data.table(areas)
  xmat <- as.data.table(xmat)
  betas <- as.data.table(betas)
  xmat_coltypes <- as.data.table(xmat_coltypes)
  if (!is.null(xmat_proj)) xmat_proj <- as.data.table(xmat_proj)
  if (!is.null(priors)) priors <- as.data.table(priors)
  if (!is.null(restrictions)) restrictions <- as.data.table(restrictions)
  
  check_dt_integrity(targets, "targets", c("from_class", "to_class", "times", "value"), err.txt)
  check_dt_integrity(areas, "areas", c("from_class", "ns", "value"), err.txt)
  check_dt_integrity(xmat, "xmat", c("ks", "ns", "value"), err.txt)
  check_dt_integrity(betas, "betas", c("ks", "from_class", "to_class", "value"), err.txt)
  check_dt_integrity(xmat_coltypes, "xmat_coltypes", c("ks", "value"), err.txt)
  
  if (!is.null(xmat_proj)) check_dt_integrity(xmat_proj, "xmat_proj", c("ns", "ks", "times", "value"), err.txt)
  if (!is.null(priors)) check_dt_integrity(priors, "priors", c("ns", "from_class", "to_class", "value"), err.txt)
  if (!is.null(restrictions)) check_dt_integrity(restrictions, "restrictions", c("ns", "from_class", "to_class", "value"), err.txt)
  
  # Value Constraints
  if (any(targets$value < 0)) stop(paste0(err.txt, "Negative targets"))
  if (any(areas$value < 0)) stop(paste0(err.txt, "Areas must be >= 0"))
  if (!is.null(priors) && any(priors$value < 0)) stop(paste0(err.txt, "Priors must be >= 0"))
  if (!is.null(restrictions) && !all(restrictions$value %in% c(0, 1))) stop(paste0(err.txt, "Restrictions must be 0 or 1"))
  max_target_sum <- if (nrow(targets) > 0) max(targets[, .(tot = sum(value, na.rm=TRUE)), by = times]$tot) else 0
  if (max_target_sum > sum(areas$value, na.rm=TRUE) + 1e-4) stop(paste0(err.txt, "Targets exceed total area"))
  
  # Target Coverage Check
  covered <- unique(betas[, .(from_class, to_class)])
  if (!is.null(priors)) covered <- unique(rbind(covered, priors[, .(from_class, to_class)]))
  
  targets_unique <- unique(targets[, .(from_class, to_class)])
  
  # Efficiently check coverage including wildcards
  wildcard <- factor(PLCHOLD_LU, levels = getOption("downscalr.global_levels")$from_class)
  
  # 1. Direct matches
  uncovered <- targets_unique[!covered, on = .(from_class, to_class)]
  
  # 2. Wildcard from_class (matches any from_class for a specific to_class)
  if (nrow(uncovered) > 0) {
    wc_from <- covered[from_class == wildcard & to_class != wildcard]
    if (nrow(wc_from) > 0) uncovered <- uncovered[!wc_from, on = .(to_class)]
  }
  
  # 3. Wildcard to_class (matches any to_class for a specific from_class)
  if (nrow(uncovered) > 0) {
    wc_to <- covered[to_class == wildcard & from_class != wildcard]
    if (nrow(wc_to) > 0) uncovered <- uncovered[!wc_to, on = .(from_class)]
  }
  
  # 4. Global wildcard (matches anything)
  if (nrow(uncovered) > 0) {
    if (nrow(covered[from_class == wildcard & to_class == wildcard]) > 0) {
      uncovered <- uncovered[0]
    }
  }
  
  if (nrow(uncovered) > 0) stop(paste0(err.txt, "Missing betas or priors for target: ", paste(uncovered[1, .(from_class, to_class)], collapse = " -> ")))
  
  # xmat Logic Validation
  valid_types <- c("static", "dynamic", "projected")
  if (!all(xmat_coltypes$value %in% valid_types)) stop(paste0(err.txt, "Invalid xmat_coltypes: must be one of ", paste(valid_types, collapse = ", ")))
  
  if ("dynamic" %in% xmat_coltypes$value && is.null(xmat_dyn_fun)) stop(paste0(err.txt, "Dynamic columns specified but xmat_dyn_fun missing"))
  if ("projected" %in% xmat_coltypes$value) {
    if (is.null(xmat_proj)) stop(paste0(err.txt, "Projected columns specified but xmat_proj missing"))
    
    proj_ks <- xmat_coltypes[value == "projected", ks]
    required_proj <- CJ(times = unique(targets$times), ks = proj_ks)
    missing_proj <- required_proj[!xmat_proj[, .(times, ks)], on = .(times, ks)]
    if (nrow(missing_proj) > 0) stop(paste0(err.txt, "Missing xmat_proj entry for ks=", missing_proj[1, ks], " at time=", missing_proj[1, times]))
  }
  
  # Completeness Checks
  if (!all(unique(xmat$ks) %in% betas$ks)) {
    missing_ks <- setdiff(unique(xmat$ks), betas$ks)
    stop(paste0(err.txt, "xmat ks not found in betas: ", paste(missing_ks, collapse = ", ")))
  }
  if (!all(unique(areas$ns) %in% xmat$ns)) stop(paste0(err.txt, "areas ns missing in xmat"))
  if (nrow(xmat) != length(unique(xmat$ks)) * length(unique(xmat$ns))) stop(paste0(err.txt, "xmat is missing ks x ns combinations"))
  
  for (dt_list in list(list(priors, "priors"), list(restrictions, "restrictions"), list(xmat_proj, "xmat_proj"))) {
    dt_obj <- dt_list[[1]]
    nm <- dt_list[[2]]
    if (!is.null(dt_obj) && !all(unique(areas$ns) %in% dt_obj$ns)) stop(paste0(err.txt, "Missing ns in ", nm, " (check reference to areas)"))
  }
  
  if (!is.null(xmat_proj) && nrow(xmat_proj) != (length(unique(xmat_proj$ks)) * length(unique(xmat_proj$ns)) * length(unique(xmat_proj$times)))) {
    stop(paste0(err.txt, "xmat_proj missing full ns x ks x times combinations"))
  }
}

#' Validate Inputs for Population Downscaling
#' @keywords internal
#' @import data.table
err_check_inputs_pop <- function(targets, xmat, betas, xmat_coltypes,
                                 xmat_proj = NULL, xmat_dyn_fun = NULL,
                                 priors = NULL, restriction = NULL,
                                 err.txt = "") {
  
  targets <- as.data.table(targets)
  xmat <- as.data.table(xmat)
  betas <- as.data.table(betas)
  xmat_coltypes <- as.data.table(xmat_coltypes)
  if (!is.null(xmat_proj)) xmat_proj <- as.data.table(xmat_proj)
  if (!is.null(priors)) priors <- as.data.table(priors)
  if (!is.null(restriction)) restriction <- as.data.table(restriction)
  
  check_dt_integrity(targets, "targets", c("pop_type", "times", "value"), err.txt)
  check_dt_integrity(xmat, "xmat", c("ks", "ns", "value"), err.txt)
  check_dt_integrity(betas, "betas", c("ks", "pop_type", "value"), err.txt)
  check_dt_integrity(xmat_coltypes, "xmat_coltypes", c("ks", "value"), err.txt)
  if (!is.null(xmat_proj)) check_dt_integrity(xmat_proj, "xmat_proj", c("ns", "ks", "times", "value"), err.txt)
  
  if (any(targets$value < 0)) stop(paste0(err.txt, "Negative targets"))
  
  covered <- unique(betas[, .(pop_type)])
  if (!is.null(priors)) covered <- unique(rbind(covered, priors[, .(pop_type)]))
  uncovered <- unique(targets[, .(pop_type)])[!covered, on = .(pop_type)]
  if (nrow(uncovered) > 0) stop(paste0(err.txt, "Missing betas or priors for target: ", paste(uncovered[1, .(pop_type)], collapse = " -> ")))
  
  valid_types <- c("static", "dynamic", "projected")
  if (!all(xmat_coltypes$value %in% valid_types)) stop(paste0(err.txt, "Invalid xmat_coltypes"))
  if ("dynamic" %in% xmat_coltypes$value && is.null(xmat_dyn_fun)) stop(paste0(err.txt, "Missing xmat_dyn_fun for dynamic xmat"))
  
  if ("projected" %in% xmat_coltypes$value) {
    if (is.null(xmat_proj)) stop(paste0(err.txt, "xmat_proj missing"))
    proj_ks <- xmat_coltypes[value == "projected", ks]
    required_proj <- CJ(times = unique(targets$times), ks = proj_ks)
    missing_proj <- required_proj[!xmat_proj[, .(times, ks)], on = .(times, ks)]
    if (nrow(missing_proj) > 0) stop(paste0(err.txt, "Missing xmat_proj for ks=", missing_proj[1, ks], " at time=", missing_proj[1, times]))
  }
  
  if (!all(unique(xmat$ks) %in% betas$ks)) stop(paste0(err.txt, "xmat ks missing in betas"))
  if (nrow(xmat) != (length(unique(xmat$ks)) * length(unique(xmat$ns)))) stop(paste0(err.txt, "xmat missing ks x ns combinations"))
  
  if (!is.null(xmat_proj)) {
    if (!all(unique(xmat$ns) %in% xmat_proj$ns)) stop(paste0(err.txt, "xmat_proj missing some ns"))
    if (nrow(xmat_proj) != (length(unique(xmat_proj$ks)) * length(unique(xmat_proj$ns)) * length(unique(xmat_proj$times)))) stop(paste0(err.txt, "Incomplete xmat_proj"))
  }
}
