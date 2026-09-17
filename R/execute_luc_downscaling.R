# ===============================================================
# execute_luc_downscaling.R: Generalized Land Use Downscaling
# ===============================================================
require(data.table)
require(cli)
require(progress)

#' Generalized Land Use Downscaling Orchestrator
#'
#' @description
#' A high-level wrapper to orchestrate the execution of the spatial downscaling solver
#' across multiple scenarios, timesteps, and spatial regions sequentially or in parallel.
#'
#' @param targets A `data.table` of aggregate thematic land use targets.
#' @param start_map A `data.table` of initial high-resolution pixel-level areas.
#' @param mapping_grid A `data.table` linking high-resolution pixels to regional definitions.
#' @param priors Optional `data.table` of spatial priors.
#' @param restrictions Optional `data.table` of transition restrictions.
#' @param xmat_full Optional `data.table` of static spatially explicit explanatory variables.
#' @param xmat_proj_full Optional `data.table` of dynamic projected explanatory variables over time.
#' @param betas_full Optional `data.table` of transition coefficients matching `xmat`.
#' @param exo_priors_all Optional `data.table` of exogenous spatial priors (e.g., policy constraints).
#' @param region_cols Character vector specifying the spatial region grouping columns.
#' @param scen_cols Character string specifying the scenario column.
#' @param time_col Character string specifying the timestep column.
#' @param ncores Integer specifying the number of cores for parallel execution.
#' @param downscaling_mode Character string defining the solver mode (`"exact"`, `"fast"`, etc.).
#'
#' @return A `data.table` containing the fully downscaled land use transition trajectories.
#' @export
# ---------------------------------------------------------------------------
# Which betas apply to one region?
# ---------------------------------------------------------------------------
# Resolution order, most specific first:
#   1. `.REGION_ID.`  the betas are already keyed at the solver's own region level  -> exact match
#   2. `country_col`  the betas are keyed by country and the grid knows the region's country
#   3. HIERARCHICAL   the betas are keyed at some COARSER level of the same code scheme, so the
#                     longest key that PREFIXES the region id wins: region "DE11" takes "DE11" if
#                     present, else "DE1", else "DE".
#   4. give up        no key matches -> every beta, which is the historical behaviour
#
# WHY 3 EXISTS. A prior fitted with country random effects is keyed by country (CAPRI 2-char),
# while the solver splits on CAPRI_NUTS regions. Neither 1 nor 2 fires -- the betas carry no
# `.REGION_ID.` and no `country_col` -- so every region silently received the coefficients of ALL
# countries at once. That is not a loud failure: the solve still runs and still produces a map.
#
# WHY IT IS SAFE FOR OTHER MODELS. This branch is only reached when 1 and 2 both miss, so a model
# with a 1:1 region keying (or a country column, as GLOBIOM has) never enters it and is bit-for-bit
# unaffected. A model whose codes are not hierarchical finds no prefix and falls through to 4,
# exactly as before. It can only ever turn "pool everything" into "use the most specific match".
#
# The longest-prefix winner is unique by construction: two DISTINCT keys of equal length cannot
# both prefix the same string, so there is no tie to break and no ordering dependence.
#
# Matching is case-sensitive and purely positional, which suits NUTS/CAPRI-style codes where the
# hierarchy is carried in the leading characters. It is NOT suitable for name-based keys
# ("Austria"), but those come in through branch 2.
.resolve_region_betas <- function(betas_full, region_id, country_col = NULL, temp_country = NULL,
                                  key_col = NULL) {
  if (is.null(betas_full) || !nrow(betas_full)) return(betas_full)

  # 1. already keyed at region level
  if (".REGION_ID." %in% names(betas_full)) return(betas_full[.REGION_ID. == region_id])

  # 2. keyed by country, and we know this region's country
  if (!is.null(country_col) && country_col %in% names(betas_full) && !is.null(temp_country)) {
    return(betas_full[get(country_col) %in% temp_country])
  }

  # 3. keyed at a coarser level of the same scheme -> longest prefix of the region id
  kc <- key_col
  if (is.null(kc)) {
    cand <- c(getOption("cascade.beta_key_col", NULL), "group", "Ns_reg", "region", "country")
    cand <- cand[!vapply(cand, is.null, logical(1))]
    kc <- intersect(cand, names(betas_full))
    kc <- if (length(kc)) kc[1] else NULL
  }
  rid <- as.character(region_id)[1]
  if (!is.null(kc) && !is.na(rid) && nzchar(rid)) {
    keys <- unique(as.character(betas_full[[kc]]))
    keys <- keys[!is.na(keys) & nzchar(keys)]
    hits <- keys[startsWith(rid, keys)]
    if (length(hits)) {
      best <- hits[which.max(nchar(hits))]
      if (!isTRUE(.cascade_beta_state$announced)) {
        .cascade_beta_state$announced <- TRUE
        message(sprintf(">>> betas: no region/country key on the betas; resolving hierarchically on `%s` (e.g. region %s -> %s). Reported once per session.",
                        kc, rid, best))
      }
      return(betas_full[get(kc) == best])
    }
  }

  # 4. nothing matched
  if (!isTRUE(.cascade_beta_state$warned)) {
    .cascade_beta_state$warned <- TRUE
    warning(sprintf("betas carry no key matching region '%s' (looked for `.REGION_ID.`, `%s`%s) -- EVERY beta is being applied to EVERY region. Reported once per session.",
                    rid, if (is.null(country_col)) "country_col" else country_col,
                    if (is.null(kc)) "" else sprintf(", and prefixes of `%s`", kc)),
            call. = FALSE, immediate. = TRUE)
  }
  betas_full
}

# Once-per-session reporting flags, so the region loop does not emit one line per region.
.cascade_beta_state <- new.env(parent = emptyenv())


execute_luc_downscaling <- function(
  targets,
  start_map,
  mapping_grid,
  priors = NULL,
  restrictions = NULL,
  xmat_full = NULL,
  xmat_proj_full = NULL,
  betas_full = NULL,
  exo_priors_all = NULL,
  region_cols = "Ns_reg",
  scen_cols = "scenario",
  time_col = "times",
  country_col = "Ns_reg",
  pixel_id_col = "ns",
  from_col = "from_class",
  to_col = "to_class",
  value_col = "value",
  files = list(
    lu_output = "outputs/downscaled_results_lu.csv",
    luc_output = "outputs/downscaled_results_luc.csv",
    diag_output = "outputs/downscaled_results_diag.csv"
  ),
  baseline_scen = "BASE",
  tol_rel = 1e-6,
  scen_filter = NULL,
  year_filter = NULL,
  country_filter = NULL,
  countries_to_run = NULL,
  downscale_resolution = "5km",
  calc_prior_impact_metric = FALSE
) {
  # Ensure inputs are data.tables
  if (!is.data.table(targets)) setDT(targets)
  if (!is.data.table(start_map)) setDT(start_map)
  if (!is.data.table(mapping_grid)) setDT(mapping_grid)
  if (!is.null(priors) && !is.data.table(priors)) setDT(priors)
  if (!is.null(restrictions) && !is.data.table(restrictions)) setDT(restrictions)

  cli_h1("1. Preparing Generalized Data for LUC Optimization")

  # Resolve combined ID mapping for region and scenario dimensions
  if (length(region_cols) > 1) {
    targets[, .REGION_ID. := do.call(paste, c(.SD, sep = "::")), .SDcols = region_cols]
    start_map[, .REGION_ID. := do.call(paste, c(.SD, sep = "::")), .SDcols = region_cols]
    mapping_grid[, .REGION_ID. := do.call(paste, c(.SD, sep = "::")), .SDcols = region_cols]
    if (!is.null(xmat_full) && all(region_cols %in% names(xmat_full))) xmat_full[, .REGION_ID. := do.call(paste, c(.SD, sep = "::")), .SDcols = region_cols]
    if (!is.null(xmat_proj_full) && all(region_cols %in% names(xmat_proj_full))) xmat_proj_full[, .REGION_ID. := do.call(paste, c(.SD, sep = "::")), .SDcols = region_cols]
    if (!is.null(betas_full) && all(region_cols %in% names(betas_full))) betas_full[, .REGION_ID. := do.call(paste, c(.SD, sep = "::")), .SDcols = region_cols]
  } else {
    if (region_cols %in% names(targets)) setnames(targets, region_cols, ".REGION_ID.")
    if (region_cols %in% names(start_map)) setnames(start_map, region_cols, ".REGION_ID.")
    if (region_cols %in% names(mapping_grid)) setnames(mapping_grid, region_cols, ".REGION_ID.")
    if (!is.null(xmat_full) && region_cols %in% names(xmat_full)) setnames(xmat_full, region_cols, ".REGION_ID.")
    if (!is.null(xmat_proj_full) && region_cols %in% names(xmat_proj_full)) setnames(xmat_proj_full, region_cols, ".REGION_ID.")
    if (!is.null(betas_full) && region_cols %in% names(betas_full)) setnames(betas_full, region_cols, ".REGION_ID.")
  }

  if (length(scen_cols) > 1) {
    targets[, .SCEN_ID. := do.call(paste, c(.SD, sep = "::")), .SDcols = scen_cols]
    if (!is.null(priors) && all(scen_cols %in% names(priors))) priors[, .SCEN_ID. := do.call(paste, c(.SD, sep = "::")), .SDcols = scen_cols]
    if (!is.null(restrictions) && all(scen_cols %in% names(restrictions))) restrictions[, .SCEN_ID. := do.call(paste, c(.SD, sep = "::")), .SDcols = scen_cols]
  } else {
    if (scen_cols %in% names(targets)) setnames(targets, scen_cols, ".SCEN_ID.")
    if (!is.null(priors) && scen_cols %in% names(priors)) setnames(priors, scen_cols, ".SCEN_ID.")
    if (!is.null(restrictions) && scen_cols %in% names(restrictions)) setnames(restrictions, scen_cols, ".SCEN_ID.")
  }

  # Save mapping context for restoration
  region_map <- unique(mapping_grid[, c(".REGION_ID.", if (length(region_cols)>1) region_cols else NULL), with=FALSE])
  scen_map <- unique(targets[, c(".SCEN_ID.", if (length(scen_cols)>1) scen_cols else NULL), with=FALSE])

  # Standardize core column names
  if (time_col %in% names(targets)) setnames(targets, time_col, "times")
  if (from_col %in% names(targets)) setnames(targets, from_col, "from_class")
  if (to_col %in% names(targets)) setnames(targets, to_col, "to_class")
  if (value_col %in% names(targets)) setnames(targets, value_col, "value")

  if (!is.null(priors)) {
    if (time_col %in% names(priors)) setnames(priors, time_col, "times")
    if (from_col %in% names(priors)) setnames(priors, from_col, "from_class")
    if (to_col %in% names(priors)) setnames(priors, to_col, "to_class")
  }

  if (!is.null(restrictions)) {
    if (time_col %in% names(restrictions)) setnames(restrictions, time_col, "times")
    if (to_col %in% names(restrictions)) setnames(restrictions, to_col, "to_class")
  }

  if (from_col %in% names(start_map)) setnames(start_map, from_col, "from_class")
  if (value_col %in% names(start_map)) setnames(start_map, value_col, "value")

  if (pixel_id_col %in% names(mapping_grid)) setnames(mapping_grid, pixel_id_col, "ns")
  if (pixel_id_col %in% names(start_map)) setnames(start_map, pixel_id_col, "ns")
  if (!is.null(priors) && pixel_id_col %in% names(priors)) setnames(priors, pixel_id_col, "ns")
  if (!is.null(exo_priors_all) && pixel_id_col %in% names(exo_priors_all)) setnames(exo_priors_all, pixel_id_col, "ns")

  # Mock mappings expected by engine
  mappings <- list(mapping_grid = mapping_grid)

  # Base-year zero-transitions are injected once, upstream, by the caller (see
  # process_downscaling_v2.R's adapt_targets()) — not repeated here to avoid
  # double-injecting duplicate rows.
  base_year_val <- if (exists("BASE_YEAR")) BASE_YEAR else 2000

  # Ensure character typing
  targets[, times := as.character(times)]
  if (!is.null(priors) && "times" %in% names(priors)) priors[, times := as.character(times)]
  if (!is.null(restrictions) && "times" %in% names(restrictions)) restrictions[, times := as.character(times)]
  if (!is.null(priors)) {
    if (!".SCEN_ID." %in% names(priors)) stop(sprintf("execute_luc_downscaling: scen_cols '%s' not found in `priors` — cannot resolve .SCEN_ID.", paste(scen_cols, collapse=", ")))
    keep_cols_priors <- intersect(names(priors), c("ns", "from_class", "to_class", "times", ".SCEN_ID.", country_col, "GLOB_LUID", "GLOB_AEZ", "GLOB_altitude", "GLOB_slope", "GLOB_soil", "value"))
    priors <- priors[, ..keep_cols_priors]
  } else {
    priors <- data.table(ns=character(), from_class=character(), to_class=character(), times=character(), .SCEN_ID.=character(), value=numeric())
    if (!is.null(country_col)) priors[, (country_col) := character()]
  }

  if (!is.null(restrictions)) {
    if (!".SCEN_ID." %in% names(restrictions)) stop(sprintf("execute_luc_downscaling: scen_cols '%s' not found in `restrictions` — cannot resolve .SCEN_ID.", paste(scen_cols, collapse=", ")))
    keep_cols_restr <- intersect(names(restrictions), c("ns", "to_class", "times", ".SCEN_ID.", country_col, "GLOB_LUID", "GLOB_AEZ", "GLOB_altitude", "GLOB_slope", "GLOB_soil", "value"))
    restrictions <- restrictions[, ..keep_cols_restr]
  } else {
    restrictions <- data.table(ns=character(), to_class=character(), times=character(), .SCEN_ID.=character(), value=numeric())
    if (!is.null(country_col)) restrictions[, (country_col) := character()]
  }

# ---------------------------------------------------------------
# 2. Applying Scoping & Pre-splitting
# ---------------------------------------------------------------
if (!is.null(scen_filter))    targets <- targets[.SCEN_ID. %in% scen_filter]
if (!is.null(year_filter))    targets <- targets[times %in% as.character(c(base_year_val, year_filter))]
if (!is.null(country_filter) && !is.null(country_col) && country_col %in% names(targets)) targets <- targets[get(country_col) %in% country_filter]

# Protect against data.table NA-split crashes
if (!is.null(xmat_full)) xmat_full <- xmat_full[!is.na(.REGION_ID.)]
if (!is.null(xmat_proj_full)) xmat_proj_full <- xmat_proj_full[!is.na(.REGION_ID.)]
if (!is.null(targets)) targets <- targets[!is.na(.REGION_ID.) & !is.na(.SCEN_ID.)]

# Pre-split data for performance in regional loops
cli_alert_info("Pre-splitting data for regional loops...")

# Pre-split large inputs by region or scenario
# We drop the original objects immediately to free memory
if(!is.null(xmat_full)) {
  xmat_split <- split(xmat_full, by = ".REGION_ID.")
  rm(xmat_full)
} else {
  xmat_split <- NULL
}

if(!is.null(xmat_proj_full)) {
  xmat_proj_split <- split(xmat_proj_full, by = ".REGION_ID.")
  rm(xmat_proj_full)
} else {
  xmat_proj_split <- NULL
}

# The prior and restriction matrices can be huge, subset inside loop instead of full split
# We don't use split() here to avoid N copies

# Filter out targets for regions that have no start map data
start_map_sums <- start_map[, .(total_area = sum(value, na.rm = TRUE)), by = .REGION_ID.]
valid_regions <- start_map_sums[total_area > 1e-6, .REGION_ID.]
n_filtered <- length(unique(targets$.REGION_ID.)) - length(intersect(unique(targets$.REGION_ID.), valid_regions))
if (n_filtered > 0) {
  cli_alert_warning("Filtering out {n_filtered} regions from targets that have no start map data (total area = 0).")
  targets <- targets[.REGION_ID. %in% valid_regions]
}

targets_split      <- split(targets, by = c(".REGION_ID.", ".SCEN_ID."), keep.by = TRUE)
start_map_split    <- split(start_map, by = ".REGION_ID.", keep.by = TRUE)

rm(targets, start_map)
gc()

# ---------------------------------------------------------------
# 3. Optimization Loop (Regional & Scenario-based)
# ---------------------------------------------------------------
ds_grid <- unique(rbindlist(targets_split)[, .(.SCEN_ID., .REGION_ID.)])
ds_grid <- ds_grid[order(.SCEN_ID. != baseline_scen, .REGION_ID., .SCEN_ID.)]

cli_h1("2. Starting Land Use Downscaling (LUC)")
pb <- cli_progress_bar(total = nrow(ds_grid))

# Initializing output files (clear existing)
unlink(c(files$luc_output, files$lu_output, files$diag_output))
for (iter in seq_len(nrow(ds_grid))) {
  temp_Ns   <- ds_grid[iter, .REGION_ID.]
  temp_scen <- ds_grid[iter, .SCEN_ID.]
  
  cli_alert_info("Processing Scenario: {temp_scen} | Region: {temp_Ns}")
  
  # 3.1 Subset data for current loop
  temp_grid <- mappings$mapping_grid[.REGION_ID. == temp_Ns]
  temp_start_map <- start_map_split[[as.character(temp_Ns)]]
  temp_targets   <- targets_split[[paste(temp_Ns, temp_scen, sep = ".")]]
  temp_priors    <- priors[.SCEN_ID. == temp_scen]
  temp_restrict  <- restrictions[.SCEN_ID. == temp_scen]
  
  if (is.null(temp_targets) || nrow(temp_targets) == 0) {
    cli_alert_warning("No targets found for this region/scenario subset. Skipping.")
    cli_progress_update(id = pb)
    next
  }
  
  if (is.null(temp_start_map) || nrow(temp_start_map) == 0) {
    cli_alert_warning("No start map area found for region {temp_Ns}. Skipping.")
    cli_progress_update(id = pb)
    next
  }
  # We keep all targets even if their from_class is missing from the current state.
  # The solver will natively constrain them to 0 area and bubble up the unallocated residuals.
  # Chronological check: a from_class is valid if it exists in the base map OR was created in an earlier timestep
  valid_classes <- unique(temp_start_map$from_class)
  invalid_transitions <- data.table()
  
  ordered_times <- sort(unique(temp_targets$times))
  for (t in ordered_times) {
    t_targets <- temp_targets[times == t]
    missing_from <- t_targets[!(from_class %in% valid_classes)]
    
    if (nrow(missing_from) > 0) {
      invalid_transitions <- rbind(invalid_transitions, missing_from)
    }
    # Any crop generated as a to_class in this timestep becomes a valid from_class for future timesteps
    valid_classes <- unique(c(valid_classes, t_targets$to_class))
  }
  
  if (nrow(invalid_transitions) > 0) {
    cli_alert_warning("Found {nrow(invalid_transitions)} target transitions for classes that mathematically do not exist yet (e.g. {paste(head(unique(invalid_transitions$from_class), 3), collapse=', ')}). These will result in true unallocated residuals.")
  }
  temp_targets <- temp_targets[from_class != "NODATA" & to_class != "NODATA"]
  
  if (nrow(temp_targets) == 0) {
    cli_alert_warning("All targets were scrubbed as unfeasible. Skipping.")
    cli_progress_update(id = pb)
    next
  }



  # 3.2 Prepare Priors for Region
  temp_country <- if (!is.null(country_col) && country_col %in% names(temp_grid)) unique(temp_grid[[country_col]]) else NULL
  temp_grid_cols <- intersect(c("ns", country_col, "GLOB_LUID", "GLOB_AEZ", "GLOB_altitude", "GLOB_slope", "GLOB_soil"), names(temp_grid))
  temp_grid_ns <- unique(temp_grid[, ..temp_grid_cols])
  # Integrity check before joining priors
  if (anyNA(temp_grid_ns)) {
    cli_alert_danger("NAs found in mapping grid for region {temp_Ns}. Missing spatial IDs?")
  }


  if (!is.null(country_col) && country_col %in% names(temp_priors)) {
    priors_region <- temp_priors[get(country_col) %in% temp_country]
  } else {
    priors_region <- temp_priors
  }
  # Pre-filter priors by active transition triplets in targets to prevent many-to-many Cartesian join explosion
  active_triplets <- unique(temp_targets[, .(times, from_class, to_class)])
  priors_region <- priors_region[active_triplets, on = .(times, from_class, to_class), nomatch = 0L]
  # Determine the spatial keys available in priors
  spatial_keys <- intersect(c("ns", country_col, "GLOB_LUID", "GLOB_AEZ", "GLOB_altitude", "GLOB_slope", "GLOB_soil"), names(priors_region))
  spatial_keys <- intersect(spatial_keys, names(temp_grid_ns))
  priors_region <- priors_region[temp_grid_ns, on = spatial_keys, allow.cartesian = TRUE, nomatch = 0]
  priors_region[, value := value + runif(.N, 1e-7, 1e-6)] # Add jitter
  priors_region <- priors_region[, .(value = suppressWarnings(max(value))), by = .(ns, from_class, to_class, times)]

  if (!is.null(exo_priors_all)) {
    exo_priors_region <- exo_priors_all[ns %in% temp_grid_ns$ns]
    if (nrow(exo_priors_region) > 0) {
      priors_region <- rbind(priors_region, exo_priors_region, use.names = TRUE, fill = TRUE)
      priors_region <- priors_region[, .(value = suppressWarnings(max(value, na.rm = TRUE))), by = .(ns, from_class, to_class, times)]
    }
  }

  # 3.3 Prepare Restrictions for Region
  if (!is.null(country_col) && country_col %in% names(temp_restrict)) {
    restrict_region <- temp_restrict[get(country_col) %in% temp_country]
  } else {
    restrict_region <- temp_restrict
  }
  # Pre-filter restrictions by active target to_classes to prevent Cartesian join explosion
  restrict_region <- restrict_region[to_class %in% unique(temp_targets$to_class)]
  spatial_keys_restrict <- intersect(c("ns", country_col, "GLOB_LUID", "GLOB_AEZ", "GLOB_altitude", "GLOB_slope", "GLOB_soil"), names(restrict_region))
  spatial_keys_restrict <- intersect(spatial_keys_restrict, names(temp_grid_ns))
  restrict_region <- restrict_region[temp_grid_ns, on = spatial_keys_restrict, allow.cartesian = TRUE, nomatch = 0]
  
  by_cols <- intersect(names(restrict_region), c("ns", "from_class", "to_class", "times"))
  restrict_region <- restrict_region[, .(value = suppressWarnings(max(value, na.rm = TRUE))), by = by_cols]
  
  # 3.3.5 Prepare xmat and betas for Region
  if (!is.null(xmat_split) && !is.null(betas_full)) {
    # Subset xmat to only the 'ns' (spatial IDs) present in the current region
    xmat_region <- xmat_split[[as.character(temp_Ns)]]
    # --- DYNAMIC TYPE ASSIGNMENT (as data.table) ---
    # 1. Initialize all as "static"
    xmat_coltype <- data.table(
      ks = unique(xmat_region$ks),
      value = "static"
    )
    
    # 2. Identify LUM classes present in the current start map
    lum_classes <- unique(temp_start_map$from_class)
    
    # 3. Update to "dynamic" where ks matches a LUM class
    xmat_coltype[ks %in% lum_classes, value := "dynamic"]
    
    # message(">>> Covariate types assigned: ", paste(capture.output(xmat_coltype[, .N, by=type]), collapse="; "))
    
    if (length(lum_classes) > 0) {
      # This ensures that any covariate belonging to a LUM class is marked dynamic
      # Note: Adjust the logic if xmat_coltype is a named vector or a list
      message(">>> Setting ", length(lum_classes), " LUM covariates to type 'dynamic'.")
    }
    
    xmat_proj_region <- if (!is.null(xmat_proj_split)) {
      xmat_proj_split[[as.character(temp_Ns)]]
    } else {
      NULL
    }
    
    # Subset betas if they are region-specific. Exact keys first; then the hierarchical
    # fallback, which is what lets a country-keyed prior serve finer regions (see
    # .resolve_region_betas). Only if nothing matches at all do we hand over every beta.
    betas_region <- .resolve_region_betas(betas_full, temp_Ns, country_col, temp_country)
    
    # 0.5 Taxonomy Aliasing / Inheritance for unmodeled target classes
    # If targets demand classes that have mapped parents in betas (e.g. Forests_managed -> Forests,
    # Natural_unmanaged -> Natural_unmanaged_{fallow, set_aside, other}), clone parent model parameters
    if (nrow(betas_region) > 0 && "from_class" %in% names(betas_region) && "to_class" %in% names(betas_region)) {
      target_classes <- unique(c(as.character(temp_targets$from_class), as.character(temp_targets$to_class)))
      beta_classes <- unique(c(as.character(betas_region$from_class), as.character(betas_region$to_class)))
      
      alias_map <- c(
        Forests                     = "Forests_managed",
        Natural_unmanaged_fallow    = "Natural_unmanaged",
        Natural_unmanaged_set_aside = "Natural_unmanaged",
        Natural_unmanaged_other     = "Natural_unmanaged"
      )
      needed_aliases <- alias_map[names(alias_map) %in% target_classes & !names(alias_map) %in% beta_classes & alias_map %in% beta_classes]
      
      if (length(needed_aliases) > 0) {
        new_chunks <- list()
        for (alias_cls in names(needed_aliases)) {
          parent_cls <- needed_aliases[[alias_cls]]
          m_to <- betas_region[to_class == parent_cls]
          if (nrow(m_to) > 0) {
            c_to <- copy(m_to)
            c_to[, to_class := alias_cls]
            new_chunks[[paste0("to_", alias_cls)]] <- c_to
          }
          m_from <- betas_region[from_class == parent_cls]
          if (nrow(m_from) > 0) {
            c_from <- copy(m_from)
            c_from[, from_class := alias_cls]
            new_chunks[[paste0("from_", alias_cls)]] <- c_from
          }
        }
        for (from_alias in names(needed_aliases)) {
          from_parent <- needed_aliases[[from_alias]]
          for (to_alias in names(needed_aliases)) {
            to_parent <- needed_aliases[[to_alias]]
            m_both <- betas_region[from_class == from_parent & to_class == to_parent]
            if (nrow(m_both) > 0) {
              c_both <- copy(m_both)
              c_both[, from_class := from_alias]
              c_both[, to_class := to_alias]
              new_chunks[[paste0("both_", from_alias, "_", to_alias)]] <- c_both
            }
          }
        }
        if (length(new_chunks) > 0) {
          betas_region <- rbindlist(c(list(betas_region), new_chunks), use.names = TRUE, fill = TRUE)
          betas_region <- unique(betas_region, by = intersect(names(betas_region), c("from_class", "to_class", "ks")))
          message(">>> Inherited econometric betas for ", length(needed_aliases), " target classes (", paste(names(needed_aliases), collapse = ", "), ") in region ", temp_Ns)
        }
      }
    }
    
    # 1. Filter betas to active triplets BEFORE injecting KS
    active_beta_pairs <- unique(temp_targets[from_class != to_class & value > 0, .(from_class, to_class)])
    if (nrow(active_beta_pairs) > 0 && nrow(betas_region) > 0) {
      if ("from_class" %in% names(betas_region)) {
        betas_region <- betas_region[active_beta_pairs, on = .(from_class, to_class), nomatch = 0L]
      } else {
        betas_region <- betas_region[to_class %in% unique(active_beta_pairs$to_class)]
      }
    } else {
      betas_region <- betas_region[0]
    }
    
    # 1.5 Inject neutral (zero) betas for active demanded target transitions missing from betas & priors
    covered_pairs <- if (nrow(betas_region) > 0) unique(betas_region[, .(from_class, to_class)]) else data.table(from_class = character(), to_class = character())
    if (!is.null(priors_region) && nrow(priors_region) > 0) {
      covered_pairs <- unique(rbind(covered_pairs, unique(priors_region[, .(from_class, to_class)])))
    }
    missing_target_pairs <- active_beta_pairs[!covered_pairs, on = .(from_class, to_class)]
    if (nrow(missing_target_pairs) > 0) {
      all_ks <- unique(xmat_region$ks)
      if (length(all_ks) == 0) all_ks <- "intercept"
      neutral_betas <- missing_target_pairs[CJ(from_class = missing_target_pairs$from_class, ks = all_ks), on = .(from_class), allow.cartesian = TRUE]
      neutral_betas <- unique(neutral_betas)
      neutral_betas[, value := 0]
      if ("value_q025" %in% names(betas_region)) neutral_betas[, value_q025 := 0]
      if ("value_q975" %in% names(betas_region)) neutral_betas[, value_q975 := 0]
      if ("pooled" %in% names(betas_region)) neutral_betas[, pooled := FALSE]
      if (!is.null(country_col) && country_col %in% names(betas_region)) neutral_betas[, (country_col) := temp_country[1]]
      betas_region <- rbind(betas_region, neutral_betas, fill = TRUE)
      message(">>> Injected neutral (zero) betas for ", nrow(missing_target_pairs), " target transitions lacking econometric models in region ", temp_Ns)
    }
    
    # 2. Add missing KS for this region safely
    missing_region_ks <- setdiff(unique(xmat_region$ks), unique(betas_region$ks))
    if (length(missing_region_ks) > 0 && nrow(betas_region) > 0) {
       key_cols <- intersect(names(betas_region), c("from_class", "to_class", country_col))
       valid_keys_reg <- unique(betas_region[, ..key_cols])
       new_betas_reg <- valid_keys_reg[, .(ks = missing_region_ks), by = names(valid_keys_reg)]
       new_betas_reg[, value := 0]
       if ("value_q025" %in% names(betas_region)) new_betas_reg[, value_q025 := 0]
       if ("value_q975" %in% names(betas_region)) new_betas_reg[, value_q975 := 0]
       if ("pooled" %in% names(betas_region)) new_betas_reg[, pooled := FALSE]
       betas_region <- rbind(betas_region, new_betas_reg, fill=TRUE)
       message(">>> Injected ", length(missing_region_ks), " missing covariates locally for region ", temp_Ns)
    }
    
    # 3. Dynamic Prior Weighting
    priors_region[, weight := 1.0]
    if (nrow(betas_region) > 0) {
      if ("from_class" %in% names(betas_region)) {
        beta_avail <- unique(betas_region[, .(from_class, to_class)])
        beta_avail[, has_beta := TRUE]
        priors_region[beta_avail, on = .(from_class, to_class), weight := 0.5]
        
        to_class_with_betas <- unique(betas_region$to_class)
        priors_region[from_class == "__NA_LU__" & to_class %in% to_class_with_betas, weight := 0.5]
        
        from_class_with_betas <- unique(betas_region$from_class)
        priors_region[to_class == "__NA_LU__" & from_class %in% from_class_with_betas, weight := 0.5]
      } else {
        beta_avail <- unique(betas_region[, .(to_class)])
        beta_avail[, has_beta := TRUE]
        priors_region[beta_avail, on = .(to_class), weight := 0.5]
        
        to_class_with_betas <- unique(betas_region$to_class)
        priors_region[from_class == "__NA_LU__" & to_class %in% to_class_with_betas, weight := 0.5]
      }
    }
  } else {
    xmat_region <- NULL
    xmat_proj_region <- NULL
    betas_region <- NULL
    xmat_coltype <- NULL
  }
  
  # 3.4 Call Optimization Engine (downscalr)
  temp_res <- downscale(
    targets       = temp_targets[from_class != to_class & value > 0],
    start_areas   = temp_start_map,
    restrictions  = restrict_region,
    priors        = priors_region,
    xmat          = xmat_region,
    xmat_coltypes = xmat_coltype,
    xmat_dyn_fun = xmat_share_to,
    #xmat_proj    = xmat_proj_region,
    betas = if (is.null(betas_region) || nrow(betas_region) == 0) NULL else betas_region,
    options       = downscale_control(
      mode = if (exists("DOWNSCALE_MODE")) DOWNSCALE_MODE else "exact",
      track_performance = TRUE,
      verbose = FALSE,
      rebalance_targets = TRUE
    )
  )
  
  # -----------------------------------------------------------
  # 3.5 Optional: Prior Impact Diagnostic (Counterfactual Run)
  # -----------------------------------------------------------
  if (calc_prior_impact_metric && exists("files") && !is.null(files$prior_impact_output)) {
    cli_alert_info("Running flat-prior counterfactual for prior impact diagnostic...")
    
    # Run the exact same optimization but with no priors and no xmat
    temp_res_flat <- downscale(
      targets       = temp_targets[from_class != to_class & value > 0],
      start_areas   = temp_start_map,
      restrictions  = restrict_region,
      priors        = NULL,
      xmat          = NULL,
      xmat_coltypes = NULL,
      xmat_dyn_fun  = xmat_share_to,
      betas         = NULL,
      options       = downscale_control(
        mode = if (exists("DOWNSCALE_MODE")) DOWNSCALE_MODE else "exact",
        track_performance = FALSE,
        verbose = FALSE,
        rebalance_targets = TRUE
      )
    )
    
    tt_res_flat <- temp_res_flat$out_res[value != 0]
    setnames(tt_res_flat, "value", "flat_value")
    
    # Join full model results with flat prior results
    tt_res_model <- temp_res$out_res[value != 0]
    tt_res_impact <- merge(tt_res_model, tt_res_flat, by = c("times", "ns", "from_class", "to_class"), all = TRUE)
    tt_res_impact[is.na(value), value := 0]
    tt_res_impact[is.na(flat_value), flat_value := 0]
    
    # Pixel-level metrics
    tt_res_impact[, shift_km2 := value - flat_value]
    
    # Regional KL divergence (Information Gain)
    # KL(P || Q) = sum( P * log2(P / Q) ) where P and Q are probabilities (shares of target)
    # Group by transition to normalize shares
    impact_reg <- tt_res_impact[, .(
      target_area = sum(value),
      total_relocated_area_km2 = sum(abs(shift_km2)) / 2, # divide by 2 to avoid double counting
      max_pixel_shift_km2 = max(abs(shift_km2)),
      KL_divergence = {
        P <- value / sum(value)
        Q <- flat_value / sum(flat_value)
        # Avoid log(0) and div-by-zero
        P_safe <- pmax(P, 1e-12)
        Q_safe <- pmax(Q, 1e-12)
        sum(P_safe * log2(P_safe / Q_safe))
      }
    ), by = .(times, from_class, to_class)]
    
    impact_reg[, `:=`(.SCEN_ID. = temp_scen, .REGION_ID. = temp_Ns)]
    
    fwrite(nThread = 1, impact_reg, file = files$prior_impact_output, append = file.exists(files$prior_impact_output), col.names = !file.exists(files$prior_impact_output))
    
    if (!is.null(files$prior_impact_pixels_output)) {
      tt_res_impact_pixels <- tt_res_impact[abs(shift_km2) > 1e-4]
      if (nrow(tt_res_impact_pixels) > 0) {
        tt_res_impact_pixels[, `:=`(.SCEN_ID. = temp_scen, .REGION_ID. = temp_Ns)]
        fwrite(nThread = 1, tt_res_impact_pixels, file = files$prior_impact_pixels_output, append = file.exists(files$prior_impact_pixels_output), col.names = !file.exists(files$prior_impact_pixels_output))
      }
    }
    rm(temp_res_flat, tt_res_flat, tt_res_impact, tt_res_impact_pixels, impact_reg)
  }

  # -----------------------------------------------------------
  # 4. Results Consolidation & Saving
  # -----------------------------------------------------------
  # 4.1 Post-process transitioning areas (LUC)
  tt_res <- temp_res$out_res[value != 0]
  tt_res[, `:=`(.SCEN_ID. = temp_scen, .REGION_ID. = temp_Ns, downscaling_area_flag = 0L)]
  
  # 4.2 Post-process aggregate land use (LU)
  tt_agg <- temp_res$out_res[, .(value = sum(value)), by = .(to_class, times, ns)]
  tt_agg[, `:=`(.SCEN_ID. = temp_scen, .REGION_ID. = temp_Ns, downscaling_area_flag = 0L)]
  setnames(tt_agg, "to_class", "lu")


  # 4.3 Diagnostics & Consistency Checks
  temp_targets_reg <- temp_targets[from_class != to_class, .(value = sum(value)), by = .(.SCEN_ID., .REGION_ID., from_class, to_class, times)]
  
  tt_res_trans <- tt_res[from_class != to_class, .(ds_val = sum(value)), by = .(from_class, to_class, times)]
  
  chck <- merge(temp_targets_reg, tt_res_trans, by = c("from_class", "to_class", "times"), all = TRUE)
  chck[is.na(value), value := 0]
  chck[is.na(ds_val), ds_val := 0]
  chck[, `:=`(diff = abs(value - ds_val), rel_diff = abs(value - ds_val) / pmax(value, 1e-12))]
  
  if (any(chck$diff > tol_rel)) {
    tot_demanded   <- sum(chck$value, na.rm = TRUE)
    tot_downscaled <- sum(chck$ds_val, na.rm = TRUE)
    pct_achieved   <- if (tot_demanded > 0) (tot_downscaled / tot_demanded) * 100 else 100
    
    cli_alert_warning("Target mismatch in {temp_Ns}: {round(sqrt(mean(chck$diff^2)), 4)} RMSD | Downscaled: {round(tot_downscaled, 2)} / Demanded: {round(tot_demanded, 2)} ({round(pct_achieved, 2)}%)")
    
    # Print top deviations for debugging
    message(">>> Top Deviations:")
    print_dt <- copy(head(chck[diff > 1e-4][order(-diff, -value)], 20))
    if (nrow(print_dt) == 0) print_dt <- data.table(Message = "No significant deviations found") else {
      print_dt[, diff := format(round(diff, 2), scientific = FALSE)]
      print_dt[, rel_diff := paste0(format(round(rel_diff * 100, 2), scientific = FALSE), "%")]
      print_dt[, value := format(round(value, 2), scientific = FALSE)]
      print_dt[, ds_val := format(round(ds_val, 2), scientific = FALSE)]
    }
    print(print_dt)
    
    faulty <- unique(chck[diff >= 1 | rel_diff > 0.01, .(times, from_class, to_class)])
    if (nrow(faulty) > 0) {
      tt_res[faulty, downscaling_area_flag := 1L, on = .(times, from_class, to_class)]
      
      faulty_agg <- unique(faulty[, .(times, lu = to_class)])
      tt_agg[faulty_agg, downscaling_area_flag := 1L, on = .(times, lu)]
    }
  }

  # 4.3.5 Area Conservation Check (per-pixel and per-region)
  # The total area per pixel should be conserved across the downscaling solve.
  # Start map: sum(value) by (ns, times=base_year) = pixel total area
  # Downscaled: sum(value) by (ns, times) over all to_class = pixel total area
  if (exists("files") && !is.null(files$area_conservation_output)) {
    start_area <- temp_start_map[, .(start_total = sum(value, na.rm = TRUE)), by = .(ns)]
    
    ds_area <- temp_res$out_res[, .(ds_total = sum(value, na.rm = TRUE)), by = .(ns, times)]
    
    conservation <- merge(ds_area, start_area, by = "ns", all.x = TRUE)
    conservation[is.na(start_total), start_total := 0]
    conservation[, pixel_deviation := ds_total - start_total]
    conservation[, pct_deviation := round(pixel_deviation / pmax(start_total, 1e-12) * 100, 4)]
    conservation[, .SCEN_ID. := temp_scen]
    conservation[, .REGION_ID. := temp_Ns]
    
    # Check against theoretical maximum area
    res_km <- suppressWarnings(as.numeric(gsub("km", "", downscale_resolution)))
    if (!is.na(res_km)) {
      max_pixel_area <- res_km^2
      conservation[, area_exceeds_max := ds_total > (max_pixel_area + 0.01)]
      conservation[, area_lower_max := ds_total < (max_pixel_area - 0.01)]
      exceeded_dt <- conservation[area_exceeds_max == TRUE]
      n_exceeded <- nrow(exceeded_dt)
      n_lower <- sum(conservation$area_lower_max, na.rm = TRUE)
      if (n_exceeded > 0) {
        exceed_diff <- exceeded_dt$ds_total - max_pixel_area
        max_ex <- max(exceed_diff, na.rm = TRUE)
        mean_ex <- mean(exceed_diff, na.rm = TRUE)
        severity <- if (max_ex < 0.1) "minor numerical inaccuracy" else "substantial excess"
        cli_alert_danger(sprintf("CRITICAL: %d pixels in %s exceed theoretical max pixel area of %.1f km²! Exceedance: Max +%.4f km², Mean +%.4f km² (%s).", n_exceeded, temp_Ns, max_pixel_area, max_ex, mean_ex, severity))
      }
      if (n_lower > 0) {
        cli_alert_info(sprintf("Info: %d pixels in %s have area lower than %.1f km² (likely coastal/border pixels).", n_lower, temp_Ns, max_pixel_area))
      }
    }
    
    # General check: total area consistency across projection timesteps
    region_area_by_time <- ds_area[, .(total_area = sum(ds_total, na.rm = TRUE)), by = times]
    if (nrow(region_area_by_time) > 1) {
       area_variance <- max(region_area_by_time$total_area) - min(region_area_by_time$total_area)
       if (area_variance > 0.01) {
           cli_alert_danger(sprintf("CRITICAL: Total area for region %s fluctuates by %.2f km² across projection timesteps! Area should remain constant.", temp_Ns, area_variance))
       } else {
           cli_alert_success(sprintf("Area conservation: Total area for %s is perfectly stable across all timesteps.", temp_Ns))
       }
    }
    
    # Summarize at region level
    region_conservation <- conservation[, .(
      start_total = sum(start_total, na.rm = TRUE),
      ds_total = sum(ds_total, na.rm = TRUE),
      max_pixel_deviation = suppressWarnings(max(abs(pixel_deviation), na.rm = TRUE)),
      mean_pixel_deviation = mean(pixel_deviation, na.rm = TRUE),
      n_pixels_with_deviation = sum(abs(pixel_deviation) > 0.01)
    ), by = .(.SCEN_ID., .REGION_ID., times)]
    region_conservation[, region_deviation := ds_total - start_total]
    
    fwrite(nThread = 1, region_conservation, file = files$area_conservation_output,
           append = file.exists(files$area_conservation_output),
           col.names = !file.exists(files$area_conservation_output))
    
    n_deviations <- sum(abs(conservation$pixel_deviation) > 0.01)
    if (n_deviations > 0) {
      cli_alert_warning("Area conservation: {n_deviations} pixels with >0.01 km² deviation in {temp_Ns}")
    }
  }
  # 4.4 Incremental Saving to Disk
  if (length(scen_cols) > 1) tt_res <- merge(tt_res, scen_map, by = ".SCEN_ID.", all.x = TRUE)
  else setnames(tt_res, ".SCEN_ID.", scen_cols)
  if (length(region_cols) > 1) tt_res <- merge(tt_res, region_map, by = ".REGION_ID.", all.x = TRUE)
  else setnames(tt_res, ".REGION_ID.", region_cols)
  cols_to_rm <- intersect(names(tt_res), c(".SCEN_ID.", ".REGION_ID."))
  if (length(cols_to_rm) > 0) set(tt_res, j = cols_to_rm, value = NULL)
  setnames(tt_res, c("times", "from_class", "to_class", "ns", "value"), c(time_col, from_col, to_col, pixel_id_col, value_col), skip_absent=TRUE)
  setcolorder(tt_res, c(scen_cols, region_cols, time_col, from_col, to_col, pixel_id_col, value_col))
  
  if (length(scen_cols) > 1) tt_agg <- merge(tt_agg, scen_map, by = ".SCEN_ID.", all.x = TRUE)
  else setnames(tt_agg, ".SCEN_ID.", scen_cols)
  if (length(region_cols) > 1) tt_agg <- merge(tt_agg, region_map, by = ".REGION_ID.", all.x = TRUE)
  else setnames(tt_agg, ".REGION_ID.", region_cols)
  cols_to_rm_agg <- intersect(names(tt_agg), c(".SCEN_ID.", ".REGION_ID."))
  if (length(cols_to_rm_agg) > 0) set(tt_agg, j = cols_to_rm_agg, value = NULL)
  setnames(tt_agg, c("times", "lu", "ns", "value"), c(time_col, "lu", pixel_id_col, value_col), skip_absent=TRUE)
  setcolorder(tt_agg, c(scen_cols, region_cols, time_col, "lu", pixel_id_col, value_col))
  
  fwrite(nThread = 1, tt_res, file = files$luc_output, append = file.exists(files$luc_output), col.names = !file.exists(files$luc_output))
  fwrite(nThread = 1, tt_agg, file = files$lu_output,  append = file.exists(files$lu_output),  col.names = !file.exists(files$lu_output))
  
  if (!is.null(temp_res$out_diag)) {
    tt_diag <- temp_res$out_diag[, `:=`(.SCEN_ID. = temp_scen, .REGION_ID. = temp_Ns)]
    tt_diag_final <- copy(tt_diag)
    if (length(scen_cols) > 1) tt_diag_final <- merge(tt_diag_final, scen_map, by = ".SCEN_ID.", all.x = TRUE)
    else setnames(tt_diag_final, ".SCEN_ID.", scen_cols)
    if (length(region_cols) > 1) tt_diag_final <- merge(tt_diag_final, region_map, by = ".REGION_ID.", all.x = TRUE)
    else setnames(tt_diag_final, ".REGION_ID.", region_cols)
    cols_to_rm_diag <- intersect(names(tt_diag_final), c(".SCEN_ID.", ".REGION_ID."))
    if (length(cols_to_rm_diag) > 0) set(tt_diag_final, j = cols_to_rm_diag, value = NULL)
    fwrite(nThread = 1, tt_diag_final, file = files$diag_output, append = file.exists(files$diag_output), col.names = !file.exists(files$diag_output))
    
    if (exists("files") && !is.null(files$report_output)) {
      chck_dt <- as.data.table(chck)
      # Ensure consistent types for merge
      chck_dt[, times := as.integer(times)]
      tt_diag[, times := as.integer(times)]
      
      perf_report <- merge(chck_dt, tt_diag, by = c(".SCEN_ID.", ".REGION_ID.", "times", "from_class", "to_class"), all.x = TRUE)
      perf_report <- perf_report[from_class != to_class]
      setcolorder(perf_report, c(".SCEN_ID.", ".REGION_ID.", "times", "from_class", "to_class", "value", "ds_val", "diff", "rel_diff"))
      fwrite(nThread = 1, perf_report, file = files$report_output, append = file.exists(files$report_output), col.names = !file.exists(files$report_output))
    }
  }
  
  # Final loop cleanup
  rm(tt_res, tt_agg, temp_res, priors_region, restrict_region, xmat_region, xmat_proj_region, betas_region); gc(verbose = FALSE)
  cli_progress_update(id = pb)
}

cli_progress_done(id = pb)
cli_alert_success("Downscaling completed successfully. Results saved to: {files$lu_output}")


  if (length(region_cols) == 1 && ".REGION_ID." %in% names(mapping_grid)) {
    setnames(mapping_grid, ".REGION_ID.", region_cols)
  }
  if (length(scen_cols) == 1 && ".SCEN_ID." %in% names(mapping_grid)) {
    setnames(mapping_grid, ".SCEN_ID.", scen_cols)
  }

  invisible(TRUE)
}
