# ===============================================================
# zonal_swapper.R
# ===============================================================
# Implements Zonal Bilateral Swapping (Hole & Fill)
#
# mode = "fragment" (default)
#   Removes entire sub-threshold values. Pixel and zone area
#   are conserved exactly.
#
# mode = "decimal"
#   Removes the sub-threshold fractional tail of every value so
#   all outputs are exact multiples of threshold. Coin allocation
#   uses the largest-remainder method to recover as much of the
#   pooled remainder as possible. Zone-level area is conserved
#   to within one coin (threshold) per zone; pixel-level totals
#   may shift by at most ±threshold per class.
#
# preserve_pixels = TRUE (default)
#   After reconstruction, any pixel that was non-zero in the input
#   but disappeared entirely is rescued by transferring min_pixel_area
#   from the zone's richest donor. Zone area is conserved exactly.
#   Truly unrecoverable pixels (no donor left) are reported as warnings.
#
# min_pixel_area = threshold (default)
#   Minimum area assigned to a rescued pixel. Defaults to threshold so
#   the rescued pixel sits exactly at the noise floor. Set higher to
#   guarantee a larger minimum presence. In decimal mode this is
#   automatically rounded up to the nearest threshold multiple.
# ===============================================================
library(data.table)
library(cli)
library(Rcpp)

# Rcpp functions are now compiled into the package namespace natively.

#' Spatial Constraint Enforcement via Zonal Swapping
#'
#' @description
#' Enforces zonal area constraints by swapping pixel-level allocations across zone boundaries
#' iteratively until zonal budget limits are strictly satisfied, preventing overshoot.
#'
#' @param dt A `data.table` of pixel-level thematic areas.
#' @param zone_cols Character vector of column names defining spatial zones.
#' @param pixel_id_col Character string of the pixel ID column name.
#' @param class_col Character string of the thematic class column name.
#' @param area_col Character string of the area value column name.
#' @param coord_cols Character vector of X and Y coordinate column names.
#' @param class_aggregate Optional grouping to conserve overall pixel area totals.
#' @param budget_table A `data.table` specifying maximum allowed areas per zone and class.
#' @param budget_zone_cols Character vector defining the zones in `budget_table`.
#' @param budget_class_col Character string of the class column in `budget_table`.
#' @param max_iters Integer maximum number of swapping iterations.
#'
#' @return A `data.table` with compliant pixel-level areas.
#' @export
execute_zonal_swaps <- function(dt,
                                zone_cols       = c("GLOB_country", "GLOB_LUID"),
                                pixel_id_col    = "ns",
                                class_col       = "GLOB_lu_class",
                                area_col        = "LUM_area_km2",
                                threshold       = 0.01,
                                allowed_classes = NULL,
                                pixel_meta_cols = c("GLOB_country", "GLOB_LUID", "GLOB_AEZ", "GLOB_REGION59"),
                                mode            = c("fragment", "decimal"),
                                preserve_pixels = TRUE,
                                min_pixel_area  = threshold) {

  mode <- match.arg(mode)

  # In decimal mode, min_pixel_area must be a coin multiple so the threshold invariant holds
  effective_min <- if (mode == "decimal" && (min_pixel_area %% threshold) > 1e-12) {
    snapped <- ceiling(min_pixel_area / threshold) * threshold
    cli_alert_info("min_pixel_area {min_pixel_area} rounded up to nearest threshold multiple: {snapped}")
    snapped
  } else {
    min_pixel_area
  }

  cli_alert_info("Zonal Bilateral Swaps [{mode} mode] — Threshold: {threshold}, Min pixel area: {effective_min}")

  keep_cols   <- unique(c(zone_cols, pixel_id_col, class_col, area_col))
  dt_internal <- dt[, ..keep_cols]
  setnames(dt_internal, c(pixel_id_col, class_col, area_col), c("pixel_id", "class_id", "val"))

  # Compute noise amount and clean base for each row
  if (mode == "fragment") {
    is_noisy <- if (!is.null(allowed_classes)) {
      dt_internal$val > 0 & dt_internal$val < threshold & dt_internal$class_id %in% allowed_classes
    } else {
      dt_internal$val > 0 & dt_internal$val < threshold
    }
    dt_internal[, noise_rem := fifelse(is_noisy, val, 0)]
  } else {
    calc_noise <- function(v, th) {
      rem <- v %% th
      # handle floating point precision (e.g. 0.02 %% 0.01 = 0.009999999)
      fifelse(rem > th - 1e-9 | rem < 1e-9, 0, rem)
    }
    dt_internal[, noise_rem := if (!is.null(allowed_classes)) {
      fifelse(class_id %in% allowed_classes, calc_noise(val, threshold), 0)
    } else {
      calc_noise(val, threshold)
    }]
  }
  dt_internal[, clean_val := val - noise_rem]

  dt_internal[, zone_id := .GRP, by = zone_cols]
  setkey(dt_internal, zone_id)
  zones <- unique(dt_internal$zone_id)

  # Class Extinction Protection (The "Ark" Protocol)
  # If a crop is so microscopic that ALL of its pixels are < threshold, it has 0 valid edges.
  # The C++ swapper would drop it, and the Mop-up step would overwrite it, causing extinction.
  # We artificially grant amnesty to the LARGEST pixel of that crop to serve as an "Ark".
  # The C++ swapper will then dump all the tiny fractions of that crop perfectly into the Ark!
  extinct_classes <- dt_internal[val > 0, .(total_clean = sum(clean_val)), by = .(zone_id, class_id)][total_clean == 0]
  if (nrow(extinct_classes) > 0) {
    arks <- dt_internal[extinct_classes, on = .(zone_id, class_id), nomatch = 0L]
    setorder(arks, zone_id, class_id, -val)
    arks <- arks[, .SD[1], by = .(zone_id, class_id)]
    
    dt_internal[arks, on = .(zone_id, class_id, pixel_id), `:=`(noise_rem = 0, clean_val = val)]
    cli_alert_info("Ark Protocol: Protected {nrow(extinct_classes)} microscopic crop(s) from extinction by designating host pixels.")
  }

  total_area_in  <- sum(dt_internal$val)
  new_fills_list <- vector("list", length(zones))

  cli_progress_bar("Swapping Zones", total = length(zones))

  for (i in seq_along(zones)) {
    z    <- zones[i]
    z_dt <- dt_internal[.(z)]

    if (any(z_dt$noise_rem > 0)) {

      # 1. Pool noise per pixel (holes) and per class
      hole_dt <- z_dt[noise_rem > 0, .(hole_size = sum(noise_rem)), by = pixel_id]
      pool_dt <- z_dt[noise_rem > 0, .(pool_size = sum(noise_rem)), by = class_id]

      # 2. Determine each class's allocatable pool
      if (mode == "decimal") {
        # Round 1: base coins per class (floor)
        pool_dt[, base_coins := floor(pool_size / threshold)]
        pool_dt[, frac_rem   := pool_size - base_coins * threshold]
        # Round 2: bonus coins from the zone-level remainder, largest-remainder method
        zone_rem <- sum(pool_dt$frac_rem)
        n_bonus  <- floor(zone_rem / threshold)
        if (n_bonus > 0) {
          setorder(pool_dt, -frac_rem)
          pool_dt[seq_len(min(n_bonus, .N)), base_coins := base_coins + 1L]
        }
        pool_dt[, coin_pool  := base_coins * threshold]
        # Hole targets: round up to nearest coin so every pixel aims for full coverage
        hole_targets <- ceiling(hole_dt$hole_size / threshold) * threshold
      } else {
        pool_dt[, coin_pool  := pool_size]
        hole_targets <- hole_dt$hole_size
      }

      # 3. Valid edges: pixels with holes, class has pool, class already present (healthy base)
      edges <- z_dt[pixel_id %in% hole_dt$pixel_id &
                      class_id %in% pool_dt$class_id &
                      clean_val > 0,
                    .(pixel_id, class_id, exist_area = clean_val)]

      filled_dt <- data.table()

      if (nrow(edges) > 0 && any(pool_dt$coin_pool > 0)) {
        setorder(edges, -exist_area)
        edges[, p_idx := match(pixel_id, hole_dt$pixel_id)]
        edges[, c_idx := match(class_id, pool_dt$class_id)]

        # C++ greedy allocator: flow = min(class_pool, pixel_hole) in dominance order
        res <- allocate_dominance_flows(pool_size = pool_dt$coin_pool,
                                        hole_size = hole_targets,
                                        edge_p    = edges$p_idx,
                                        edge_c    = edges$c_idx)

        # In decimal mode, snap fills to coin multiples to maintain the threshold invariant
        edges[, fill_area := if (mode == "decimal") {
          floor(res$edge_fill / threshold) * threshold
        } else {
          res$edge_fill
        }]
        filled_dt       <- edges[fill_area > 0, .(pixel_id, class_id, fill_area)]
        if (mode == "decimal") {
          alloc_per_pixel <- edges[, .(alloc = sum(fill_area)), by = pixel_id]
          hole_dt[, alloc := alloc_per_pixel$alloc[match(pixel_id, alloc_per_pixel$pixel_id)]]
          hole_dt[is.na(alloc), alloc := 0]
          hole_dt$hole_size <- hole_dt$hole_size - hole_dt$alloc
          hole_dt[, alloc := NULL]
        } else {
          hole_dt$hole_size <- res$hole_remainder
        }
      }

      # 4. Mop-up for unallocated holes
      #    fragment: all remaining holes → largest existing class in pixel
      #    decimal:  only for pixels that are 100% noise (clean_val == 0 everywhere);
      #              sub-threshold remainders in other pixels are truly irrecoverable
      mop_dt <- data.table()
      hole_dt <- hole_dt[hole_size > 1e-12]

      if (nrow(hole_dt) > 0) {
        largest_class <- z_dt[pixel_id %in% hole_dt$pixel_id & clean_val > 0]
        setorder(largest_class, pixel_id, -val)
        largest_class <- largest_class[, .SD[1], by = pixel_id]

        mop_dt <- merge(hole_dt, largest_class[, .(pixel_id, class_id)], by = "pixel_id", all.x = TRUE)

        na_p <- mop_dt[is.na(class_id), pixel_id]
        if (length(na_p) > 0) {
          cli_alert_warning("Zone {z}: {length(na_p)} pixel(s) are 100% noise — assigning remainder to largest noise class.")
          fallback <- z_dt[pixel_id %in% na_p & noise_rem > 0, .SD[1], by = pixel_id]
          mop_dt[is.na(class_id), class_id := fallback$class_id[match(pixel_id, fallback$pixel_id)]]
        }

        # Decimal mop-up: running remainder stochastic allocation for perfect mass conservation
        mop_dt <- if (mode == "decimal") {
          total_in <- sum(z_dt$noise_rem) - if (exists("filled_dt") && nrow(filled_dt) > 0) sum(filled_dt$fill_area) else 0
          
          total_coins <- floor(total_in / threshold)
          exact_rem <- total_in - (total_coins * threshold)
          
          mop_dt[, fill_area := 0]
          if (total_coins > 0) {
            while(total_coins > 0 && nrow(mop_dt) > 0) {
              setorder(mop_dt, -hole_size)
              idx <- min(total_coins, nrow(mop_dt))
              mop_dt[seq_len(idx), fill_area := fill_area + threshold]
              mop_dt[seq_len(idx), hole_size := hole_size - threshold]
              total_coins <- total_coins - idx
            }
          }
          
          if (exact_rem > 1e-12 && nrow(mop_dt) > 0) {
            setorder(mop_dt, -hole_size)
            mop_dt[1, fill_area := fill_area + exact_rem]
          }
          
          mop_dt[fill_area > 0, .(pixel_id, class_id, fill_area)]
        } else {
          mop_dt[, .(pixel_id, class_id, fill_area = hole_size)]
        }
      }

      # 5. Reconstruct zone: base is clean_val, then add fills
      z_dt[, val := clean_val]
      healthy_dt <- z_dt[val > 0, .(pixel_id, class_id, fill_area = val)]

      z_out <- rbindlist(list(healthy_dt, filled_dt, mop_dt), use.names = TRUE, fill = TRUE)
      z_out <- if (nrow(z_out) > 0) {
        z_out[, .(val = sum(fill_area, na.rm = TRUE)), by = .(pixel_id, class_id)][val > 0]
      } else {
        data.table(pixel_id = integer(), class_id = character(), val = numeric())
      }

      # 6. Pixel survival: ensure no input pixel disappears entirely.
      #    Strategy: for each missing pixel, transfer effective_min from the
      #    zone's richest donor to the pixel's original dominant class. This is a
      #    strict zone-internal transfer — area is conserved exactly.
      if (preserve_pixels && nrow(z_out) > 0) {
        missing_pix <- setdiff(unique(z_dt$pixel_id), z_out$pixel_id)

        if (length(missing_pix) > 0) {
          # Reconstruct original val (z_dt$val is now clean_val, so add noise_rem back)
          orig_dom <- z_dt[pixel_id %in% missing_pix,
                           .(pixel_id, class_id, orig_val = val + noise_rem)]
          setorder(orig_dom, pixel_id, -orig_val)
          orig_dom <- orig_dom[, .SD[1], by = pixel_id][, .(pixel_id, recv_class = class_id)]

          setorder(z_out, -val)
          n_rescued <- 0L
          n_lost    <- 0L

          for (j in seq_len(nrow(orig_dom))) {
            pix  <- orig_dom$pixel_id[j]
            rcls <- orig_dom$recv_class[j]

            # Donor must have strictly more than effective_min so it stays above zero after transfer
            donor_idx <- which(z_out$val > effective_min)[1]
            if (is.na(donor_idx)) { n_lost <- n_lost + 1L; next }

            # Transfer effective_min from donor to rescued pixel
            z_out[donor_idx, val := val - effective_min]
            existing_idx <- which(z_out$pixel_id == pix & z_out$class_id == rcls)
            if (length(existing_idx) > 0) {
              z_out[existing_idx, val := val + effective_min]
            } else {
              z_out <- rbindlist(list(z_out,
                                      data.table(pixel_id = pix, class_id = rcls, val = effective_min)))
            }
            z_out <- z_out[val > 0]
            setorder(z_out, -val)
            n_rescued <- n_rescued + 1L
          }

          if (n_lost    > 0) cli_alert_warning("Zone {z}: {n_lost} pixel(s) unrecoverable — no donor with val > {effective_min} available.")
          if (n_rescued > 0) cli_alert_info("Zone {z}: {n_rescued} pixel(s) rescued with {effective_min} from zone donor.")
        }
      }

    } else {
      z_out <- z_dt[val > 0, .(pixel_id, class_id, val)]
    }

    if (nrow(z_out) > 0) {
      for (zc in zone_cols) z_out[, (zc) := z_dt[[zc]][1]]
      new_fills_list[[i]] <- z_out
    }

    cli_progress_update()
  }

  cli_alert_success("Zonal sweeps completed. Reconstructing output...")

  if (all(sapply(new_fills_list, is.null))) return(dt[0])

  all_fills <- rbindlist(new_fills_list)

  # Conservation check: warn if zone-level area drift exceeds floating-point tolerance
  total_area_out <- sum(all_fills$val)
  if (abs(total_area_out - total_area_in) > 1e-6 * total_area_in) {
    cli_alert_warning(
      "Area conservation violation: in={round(total_area_in,4)}, out={round(total_area_out,4)}, diff={round(total_area_out-total_area_in,6)}"
    )
  }

  # ── Diagnostics: quantify what actually moved ──────────────────────────────
  # dt_internal still holds the original val (zone subsets were copies), so we
  # can compare input vs output cell by cell.
  diag_cmp <- merge(
    dt_internal[, .(val_in = sum(val)),  by = c(zone_cols, "pixel_id", "class_id")],
    all_fills[,   .(val_out = sum(val)), by = c(zone_cols, "pixel_id", "class_id")],
    by = c(zone_cols, "pixel_id", "class_id"), all = TRUE
  )
  diag_cmp[is.na(val_in),  val_in  := 0]
  diag_cmp[is.na(val_out), val_out := 0]

  moved_area <- sum(abs(diag_cmp$val_out - diag_cmp$val_in)) / 2
  n_changed  <- diag_cmp[abs(val_out - val_in) > 1e-9, .N]
  n_zeroed   <- diag_cmp[val_in > 0 & val_out <= 1e-12, .N]
  n_created  <- diag_cmp[val_in <= 1e-12 & val_out > 0, .N]
  n_pix_chg  <- diag_cmp[abs(val_out - val_in) > 1e-9, uniqueN(pixel_id)]

  cli_alert_info(paste0(
    "Swap diagnostics: {round(moved_area, 2)} km² shifted ",
    "({signif(100 * moved_area / total_area_in, 3)}% of total area) | ",
    "{n_changed} pixel-class cells changed in {n_pix_chg} pixels | ",
    "{n_zeroed} cells zeroed, {n_created} cells created"
  ))

  setnames(all_fills, c("pixel_id", "class_id", "val"), c(pixel_id_col, class_col, area_col))

  # Re-attach pixel-level metadata using composite key (pixel + zone) to handle AEZ splits
  present_meta_cols <- pixel_meta_cols[pixel_meta_cols %in% names(dt)]
  meta_merge_keys   <- unique(c(pixel_id_col, zone_cols))

  if (length(present_meta_cols) > 0) {
    pixel_meta_raw <- unique(dt[, unique(c(meta_merge_keys, present_meta_cols)), with = FALSE])
    n_splits       <- uniqueN(dt[, ..meta_merge_keys])
    if (nrow(pixel_meta_raw) > n_splits) {
      cli_alert_warning("{nrow(pixel_meta_raw) - n_splits} sub-pixel zone(s) have ambiguous metadata — taking first row.")
    }
    meta_dt   <- unique(pixel_meta_raw, by = meta_merge_keys)
    all_fills <- merge(all_fills, meta_dt, by = meta_merge_keys, all.x = TRUE)
  }

  return(all_fills)
}
