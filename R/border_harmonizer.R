# ===============================================================
# border_harmonizer.R
# ===============================================================
# Within-zone gravity redistribution of land-use class area.
#
# WHAT GETS TAKEN:
#   Thin pixels — those with class area below thin_threshold —
#   are the explicit donors. Each gives blend × s_p of its area.
#   Thick pixels are left untouched regardless of position.
#
# WHERE IT GOES:
#   The pooled area is distributed back into the zone, weighted
#   by g_p × s_p:  g_p is the gravitational pull from all
#   neighboring zones where the class is strongly present
#   (magnets), s_p is the pixel's existing class area (so we
#   never create presence where there was none). The result is
#   that thin interior scatter is drawn toward border pixels
#   facing strong neighboring presence.
#
# PIXEL BALANCE:
#   Every change to class X in a pixel is inverted across the
#   complement classes in the same aggregate group, keeping
#   pixel-level area totals constant. Requires class_aggregate.
#
# Parameters
# ----------
# thin_threshold   area below which a pixel is a donor (same
#                  units as area_col). Set to the noise threshold
#                  used in zonal_swapper for consistency.
# gravity_decay    km half-life of the gravitational pull.
# blend            fraction of each thin pixel's area to move
#                  (0 = nothing, 1 = drain thin pixels fully).
# magnet_threshold minimum mean class density in a neighbor's
#                  border pixels to treat it as a magnet.
# class_aggregate  2-col data.table (class_col, agg_group) or
#                  named character vector. Required for pixel-
#                  total conservation. NULL = no balancing.
# ===============================================================
library(data.table)
library(cli)

#' Focal Border Harmonization for Zonal Aggregation
#'
#' @description
#' Smoothes rigid administrative zone boundaries by reallocating target areas across
#' zone borders according to a spatial distance decay function. This prevents hard
#' thematic discontinuities at administrative borders.
#'
#' @param dt A `data.table` of pixel-level thematic areas.
#' @param zone_cols Character vector of column names defining spatial zones.
#' @param pixel_id_col Character string of the pixel ID column name.
#' @param class_col Character string of the thematic class column name.
#' @param area_col Character string of the area value column name.
#' @param coord_cols Character vector of X and Y coordinate column names.
#' @param border_depth_km Numeric depth in kilometers to define the active border region.
#' @param gravity_decay Numeric half-life decay parameter for the gravitational pull.
#' @param blend Numeric fraction (0 to 1) defining the blend intensity across borders.
#' @param magnet_threshold Numeric threshold to qualify a neighboring zone as a pull magnet.
#' @param class_aggregate Optional grouping to conserve overall pixel area totals.
#'
#' @return A `data.table` with smoothed and reallocated areas.
#' @export
harmonize_zone_borders <- function(dt,
                                   zone_cols            = c("GLOB_country", "GLOB_LUID"),
                                   pixel_id_col         = "ns",
                                   class_col            = "GLOB_lu_class",
                                   area_col             = "LUM_area_km2",
                                   coord_cols           = c("x", "y"),
                                   pixel_size           = 10,
                                   thin_threshold       = 0.1,
                                   gravity_decay        = 30,
                                   blend                = 0.8,
                                   magnet_threshold     = 0,
                                   class_aggregate      = NULL,
                                   classes_to_harmonize = NULL) {

  cli_alert_info(
    "Border gravity redistribution — thin<{thin_threshold}, decay={gravity_decay}km, blend={blend}"
  )
  if (is.null(class_aggregate))
    cli_alert_warning("class_aggregate = NULL: pixel totals will not be conserved.")

  decay_steps <- gravity_decay / pixel_size

  # ── 0. Working copy ─────────────────────────────────────────────────────────
  need_cols <- unique(c(zone_cols, pixel_id_col, coord_cols, class_col, area_col))
  dt_w      <- copy(dt[, intersect(need_cols, names(dt)), with = FALSE])
  setnames(dt_w,
           c(pixel_id_col, coord_cols[1], coord_cols[2], class_col, area_col),
           c("pid",        "px",          "py",          "cls",      "val"))

  if (!is.null(classes_to_harmonize)) dt_w <- dt_w[cls %in% classes_to_harmonize]

  dt_w[, zone_id := .GRP, by = zone_cols]
  dt_w[, val_new := val]
  total_area_in <- sum(dt_w$val)

  # ── Calibration diagnostic: how much of the data is donor-eligible? ────────
  pos_vals <- dt_w[val > 0, val]
  qv <- quantile(pos_vals, probs = c(0.10, 0.25, 0.50, 0.75, 0.90))
  pct_eligible  <- 100 * mean(pos_vals < thin_threshold)
  area_eligible <- sum(pos_vals[pos_vals < thin_threshold])
  cli_alert_info(paste0(
    "Cell-area quantiles (>0): p10={signif(qv[1],3)} p25={signif(qv[2],3)} ",
    "p50={signif(qv[3],3)} p75={signif(qv[4],3)} p90={signif(qv[5],3)} km²"
  ))
  cli_alert_info(paste0(
    "thin_threshold={thin_threshold}: {round(pct_eligible,1)}% of cells donor-eligible, ",
    "holding {round(area_eligible,1)} km² ",
    "({signif(100*area_eligible/total_area_in,3)}% of total area) — ",
    "max movable = blend x that = {round(blend*area_eligible,1)} km²"
  ))

  agg_map <- NULL
  if (!is.null(class_aggregate)) {
    agg_map <- if (is.data.table(class_aggregate)) {
      setNames(class_aggregate[[2]], class_aggregate[[1]])
    } else {
      class_aggregate
    }
  }

  pix_meta <- unique(dt_w[val > 0, .(pid, px, py, zone_id)], by = "pid")
  setkey(pix_meta, px, py)

  # ── 1. 8-connected cross-zone adjacency ─────────────────────────────────────
  offs <- CJ(dx = c(-pixel_size, 0L, pixel_size),
             dy = c(-pixel_size, 0L, pixel_size))[!(dx == 0L & dy == 0L)]

  cross_adj <- rbindlist(lapply(seq_len(nrow(offs)), function(k) {
    nb <- pix_meta[, .(pA = pid, zA = zone_id,
                        nx = px + offs$dx[k], ny = py + offs$dy[k])]
    m  <- pix_meta[nb, on = .(px = nx, py = ny), nomatch = NULL]
    m[zone_id != zA, .(pA, zA, pB = pid, zB = zone_id)]
  }))
  cross_adj <- unique(cross_adj)

  if (nrow(cross_adj) == 0) {
    cli_alert_warning("No cross-zone adjacencies found — check coord_cols and pixel_size.")
    return(invisible(dt))
  }

  # Border index: for each (zone, neighbor), which zone pixels face the neighbor?
  border_idx <- rbindlist(list(
    cross_adj[, .(zone = zA, neigh = zB, bp_zone = pA, bp_neigh = pB)],
    cross_adj[, .(zone = zB, neigh = zA, bp_zone = pB, bp_neigh = pA)]
  ))
  border_idx <- unique(border_idx)

  all_zones <- unique(dt_w$zone_id)
  cli_alert_info(
    "{length(all_zones)} zones, {nrow(unique(border_idx[, .(zone, neigh)]))} directed border pairs"
  )

  # ── 2. Per-zone, per-class redistribution ───────────────────────────────────
  setkeyv(dt_w, c("pid", "cls"))

  n_moved       <- 0L
  n_pool_pix    <- 0L   # donor pixel-class instances drained
  n_skip_no_g   <- 0L   # class×zone skipped: no gravity signal
  n_skip_nopool <- 0L   # class×zone skipped: no donor-eligible pixels
  cli_progress_bar("Zones", total = length(all_zones))

  for (z in all_zones) {
    pix_z  <- pix_meta[zone_id == z, .(pid, px, py)]
    if (nrow(pix_z) == 0) { cli_progress_update(); next }

    neighs <- unique(border_idx[zone == z, neigh])
    if (length(neighs) == 0) { cli_progress_update(); next }

    classes_z <- unique(dt_w[zone_id == z & val_new > 0, cls])
    if (!is.null(classes_to_harmonize)) classes_z <- intersect(classes_z, classes_to_harmonize)

    for (cls_j in classes_z) {

      # Current class values for this zone
      s_dt <- dt_w[zone_id == z & cls == cls_j, .(pid, s = val_new)]
      if (nrow(s_dt) == 0 || sum(s_dt$s) == 0) next

      # ── 2a. Gravity: cross-zone + internal, additive ────────────────────────
      # Both sources contribute to the same field. A pixel near a strong
      # neighbor border AND near an internal dense pixel gets double gravity,
      # which is correct — it is the most reinforced location for the class.
      g_dt <- data.table(pid = pix_z$pid, g = 0.0)

      # Source 1: cross-zone magnets (neighbor border pixels with class present)
      for (znb in neighs) {
        bp_neigh  <- unique(border_idx[zone == z & neigh == znb, bp_neigh])
        strength  <- dt_w[pid %in% bp_neigh & cls == cls_j, mean(val_new)]
        if (is.na(strength) || strength <= magnet_threshold) next

        bp_zone_side <- unique(border_idx[zone == z & neigh == znb, bp_zone])
        b_coords     <- pix_meta[pid %in% bp_zone_side, .(bx = px, by = py)]
        if (nrow(b_coords) == 0) next

        d_dt <- pix_z[b_coords, on = character(), allow.cartesian = TRUE][
          , .(d = min(pmax(abs(px - bx), abs(py - by)) / pixel_size)), by = pid]

        g_dt[d_dt, on = "pid", g := g + strength * exp(-d / decay_steps)]
      }

      # Source 2: internal concentration magnets (pixels with s above zone mean)
      # Handles classes with no cross-zone signal by pulling scattered thin
      # pixels toward the zone's own existing dense pixels.
      s_mean_internal <- s_dt[s > 0, mean(s)]
      centers         <- s_dt[s > s_mean_internal]
      if (nrow(centers) > 0) {
        strength_int <- centers[, mean(s)]
        if (strength_int > magnet_threshold) {
          c_coords <- pix_meta[pid %in% centers$pid, .(bx = px, by = py)]
          d_dt_int <- pix_z[c_coords, on = character(), allow.cartesian = TRUE][
            , .(d = min(pmax(abs(px - bx), abs(py - by)) / pixel_size)), by = pid]
          g_dt[d_dt_int, on = "pid", g := g + strength_int * exp(-d / decay_steps)]
        }
      }

      if (all(g_dt$g == 0)) { n_skip_no_g <- n_skip_no_g + 1L; next }

      # ── 2b. Split thin pixels by gravity ────────────────────────────────────
      # Pool  = thin pixels with g < g_mean: scattered, far from any magnet.
      # Thin pixels with g >= g_mean sit near a border or internal dense core
      # — gravity reinforces them, so they are receivers, not donors.
      mix <- s_dt[g_dt, on = "pid"]
      mix[is.na(s), s := 0]

      g_mean <- mix[s > 0, weighted.mean(g, w = s)]
      if (is.na(g_mean) || g_mean == 0) next

      pool_pix <- mix[s > 0 & s < thin_threshold & g < g_mean]
      if (nrow(pool_pix) == 0) { n_skip_nopool <- n_skip_nopool + 1L; next }

      pool_pix[, donation := blend * s]
      T <- sum(pool_pix$donation)
      if (T < 1e-12) next

      # ── 2c. Receiver weights: gravity × existing presence ───────────────────
      # All pixels with g > 0 are eligible; weighting by s ensures we only
      # reinforce existing presence (never seed a class from zero).
      mix[, recv_w := g * s]
      total_rw <- sum(mix$recv_w)
      if (total_rw == 0) next

      mix[, addition := T * recv_w / total_rw]

      # ── 2d. Net per-pixel delta ─────────────────────────────────────────────
      # Pool pixels: lose donation, may also receive if they have g > 0
      # (net could be positive if they're slightly thin but still above g_mean —
      # but those are excluded from the pool, so pool pixels always have g < g_mean
      # and thus low recv_w → net is negative for them as intended).
      delta_dt <- merge(
        pool_pix[, .(pid, delta = -donation)],
        mix[addition > 0, .(pid, delta = addition)],
        by = "pid", all = TRUE, suffixes = c("_d", "_r")
      )
      delta_dt[is.na(delta_d), delta_d := 0]
      delta_dt[is.na(delta_r), delta_r := 0]
      delta_dt[, delta := delta_d + delta_r]

      # Compute final values and floor at 0
      final_vals <- s_dt[delta_dt[, .(pid, delta)], on = "pid"]
      final_vals[is.na(s), s := 0]
      final_vals[, val_new := pmax(0, s + delta)]

      # Renormalize to guarantee exact zone-total conservation
      # (floating point and pmax(0) may introduce tiny drift)
      scale <- sum(s_dt$s) / sum(final_vals$val_new)
      final_vals[, val_new := val_new * scale]
      final_vals[, delta   := val_new - s]

      if (max(abs(final_vals$delta)) < 1e-12) next

      # Write primary class
      dt_w[final_vals[, .(pid, cls = cls_j, val_new)], val_new := i.val_new]

      # ── 2e. Complement compensation ─────────────────────────────────────────
      if (!is.null(agg_map) && cls_j %in% names(agg_map)) {
        grp      <- agg_map[[cls_j]]
        comp_cls <- names(agg_map)[agg_map == grp & names(agg_map) != cls_j]

        if (length(comp_cls) > 0) {
          comp <- dt_w[zone_id == z & cls %in% comp_cls, .(pid, cls, val_new)]
          comp[final_vals[, .(pid, delta)], on = "pid", delta := i.delta]
          comp[is.na(delta), delta := 0]
          comp[, grp_tot := sum(val_new), by = pid]
          # Invert: where cls_j grew (delta > 0), complement shrinks, and vice versa
          comp[, val_new := pmax(0,
                                  val_new - fifelse(grp_tot > 0,
                                                    val_new / grp_tot * delta,
                                                    0))]
          dt_w[comp[, .(pid, cls, val_new)], val_new := i.val_new]
        }
      }

      n_moved    <- n_moved + 1L
      n_pool_pix <- n_pool_pix + nrow(pool_pix)
    }

    cli_progress_update()
  }

  cli_alert_success("{n_moved} class redistribution(s) applied")
  cli_alert_info(paste0(
    "Skipped class x zone combos: {n_skip_no_g} without gravity signal, ",
    "{n_skip_nopool} without donor-eligible pixels"
  ))

  total_area_out <- sum(dt_w$val_new)
  if (abs(total_area_out - total_area_in) > 1e-6 * total_area_in)
    cli_alert_warning("Area drift: {round(total_area_out - total_area_in, 6)}")

  # ── Movement diagnostics: input vs output cell by cell ──────────────────────
  moved_area <- dt_w[, sum(abs(val_new - val))] / 2
  n_cells    <- dt_w[abs(val_new - val) > 1e-9, .N]
  n_pix      <- dt_w[abs(val_new - val) > 1e-9, uniqueN(pid)]
  cli_alert_info(paste0(
    "Harmonizer diagnostics: {round(moved_area, 2)} km² shifted ",
    "({signif(100 * moved_area / total_area_in, 3)}% of total area) | ",
    "{n_cells} pixel-class cells changed in {n_pix} pixels | ",
    "{n_pool_pix} donor drain events"
  ))

  diag_cls <- dt_w[, .(moved_km2 = round(sum(abs(val_new - val)) / 2, 2),
                       n_cells   = sum(abs(val_new - val) > 1e-9)),
                   by = cls][moved_km2 > 0][order(-moved_km2)]
  if (nrow(diag_cls) > 0) {
    cli_alert_info("Top classes by shifted area:")
    print(head(diag_cls, 10))
  }

  # ── 3. Write back ────────────────────────────────────────────────────────────
  val_upd <- dt_w[, .(pid, cls, val_new)]
  setnames(val_upd, c("pid", "cls"), c(pixel_id_col, class_col))
  setkeyv(val_upd, c(pixel_id_col, class_col))

  dt_out <- copy(dt)
  setkeyv(dt_out, c(pixel_id_col, class_col))
  dt_out[val_upd, (area_col) := i.val_new]

  return(dt_out)
}
