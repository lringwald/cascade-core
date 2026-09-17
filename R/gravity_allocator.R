# ===============================================================
# gravity_allocator.R
# ===============================================================
# Spatially weighted allocation of zone-level land-use transition
# flows to pixels. Replaces uniform share application
# (area_p x share) with a transportation problem solved by IPF
# for every (zone, from_class):
#
#   row margins  r_p = pixel area of from_class
#                      -> pixel totals conserved exactly
#   col margins  c_B = share(A->B) x sum_p(r_p)
#                      -> zone flow budgets hit by construction
#   IPF seed     M0[p,B] = gravity_p(B) + gravity_floor
#
# gravity_p(B) is the neighborhood density of target class B
# around pixel p on the current (pre-transition) map, computed on
# coarse blocks with exponential ring decay.
# WITH TARGET SMOOTHING:
# The IPF seed is now additionally weighted by a smoothed field of 
# the administrative zone shares (A->B target share). This prevents
# allocations from clustering right against the border of neighboring
# zones that have zero target for that transition.
# ===============================================================
library(data.table)
library(cli)

#' Gravity-Weighted Spatial Allocation via Iterative Proportional Fitting
#'
#' @description
#' Spatially allocates zone-level transition flows to pixels using Iterative Proportional Fitting (IPF).
#' This prevents hard discontinuities at zone borders by using a continuous Gaussian gravity field
#' based on the spatial neighborhood density of the target class.
#'
#' @param pixel_dt A `data.table` of initial pixel-level areas containing coordinates.
#' @param shares_dt A `data.table` of zone-level transition shares.
#' @param zone_cols Character vector of column names defining the spatial zones.
#' @param pixel_id_col Character string of the pixel ID column name.
#' @param from_col Character string of the starting class column name.
#' @param to_col Character string of the target class column name.
#' @param share_col Character string of the zone share value column name.
#' @param area_col Character string of the pixel area value column name.
#' @param coord_cols Character vector of the X and Y coordinate column names.
#' @param pixel_size Numeric pixel resolution in kilometers.
#' @param block_pixels Integer block size for rasterization speedup.
#' @param gravity_reach Numeric reach of the Gaussian gravity field.
#' @param gravity_floor Numeric baseline probability to allow transitions in empty regions.
#' @param ipf_iter Integer maximum number of IPF iterations.
#' @param ipf_tol Numeric convergence tolerance for IPF.
#'
#' @return A `data.table` of pixel-level transition flows.
#' @export
allocate_transitions_spatial <- function(pixel_dt,
                                         shares_dt,
                                         zone_cols     = c("GLOB_country", "GLOB_LUID", "GLOB_AEZ"),
                                         pixel_id_col  = "ns",
                                         from_col      = "from_class",
                                         to_col        = "to_class",
                                         share_col     = "value",
                                         area_col      = "LUM_area_km2",
                                         coord_cols    = c("x_ns", "y_ns"),
                                         pixel_size    = 10,
                                         block_pixels  = 4,
                                         gravity_reach = 2,
                                         gravity_floor = 0.05,
                                         rare_area_km2 = 0,
                                         ipf_iter      = 50,
                                         ipf_tol       = 1e-9) {

  grid_res_meters <- pixel_size * 1000
  radius_meters <- (gravity_reach + 0.5) * block_pixels * grid_res_meters

  cli_alert_info(paste0(
    "Continuous Gravity-weighted allocation WITH TARGET SMOOTHING — {grid_res_meters/1000} km resolution, ",
    "Gaussian reach ~{radius_meters/1000} km, floor {gravity_floor}"
  ))

  # ── 0. Internal working copies ──────────────────────────────────────────────
  px <- pixel_dt[, unique(c(zone_cols, pixel_id_col, coord_cols, from_col, area_col)), with = FALSE]
  setnames(px,
           c(pixel_id_col, coord_cols[1], coord_cols[2], from_col, area_col),
           c("pid",        "cx",          "cy",          "fcls",   "val"))
  px <- px[val > 0, .(val = sum(val)), by = c(zone_cols, "pid", "cx", "cy", "fcls")]

  n_nocoord <- px[is.na(cx) | is.na(cy), uniqueN(pid)]
  if (n_nocoord > 0)
    cli_alert_warning("{n_nocoord} pixel(s) without coordinates — they get uniform (floor-only) gravity.")

  sh <- shares_dt[, unique(c(zone_cols, from_col, to_col, share_col)), with = FALSE]
  setnames(sh, c(from_col, to_col, share_col), c("fcls", "tcls", "shr"))
  sh <- sh[!is.na(shr) & shr > 0]
  sh[, shr := shr / sum(shr), by = c(zone_cols, "fcls")]

  total_in <- sum(px$val)

  # ── 1. Gravity field: continuous focal class densities ──────────────────────
  px[, `:=`(x_snap = floor(cx / grid_res_meters) * grid_res_meters + grid_res_meters / 2,
            y_snap = floor(cy / grid_res_meters) * grid_res_meters + grid_res_meters / 2)]

  blk_d <- px[!is.na(x_snap), .(dens = sum(val)), by = .(x_snap, y_snap, cls = fcls)]

  x_rng <- range(px$x_snap, na.rm = TRUE)
  y_rng <- range(px$y_snap, na.rm = TRUE)
  full_grid <- CJ(x_snap = seq(x_rng[1], x_rng[2], by = grid_res_meters), 
                  y_snap = seq(y_rng[1], y_rng[2], by = grid_res_meters))

  cls_list <- unique(blk_d$cls)
  res_list <- vector("list", length(cls_list))
  
  r_ext <- terra::rast(terra::ext(x_rng[1]-grid_res_meters/2, x_rng[2]+grid_res_meters/2,
                                  y_rng[1]-grid_res_meters/2, y_rng[2]+grid_res_meters/2),
                       res=grid_res_meters)
  # focalMat Gauss semantics: d = c(sigma, window_radius). Sigma must be the
  # SMALLER value or the kernel degenerates to a flat truncated disk.
  # terra::focal requires the window to be at most twice the raster dimension;
  # cap the radius for small extents.
  min_dim   <- min(terra::nrow(r_ext), terra::ncol(r_ext))
  max_win   <- 2L * min_dim - 1L
  des_win   <- 2L * ceiling(radius_meters / grid_res_meters) + 1L
  radius_eff <- if (des_win > max_win) {
    capped <- (max_win - 1L) %/% 2L * grid_res_meters
    cli_alert_warning("Gravity radius {radius_meters/1000} km exceeds raster extent — capped to {capped/1000} km.")
    capped
  } else {
    radius_meters
  }
  w <- terra::focalMat(r_ext, d = c(radius_eff / 3, radius_eff), type = "Gauss")

  # Coverage denominator: weighted share of the kernel that lies on land cells.
  # Dividing by it turns focal sums into focal means over land, so coastal and
  # map-edge pixels are not systematically attenuated versus interior pixels.
  cov_dt <- unique(px[!is.na(x_snap), .(x_snap, y_snap)])[, z := 1]
  z_cov  <- merge(full_grid, cov_dt, by = c("x_snap", "y_snap"), all.x = TRUE)
  z_cov[is.na(z), z := 0]
  r_cov   <- terra::rast(z_cov, type = "xyz", digits = 1)
  cov_s   <- terra::focal(r_cov, w, fun = "sum", na.rm = TRUE)
  cov_div <- terra::ifel(cov_s > 0, cov_s, NA)

  for (j in seq_along(cls_list)) {
    c_name <- cls_list[j]
    z_dt <- blk_d[cls == c_name, .(x_snap, y_snap, z = dens)]
    z_full <- merge(full_grid, z_dt, by = c("x_snap", "y_snap"), all.x = TRUE)
    z_full[is.na(z), z := 0]

    r <- terra::rast(z_full, type = "xyz", digits = 1)
    r_smooth <- terra::focal(r, w, fun = "sum", na.rm = TRUE) / cov_div

    s_dt <- as.data.table(r_smooth, xy = TRUE)
    setnames(s_dt, c("x", "y", "g"))
    s_dt[, `:=`(x = round(x), y = round(y))]
    s_dt <- s_dt[g > 0]
    s_dt[, cls := c_name]
    res_list[[j]] <- s_dt
  }
  
  blk_g <- rbindlist(res_list)
  blk_g[, g := g / mean(g[g > 0]), by = cls]
  setkey(blk_g, cls, x, y)

  # ── 1.5 Gravity field of target shares ──────────────────────────────────────
  # Smooth the zone-level share of target transitions so they drop off near
  # borders of zones that have 0 transition target.
  px[, zid := .GRP, by = zone_cols]
  zmap <- unique(px[, c("zid", zone_cols), with = FALSE])
  sh   <- merge(sh, zmap, by = zone_cols)
  
  ft_pairs <- unique(sh[, .(fcls, tcls)])
  res_shr_list <- vector("list", nrow(ft_pairs))
  
  for (j in seq_len(nrow(ft_pairs))) {
    f_c <- ft_pairs$fcls[j]
    t_c <- ft_pairs$tcls[j]
    
    # We only care about smoothing where the from_class exists
    b_px <- px[fcls == f_c & !is.na(x_snap), .(x_snap, y_snap, zid)]
    s_sub <- sh[fcls == f_c & tcls == t_c, .(zid, shr)]
    b_px <- merge(b_px, s_sub, by = "zid", all.x = TRUE)
    b_px[is.na(shr), shr := 0]
    
    z_dt <- b_px[, .(z = mean(shr)), by = .(x_snap, y_snap)]
    z_full <- merge(full_grid, z_dt, by = c("x_snap", "y_snap"), all.x = TRUE)
    z_full[is.na(z), z := 0]

    r <- terra::rast(z_full, type = "xyz", digits = 1)
    r_smooth <- terra::focal(r, w, fun = "sum", na.rm = TRUE) / cov_div

    s_dt <- as.data.table(r_smooth, xy = TRUE)
    setnames(s_dt, c("x", "y", "shr_s"))
    s_dt[, `:=`(x = round(x), y = round(y), fcls = f_c, tcls = t_c)]
    res_shr_list[[j]] <- s_dt[shr_s > 0]
  }
  
  blk_shr <- rbindlist(res_shr_list)
  setkey(blk_shr, fcls, tcls, x, y)

  # ── 2. Per-(zone, from_class) IPF allocation ────────────────────────────────
  setkey(px, zid)
  setkey(sh, zid)

  zone_ids <- unique(px$zid)
  out_list <- vector("list", length(zone_ids))

  tot_offdiag <- 0.0
  max_col_err <- 0.0
  n_ipf       <- 0L

  cli_progress_bar("Allocating zones", total = length(zone_ids))

  for (i in seq_along(zone_ids)) {
    z    <- zone_ids[i]
    z_px <- px[.(z)]
    z_sh <- if (z %in% sh$zid) sh[.(z)] else sh[0]

    fcls_shared <- unique(z_sh$fcls)

    res_z <- list(z_px[!fcls %in% fcls_shared, .(pid, tcls = fcls, val)])

    for (f in fcls_shared) {
      rows <- z_px[fcls == f]
      if (nrow(rows) == 0) next

      sv   <- z_sh[fcls == f, .(tcls, shr)]
      r    <- rows$val
      ctot <- sum(r)
      cvec <- sv$shr * ctot
      npix <- nrow(rows)
      ntgt <- nrow(sv)

      if (ntgt == 1) {
        res_z[[length(res_z) + 1L]] <- data.table(pid = rows$pid, tcls = sv$tcls, val = r)
        if (sv$tcls != f) tot_offdiag <- tot_offdiag + ctot
        next
      }

      G <- vapply(sv$tcls, function(tc) {
        gv <- blk_g[.(tc, round(rows$x_snap), round(rows$y_snap)), g]
        if (all(is.na(gv))) {
          gv_parent <- blk_g[.(f, round(rows$x_snap), round(rows$y_snap)), g]
          if (!all(is.na(gv_parent))) gv <- gv_parent
        }
        fifelse(is.na(gv), 0, gv)
      }, numeric(npix))
      G <- matrix(G, nrow = npix)
      
      S_shr <- vapply(sv$tcls, function(tc) {
        sv_s <- blk_shr[.(f, tc, round(rows$x_snap), round(rows$y_snap)), shr_s]
        fifelse(is.na(sv_s), 0, sv_s)
      }, numeric(npix))
      S_shr <- matrix(S_shr, nrow = npix)

      # For rare classes (tiny area in this zone), the target-share smoothing
      # field creates artificial stripes at zone borders because the Gaussian
      # gradient dominates the near-zero gravity signal. Skip S_shr modulation
      # and fall back to pure gravity + floor for these cases.
      if (rare_area_km2 > 0 && ctot < rare_area_km2) {
        M <- G + gravity_floor
      } else {
        M <- (G + gravity_floor) * (S_shr + 1e-6)
      }

      for (it in seq_len(ipf_iter)) {
        M  <- M * (r / rowSums(M))
        cs <- colSums(M)
        M  <- M * rep(cvec / cs, each = npix)
        if (max(abs(rowSums(M) - r) / pmax(r, 1e-12)) < ipf_tol) break
      }
      
      M <- M * (r / rowSums(M))
      max_col_err <- max(max_col_err, max(abs(colSums(M) - cvec)) / max(ctot, 1e-12))
      n_ipf <- n_ipf + 1L

      tot_offdiag <- tot_offdiag + sum(M[, sv$tcls != f, drop = FALSE])

      res_z[[length(res_z) + 1L]] <- data.table(
        pid  = rep(rows$pid, times = ntgt),
        tcls = rep(sv$tcls,  each  = npix),
        val  = as.vector(M)
      )[val > 1e-12]
    }

    z_out <- rbindlist(res_z, use.names = TRUE)
    if (nrow(z_out) > 0) {
      z_out <- z_out[, .(val = sum(val)), by = .(pid, tcls)]
      z_out[, zid := z]
      out_list[[i]] <- z_out
    }

    cli_progress_update()
  }

  out <- rbindlist(out_list, use.names = TRUE)
  out <- merge(out, zmap, by = "zid")
  out[, zid := NULL]

  # ── 3. Diagnostics & conservation ───────────────────────────────────────────
  total_out <- sum(out$val)
  if (abs(total_out - total_in) > 1e-6 * total_in)
    cli_alert_warning("Allocation area drift: {round(total_out - total_in, 6)} km²")

  cli_alert_info(paste0(
    "Allocation diagnostics: {n_ipf} IPF problems solved | ",
    "{round(tot_offdiag, 1)} km² transitioned off-diagonal ",
    "({signif(100 * tot_offdiag / total_in, 3)}% of total area) | ",
    "worst flow-budget deviation {signif(max_col_err, 3)} (relative)"
  ))

  setnames(out, c("pid", "tcls", "val"), c(pixel_id_col, to_col, area_col))
  setcolorder(out, c(zone_cols, pixel_id_col, to_col, area_col))
  return(out[])
}
