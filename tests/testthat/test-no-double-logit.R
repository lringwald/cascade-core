# Regression test for single MNL logit normalization in downscale() / solve_biascorr_lu()
# Before fixing R/solve_biascorr_lu.R, curr_priors was normalized by (1 + rowSums(exp(U)))
# BEFORE passing into allocate_mnl(), which divides by (1 + rowSums(...)) a second time,
# flattening large utility differences across pixels.

analytical_mnl_single_logit <- function(u, areas, target) {
  fn <- function(log_x) {
    odds <- exp(u + log_x)
    sum(areas * odds / (1 + odds)) - target
  }
  root <- uniroot(fn, interval = c(-50, 50), tol = 1e-12)$root
  odds <- exp(u + root)
  areas * odds / (1 + odds)
}

run_two_pixel_mnl <- function(optimize_log_scale) {
  start_areas <- data.table::data.table(
    ns = c("p1", "p2"),
    from_class = "low",
    value = c(10, 10)
  )
  targets <- data.table::data.table(
    times = "2030",
    from_class = "low",
    to_class = "high",
    value = 10
  )
  xmat <- data.table::data.table(
    ns = c("p1", "p2"),
    ks = "k1",
    value = c(5, 1)
  )
  betas <- data.table::data.table(
    from_class = "low",
    to_class = "high",
    ks = "k1",
    value = 1
  )
  opts <- downscale_control(
    optimize_log_scale = optimize_log_scale,
    verbose = FALSE
  )
  res <- suppressMessages(suppressWarnings(utils::capture.output(
    out <- downscale(
      targets = targets,
      start_areas = start_areas,
      times = "2030",
      xmat = xmat,
      betas = betas,
      options = opts
    )
  )))
  dt <- out$out_res[from_class == "low" & to_class == "high"]
  data.table::setorder(dt, ns)
  dt$value
}

test_that("downscale matches exact single-logit MNL shares (no double logit)", {
  expected <- analytical_mnl_single_logit(u = c(5, 1), areas = c(10, 10), target = 10)

  got_log <- run_two_pixel_mnl(optimize_log_scale = TRUE)
  expect_equal(got_log, expected, tolerance = 1e-4)

  got_lin <- run_two_pixel_mnl(optimize_log_scale = FALSE)
  expect_equal(got_lin, expected, tolerance = 1e-4)
})

