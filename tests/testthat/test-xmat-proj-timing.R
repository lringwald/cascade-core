# Projected (exogenous) covariates are contemporaneous with the period being solved.
# The prior models are fit on covariates of the outcome year (gamble-core DRIVER_COV_YEARS), so the
# solve for period t must use xmat_proj[times == t]. Before the fix, downscale() swapped them in
# after solving t, so period t ran on t-1's covariates (the first period on the base xmat).

analytical_single_logit <- function(u, areas, target) {
  fn <- function(log_x) sum(areas * exp(u + log_x) / (1 + exp(u + log_x))) - target
  r <- uniroot(fn, interval = c(-50, 50), tol = 1e-12)$root
  areas * exp(u + r) / (1 + exp(u + r))
}

run_two_period <- function() {
  start_areas <- data.table::data.table(ns = c("p1", "p2"), from_class = "low", value = c(10, 10))
  targets <- data.table::data.table(times = c("2020", "2030"), from_class = "low", to_class = "high", value = c(2, 2))
  xmat <- data.table::data.table(ns = c("p1", "p2"), ks = "k1", value = c(0, 0))          # neutral base
  xmat_proj <- data.table::data.table(ns = c("p1", "p2", "p1", "p2"), ks = "k1",
                                      times = c("2020", "2020", "2030", "2030"),
                                      value = c(5, 1, 1, 5))                              # p1 first, then p2
  betas <- data.table::data.table(from_class = "low", to_class = "high", ks = "k1", value = 1)
  utils::capture.output(out <- suppressMessages(suppressWarnings(downscale(
    targets = targets, start_areas = start_areas, times = c("2020", "2030"),
    xmat = xmat, betas = betas,
    xmat_coltypes = data.table::data.table(ks = "k1", value = "projected"),
    xmat_proj = xmat_proj,
    options = downscale_control(verbose = FALSE)
  ))))
  r <- out$out_res[from_class == "low" & to_class == "high"]
  r[, times := as.character(times)]
  data.table::setorder(r, times, ns)
  r
}

test_that("each period is solved with its own projected covariates", {
  r <- run_two_period()
  a2020 <- r[times == "2020", value]
  expect_equal(a2020, analytical_single_logit(c(5, 1), c(10, 10), 2), tolerance = 1e-4)
  expect_gt(a2020[1], a2020[2])                       # the old timing split 2020 evenly (base xmat)
  a2030 <- r[times == "2030", value]
  expect_equal(a2030, analytical_single_logit(c(1, 5), c(10, 10) - a2020, 2), tolerance = 1e-4)
  expect_gt(a2030[2], a2030[1])
})

test_that("a period without projection rows stops instead of keeping old covariates", {
  start_areas <- data.table::data.table(ns = c("p1", "p2"), from_class = "low", value = c(10, 10))
  targets <- data.table::data.table(times = c("2020", "2030"), from_class = "low", to_class = "high", value = c(2, 2))
  xmat <- data.table::data.table(ns = c("p1", "p2"), ks = "k1", value = c(0, 0))
  xmat_proj <- data.table::data.table(ns = c("p1", "p2"), ks = "k1", times = "2020", value = c(5, 1))
  betas <- data.table::data.table(from_class = "low", to_class = "high", ks = "k1", value = 1)
  expect_error(suppressMessages(suppressWarnings(utils::capture.output(downscale(
    targets = targets, start_areas = start_areas, times = c("2020", "2030"),
    xmat = xmat, betas = betas,
    xmat_coltypes = data.table::data.table(ks = "k1", value = "projected"),
    xmat_proj = xmat_proj, options = downscale_control(verbose = FALSE))))),
    "xmat_proj")
})
