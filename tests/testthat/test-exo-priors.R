# Exogenous priors that do not name a source class, destination class or period apply to
# all of them (wildcards), and named classes that occur in no target are reported.

test_that("missing or empty from_class, to_class and times become wildcards", {
  exo <- data.table::data.table(ns = c("p1", "p2", "p3"), value = c(1, 0.5, 0.2),
                                to_class = c("high", NA, ""))
  out <- expect_message(.complete_exo_priors(exo), "applied to all")
  expect_identical(out$from_class, rep("__NA_LU__", 3))
  expect_identical(out$times, rep("__NA_TIME__", 3))
  expect_identical(out$to_class, c("high", "__NA_LU__", "__NA_LU__"))
  expect_false("from_class" %in% names(exo))                       # the caller's table is unchanged
  full <- data.table::data.table(ns = "p1", value = 1, from_class = "low", to_class = "high", times = "2030")
  expect_silent(.complete_exo_priors(full))
  expect_identical(.complete_exo_priors(full)$times, "2030")
})

test_that("classes that occur in no target are reported; ns and value are required", {
  targets <- data.table::data.table(from_class = "low", to_class = "high")
  exo <- data.table::data.table(ns = c("p1", "p2"), value = 1, to_class = c("high", "Forest_protected"))
  expect_warning(suppressMessages(.complete_exo_priors(exo, targets)), "1 of 2 rows can never take effect: Forest_protected")
  expect_error(.complete_exo_priors(data.table::data.table(ns = "p1")), "needs the columns")
})

run_exo <- function(exo = NULL) {
  out <- withr::local_tempdir()
  files <- list(lu_output = file.path(out, "lu.csv"), luc_output = file.path(out, "luc.csv"),
                diag_output = file.path(out, "diag.csv"))
  targets <- data.table::data.table(Ns_reg = "R1", scenario = "S1", times = "2030",
                                    from_class = "low", to_class = "high", value = 4)
  start_map <- data.table::data.table(ns = c("p1", "p2"), Ns_reg = "R1", from_class = "low", value = c(10, 10))
  mapping_grid <- data.table::data.table(ns = c("p1", "p2"), Ns_reg = "R1")
  suppressMessages(suppressWarnings(utils::capture.output(
    execute_luc_downscaling(targets = targets, start_map = start_map, mapping_grid = mapping_grid,
                            exo_priors_all = exo, files = files, base_year = 2020))))
  lu <- data.table::fread(files$lu_output)
  lu[times == 2030 & lu == "high", .(area = sum(value)), keyby = ns]
}

test_that("an exogenous prior without source class or period steers the allocation", {
  base <- run_exo()
  # only the destination is named: the prior applies to every source class and period
  steered <- run_exo(data.table::data.table(ns = c("p1", "p2"), to_class = "high", value = c(1e-6, 1)))
  expect_equal(sum(steered$area), 4, tolerance = 1e-6)               # the target is still met
  expect_gt(steered[ns == "p2", area], base[ns == "p2", area] + 0.5)  # and goes where the prior points
})
