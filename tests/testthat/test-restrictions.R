# Restrictions in execute_luc_downscaling(): pair rules, scenario/period keys,
# and unmet-demand diagnostics. One region, two pixels, sources `low` and `mid`
# (10 km2 each per pixel), destination `high`.

run_case <- function(restrictions, target_low = 4, target_mid = 4) {
  out <- withr::local_tempdir(.local_envir = parent.frame())
  files <- list(lu_output   = file.path(out, "lu.csv"),
                luc_output  = file.path(out, "luc.csv"),
                diag_output = file.path(out, "diag.csv"))
  mapping_grid <- data.table::data.table(ns = c("p1", "p2"), Ns_reg = "R1")
  start_map <- data.table::data.table(ns = c("p1", "p1", "p2", "p2"), Ns_reg = "R1",
                                      from_class = c("low", "mid", "low", "mid"),
                                      value = 10)
  targets <- data.table::data.table(Ns_reg = "R1", scenario = "S1", times = "2030",
                                    from_class = c("low", "mid"), to_class = "high",
                                    value = c(target_low, target_mid))
  suppressMessages(suppressWarnings(utils::capture.output(
    execute_luc_downscaling(targets = targets, start_map = start_map,
                            mapping_grid = mapping_grid, restrictions = restrictions,
                            region_cols = "Ns_reg", scen_cols = "scenario",
                            time_col = "times", country_col = "Ns_reg",
                            files = files, base_year = 2020))))
  luc  <- data.table::fread(files$luc_output)
  diag <- data.table::fread(files$diag_output)
  flow <- function(f) luc[from_class == f & to_class == "high" & times == 2030, sum(value)]
  list(low = flow("low"), mid = flow("mid"), diag = diag)
}

pair_rule <- function(...) {
  data.table::data.table(ns = c("p1", "p2"), from_class = "low", to_class = "high",
                         value = 1, ...)
}

test_that("a pair rule blocks only its own source", {
  r <- run_case(pair_rule(scenario = "S1"))
  expect_equal(r$low, 0, tolerance = 1e-6)
  expect_equal(r$mid, 4, tolerance = 1e-6)
})

test_that("the scenario wildcard applies the rule to every scenario", {
  r <- run_case(pair_rule(scenario = PLCHOLD_SCEN))
  expect_equal(r$low, 0, tolerance = 1e-6)
  expect_equal(r$mid, 4, tolerance = 1e-6)
})

test_that("explicit periods work like an omitted times column", {
  r <- run_case(pair_rule(scenario = "S1", times = "2030"))
  expect_equal(r$low, 0, tolerance = 1e-6)
  expect_equal(r$mid, 4, tolerance = 1e-6)
})

test_that("rows without from_class still restrict every source", {
  dest_only <- data.table::data.table(ns = c("p1", "p2"), to_class = "high",
                                      value = 1, scenario = "S1")
  r <- run_case(dest_only)
  expect_equal(r$low, 0, tolerance = 1e-6)
  expect_equal(r$mid, 0, tolerance = 1e-6)
})

test_that("unmatched scenario or period tags fail instead of dropping rows", {
  expect_error(run_case(pair_rule(scenario = "ALL")), "match no target scenario")
  expect_error(run_case(pair_rule(scenario = "S1", times = "ALL")), "match no target period")
})

test_that("non-binary restriction values are rejected before solving", {
  frac <- pair_rule(scenario = "S1"); frac[, value := 0.5]
  expect_error(run_case(frac), "must be 0 or 1")
})

test_that("diagnostics report requested, capped, achieved and unmet area", {
  # low -> high is restricted in p1 only, so 10 km2 (p2) is allowed.
  one_pixel <- pair_rule(scenario = "S1"); one_pixel[ns == "p2", value := 0]
  r <- run_case(one_pixel, target_low = 15, target_mid = 4)
  d_low <- r$diag[from_class == "low" & to_class == "high"]
  d_mid <- r$diag[from_class == "mid" & to_class == "high"]
  expect_equal(d_low$target_requested, 15)
  expect_equal(d_low$target_capped, 10, tolerance = 1e-6)
  expect_equal(d_low$achieved, 10, tolerance = 1e-4)
  expect_equal(d_low$unmet, 5, tolerance = 1e-4)
  expect_true(d_low$restriction_blocking)
  expect_equal(d_mid$unmet, 0, tolerance = 1e-4)
  expect_equal(r$low, 10, tolerance = 1e-4)
})
