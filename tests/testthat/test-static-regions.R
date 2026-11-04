# Start-map regions without targets (e.g. city states) are written unchanged to the
# outputs for every period of every scenario, flagged downscaling_area_flag = 2.

run_static <- function(targets, start_map, mapping_grid, ...) {
  out <- withr::local_tempdir(.local_envir = parent.frame())
  files <- list(lu_output = file.path(out, "lu.csv"), luc_output = file.path(out, "luc.csv"),
                diag_output = file.path(out, "diag.csv"))
  suppressMessages(suppressWarnings(utils::capture.output(
    execute_luc_downscaling(targets = targets, start_map = start_map, mapping_grid = mapping_grid,
                            files = files, base_year = 2020, ...))))
  list(lu = data.table::fread(files$lu_output), luc = data.table::fread(files$luc_output))
}

base_inputs <- function() {
  list(
    targets = data.table::data.table(Ns_reg = "R1", scenario = rep(c("S1", "S2"), each = 2),
                                     times = rep(c("2030", "2040"), 2),
                                     from_class = "low", to_class = "high", value = c(4, 2, 3, 1)),
    start_map = data.table::data.table(ns = c("p1", "p2", "p3", "p3"), Ns_reg = c("R1", "R1", "CITY", "CITY"),
                                       from_class = c("low", "low", "Cropland_other", "Urban"),
                                       value = c(10, 10, 6, 19)),
    mapping_grid = data.table::data.table(ns = c("p1", "p2", "p3"), Ns_reg = c("R1", "R1", "CITY")))
}

test_that("a region without targets appears in the LU output for every scenario and period", {
  x <- base_inputs()
  r <- run_static(x$targets, x$start_map, x$mapping_grid)
  city <- r$lu[Ns_reg == "CITY"]
  expect_setequal(unique(city$scenario), c("S1", "S2"))
  expect_setequal(unique(city$times), c(2020, 2030, 2040))
  expect_true(all(city$downscaling_area_flag == 2L))
  per <- city[, .(area = sum(value)), by = .(scenario, times, lu)]
  expect_true(all(per[lu == "Cropland_other", area] == 6))
  expect_true(all(per[lu == "Urban", area] == 19))
})

test_that("static regions get stay rows in the LUC output and solved regions are unchanged", {
  x <- base_inputs()
  r <- run_static(x$targets, x$start_map, x$mapping_grid)
  city <- r$luc[Ns_reg == "CITY"]
  expect_true(nrow(city) > 0 && all(city$from_class == city$to_class))
  expect_false(any(city$times == 2020))
  expect_true(all(r$luc[Ns_reg == "R1"]$downscaling_area_flag != 2L))
  expect_identical(names(r$luc), c("scenario", "Ns_reg", "times", "from_class", "to_class", "ns", "value",
                                   "rebalancing_applied", "downscaling_area_flag"))
})

test_that("static regions can be switched off", {
  x <- base_inputs()
  r <- run_static(x$targets, x$start_map, x$mapping_grid, write_static_regions = FALSE)
  expect_false("CITY" %in% r$lu$Ns_reg)
})

test_that("a region with targets only in another scenario is static in this one", {
  x <- base_inputs()
  x$targets <- rbind(x$targets, data.table::data.table(Ns_reg = "CITY", scenario = "S1", times = "2030",
                                                       from_class = "Urban", to_class = "Cropland_other", value = 1))
  r <- run_static(x$targets, x$start_map, x$mapping_grid)
  expect_true(all(r$lu[Ns_reg == "CITY" & scenario == "S1"]$downscaling_area_flag != 2L))
  expect_true(all(r$lu[Ns_reg == "CITY" & scenario == "S2"]$downscaling_area_flag == 2L))
})

test_that("the year filter limits static periods but does not make targeted regions static", {
  x <- base_inputs()
  r <- run_static(x$targets, x$start_map, x$mapping_grid, year_filter = 2030)
  expect_setequal(unique(r$lu[Ns_reg == "CITY"]$times), c(2020, 2030))
  expect_false(any(r$lu[Ns_reg == "R1"]$downscaling_area_flag == 2L))
})
