# Allowed-transition matrix, transport feasibility and fit diagnostics.
# Without a matrix, solve_gross_transitions() and execute_luc_downscaling() behave as before.

ft <- function(lu, v) data.table::data.table(lu = lu, value = v)

test_that("without a matrix the flows are unchanged", {
  from <- ft(c("A", "B", "C"), c(10, 5, 5)); to <- ft(c("A", "B", "C"), c(6, 8, 6))
  ref <- suppressWarnings(solve_gross_transitions(from, to))
  none <- data.table::data.table(from_class = "A", to_class = "B", allowed = 1)
  same <- suppressWarnings(solve_gross_transitions(from, to, allowed = none))
  data.table::setorder(ref, from_class, to_class); data.table::setorder(same, from_class, to_class)
  expect_equal(ref$value, same$value)
  expect_true(attr(ref, "ipf_diagnostics")$converged)
})

test_that("disallowed pairs get no flow and the margins are still met", {
  from <- ft(c("A", "B", "C"), c(10, 5, 5)); to <- ft(c("A", "B", "C"), c(6, 8, 6))
  z <- data.table::data.table(from_class = c("A", "C"), to_class = c("C", "A"), allowed = 0)
  out <- solve_gross_transitions(from, to, allowed = z)
  expect_false(any(out$from_class == "A" & out$to_class == "C"))
  expect_false(any(out$from_class == "C" & out$to_class == "A"))
  expect_equal(out[, sum(value), by = from_class][order(from_class)]$V1, c(10, 5, 5), tolerance = 1e-6)
  expect_equal(out[, sum(value), by = to_class][order(to_class)]$V1, c(6, 8, 6), tolerance = 1e-6)
})

test_that("infeasible margins stop with a matrix, or warn on request", {
  # B gains 5, but only A loses area and A -> B is disallowed
  from <- ft(c("A", "B"), c(10, 0)); to <- ft(c("A", "B"), c(5, 5))
  z <- data.table::data.table(from_class = "A", to_class = "B", allowed = 0)
  expect_error(solve_gross_transitions(from, to, allowed = z), "cannot be met")
  out <- suppressWarnings(solve_gross_transitions(from, to, allowed = z, on_infeasible = "warn"))
  expect_warning(solve_gross_transitions(from, to, allowed = z, on_infeasible = "warn"), "cannot be met")
  expect_false(attr(out, "ipf_diagnostics")$feasible)
})

test_that("exempt balancing classes keep all transitions", {
  from <- ft(c("A", "NODATA"), c(5, 5)); to <- ft(c("A", "B"), c(5, 5))
  z <- data.table::data.table(from_class = "NODATA", to_class = "B", allowed = 0)
  out <- solve_gross_transitions(from, to, allowed = z)
  expect_true(out[from_class == "NODATA" & to_class == "B", sum(value)] > 0)
})

test_that("fienberg reports the remaining margin error and convergence", {
  m <- matrix(c(0.9, 0.1, 0.1, 0.9), 2); m <- m / rowSums(m)
  ok <- fienberg(m, c(10, 10), c(4, 16))
  expect_true(ok$converged); expect_lt(ok$margin_error, 1e-6)
  short <- fienberg(m, c(10, 10), c(4, 16), max.iter = 1, min.iter = 1)
  expect_false(short$converged); expect_gt(short$margin_error, 1e-6)
})

test_that("fits are recorded, and missed margins warn", {
  ipf_diagnostics(reset = TRUE)
  invisible(suppressWarnings(solve_gross_transitions(ft(c("A", "B"), c(10, 5)), ft(c("A", "B"), c(7, 8)))))
  d <- ipf_diagnostics()
  expect_equal(nrow(d), 1L)
  expect_true(all(c("feasible", "iter", "stop_reason", "margin_error", "converged") %in% names(d)))
  # a zero prior cell that leaves the margins infeasible: recorded, warned, not converged
  prior <- data.table::data.table(from_class = "A", to_class = "B", value = 0)
  expect_warning(solve_gross_transitions(ft(c("A", "B"), c(10, 0)), ft(c("A", "B"), c(5, 5)), prior = prior),
                 "cannot be met|misses its margins")
  d2 <- ipf_diagnostics(reset = TRUE)
  expect_false(d2$converged[2])
  expect_equal(nrow(ipf_diagnostics()), 0L)
})

test_that("check_transport_feasibility finds the max flow", {
  s <- matrix(c(1, 0, 1, 1), 2)  # A -> A, A -> B, B -> B allowed; B -> A not
  expect_true(check_transport_feasibility(s, c(5, 5), c(5, 5))$feasible)
  s2 <- matrix(c(1, 0, 0, 1), 2)  # diagonal only
  r <- check_transport_feasibility(s2, c(8, 2), c(5, 5))
  expect_false(r$feasible); expect_equal(r$shortfall, 3, tolerance = 1e-9)
})

test_that("normalize_allowed_transitions validates its input", {
  expect_null(normalize_allowed_transitions(NULL))
  expect_error(normalize_allowed_transitions(data.frame(from_class = "A", to_class = "B")), "need columns")
  expect_error(normalize_allowed_transitions(data.frame(from_class = "A", to_class = "B", allowed = 2)), "0/1")
  expect_error(normalize_allowed_transitions(data.frame(from_class = c("A", "A"), to_class = "B", allowed = 0)), "duplicated")
  expect_message(d <- normalize_allowed_transitions(data.frame(from_class = "A", to_class = "A", allowed = 0)), "always allowed")
  expect_equal(nrow(d), 0L)
})

# execute_luc_downscaling --------------------------------------------------

run_ds <- function(targets, ...) {
  out <- withr::local_tempdir(.local_envir = parent.frame())
  files <- list(lu_output = file.path(out, "lu.csv"), luc_output = file.path(out, "luc.csv"),
                diag_output = file.path(out, "diag.csv"))
  start_map <- data.table::data.table(ns = c("p1", "p2"), Ns_reg = "R1", from_class = "low", value = c(10, 10))
  mapping_grid <- data.table::data.table(ns = c("p1", "p2"), Ns_reg = "R1")
  suppressMessages(suppressWarnings(utils::capture.output(
    execute_luc_downscaling(targets = targets, start_map = start_map, mapping_grid = mapping_grid,
                            files = files, base_year = 2020, ...))))
  data.table::fread(files$luc_output)
}
tg <- function() data.table::data.table(Ns_reg = "R1", scenario = "S1", times = "2030",
                                        from_class = c("low", "low"), to_class = c("high", "mid"), value = c(4, 2))

test_that("execute_luc_downscaling is unchanged without a matrix", {
  luc <- run_ds(tg())
  expect_true(luc[from_class == "low" & to_class == "high", sum(value)] > 0)
})

test_that("targets on disallowed transitions stop, or are dropped on request", {
  z <- data.table::data.table(from_class = "low", to_class = "mid", allowed = 0)
  expect_error(run_ds(tg(), allowed_transitions = z), "disallowed transitions")
  luc <- run_ds(tg(), allowed_transitions = z, on_disallowed = "drop")
  expect_equal(luc[from_class == "low" & to_class == "mid", sum(value)], 0)
  expect_equal(luc[from_class == "low" & to_class == "high", sum(value)], 4, tolerance = 1e-6)
})
