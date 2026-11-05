# Prior matching in solve_gross_transitions(): exact and hierarchical (prefix) rows,
# applied in order so that a later row wins; duplicated margin rows are summed.

ft <- function(lu, v) data.table::data.table(lu = lu, value = v)

test_that("a later exact row overrides an earlier prefix row, and vice versa", {
  from <- ft(c("Crop.a", "Crop.b", "Forest.x"), c(5, 5, 10)); to <- ft(c("Crop.a", "Crop.b", "Forest.x"), c(3, 3, 14))
  pre <- data.table::data.table(from_class = "Crop", to_class = "Forest.x", value = 5)
  ex  <- data.table::data.table(from_class = "Crop.a", to_class = "Forest.x", value = 0.01)
  f <- function(pr) solve_gross_transitions(from, to, prior = pr)[from_class == "Crop.a" & to_class == "Forest.x", value]
  exact_last <- f(rbind(pre, ex)); prefix_last <- f(rbind(ex, pre))
  expect_lt(exact_last, prefix_last)          # the small exact weight wins only when it comes last
  expect_equal(prefix_last, f(pre), tolerance = 1e-9)
})

test_that("without dotted classes, the last of duplicated prior rows wins", {
  from <- ft(c("A", "B"), c(10, 5)); to <- ft(c("A", "B"), c(7, 8))
  p1 <- data.table::data.table(from_class = "A", to_class = "B", value = c(0.001, 5))
  p2 <- data.table::data.table(from_class = "A", to_class = "B", value = 5)
  expect_equal(solve_gross_transitions(from, to, prior = p1)$value, solve_gross_transitions(from, to, prior = p2)$value)
})

test_that("duplicated margin rows are summed", {
  a <- solve_gross_transitions(ft(c("A", "A", "B"), c(4, 6, 5)), ft(c("A", "B"), c(7, 8)))
  b <- solve_gross_transitions(ft(c("A", "B"), c(10, 5)), ft(c("A", "B"), c(7, 8)))
  data.table::setorder(a, from_class, to_class); data.table::setorder(b, from_class, to_class)
  expect_equal(a$value, b$value)
})
