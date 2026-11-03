# Fractional restrictions: a restricted share r caps the allocation of the source
# class to a target at (1 - r) * source area in every pixel, in every solver stage.
# Pixel p1 is strongly preferred, so stages after the first see p1 already at or
# near its cap and must not add more there.

stage_case <- function(...) {
  px <- c("p1", "p2")
  opt <- downscale_control(..., max_diff = 1e-8)
  res <- suppressMessages(suppressWarnings(cascadecore:::optim_allocation(
    targets = c(high = 9.9),
    priors = matrix(c(50, 1e-4), ncol = 1, dimnames = list(px, "high")),
    restrictions = matrix(0.5, nrow = 2, ncol = 1, dimnames = list(px, "high")),
    pixel_avail = px, type = "lu", areas = setNames(c(10, 10), px),
    current_lu_from = "low", options = opt)))
  res$out_mu[, "high"]
}

test_that("mop-up respects fractional caps after a partial first stage", {
  a <- stage_case(use_ipf = TRUE, use_multi_nloptr = FALSE, use_single_nloptr = FALSE,
                  use_mopup = TRUE, min_eval = 1, max_eval = 1)
  expect_true(all(a <= 5 + 1e-6))
  expect_equal(sum(a), 9.9, tolerance = 1e-6)
})

test_that("the single-target stage respects fractional caps", {
  a <- stage_case(use_ipf = TRUE, use_multi_nloptr = TRUE, use_single_nloptr = TRUE,
                  use_mopup = FALSE, min_eval = 1, max_eval = 2)
  expect_true(all(a <= 5 + 1e-6))
})

test_that("single-target followed by mop-up respects fractional caps", {
  a <- stage_case(use_ipf = TRUE, use_multi_nloptr = TRUE, use_single_nloptr = TRUE,
                  use_mopup = TRUE, min_eval = 1, max_eval = 2)
  expect_true(all(a <= 5 + 1e-6))
  expect_equal(sum(a), 9.9, tolerance = 1e-6)
})

test_that("residual_restrictions leaves binary restrictions unchanged", {
  px <- c("p1", "p2", "p3")
  r <- matrix(c(1, 0, 0, 0, 1, 0), nrow = 3, dimnames = list(px, c("a", "b")))
  areas <- setNames(c(10, 10, 10), px)
  allocated <- matrix(c(0, 3, 1, 2, 0, 4), nrow = 3, dimnames = list(px, c("a", "b")))
  active <- areas - rowSums(allocated)
  r_eff <- cascadecore:::residual_restrictions(r, areas, allocated, active, px)
  expect_equal(r_eff, r)
})

test_that("residual_restrictions converts a fractional cap into remaining capacity", {
  px <- c("p1", "p2")
  r <- matrix(0.5, nrow = 2, ncol = 1, dimnames = list(px, "a"))
  areas <- setNames(c(10, 10), px)
  allocated <- matrix(c(4, 5), ncol = 1, dimnames = list(px, "a"))   # cap is 5 in each pixel
  active <- areas - rowSums(allocated)                                # 6 and 5 left
  r_eff <- cascadecore:::residual_restrictions(r, areas, allocated, active, px)
  expect_equal(as.vector(active * (1 - r_eff)), c(1, 0))              # 1 km2 left in p1, none in p2
})
