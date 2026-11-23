# Splits where a label is larger than what stays, above all complete splits where
# nothing stays. Measured against the implicit "stay" option, such a split needs
# very large multipliers. Before the fix, IPF and the multi-target stage did not
# converge, the tempering loop raised T up to 2, and the later stages hit the
# targets with every odds raised to 1/T (coefficients halved at T = 2).

set.seed(42)
n_px  <- 300
ns    <- sprintf("p%03d", seq_len(n_px))
x     <- rnorm(n_px)
areas <- stats::setNames(rgamma(n_px, 2) * 10, ns)
betas <- c(L1 = 0, L2 = 1, L3 = 2)
odds  <- exp(outer(x, betas)); dimnames(odds) <- list(ns, names(betas))

# Exact allocation: a multinomial logit over the labels plus "stay" (odds 1), with
# one multiplier per option so that every option meets its target. It is solved
# relative to the largest option, which keeps the fixed point well conditioned when
# "stay" is tiny.
exact_allocation <- function(targets, stay) {
  stay_target <- sum(areas) - sum(targets)
  w <- cbind(.stay = rep(1, n_px), odds[, names(targets), drop = FALSE])
  tg <- c(.stay = stay_target, targets)
  if (stay == 0) { w <- w[, -1, drop = FALSE]; tg <- tg[-1] }
  ref <- which.max(tg)
  d <- rep(0, length(tg))
  for (i in 1:50000) {
    e <- sweep(w, 2, exp(d), "*")
    q <- e / rowSums(e)
    step <- log(tg / colSums(areas * q)); step[ref] <- 0
    d <- d + step
    if (max(abs(step)) < 1e-13) break
  }
  (areas * q)[, names(targets), drop = FALSE]
}

solve_split <- function(targets, restrictions = NULL, ...) {
  suppressMessages(cascadecore:::optim_allocation(
    targets = targets, priors = odds[, names(targets), drop = FALSE],
    restrictions = restrictions, pixel_avail = ns, type = "lu", areas = areas,
    current_lu_from = "Arable", options = downscale_control(...)))
}

tv <- function(a, b) sum(abs(a - b)) / sum(areas)

test_that("a complete split is the multinomial logit over the labels, untempered", {
  targets <- c(L1 = 0.5, L2 = 0.3, L3 = 0.2) * sum(areas)
  res <- solve_split(targets)
  expect_equal(res$solver_runs$final_tempering, 1)
  expect_lt(tv(res$out_mu[, names(targets)], exact_allocation(targets, stay = 0)), 1e-6)
  expect_equal(unname(colSums(res$out_mu)), unname(targets), tolerance = 1e-8)
  expect_true(all(abs(rowSums(res$out_mu) - areas) < 1e-8 * areas))
})

test_that("an almost complete split keeps its small stay share", {
  for (stay_share in c(1e-9, 1e-6, 1e-4)) {
    targets <- c(L1 = 0.5, L2 = 0.3, L3 = 0.2) * sum(areas) * (1 - stay_share)
    res <- solve_split(targets)
    expect_equal(res$solver_runs$final_tempering, 1)
    expect_lt(tv(res$out_mu[, names(targets)], exact_allocation(targets, stay = 1)), 1e-6)
    expect_equal(sum(areas) - sum(res$out_mu), stay_share * sum(areas),
                 tolerance = 1e-6 * sum(areas))
  }
})

test_that("no tempering at any stay share, also with strong covariates", {
  for (scale in c(1, 3)) {
    odds <<- exp(outer(scale * x, betas)); dimnames(odds) <<- list(ns, names(betas))
    for (stay_share in c(0.45, 0.1, 0.01, 0.005, 0.002, 0.0011)) {
      targets <- c(L1 = 0.5, L2 = 0.3, L3 = 0.2) * sum(areas) * (1 - stay_share)
      res <- solve_split(targets)
      expect_equal(res$solver_runs$final_tempering, 1)
      expect_lt(tv(res$out_mu[, names(targets)], exact_allocation(targets, stay = 1)), 1e-6)
    }
  }
  odds <<- exp(outer(x, betas)); dimnames(odds) <<- list(ns, names(betas))
})

test_that("a split where staying is the largest option uses the general solver", {
  targets <- c(L1 = 0.2, L2 = 0.15, L3 = 0.1) * sum(areas)
  res <- solve_split(targets)
  expect_equal(res$solver_runs$final_tempering, 1)
  expect_null(res$solver_runs$reference_label)
  expect_lt(tv(res$out_mu[, names(targets)], exact_allocation(targets, stay = 1)), 1e-6)
})

test_that("a single label taking the whole area receives every pixel's area", {
  targets <- c(L2 = sum(areas))
  res <- solve_split(targets)
  expect_equal(unname(res$out_mu[, "L2"]), unname(areas), tolerance = 1e-10)
})

test_that("the reference label is one without restrictions", {
  targets <- c(L1 = 0.5, L2 = 0.3, L3 = 0.2) * sum(areas)
  r <- matrix(0, n_px, 3, dimnames = list(ns, names(targets)))
  r[x > 1.5, "L1"] <- 1                     # L1 is the largest target but restricted
  res <- solve_split(targets, restrictions = r)
  expect_equal(res$solver_runs$reference_label, "L2")
  expect_true(all(res$out_mu[x > 1.5, "L1"] < 1e-8))
  expect_equal(unname(colSums(res$out_mu)), unname(targets), tolerance = 1e-6)
})

test_that("with every label restricted it falls back to the general solver", {
  targets <- c(L1 = 0.5, L2 = 0.3, L3 = 0.2) * sum(areas)
  r <- matrix(0, n_px, 3, dimnames = list(ns, names(targets)))
  r[1, ] <- 0.5
  res <- solve_split(targets, restrictions = r)
  expect_null(res$solver_runs$reference_label)
  expect_true(all(res$out_mu[1, ] <= 0.5 * areas[1] + 1e-8))
})

test_that("tempering above 1 is reported, not silent", {
  targets <- c(L1 = 0.5, L2 = 0.3, L3 = 0.2) * sum(areas)
  r <- matrix(0, n_px, 3, dimnames = list(ns, names(targets)))
  r[1, ] <- 0.5                             # forces the general solver on a full split
  expect_message(
    cascadecore:::optim_allocation(
      targets = targets, priors = odds, restrictions = r, pixel_avail = ns, type = "lu",
      areas = areas, current_lu_from = "Arable", options = downscale_control()),
    "tempering")
})
