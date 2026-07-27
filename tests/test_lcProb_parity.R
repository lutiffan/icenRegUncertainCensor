if (!requireNamespace("testthat", quietly = TRUE)) {
  stop("testthat required")
}

test_that("lcProb NULL matches CRAN icenReg fit", {
  skip_if_not_installed("icenReg")
  skip_if_not(requireNamespace("pkgload", quietly = TRUE))
  skip_if_not(requireNamespace("survival", quietly = TRUE))
  library(survival)

  set.seed(42)
  sim_data <- icenReg::simIC_weib(n = 80, inspections = 4, inspectLength = 1)
  form <- Surv(l, u, type = "interval2") ~ x1 + x2

  ref <- icenReg::ic_sp(form, data = sim_data, bs_samples = 0)
  pkgload::load_all(".", quiet = TRUE)
  dev <- icenRegDev::ic_sp(form, data = sim_data, lcProb = NULL, bs_samples = 0)

  expect_equal(as.numeric(dev$coefficients), as.numeric(ref$coefficients), tolerance = 1e-6)
  expect_equal(dev$llk, ref$llk, tolerance = 1e-8)
  expect_equal(as.numeric(dev$p_hat), as.numeric(ref$p_hat), tolerance = 1e-8)
})

test_that("checkLcProb validation and mixture behavior", {
  skip_if_not(requireNamespace("pkgload", quietly = TRUE))
  skip_if_not(requireNamespace("survival", quietly = TRUE))
  library(survival)
  pkgload::load_all(".", quiet = TRUE)

  y <- matrix(c(1, 1, 2, Inf, 1, 3), ncol = 2, byrow = TRUE)
  expect_error(icenRegDev:::checkLcProb(c(0, 0, 0.5), y), "interval-censored")
  expect_error(icenRegDev:::checkLcProb(c(0, 0.5, 0), y), "right-censored")
  expect_error(icenRegDev:::checkLcProb(c(1.1, 0, 0), y), "lcProb must be in")

  y_lc <- matrix(c(1, 1, 2, Inf, 0, 3), ncol = 2, byrow = TRUE)
  resolved <- icenRegDev:::checkLcProb(c(0.5, NA, NA), y_lc)
  expect_equal(resolved, c(0.5, 0, 1))

  set.seed(1)
  n <- 40
  x1 <- rnorm(n)
  t <- exp(0.3 * x1 + rnorm(n, sd = 0.5))
  yMat <- cbind(t, t)
  yMat[sample(n, 10), 2] <- Inf
  dat <- data.frame(l = yMat[, 1], u = yMat[, 2], x1 = x1)
  lc0 <- rep(0, n)
  lc0[yMat[, 1] == yMat[, 2]] <- 0

  # Explicit lcProb = 0 on exact obs should match omitting lcProb (all implied 0).
  fit0 <- icenRegDev::ic_sp(Surv(l, u, type = "interval2") ~ x1, data = dat, lcProb = lc0, bs_samples = 0)
  fitNull <- icenRegDev::ic_sp(Surv(l, u, type = "interval2") ~ x1, data = dat, bs_samples = 0)
  expect_equal(as.numeric(fit0$coefficients), as.numeric(fitNull$coefficients), tolerance = 1e-5)

  # Interior mixture (pi = 0.5 on exact obs) changes the likelihood vs pi = 0.
  lcHalf <- lc0
  exact_idx <- which(yMat[, 1] == yMat[, 2])
  lcHalf[exact_idx] <- 0.5
  fitHalf <- icenRegDev::ic_sp(Surv(l, u, type = "interval2") ~ x1, data = dat, lcProb = lcHalf, bs_samples = 0)
  expect_false(isTRUE(all.equal(fitHalf$llk, fit0$llk, tolerance = 1e-8)))

  # pi = 1 on one exact obs (left-censored at R) also differs from pi = 0.
  lcOne <- lc0
  lcOne[exact_idx[1]] <- 1
  fitOne <- icenRegDev::ic_sp(Surv(l, u, type = "interval2") ~ x1, data = dat, lcProb = lcOne, bs_samples = 0)
  expect_false(isTRUE(all.equal(fitOne$llk, fit0$llk, tolerance = 1e-8)))
})

test_that("PH interior-pi uses analytic regression derivatives", {
  skip_if_not(requireNamespace("pkgload", quietly = TRUE))
  skip_if_not(requireNamespace("survival", quietly = TRUE))
  library(survival)
  pkgload::load_all(".", quiet = TRUE)

  set.seed(31415)
  sim_data <- icenRegDev::simIC_weib(n = 80, inspections = 4, inspectLength = 1, b1 = 0.3)
  form <- Surv(l, u, type = "interval2") ~ x1 + x2
  lc <- rep(0, nrow(sim_data))
  exact <- sim_data$l == sim_data$u & is.finite(sim_data$u)
  lc[exact] <- 0.5

  fit <- icenRegDev::ic_sp(form, data = sim_data, lcProb = lc, model = "ph", bs_samples = 0)

  expect_true(is.finite(fit$llk))
  expect_true(all(is.finite(fit$coefficients)))
  expect_gt(fit$iterations, 0L)
  expect_equal(as.numeric(fit$coefficients), c(0.6137077738, -0.7102174532), tolerance = 1e-5)
  expect_equal(fit$llk, -84.12603, tolerance = 1e-5)
})

test_that("PO interior-pi still fits with numeric derivatives", {
  skip_if_not(requireNamespace("pkgload", quietly = TRUE))
  skip_if_not(requireNamespace("survival", quietly = TRUE))
  library(survival)
  pkgload::load_all(".", quiet = TRUE)

  set.seed(31415)
  sim_data <- icenRegDev::simIC_weib(n = 80, inspections = 4, inspectLength = 1, b1 = 0.3)
  form <- Surv(l, u, type = "interval2") ~ x1 + x2
  lc <- rep(0, nrow(sim_data))
  exact <- sim_data$l == sim_data$u & is.finite(sim_data$u)
  lc[exact] <- 0.5

  fit <- icenRegDev::ic_sp(form, data = sim_data, lcProb = lc, model = "po", bs_samples = 0)

  expect_true(is.finite(fit$llk))
  expect_true(all(is.finite(fit$coefficients)))
  expect_gt(fit$iterations, 0L)
})
