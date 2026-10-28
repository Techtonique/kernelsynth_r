x <- seq(0, 1, length.out = 50)

test_that("base kernels match scikit-learn definitions", {
  d <- abs(outer(x, x, "-"))
  expect_equal(kernel_matrix(rbf(0.3), x), exp(-0.5 * (d / 0.3)^2))
  expect_equal(kernel_matrix(exp_sine_squared(0.2, 1.5), x),
               exp(-2 * sin(pi * d / 0.2)^2 / 1.5^2))
  expect_equal(kernel_matrix(rational_quadratic(0.5, 2), x),
               (1 + d^2 / (2 * 2 * 0.25))^(-2))
  expect_equal(kernel_matrix(dot_product(1), x), 1 + outer(x, x))
  expect_equal(kernel_matrix(white_kernel(0.1), x), diag(0.1, 50))
  expect_equal(kernel_matrix(white_kernel(0.1), x, x), matrix(0, 50, 50))
  expect_equal(kernel_matrix(constant_kernel(2), x), matrix(2, 50, 50))
})

test_that("kernel algebra works and keeps labels", {
  k <- rbf(0.1) + exp_sine_squared(0.2)
  p <- k * dot_product(0)
  expect_s3_class(p, "ks_kernel")
  expect_equal(kernel_matrix(p, x),
               (kernel_matrix(rbf(0.1), x) + kernel_matrix(exp_sine_squared(0.2), x)) *
                 kernel_matrix(dot_product(0), x))
  expect_match(p$label, "^\\(RBF")
  expect_equal(kernel_matrix(2 * rbf(1), x), 2 * kernel_matrix(rbf(1), x))
  expect_error(rbf(1) - rbf(2))
})

test_that("all kernels give symmetric PSD matrices", {
  ks <- list(matern(0.2, 0.5), matern(0.2, 1.5), matern(0.2, 2.5), brownian(),
             polynomial(3), cosine_kernel(0.1), spectral_mixture(0.1, 0.3),
             locally_periodic(0.1), changepoint(rbf(0.1), exp_sine_squared(0.1)),
             warp_kernel(exp_sine_squared(0.1), 2, 0.5))
  for (k in ks) {
    K <- kernel_matrix(k, x)
    expect_equal(K, t(K), tolerance = 1e-10)
    expect_gt(min(eigen(K, symmetric = TRUE, only.values = TRUE)$values), -1e-6)
  }
})

test_that("simulate() returns ts objects", {
  y <- simulate(rbf(0.2), length = 60, frequency = 12, seed = 1)
  expect_s3_class(y, "ts")
  expect_equal(frequency(y), 12)
  expect_length(y, 60)
  Y <- simulate(rbf(0.2), nsim = 3, length = 60, seed = 1)
  expect_true(is.mts(Y))
  expect_identical(simulate(rbf(0.2), length = 60, seed = 9),
                   simulate(rbf(0.2), length = 60, seed = 9))
})

test_that("GP sampling handles degenerate kernels", {
  Y <- sample_gp_prior(dot_product(0) * dot_product(0), x, 2)
  expect_equal(dim(Y), c(50, 2))
  expect_true(all(is.finite(Y)))
  Y2 <- sample_gp_prior(rbf(10), x, method = "eigen")
  expect_true(all(is.finite(Y2)))
})
