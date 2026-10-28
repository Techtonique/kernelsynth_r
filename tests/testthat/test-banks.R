test_that("chronos bank reproduces kernel-synth.py", {
  b <- kernel_bank_chronos()
  expect_length(b, 33)
  expect_equal(b[[1]]$params$periodicity, 24 / 1024)
  expect_equal(sum(vapply(b, function(k) k$type == "periodic", logical(1))), 21)
})

test_that("extended bank is frequency-aware and grouped", {
  b <- kernel_bank_extended(120, frequency = 12)
  g <- attr(b, "groups")
  expect_setequal(unique(g), c("periodic", "quasi_periodic", "smooth", "trend", "noise"))
  per <- vapply(b[g == "periodic"], function(k) k$params$periodicity, numeric(1)) * 119
  expect_equal(sort(per), seasonal_periods(12))
  set.seed(1); b1 <- kernel_bank_extended(100, 12, randomize = FALSE)
  set.seed(2); b2 <- kernel_bank_extended(100, 12, randomize = FALSE)
  expect_identical(vapply(b1, `[[`, "", "label"), vapply(b2, `[[`, "", "label"))
})

test_that("composite kernels respect max_kernels", {
  set.seed(1)
  for (i in 1:20) {
    k <- sample_composite_kernel(kernel_bank_chronos(), max_kernels = 1)
    expect_true(k$type %in% c("periodic", "smooth", "trend", "noise", "constant"))
  }
  w <- group_weights(periodic = 1, quasi_periodic = 0, smooth = 0, trend = 0, noise = 0)
  k <- sample_composite_kernel(kernel_bank_extended(100, 12), 1, weights = w)
  expect_equal(k$type, "periodic")
})
