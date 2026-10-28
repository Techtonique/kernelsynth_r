test_that("generate_time_series() is a faithful, working port", {
  set.seed(1)
  y <- generate_time_series(length = 128)
  expect_s3_class(y, "ts")
  expect_length(y, 128)
  expect_true(all(is.finite(y)))
  expect_type(attr(y, "kernel"), "character")
})

test_that("kernel_synth() returns ts objects with metadata", {
  ys <- kernel_synth(n = 5, length = 96, frequency = 12, seed = 1)
  expect_type(ys, "list")
  expect_length(ys, 5)
  for (y in ys) {
    expect_s3_class(y, "ts")
    expect_equal(frequency(y), 12)
    expect_length(y, 96)
    expect_false(is.null(attr(y, "kernel")))
  }
  y <- kernel_synth(length = 50, seed = 2)
  expect_s3_class(y, "ts")
  info <- ks_info(y)
  expect_named(info, c("kernel", "augmentations", "bank", "frequency", "length", "n_vars"))
})

test_that("kernel_synth() is reproducible, also in parallel", {
  a <- kernel_synth(n = 4, length = 64, seed = 10)
  b <- kernel_synth(n = 4, length = 64, seed = 10)
  expect_identical(a, b)
  skip_on_os("windows")
  skip_on_cran()
  c2 <- kernel_synth(n = 4, length = 64, seed = 10, n_cores = 2)
  expect_identical(lapply(a, as.numeric), lapply(c2, as.numeric))
})

test_that("variable lengths, chronos bank, mts and multivariate", {
  ys <- kernel_synth(n = 6, length = c(40, 80), seed = 3)
  lens <- vapply(ys, length, numeric(1))
  expect_true(all(lens >= 40 & lens <= 80))
  y <- kernel_synth(length = 64, bank = "chronos", augment = FALSE, seed = 4)
  expect_equal(attr(y, "bank"), "chronos")
  expect_length(attr(y, "augmentations"), 0)
  m <- kernel_synth(n = 3, length = 48, frequency = 4, as_mts = TRUE, seed = 5)
  expect_true(is.mts(m))
  expect_equal(dim(m), c(48, 3))
  mv <- kernel_synth(length = 72, frequency = 24, n_vars = 3, seed = 6)
  expect_true(is.mts(mv))
  expect_equal(ncol(mv), 3)
})
