test_that("augmentations preserve ts attributes", {
  set.seed(1)
  y <- simulate(rbf(0.2), length = 100, frequency = 12, start = c(2000, 1))
  fs <- list(aug_integrate, aug_trend, aug_level_shift, aug_positive,
             aug_intermittent, aug_counts, aug_heteroscedastic, aug_noise,
             aug_outliers, aug_scale, aug_missing)
  for (f in fs) {
    z <- f(y)
    expect_s3_class(z, "ts")
    expect_equal(tsp(z), tsp(y))
    expect_equal(attr(z, "kernel"), attr(y, "kernel"))
  }
  expect_true(all(aug_positive(y) > 0))
  z <- aug_counts(y); expect_true(all(z == round(z)))
  expect_gt(mean(aug_intermittent(y, p_zero = 0.7) == 0), 0.5)
  expect_true(anyNA(aug_missing(y, 0.1)))
})

test_that("augment_series() records the pipeline", {
  set.seed(2)
  y <- simulate(rbf(0.2), length = 100)
  ctrl <- augmentation_control(trend = 1, outliers = 1, positive = 0,
                               counts = 0, intermittent = 0, integrate = 0)
  z <- augment_series(y, ctrl)
  expect_true(all(c("trend", "outliers") %in% attr(z, "augmentations")))
  expect_error(augmentation_control(trend = 2))
})

test_that("ts_mixup() mixes series", {
  ys <- kernel_synth(n = 6, length = c(60, 90), frequency = 12, seed = 1)
  m <- ts_mixup(ys, n = 4, k_max = 3, seed = 2)
  expect_length(m, 4)
  for (z in m) {
    expect_s3_class(z, "ts")
    expect_equal(sum(attr(z, "mixup")$weight), 1)
  }
})

test_that("make_tasks() splits consistently", {
  ys <- kernel_synth(n = 3, length = 60, frequency = 12, seed = 1)
  tk <- make_tasks(ys, h = function(f) f)
  expect_length(tk, 3)
  t1 <- tk[[1]]
  expect_equal(length(t1$x), 48)
  expect_equal(length(t1$xx), 12)
  expect_equal(tsp(t1$xx)[1], tsp(t1$x)[2] + 1 / 12)
  expect_equal(c(as.numeric(t1$x), as.numeric(t1$xx)), as.numeric(ys[[1]]))
})

test_that("meta_features() returns one row per series", {
  ys <- kernel_synth(n = 5, length = 80, frequency = 12, seed = 1)
  mf <- meta_features(ys)
  expect_s3_class(mf, "data.frame")
  expect_equal(nrow(mf), 5)
  expect_true(all(c("trend_strength", "seasonal_strength", "spectral_entropy") %in% names(mf)))
  expect_equal(nrow(meta_features(ts(rep(1, 30)))), 1)
})

test_that("R6 interface works", {
  gen <- KernelSynth$new(length = 60, frequency = 4, seed = 1)
  a <- gen$generate(2)
  b <- gen$generate(2)
  expect_false(identical(a, b))
  gen2 <- KernelSynth$new(length = 60, frequency = 4, seed = 1)
  expect_identical(a, gen2$generate(2))
  expect_s3_class(gen$sample_kernel(), "ks_kernel")
  expect_length(gen$tasks(3, h = 4), 3)
  expect_output(print(gen), "KernelSynth")
})
