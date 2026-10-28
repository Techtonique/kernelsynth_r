#' Augmentation control
#'
#' Probabilities with which each augmentation is applied by
#' [augment_series()] (and by [kernel_synth()] when `augment = TRUE`).
#' `positive`, `intermittent` and `counts` change the support of the series
#' and are mutually exclusive (at most one is applied, checked in the order
#' counts, intermittent, positive).
#'
#' @param time_warp Input warping of the GP kernel (applied at generation time).
#' @param integrate Cumulative sum (unit-root-like behaviour), see [aug_integrate()].
#' @param trend Deterministic trend, see [aug_trend()].
#' @param level_shift Structural breaks in level, see [aug_level_shift()].
#' @param positive Positive support (multiplicative effects), see [aug_positive()].
#' @param intermittent Intermittent counts, see [aug_intermittent()].
#' @param counts Poisson counts, see [aug_counts()].
#' @param heteroscedastic Time-varying noise variance, see [aug_heteroscedastic()].
#' @param noise Additive white noise, see [aug_noise()].
#' @param outliers Additive spikes, see [aug_outliers()].
#' @param scale Random affine rescaling, see [aug_scale()].
#' @param missing Missing values, see [aug_missing()].
#'
#' @return A list of class `ks_augment_control`.
#' @export
#' @examples
#' augmentation_control(missing = 0.1, counts = 0)
augmentation_control <- function(time_warp = 0.10, integrate = 0.10, trend = 0.20,
                                 level_shift = 0.10, positive = 0.15,
                                 intermittent = 0.05, counts = 0.05,
                                 heteroscedastic = 0.10, noise = 0.20,
                                 outliers = 0.10, scale = 0.50, missing = 0) {
  ctrl <- list(time_warp = time_warp, integrate = integrate, trend = trend,
               level_shift = level_shift, positive = positive,
               intermittent = intermittent, counts = counts,
               heteroscedastic = heteroscedastic, noise = noise,
               outliers = outliers, scale = scale, missing = missing)
  ok <- vapply(ctrl, function(p) is.numeric(p) && length(p) == 1L && p >= 0 && p <= 1,
               logical(1))
  if (!all(ok)) stop("all probabilities must be single numbers in [0, 1]: ",
                     paste(names(ctrl)[!ok], collapse = ", "))
  structure(ctrl, class = "ks_augment_control")
}

#' @export
print.ks_augment_control <- function(x, ...) {
  cat("<ks_augment_control>\n")
  print(unlist(unclass(x)))
  invisible(x)
}

#' Augment a time series
#'
#' Applies a random pipeline of augmentations, each with the probability
#' given in `control`, in the order: integrate, trend, level shift,
#' support transform (counts / intermittent / positive), heteroscedastic
#' noise, white noise, outliers, scaling, missing values. Applied steps are
#' recorded in `attr(y, "augmentations")`. Multivariate series are augmented
#' column by column.
#'
#' @param y A numeric vector, `ts` or `mts`.
#' @param control See [augmentation_control()].
#'
#' @return The augmented series, of the same kind as `y`.
#' @export
#' @examples
#' set.seed(3)
#' y <- simulate(rbf(0.2), length = 200)
#' z <- augment_series(y, augmentation_control(trend = 1, outliers = 1))
#' attr(z, "augmentations")
augment_series <- function(y, control = augmentation_control()) {
  if (is.matrix(y)) {
    cols <- lapply(seq_len(ncol(y)), function(j) {
      yj <- if (stats::is.ts(y)) y[, j] else stats::ts(y[, j])
      augment_series(yj, control)
    })
    out <- do.call(cbind, lapply(cols, as.numeric))
    colnames(out) <- colnames(y)
    if (stats::is.ts(y)) out <- stats::ts(out, start = stats::start(y),
                                          frequency = stats::frequency(y))
    out <- .copy_attrs(out, y)
    attr(out, "augmentations") <- lapply(cols, attr, "augmentations")
    return(out)
  }
  applied <- character()
  hit <- function(name) control[[name]] > 0 && stats::runif(1) < control[[name]]
  step <- function(name, f) {
    if (hit(name)) { y <<- f(y); applied <<- c(applied, name); TRUE } else FALSE
  }
  step("integrate", aug_integrate)
  step("trend", aug_trend)
  step("level_shift", aug_level_shift)
  support <- if (step("counts", aug_counts)) "counts"
             else if (step("intermittent", aug_intermittent)) "intermittent"
             else if (step("positive", aug_positive)) "positive"
             else "real"
  if (support == "real") {
    step("heteroscedastic", aug_heteroscedastic)
    step("noise", aug_noise)
    step("outliers", aug_outliers)
    step("scale", aug_scale)
  } else if (support == "positive") {
    step("outliers", function(v) aug_outliers(v, positive = TRUE))
    step("scale", function(v) aug_scale(v, location = 0))
  }
  step("missing", aug_missing)
  attr(y, "augmentations") <- c(attr(y, "augmentations"), applied)
  y
}

#' Individual augmentations
#'
#' Building blocks used by [augment_series()]. Each takes a numeric vector
#' or univariate `ts` and returns an object of the same kind (time
#' attributes preserved). Unspecified parameters are drawn at random.
#'
#' * `aug_integrate()`: cumulative sum, rescaled to the input's mean/sd.
#' * `aug_trend()`: adds a linear, quadratic, exponential, logistic or
#'   piecewise-linear trend of random strength.
#' * `aug_level_shift()`: adds 1-3 permanent level shifts.
#' * `aug_positive()`: maps to positive values (`exp` = multiplicative
#'   dynamics, `softplus`, or `shift`).
#' * `aug_intermittent()`: intermittent (sparse) Poisson counts.
#' * `aug_counts()`: Poisson counts with a softplus link.
#' * `aug_heteroscedastic()`: noise with a smoothly time-varying variance.
#' * `aug_noise()`: additive Gaussian noise at a random signal-to-noise ratio.
#' * `aug_outliers()`: additive spikes.
#' * `aug_scale()`: affine map `location + scale * y`.
#' * `aug_missing()`: scattered and block-missing values (`NA`).
#'
#' @param y Numeric vector or univariate `ts`.
#' @param type Trend type: one of `"linear"`, `"quadratic"`, `"exponential"`,
#'   `"logistic"`, `"piecewise"`.
#' @param strength Trend amplitude in units of the series' standard deviation.
#' @param n_shifts Number of level shifts.
#' @param method Positivity transform: `"exp"`, `"softplus"` or `"shift"`.
#' @param p_zero Target proportion of zeros.
#' @param rate Counts scale (`aug_counts`), outlier rate (`aug_outliers`) or
#'   missing rate (`aug_missing`).
#' @param snr Signal-to-noise (variance) ratio.
#' @param positive Only upward spikes (for positive series).
#' @param location,scale Affine parameters.
#'
#' @return Augmented series.
#' @name augmentations
#' @examples
#' set.seed(4)
#' y <- simulate(exp_sine_squared(12 / 119), length = 120, frequency = 12)
#' plot(aug_positive(aug_trend(y, "logistic"), "exp"))
NULL

#' @rdname augmentations
#' @export
aug_integrate <- function(y) {
  .with_ts(y, function(v) .rescale_like(cumsum(.std(v)), v))
}

#' @rdname augmentations
#' @export
aug_trend <- function(y, type = NULL, strength = NULL) {
  types <- c("linear", "quadratic", "exponential", "logistic", "piecewise")
  type <- if (is.null(type)) .pick(types) else match.arg(type, types)
  if (is.null(strength)) strength <- stats::runif(1, 0.5, 4) * .pick(c(-1, 1))
  .with_ts(y, function(v) {
    tt <- .linspace(length(v))
    tr <- switch(type,
      linear = tt,
      quadratic = (tt - stats::runif(1, -0.2, 1.2))^2,
      exponential = { r <- stats::runif(1, 1, 6); (exp(r * tt) - 1) / (exp(r) - 1) },
      logistic = stats::plogis(stats::runif(1, 5, 40) * (tt - stats::runif(1, 0.2, 0.8))),
      piecewise = {
        cps <- sort(stats::runif(sample.int(3L, 1L), 0.1, 0.9))
        slopes <- stats::rnorm(length(cps) + 1L)
        slopes[1] * tt + colSums(slopes[-1] * t(outer(tt, cps, function(a, b) pmax(a - b, 0))))
      })
    rng <- diff(range(tr))
    tr <- if (rng > 0) (tr - min(tr)) / rng else tr * 0
    v + strength * .safe_sd(v) * tr
  })
}

#' @rdname augmentations
#' @export
aug_level_shift <- function(y, n_shifts = NULL) {
  .with_ts(y, function(v) {
    n <- length(v)
    if (n < 10) return(v)
    k <- if (is.null(n_shifts)) sample.int(3L, 1L) else n_shifts
    pos <- sample(seq.int(ceiling(0.1 * n), floor(0.9 * n)), k)
    shift <- numeric(n)
    for (p in pos) shift[p:n] <- shift[p:n] + stats::rnorm(1, 0, 1.5) * .safe_sd(v)
    v + shift
  })
}

#' @rdname augmentations
#' @export
aug_positive <- function(y, method = NULL) {
  method <- if (is.null(method)) .pick(c("exp", "exp", "softplus", "shift"))
            else match.arg(method, c("exp", "softplus", "shift"))
  .with_ts(y, function(v) {
    z <- .std(v)
    switch(method,
      exp = exp(stats::runif(1, 0.1, 1) * z),
      softplus = .softplus(stats::runif(1, 0.5, 3) * z),
      shift = z - min(z, na.rm = TRUE) + stats::runif(1, 0.01, 1))
  })
}

#' @rdname augmentations
#' @export
aug_intermittent <- function(y, p_zero = NULL) {
  if (is.null(p_zero)) p_zero <- stats::runif(1, 0.3, 0.9)
  .with_ts(y, function(v) {
    z <- .std(v)
    thr <- stats::quantile(z, p_zero, na.rm = TRUE, names = FALSE)
    lambda <- pmax(z - thr, 0) * .rlogunif(1, 1, 20)
    out <- numeric(length(z))
    pos <- which(lambda > 0)
    out[pos] <- 1 + stats::rpois(length(pos), lambda[pos])
    out
  })
}

#' @rdname augmentations
#' @export
aug_counts <- function(y, rate = NULL) {
  if (is.null(rate)) rate <- .rlogunif(1, 0.5, 200)
  .with_ts(y, function(v) {
    lambda <- .softplus(stats::runif(1, 0.5, 2) * .std(v)) * rate
    as.numeric(stats::rpois(length(lambda), lambda))
  })
}

#' @rdname augmentations
#' @export
aug_heteroscedastic <- function(y) {
  .with_ts(y, function(v) {
    tt <- .linspace(length(v))
    g <- stats::runif(1, 0.3, 1.5) * sin(2 * pi * stats::runif(1, 0.5, 4) * tt +
                                           stats::runif(1, 0, 2 * pi)) +
         stats::rnorm(1, 0, 1) * tt
    v + stats::runif(1, 0.1, 0.5) * .safe_sd(v) * exp(g) * stats::rnorm(length(v))
  })
}

#' @rdname augmentations
#' @export
aug_noise <- function(y, snr = NULL) {
  if (is.null(snr)) snr <- .rlogunif(1, 1, 100)
  .with_ts(y, function(v) v + stats::rnorm(length(v), 0, .safe_sd(v) / sqrt(snr)))
}

#' @rdname augmentations
#' @export
aug_outliers <- function(y, rate = NULL, positive = FALSE) {
  if (is.null(rate)) rate <- stats::runif(1, 0.002, 0.02)
  .with_ts(y, function(v) {
    n <- length(v)
    k <- max(1L, stats::rbinom(1, n, rate))
    idx <- sample.int(n, min(k, n))
    sgn <- if (positive) 1 else .pick(c(-1, 1), length(idx), replace = TRUE)
    v[idx] <- v[idx] + sgn * stats::runif(length(idx), 3, 8) * .safe_sd(v)
    v
  })
}

#' @rdname augmentations
#' @export
aug_scale <- function(y, location = NULL, scale = NULL) {
  if (is.null(scale)) scale <- 10^stats::runif(1, -1, 3)
  if (is.null(location)) location <- stats::rnorm(1, 0, 3) * scale
  .with_ts(y, function(v) location + scale * v)
}

#' @rdname augmentations
#' @export
aug_missing <- function(y, rate = NULL) {
  if (is.null(rate)) rate <- stats::runif(1, 0.01, 0.1)
  .with_ts(y, function(v) {
    n <- length(v)
    v[sample.int(n, floor(rate * n / 2))] <- NA
    blk <- floor(rate * n / 2)
    if (blk > 0) { s <- sample.int(max(1, n - blk), 1); v[s:(s + blk - 1)] <- NA }
    v
  })
}
