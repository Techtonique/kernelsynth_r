#' Forecasting tasks for meta-learning
#'
#' Splits each series into a training part (`x`) and a test part (`xx`) of
#' length `h`, in the list format used by the 'Mcomp', 'Tcomp' and
#' 'M4metalearning' packages (`x`, `xx`, `h`, `period`, `n`). Time attributes
#' are preserved so that `xx` starts right after `x`.
#'
#' @param series A `ts`, or a list of `ts` objects (univariate or `mts`).
#' @param h Forecasting horizon. A function of the frequency is also accepted,
#'   e.g. `function(f) 2 * f`.
#' @param min_train Minimum length of the training part.
#'
#' @return A list of tasks; each task is a list with elements `x`, `xx`, `h`,
#'   `period` (frequency), `n` (training length) and `kernel`.
#' @export
#' @examples
#' ys <- kernel_synth(n = 5, length = 120, frequency = 12, seed = 1)
#' tasks <- make_tasks(ys, h = 12)
#' str(tasks[[1]], max.level = 1)
make_tasks <- function(series, h, min_train = 10L) {
  if (stats::is.ts(series)) series <- list(series)
  lapply(series, function(y) {
    if (!stats::is.ts(y)) y <- stats::ts(y)
    f <- stats::frequency(y)
    hh <- if (is.function(h)) h(f) else h
    n <- NROW(y)
    if (n - hh < min_train) stop("series too short for h = ", hh)
    s <- stats::tsp(y)[1]
    sub <- function(i) if (is.matrix(y)) y[i, , drop = FALSE] else as.numeric(y)[i]
    list(x = stats::ts(sub(seq_len(n - hh)), start = s, frequency = f),
         xx = stats::ts(sub((n - hh + 1L):n), start = s + (n - hh) / f, frequency = f),
         h = hh, period = f, n = n - hh,
         kernel = attr(y, "kernel"))
  })
}

#' Meta-features of time series
#'
#' Cheap, dependency-free features commonly used for meta-learning and for
#' checking the diversity of a synthetic corpus: length, frequency, mean,
#' standard deviation, skewness, excess kurtosis, lag-1 autocorrelation (of
#' the series and its first differences), seasonal autocorrelation, trend and
#' seasonal strength (STL-based, Wang et al. 2006), spectral entropy,
#' proportion of zeros and missing values, and whether values are integers.
#'
#' @param x A `ts`, numeric vector, `mts`, or a list of these.
#'
#' @return A one-row `data.frame` for a univariate series, otherwise one row
#'   per series (or per column of an `mts`).
#' @export
#' @examples
#' ys <- kernel_synth(n = 20, length = 120, frequency = 12, seed = 1)
#' mf <- meta_features(ys)
#' summary(mf[, c("trend_strength", "seasonal_strength", "spectral_entropy")])
meta_features <- function(x) {
  if (is.list(x) && !stats::is.ts(x)) {
    return(do.call(rbind, lapply(x, meta_features)))
  }
  if (is.matrix(x)) {
    f <- if (stats::is.ts(x)) stats::frequency(x) else 1
    return(do.call(rbind, lapply(seq_len(ncol(x)), function(j)
      meta_features(stats::ts(x[, j], frequency = f)))))
  }
  .meta_features_one(x)
}

.meta_features_one <- function(y) {
  f <- if (stats::is.ts(y)) stats::frequency(y) else 1
  v <- as.numeric(y)
  n <- length(v)
  prop_na <- mean(is.na(v))
  if (anyNA(v)) {
    ok <- which(!is.na(v))
    v <- if (length(ok) >= 2) stats::approx(ok, v[ok], xout = seq_len(n), rule = 2)$y
         else rep(0, n)
  }
  s <- stats::sd(v)
  const <- !is.finite(s) || s < .Machine$double.eps
  z <- if (const) v * 0 else (v - mean(v)) / s
  acf_at <- function(u, lag) {
    if (length(u) <= lag + 1 || stats::sd(u) < .Machine$double.eps) return(NA_real_)
    stats::acf(u, lag.max = lag, plot = FALSE, na.action = stats::na.pass)$acf[lag + 1]
  }
  # trend / seasonal strength
  trend_strength <- seasonal_strength <- NA_real_
  if (!const && n >= 8) {
    yt <- stats::ts(v, frequency = f)
    if (f > 1 && n > 2 * f + 1 && f == round(f)) {
      dec <- stats::stl(yt, s.window = "periodic", robust = TRUE)$time.series
      Tr <- dec[, "trend"]; S <- dec[, "seasonal"]; R <- dec[, "remainder"]
      seasonal_strength <- max(0, 1 - stats::var(R) / stats::var(S + R))
    } else {
      Tr <- stats::lowess(seq_len(n), v, f = 0.3)$y; R <- v - Tr
    }
    trend_strength <- max(0, 1 - stats::var(R) / stats::var(Tr + R))
  }
  # spectral entropy
  spec_ent <- NA_real_
  if (!const && n >= 8) {
    sp <- stats::spec.pgram(v, taper = 0, detrend = TRUE, plot = FALSE,
                            fast = FALSE)$spec
    p <- sp / sum(sp)
    p <- p[p > 0]
    spec_ent <- -sum(p * log(p)) / log(length(sp))
  }
  data.frame(
    length = n, frequency = f,
    mean = mean(v), sd = if (const) 0 else s,
    skewness = if (const) NA_real_ else mean(z^3),
    kurtosis = if (const) NA_real_ else mean(z^4) - 3,
    acf1 = acf_at(v, 1L), diff_acf1 = acf_at(diff(v), 1L),
    seasonal_acf = if (f > 1 && f == round(f)) acf_at(v, as.integer(f)) else NA_real_,
    trend_strength = trend_strength, seasonal_strength = seasonal_strength,
    spectral_entropy = spec_ent,
    prop_zeros = mean(v == 0), prop_na = prop_na,
    integer_valued = all(abs(v - round(v)) < 1e-8),
    row.names = NULL)
}
