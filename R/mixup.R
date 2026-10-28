#' TSMixup: convex combinations of time series
#'
#' Data augmentation used for pretraining Chronos (Ansari et al., 2024):
#' for each new series, draw \eqn{k \sim U\{1, ..., k_{max}\}} series,
#' take a random window of common length from each, mean-scale them
#' (divide by the mean absolute value) and combine them with weights
#' drawn from a symmetric Dirichlet(\eqn{\alpha}) distribution.
#'
#' @param series A list of univariate `ts` objects (or numeric vectors), for
#'   example the output of [kernel_synth()].
#' @param n Number of mixed series to generate.
#' @param k_max Maximum number of series per mixture.
#' @param alpha Dirichlet concentration parameter.
#' @param length Length of the mixed series (default: shortest selected series).
#' @param seed Optional random seed.
#'
#' @return A `ts` (if `n = 1`) or a list of `ts`, with frequency taken from
#'   the first selected series and the mixture in `attr(, "mixup")`.
#' @export
#' @examples
#' ys <- kernel_synth(n = 10, length = 120, frequency = 12, seed = 1)
#' mixed <- ts_mixup(ys, n = 3, seed = 2)
#' attr(mixed[[1]], "mixup")
ts_mixup <- function(series, n = 1L, k_max = 3L, alpha = 1.5, length = NULL,
                     seed = NULL) {
  if (stats::is.ts(series) || is.numeric(series)) series <- list(series)
  if (any(vapply(series, function(s) NCOL(s) > 1, logical(1))))
    stop("ts_mixup() expects univariate series")
  if (!is.null(seed)) set.seed(seed)
  lens <- vapply(series, base::length, numeric(1))
  one <- function() {
    k <- sample.int(min(k_max, base::length(series)), 1L)
    idx <- sample.int(base::length(series), k)
    L <- if (is.null(length)) min(lens[idx]) else min(length, min(lens[idx]))
    parts <- lapply(idx, function(i) {
      s <- as.numeric(series[[i]])
      off <- sample.int(base::length(s) - L + 1L, 1L) - 1L
      s <- s[off + seq_len(L)]
      m <- mean(abs(s), na.rm = TRUE)
      if (is.finite(m) && m > 0) s / m else s
    })
    w <- stats::rgamma(k, shape = alpha)
    w <- w / sum(w)
    y <- Reduce(`+`, Map(`*`, parts, w))
    first <- series[[idx[1]]]
    f <- if (stats::is.ts(first)) stats::frequency(first) else 1
    out <- stats::ts(y, start = 1, frequency = f)
    attr(out, "mixup") <- data.frame(index = idx, weight = w)
    attr(out, "kernel") <- unlist(lapply(series[idx], attr, "kernel"))
    out
  }
  res <- lapply(seq_len(n), function(i) one())
  if (n == 1L) res[[1]] else res
}
