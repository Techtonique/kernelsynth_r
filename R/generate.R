# Parts of this file are an R port of KernelSynth (scripts/kernel-synth.py in
# https://github.com/amazon-science/chronos-forecasting), by the Chronos
# authors (Ansari et al., 2024), Copyright Amazon.com, Inc. or its affiliates,
# Apache-2.0. See NOTICE.

#' Generate one KernelSynth series (original algorithm)
#'
#' Faithful port of `generate_time_series()` in kernel-synth.py: draws up to
#' `max_kernels` kernels uniformly (with replacement) from
#' [kernel_bank_chronos()], folds them with random `+`/`*`, and samples one
#' path from the GP prior on an equally spaced grid on \[0, 1\]. Failed
#' factorizations are retried with a new kernel.
#'
#' @param max_kernels Maximum number of base kernels per series.
#' @param length Series length (1024 in the original script).
#' @param frequency,start Passed to [stats::ts()] (the original timestamps
#'   are arbitrary).
#' @param max_tries Maximum number of kernels to try.
#'
#' @return A univariate `ts`, with the kernel expression in `attr(, "kernel")`.
#' @export
#' @examples
#' set.seed(42)
#' y <- generate_time_series(max_kernels = 5, length = 256)
#' attr(y, "kernel")
#' plot(y)
generate_time_series <- function(max_kernels = 5L, length = 1024L,
                                 frequency = 1, start = 1, max_tries = 100L) {
  bank <- kernel_bank_chronos(length)
  x <- .linspace(length)
  for (i in seq_len(max_tries)) {
    kernel <- sample_composite_kernel(bank, max_kernels)
    y <- tryCatch(sample_gp_prior(kernel, x)[, 1], error = function(e) NULL)
    if (!is.null(y) && all(is.finite(y))) {
      out <- stats::ts(y, start = start, frequency = frequency)
      attr(out, "kernel") <- kernel$label
      return(out)
    }
  }
  stop("could not generate a valid series in ", max_tries, " tries")
}

.default_frequencies <- c(1, 4, 12, 52, 7, 24)
.default_freq_probs <- c(0.10, 0.10, 0.25, 0.15, 0.25, 0.15)

#' Generate diverse synthetic time series (enriched KernelSynth)
#'
#' Main generator of the package. For each series:
#' 1. a frequency is drawn (from `frequency`) and a length (from `length`);
#' 2. a kernel bank is built (`"extended"`: frequency-aligned periods and
#'    randomized hyperparameters; `"chronos"`: the original bank);
#' 3. a composite kernel is drawn with [sample_composite_kernel()]
#'    (sums, products and, for the extended bank, changepoints), optionally
#'    input-warped;
#' 4. a path is sampled from the GP prior (several latent paths mixed
#'    linearly if `n_vars > 1`);
#' 5. random augmentations are applied ([augment_series()]).
#'
#' Metadata are attached to each series as attributes: `"kernel"` (kernel
#' expression(s)), `"augmentations"`, and `"bank"`. Use [ks_info()] to
#' retrieve them.
#'
#' Series are generated with per-series seeds drawn from the master `seed`,
#' so results are identical whatever `n_cores`.
#'
#' @param n Number of series.
#' @param length Series length, or a range `c(min, max)` from which lengths
#'   are drawn uniformly.
#' @param max_kernels Maximum number of base kernels per latent process.
#' @param bank `"extended"` (default) or `"chronos"`.
#' @param frequency `NULL` (draw from 1, 4, 12, 52, 7, 24), a single
#'   frequency, or a vector of frequencies to draw from uniformly.
#' @param randomize Randomize kernel hyperparameters (extended bank only).
#' @param augment Apply augmentations.
#' @param control Augmentation probabilities, see [augmentation_control()].
#' @param n_vars Number of variables; `> 1` returns `mts` objects built with a
#'   linear model of coregionalization (random number of latent GPs, random
#'   mixing matrix).
#' @param p_changepoint Probability of a changepoint when combining two
#'   kernels (extended bank only).
#' @param weights Group weights for the extended bank, see [group_weights()].
#' @param start Start of the series; `NULL` uses a default depending on the
#'   frequency.
#' @param seed Optional master random seed.
#' @param n_cores Number of cores (uses [parallel::mclapply()], not on Windows).
#' @param as_mts If `TRUE` and all series are univariate with the same length
#'   and frequency, return them bound in a single `mts`.
#' @param max_tries Maximum number of kernels tried per latent process.
#'
#' @return A `ts` if `n = 1`, otherwise a list of `ts` objects (or one `mts`
#'   if `as_mts = TRUE`).
#' @export
#' @examples
#' # 4 monthly series of length 120
#' ys <- kernel_synth(n = 4, length = 120, frequency = 12, seed = 1)
#' class(ys[[1]])
#' attr(ys[[1]], "kernel")
#'
#' # original KernelSynth, no augmentation
#' y <- kernel_synth(length = 256, bank = "chronos", augment = FALSE, seed = 2)
#'
#' # bivariate hourly series
#' y2 <- kernel_synth(length = 168, frequency = 24, n_vars = 2, seed = 3)
#' plot(y2)
kernel_synth <- function(n = 1L, length = 1024L, max_kernels = 5L,
                         bank = c("extended", "chronos"), frequency = NULL,
                         randomize = TRUE, augment = TRUE,
                         control = augmentation_control(), n_vars = 1L,
                         p_changepoint = 0.1, weights = group_weights(),
                         start = NULL, seed = NULL, n_cores = 1L,
                         as_mts = FALSE, max_tries = 50L) {
  bank <- match.arg(bank)
  stopifnot(n >= 1, all(length >= 2), base::length(length) %in% 1:2,
            max_kernels >= 1, n_vars >= 1)
  if (!inherits(control, "ks_augment_control")) control <- do.call(augmentation_control, control)
  if (!is.null(seed)) set.seed(seed)
  seeds <- sample.int(.Machine$integer.max, n)

  worker <- function(i) {
    set.seed(seeds[i])
    len <- if (base::length(length) == 2L) {
      lo <- min(length); hi <- max(length); lo + sample.int(hi - lo + 1L, 1L) - 1L
    } else length
    freq <- if (is.null(frequency)) {
      if (bank == "chronos") 1 else .pick(.default_frequencies, prob = .default_freq_probs)
    } else .pick(frequency)
    .generate_one(len, max_kernels, bank, freq, randomize, augment, control,
                  n_vars, p_changepoint, weights,
                  if (is.null(start)) .default_start(freq) else start, max_tries)
  }

  res <- if (n_cores > 1L && .Platform$OS.type != "windows") {
    parallel::mclapply(seq_len(n), worker, mc.cores = n_cores)
  } else {
    lapply(seq_len(n), worker)
  }
  errs <- vapply(res, inherits, logical(1), "try-error")
  if (any(errs)) stop(as.character(res[[which(errs)[1]]]))

  if (n == 1L) return(res[[1]])
  if (as_mts) {
    lens <- vapply(res, NROW, numeric(1))
    freqs <- vapply(res, stats::frequency, numeric(1))
    if (n_vars == 1L && length(unique(lens)) == 1L && length(unique(freqs)) == 1L) {
      out <- do.call(cbind, lapply(res, as.numeric))
      colnames(out) <- paste0("Series", seq_len(n))
      out <- stats::ts(out, start = stats::start(res[[1]]), frequency = freqs[1])
      attr(out, "kernel") <- vapply(res, attr, character(1), "kernel")
      return(out)
    }
    warning("series differ in length/frequency or are multivariate; returning a list")
  }
  res
}

.sample_latent <- function(length, max_kernels, bank, frequency, randomize,
                           augment, control, p_changepoint, weights, max_tries) {
  x <- .linspace(length)
  for (i in seq_len(max_tries)) {
    b <- if (bank == "chronos") kernel_bank_chronos(length)
         else kernel_bank_extended(length, frequency, randomize)
    k <- sample_composite_kernel(b, max_kernels,
                                 p_changepoint = if (bank == "chronos") 0 else p_changepoint,
                                 weights = weights)
    warped <- augment && control$time_warp > 0 && stats::runif(1) < control$time_warp
    if (warped) k <- warp_kernel(k, a = .rlogunif(1, 0.3, 3), b = .rlogunif(1, 0.3, 3))
    y <- tryCatch(sample_gp_prior(k, x)[, 1], error = function(e) NULL)
    if (!is.null(y) && all(is.finite(y))) return(list(y = y, kernel = k$label, warped = warped))
  }
  stop("could not generate a valid series in ", max_tries, " tries")
}

.generate_one <- function(length, max_kernels, bank, frequency, randomize,
                          augment, control, n_vars, p_changepoint, weights,
                          start, max_tries) {
  draw <- function() .sample_latent(length, max_kernels, bank, frequency, randomize,
                                    augment, control, p_changepoint, weights, max_tries)
  if (n_vars == 1L) {
    lat <- draw()
    out <- stats::ts(lat$y, start = start, frequency = frequency)
    kernels <- lat$kernel
    warped <- lat$warped
  } else {
    q <- sample.int(n_vars, 1L)
    lats <- lapply(seq_len(q), function(j) draw())
    Z <- vapply(lats, `[[`, numeric(length), "y")
    Z <- matrix(Z, nrow = length)
    W <- matrix(stats::rnorm(n_vars * q), n_vars, q)
    W <- W / sqrt(rowSums(W^2))
    Y <- Z %*% t(W)
    Y <- Y + matrix(stats::rnorm(length * n_vars, 0, 0.05 * .safe_sd(Y)), length, n_vars)
    colnames(Y) <- paste0("V", seq_len(n_vars))
    out <- stats::ts(Y, start = start, frequency = frequency)
    kernels <- vapply(lats, `[[`, character(1), "kernel")
    warped <- any(vapply(lats, `[[`, logical(1), "warped"))
    attr(out, "mixing") <- W
  }
  attr(out, "kernel") <- kernels
  attr(out, "bank") <- bank
  if (augment) {
    out <- augment_series(out, control)
    if (warped) {
      a <- attr(out, "augmentations")
      attr(out, "augmentations") <- if (is.list(a)) lapply(a, function(z) c("time_warp", z))
                                    else c("time_warp", a)
    }
  } else {
    attr(out, "augmentations") <- character()
  }
  out
}

#' Metadata of a generated series
#'
#' @param x A series returned by [kernel_synth()] or [generate_time_series()].
#' @return A list with the kernel expression(s), applied augmentations,
#'   bank, frequency, length and number of variables.
#' @export
#' @examples
#' ks_info(kernel_synth(length = 100, seed = 1))
ks_info <- function(x) {
  list(kernel = attr(x, "kernel"),
       augmentations = attr(x, "augmentations"),
       bank = attr(x, "bank"),
       frequency = stats::frequency(x),
       length = NROW(x),
       n_vars = NCOL(x))
}
