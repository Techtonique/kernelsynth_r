# Parts of this file are an R port of KernelSynth (scripts/kernel-synth.py in
# https://github.com/amazon-science/chronos-forecasting), by the Chronos
# authors (Ansari et al., 2024), Copyright Amazon.com, Inc. or its affiliates,
# Apache-2.0. See NOTICE.

#' Kernel banks
#'
#' `kernel_bank_chronos()` is the exact `KERNEL_BANK` of kernel-synth.py
#' (21 periodic kernels for hourly, daily, weekly, monthly, quarterly and
#' yearly periodicities, 3 dot-product, 3 RBF, 3 rational quadratic,
#' 2 white noise and 1 constant kernel).
#'
#' `kernel_bank_extended()` is a richer, frequency-aware bank, organised in
#' groups (`"periodic"`, `"quasi_periodic"`, `"smooth"`, `"trend"`,
#' `"noise"`) stored in `attr(bank, "groups")`. Periods are expressed in
#' time steps and aligned with the seasonal periods of `frequency` (see
#' [seasonal_periods()]), so that seasonality in the generated `ts` matches
#' its `frequency`. With `randomize = TRUE`, hyperparameters (length scales,
#' alphas, offsets, noise levels, quasi-periods) are redrawn at each call
#' from log-uniform priors.
#'
#' @param length Series length the bank is built for.
#' @param frequency Seasonal frequency of the target series. `NULL` uses the
#'   union of all periods in [seasonal_periods()].
#' @param randomize Draw random hyperparameters (`TRUE`) or use fixed defaults.
#'
#' @return A list of `ks_kernel` objects.
#' @name kernel_banks
#' @examples
#' length(kernel_bank_chronos())
#' bank <- kernel_bank_extended(length = 120, frequency = 12)
#' table(attr(bank, "groups"))
NULL

#' @rdname kernel_banks
#' @export
kernel_bank_chronos <- function(length = 1024L) {
  L <- length
  periods <- c(24, 48, 96,                 # H
               24 * 7, 48 * 7, 96 * 7,     # H
               7, 14,                      # D
               30, 60,                     # D
               365, 365 * 2,               # D
               4, 26, 52,                  # W
               4, 6, 12,                   # M
               4, 4 * 10,                  # Q
               10)                         # Y
  c(lapply(periods, function(p) exp_sine_squared(periodicity = p / L)),
    list(dot_product(sigma_0 = 0), dot_product(sigma_0 = 1), dot_product(sigma_0 = 10),
         rbf(length_scale = 0.1), rbf(length_scale = 1), rbf(length_scale = 10),
         rational_quadratic(alpha = 0.1), rational_quadratic(alpha = 1),
         rational_quadratic(alpha = 10),
         white_kernel(noise_level = 0.1), white_kernel(noise_level = 1),
         constant_kernel()))
}

.seasonal_table <- list(
  "1"  = c(5, 10, 20),                                  # yearly: multi-year cycles
  "4"  = c(4, 8, 40),                                   # quarterly
  "12" = c(3, 4, 6, 12, 24),                            # monthly
  "52" = c(4, 13, 26, 52),                              # weekly
  "7"  = c(7, 14, 30.4375, 60.875, 91.3125, 365.25),    # daily
  "24" = c(12, 24, 48, 96, 168, 336, 672)               # hourly
)

#' Seasonal periods associated with a frequency
#'
#' Candidate periods (in time steps) used by the extended kernel bank:
#' frequency 1 (yearly: multi-year cycles), 4 (quarterly), 12 (monthly),
#' 52 (weekly), 7 (daily with weekly/monthly/yearly cycles), 24 (hourly with
#' daily/weekly cycles). Other frequencies get `f/2`, `f` and `2f`.
#'
#' @param frequency A positive number.
#' @return Numeric vector of periods.
#' @export
#' @examples
#' seasonal_periods(12)
seasonal_periods <- function(frequency) {
  key <- as.character(frequency)
  if (!is.null(.seasonal_table[[key]])) return(.seasonal_table[[key]])
  p <- c(frequency / 2, frequency, 2 * frequency)
  p[p >= 2]
}

#' @rdname kernel_banks
#' @export
kernel_bank_extended <- function(length = 1024L, frequency = NULL,
                                 randomize = TRUE) {
  n1 <- max(length - 1, 1)
  r <- function(lo, hi, default) if (randomize) .rlogunif(1, lo, hi) else default

  periods <- if (is.null(frequency)) sort(unique(unlist(.seasonal_table)))
             else seasonal_periods(frequency)
  periods <- periods[periods >= 2 & periods <= 2 * length]
  if (!length(periods)) periods <- max(2, length / 4)

  periodic <- lapply(periods, function(p)
    exp_sine_squared(periodicity = p / n1, length_scale = r(0.5, 2, 1)))

  quasi <- c(
    lapply(periods, function(p)
      locally_periodic(p / n1, length_scale_periodic = r(0.5, 2, 1),
                       length_scale_decay = r(0.1, 2, 0.5))),
    list(
      cosine_kernel(periodicity = (if (randomize) .pick(periods) else periods[1]) / n1),
      spectral_mixture(periodicity = r(2, max(3, length / 2), max(2, length / 8)) / n1,
                       length_scale = r(0.05, 1, 0.3)),
      spectral_mixture(periodicity = r(2, max(3, length / 2), max(2, length / 32)) / n1,
                       length_scale = r(0.05, 1, 0.1))))

  smooth <- list(
    rbf(r(0.01, 0.1, 0.1)), rbf(r(0.1, 1, 1)), rbf(r(1, 10, 10)),
    rational_quadratic(length_scale = r(0.05, 2, 1), alpha = r(0.05, 0.5, 0.1)),
    rational_quadratic(length_scale = r(0.05, 2, 1), alpha = r(0.5, 2, 1)),
    rational_quadratic(length_scale = r(0.05, 2, 1), alpha = r(2, 20, 10)),
    matern(r(0.01, 1, 0.1), nu = 0.5),
    matern(r(0.01, 1, 0.2), nu = 1.5),
    matern(r(0.01, 1, 0.3), nu = 2.5))

  trend <- list(
    dot_product(sigma_0 = r(0.01, 0.5, 0)), dot_product(sigma_0 = r(0.5, 2, 1)),
    dot_product(sigma_0 = r(2, 10, 10)),
    polynomial(degree = 2, sigma_0 = r(0.1, 2, 1)),
    polynomial(degree = 3, sigma_0 = r(0.1, 2, 1)),
    brownian(sigma_0 = r(0.01, 1, 0)))

  noise <- list(white_kernel(r(0.01, 0.1, 0.1)), white_kernel(r(0.1, 1, 1)),
                constant_kernel(r(0.1, 10, 1)))

  bank <- c(periodic, quasi, smooth, trend, noise)
  attr(bank, "groups") <- rep(c("periodic", "quasi_periodic", "smooth", "trend", "noise"),
                              c(length(periodic), length(quasi), length(smooth),
                                length(trend), length(noise)))
  bank
}

#' Default group weights for the extended bank
#'
#' Probabilities of drawing a base kernel from each group of
#' [kernel_bank_extended()].
#'
#' @param periodic,quasi_periodic,smooth,trend,noise Non-negative weights.
#' @return A named numeric vector summing to one.
#' @export
#' @examples
#' group_weights(periodic = 0.5, noise = 0)
group_weights <- function(periodic = 0.35, quasi_periodic = 0.15, smooth = 0.25,
                          trend = 0.15, noise = 0.10) {
  w <- c(periodic = periodic, quasi_periodic = quasi_periodic, smooth = smooth,
         trend = trend, noise = noise)
  stopifnot(all(w >= 0), sum(w) > 0)
  w / sum(w)
}

#' Combine two kernels with a random operator
#'
#' Port of `random_binary_map()` in kernel-synth.py: returns `a + b` or
#' `a * b` with equal probability. With `p_changepoint > 0`, the kernels are
#' instead combined with [changepoint()] (random location and steepness)
#' with probability `p_changepoint`.
#'
#' @param a,b `ks_kernel` objects.
#' @param p_changepoint Probability of a changepoint combination.
#'
#' @return A `ks_kernel`.
#' @export
#' @examples
#' set.seed(1)
#' random_binary_map(rbf(0.1), exp_sine_squared(0.1))
random_binary_map <- function(a, b, p_changepoint = 0) {
  if (p_changepoint > 0 && stats::runif(1) < p_changepoint) {
    return(changepoint(a, b, location = stats::runif(1, 0.15, 0.85),
                       steepness = .rlogunif(1, 10, 500)))
  }
  if (stats::runif(1) < 0.5) a + b else a * b
}

#' Sample a random composite kernel
#'
#' Draws between 1 and `max_kernels` base kernels (with replacement) from
#' `bank` and folds them with [random_binary_map()]. If `bank` has a
#' `"groups"` attribute (extended bank), a group is drawn first according
#' to `weights`, then a kernel uniformly within the group; otherwise kernels
#' are drawn uniformly, as in kernel-synth.py.
#'
#' @param bank A list of `ks_kernel` objects.
#' @param max_kernels Maximum number of base kernels.
#' @param p_changepoint Probability of combining two kernels with a changepoint.
#' @param weights Group weights, see [group_weights()].
#'
#' @return A `ks_kernel`.
#' @export
#' @examples
#' set.seed(2)
#' sample_composite_kernel(kernel_bank_chronos(), max_kernels = 3)
sample_composite_kernel <- function(bank, max_kernels = 5L, p_changepoint = 0,
                                    weights = group_weights()) {
  k <- sample.int(max_kernels, 1L)
  groups <- attr(bank, "groups")
  if (is.null(groups)) {
    selected <- bank[sample.int(length(bank), k, replace = TRUE)]
  } else {
    w <- weights[names(weights) %in% unique(groups)]
    w <- w[w > 0]
    g <- names(w)[sample.int(length(w), k, replace = TRUE, prob = w)]
    selected <- lapply(g, function(gi) .pick(bank[groups == gi])[[1]])
  }
  Reduce(function(a, b) random_binary_map(a, b, p_changepoint), selected)
}
