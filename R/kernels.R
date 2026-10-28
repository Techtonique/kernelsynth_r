# Parts of this file are an R port of KernelSynth (scripts/kernel-synth.py in
# https://github.com/amazon-science/chronos-forecasting), by the Chronos
# authors (Ansari et al., 2024), Copyright Amazon.com, Inc. or its affiliates,
# Apache-2.0. See NOTICE.

# Gaussian process kernels -------------------------------------------------------
#
# Base kernels follow scikit-learn's parameterization, so that the "chronos"
# bank reproduces kernel-synth.py. A kernel is an S3 object of class
# "ks_kernel" holding a function `fun(x, y = NULL)` that returns the
# covariance matrix. `y = NULL` means "training covariance K(x, x)" (this is
# what makes WhiteKernel add noise only on the diagonal, as in scikit-learn).

new_kernel <- function(fun, label, type = "base", params = list(),
                       children = NULL) {
  structure(list(fun = fun, label = label, type = type, params = params,
                 children = children),
            class = "ks_kernel")
}

.yy <- function(x, y) if (is.null(y)) x else y
.absdist <- function(x, y) abs(outer(x, .yy(x, y), "-"))

#' Gaussian process kernels
#'
#' Constructors for the base covariance kernels used by KernelSynth.
#' Kernels can be combined with `+` and `*` (numbers are promoted to
#' [constant_kernel()]), with [changepoint()] and warped with
#' [warp_kernel()].
#'
#' Definitions (with \eqn{d = |x - x'|}):
#' * `exp_sine_squared`: \eqn{\exp(-2\sin^2(\pi d / p)/\ell^2)}
#' * `rbf`: \eqn{\exp(-d^2 / (2\ell^2))}
#' * `rational_quadratic`: \eqn{(1 + d^2/(2\alpha\ell^2))^{-\alpha}}
#' * `dot_product`: \eqn{\sigma_0^2 + x x'}
#' * `white_kernel`: noise level on the diagonal of \eqn{K(x, x)}
#' * `constant_kernel`: constant \eqn{c}
#' * `matern`: Matern with \eqn{\nu \in \{1/2, 3/2, 5/2\}}
#' * `brownian`: \eqn{\sigma_0^2 + \min(x, x')} (Wiener process, \eqn{x \ge 0})
#' * `polynomial`: \eqn{(\sigma_0^2 + x x')^{degree}}
#' * `cosine_kernel`: \eqn{\cos(2\pi d / p)}
#' * `spectral_mixture`: \eqn{\exp(-d^2/(2\ell^2))\cos(2\pi d/p)}
#' * `locally_periodic`: `exp_sine_squared(p, l_p) * rbf(l_d)`
#'
#' Inputs are assumed to live on \eqn{[0, 1]}, so a period of `P` time steps
#' in a series of length `n` corresponds to `periodicity = P / (n - 1)`.
#'
#' @param periodicity Period, on the input scale.
#' @param length_scale Length scale.
#' @param alpha Scale mixture parameter of the rational quadratic kernel.
#' @param sigma_0 Inhomogeneity (offset) parameter.
#' @param noise_level Variance of the white noise.
#' @param constant_value The constant.
#' @param nu Smoothness of the Matern kernel: 0.5, 1.5 or 2.5.
#' @param degree Polynomial degree.
#' @param length_scale_periodic,length_scale_decay Length scales of the
#'   periodic and decaying (RBF) parts of the locally periodic kernel.
#'
#' @return An object of class `ks_kernel`.
#' @name kernels
#' @examples
#' k <- exp_sine_squared(periodicity = 12 / 199) * rbf(0.5) + white_kernel(0.01)
#' k
#' K <- kernel_matrix(k, seq(0, 1, length.out = 5))
#' dim(K)
NULL

#' @rdname kernels
#' @export
exp_sine_squared <- function(periodicity = 1, length_scale = 1) {
  stopifnot(periodicity > 0, length_scale > 0)
  force(periodicity); force(length_scale)
  new_kernel(
    function(x, y = NULL) {
      exp(-2 * sin(pi * .absdist(x, y) / periodicity)^2 / length_scale^2)
    },
    label = sprintf("ExpSineSquared(length_scale=%s, periodicity=%s)",
                    .fmt(length_scale), .fmt(periodicity)),
    type = "periodic",
    params = list(periodicity = periodicity, length_scale = length_scale))
}

#' @rdname kernels
#' @export
rbf <- function(length_scale = 1) {
  stopifnot(length_scale > 0)
  force(length_scale)
  new_kernel(
    function(x, y = NULL) exp(-0.5 * (.absdist(x, y) / length_scale)^2),
    label = sprintf("RBF(length_scale=%s)", .fmt(length_scale)),
    type = "smooth", params = list(length_scale = length_scale))
}

#' @rdname kernels
#' @export
rational_quadratic <- function(length_scale = 1, alpha = 1) {
  stopifnot(length_scale > 0, alpha > 0)
  force(length_scale); force(alpha)
  new_kernel(
    function(x, y = NULL) {
      (1 + .absdist(x, y)^2 / (2 * alpha * length_scale^2))^(-alpha)
    },
    label = sprintf("RationalQuadratic(alpha=%s, length_scale=%s)",
                    .fmt(alpha), .fmt(length_scale)),
    type = "smooth", params = list(length_scale = length_scale, alpha = alpha))
}

#' @rdname kernels
#' @export
dot_product <- function(sigma_0 = 1) {
  stopifnot(sigma_0 >= 0)
  force(sigma_0)
  new_kernel(
    function(x, y = NULL) sigma_0^2 + outer(x, .yy(x, y)),
    label = sprintf("DotProduct(sigma_0=%s)", .fmt(sigma_0)),
    type = "trend", params = list(sigma_0 = sigma_0))
}

#' @rdname kernels
#' @export
white_kernel <- function(noise_level = 1) {
  stopifnot(noise_level >= 0)
  force(noise_level)
  new_kernel(
    function(x, y = NULL) {
      if (is.null(y)) diag(noise_level, length(x))
      else matrix(0, length(x), length(y))
    },
    label = sprintf("WhiteKernel(noise_level=%s)", .fmt(noise_level)),
    type = "noise", params = list(noise_level = noise_level))
}

#' @rdname kernels
#' @export
constant_kernel <- function(constant_value = 1) {
  force(constant_value)
  new_kernel(
    function(x, y = NULL) matrix(constant_value, length(x), length(.yy(x, y))),
    label = sprintf("%s**2", .fmt(sqrt(abs(constant_value)))),
    type = "constant", params = list(constant_value = constant_value))
}

#' @rdname kernels
#' @export
matern <- function(length_scale = 1, nu = 1.5) {
  stopifnot(length_scale > 0, nu %in% c(0.5, 1.5, 2.5))
  force(length_scale); force(nu)
  new_kernel(
    function(x, y = NULL) {
      r <- .absdist(x, y) / length_scale
      switch(as.character(nu),
             "0.5" = exp(-r),
             "1.5" = (1 + sqrt(3) * r) * exp(-sqrt(3) * r),
             "2.5" = (1 + sqrt(5) * r + 5 * r^2 / 3) * exp(-sqrt(5) * r))
    },
    label = sprintf("Matern(length_scale=%s, nu=%s)", .fmt(length_scale), nu),
    type = "smooth", params = list(length_scale = length_scale, nu = nu))
}

#' @rdname kernels
#' @export
brownian <- function(sigma_0 = 0) {
  force(sigma_0)
  new_kernel(
    function(x, y = NULL) sigma_0^2 + outer(pmax(x, 0), pmax(.yy(x, y), 0), pmin),
    label = sprintf("Brownian(sigma_0=%s)", .fmt(sigma_0)),
    type = "trend", params = list(sigma_0 = sigma_0))
}

#' @rdname kernels
#' @export
polynomial <- function(degree = 2, sigma_0 = 1) {
  stopifnot(degree >= 1)
  force(degree); force(sigma_0)
  new_kernel(
    function(x, y = NULL) (sigma_0^2 + outer(x, .yy(x, y)))^degree,
    label = sprintf("Polynomial(degree=%s, sigma_0=%s)", degree, .fmt(sigma_0)),
    type = "trend", params = list(degree = degree, sigma_0 = sigma_0))
}

#' @rdname kernels
#' @export
cosine_kernel <- function(periodicity = 1) {
  stopifnot(periodicity > 0)
  force(periodicity)
  new_kernel(
    function(x, y = NULL) cos(2 * pi * outer(x, .yy(x, y), "-") / periodicity),
    label = sprintf("Cosine(periodicity=%s)", .fmt(periodicity)),
    type = "periodic", params = list(periodicity = periodicity))
}

#' @rdname kernels
#' @export
spectral_mixture <- function(periodicity = 1, length_scale = 1) {
  stopifnot(periodicity > 0, length_scale > 0)
  force(periodicity); force(length_scale)
  new_kernel(
    function(x, y = NULL) {
      d <- outer(x, .yy(x, y), "-")
      exp(-0.5 * (d / length_scale)^2) * cos(2 * pi * d / periodicity)
    },
    label = sprintf("SpectralMixture(periodicity=%s, length_scale=%s)",
                    .fmt(periodicity), .fmt(length_scale)),
    type = "quasi_periodic",
    params = list(periodicity = periodicity, length_scale = length_scale))
}

#' @rdname kernels
#' @export
locally_periodic <- function(periodicity = 1, length_scale_periodic = 1,
                             length_scale_decay = 0.5) {
  k <- exp_sine_squared(periodicity, length_scale_periodic) * rbf(length_scale_decay)
  k$type <- "quasi_periodic"
  k
}

#' Changepoint kernel
#'
#' Smoothly switches from kernel `k1` to kernel `k2` around `location`:
#' \eqn{K = (1-s(x))(1-s(x'))K_1 + s(x)s(x')K_2}, with
#' \eqn{s(x) = \mathrm{logistic}(steepness (x - location))}.
#'
#' @param k1,k2 `ks_kernel` objects (before and after the change).
#' @param location Changepoint location on the input scale (in \[0, 1\]).
#' @param steepness Steepness of the sigmoid transition.
#'
#' @return An object of class `ks_kernel`.
#' @export
#' @examples
#' k <- changepoint(rbf(0.5), exp_sine_squared(0.05), location = 0.6)
#' y <- simulate(k, length = 200, seed = 1)
changepoint <- function(k1, k2, location = 0.5, steepness = 50) {
  stopifnot(inherits(k1, "ks_kernel"), inherits(k2, "ks_kernel"))
  force(location); force(steepness)
  new_kernel(
    function(x, y = NULL) {
      sx <- stats::plogis(steepness * (x - location))
      sy <- if (is.null(y)) sx else stats::plogis(steepness * (y - location))
      outer(1 - sx, 1 - sy) * k1$fun(x, y) + outer(sx, sy) * k2$fun(x, y)
    },
    label = sprintf("Changepoint(location=%s, steepness=%s)[%s | %s]",
                    .fmt(location), .fmt(steepness), k1$label, k2$label),
    type = "changepoint",
    params = list(location = location, steepness = steepness),
    children = list(k1, k2))
}

#' Input-warped kernel
#'
#' Applies a monotone Kumaraswamy CDF warping
#' \eqn{w(x) = 1 - (1 - x^a)^b} to the inputs before evaluating `kernel`,
#' which turns stationary kernels into non-stationary ones (e.g. chirps,
#' accelerating/decelerating dynamics).
#'
#' @param kernel A `ks_kernel`.
#' @param a,b Positive warping parameters (`a = b = 1` is the identity).
#'
#' @return An object of class `ks_kernel`.
#' @export
#' @examples
#' k <- warp_kernel(exp_sine_squared(0.05), a = 2, b = 1)
warp_kernel <- function(kernel, a = 1, b = 1) {
  stopifnot(inherits(kernel, "ks_kernel"), a > 0, b > 0)
  force(a); force(b)
  w <- function(x) 1 - (1 - pmin(pmax(x, 0), 1)^a)^b
  new_kernel(
    function(x, y = NULL) kernel$fun(w(x), if (is.null(y)) NULL else w(y)),
    label = sprintf("Warp(a=%s, b=%s)[%s]", .fmt(a), .fmt(b), kernel$label),
    type = "warp", params = list(a = a, b = b), children = list(kernel))
}

.kernel_sum <- function(a, b) {
  new_kernel(function(x, y = NULL) a$fun(x, y) + b$fun(x, y),
             label = paste(a$label, "+", b$label),
             type = "sum", children = list(a, b))
}

.kernel_product <- function(a, b) {
  lab <- function(k) if (identical(k$type, "sum")) paste0("(", k$label, ")") else k$label
  new_kernel(function(x, y = NULL) a$fun(x, y) * b$fun(x, y),
             label = paste(lab(a), "*", lab(b)),
             type = "product", children = list(a, b))
}

#' Kernel algebra
#'
#' `+` and `*` combine kernels; numbers are promoted to [constant_kernel()].
#'
#' @param e1,e2 `ks_kernel` objects or numbers.
#' @return An object of class `ks_kernel`.
#' @export
#' @examples
#' rbf(0.1) + 2 * dot_product(1)
Ops.ks_kernel <- function(e1, e2) {
  if (!(.Generic %in% c("+", "*")))
    stop("only '+' and '*' are defined for kernels", call. = FALSE)
  if (missing(e2)) stop("unary operators are not defined for kernels", call. = FALSE)
  if (is.numeric(e1)) e1 <- constant_kernel(e1)
  if (is.numeric(e2)) e2 <- constant_kernel(e2)
  if (.Generic == "+") .kernel_sum(e1, e2) else .kernel_product(e1, e2)
}

#' Evaluate a kernel
#'
#' @param kernel A `ks_kernel`.
#' @param x Numeric input points.
#' @param y Optional second set of points. When `NULL`, the "training"
#'   covariance \eqn{K(x, x)} is returned (white noise is added on the diagonal).
#'
#' @return A `length(x)` by `length(y)` covariance matrix.
#' @export
kernel_matrix <- function(kernel, x, y = NULL) {
  stopifnot(inherits(kernel, "ks_kernel"))
  x <- as.numeric(x)
  if (!is.null(y)) y <- as.numeric(y)
  K <- kernel$fun(x, y)
  if (is.null(dim(K))) K <- matrix(K, length(x), length(.yy(x, y)))
  K
}

#' @export
format.ks_kernel <- function(x, ...) x$label

#' @export
as.character.ks_kernel <- function(x, ...) x$label

#' @export
print.ks_kernel <- function(x, ...) {
  cat("<ks_kernel>", x$label, "\n")
  invisible(x)
}

#' Simulate from a Gaussian process prior
#'
#' Draws sample paths from a zero-mean GP prior with covariance `object`
#' and returns them as a `ts` (one path) or `mts` (several paths).
#'
#' @param object A `ks_kernel`.
#' @param nsim Number of sample paths.
#' @param seed Optional random seed.
#' @param length Series length (ignored if `x` is given).
#' @param frequency,start Passed to [stats::ts()].
#' @param x Optional input points (default: equally spaced on \[0, 1\]).
#' @param method Sampling method, see [sample_gp_prior()].
#' @param ... Unused.
#'
#' @return A `ts` object.
#' @importFrom stats simulate
#' @export
#' @examples
#' y <- simulate(exp_sine_squared(12 / 143) + rbf(0.3), length = 144,
#'               frequency = 12, seed = 123)
#' plot(y)
simulate.ks_kernel <- function(object, nsim = 1, seed = NULL, length = 1024L,
                               frequency = 1, start = 1, x = NULL,
                               method = c("cholesky", "eigen"), ...) {
  if (!is.null(seed)) set.seed(seed)
  if (is.null(x)) x <- .linspace(length)
  Y <- sample_gp_prior(object, x, n_samples = nsim, method = match.arg(method))
  out <- stats::ts(if (nsim == 1) Y[, 1] else Y, start = start, frequency = frequency)
  attr(out, "kernel") <- object$label
  out
}
