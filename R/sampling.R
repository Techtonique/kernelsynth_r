# Parts of this file are an R port of KernelSynth (scripts/kernel-synth.py in
# https://github.com/amazon-science/chronos-forecasting), by the Chronos
# authors (Ansari et al., 2024), Copyright Amazon.com, Inc. or its affiliates,
# Apache-2.0. See NOTICE.

#' Draw samples from a Gaussian process prior
#'
#' Equivalent of `sample_from_gp_prior_efficient()` in kernel-synth.py.
#' The covariance matrix is factorized by Cholesky decomposition (with an
#' increasing diagonal jitter if needed), falling back to an eigendecomposition
#' with negative eigenvalues clipped to zero (`method = "eigh"` in numpy).
#'
#' @param kernel A `ks_kernel`.
#' @param x Numeric input points.
#' @param n_samples Number of independent sample paths.
#' @param method `"cholesky"` (default, fastest) or `"eigen"`.
#'
#' @return A `length(x)` by `n_samples` matrix.
#' @export
#' @examples
#' Y <- sample_gp_prior(rbf(0.1), seq(0, 1, length.out = 100), n_samples = 3)
#' matplot(Y, type = "l")
sample_gp_prior <- function(kernel, x, n_samples = 1L,
                            method = c("cholesky", "eigen")) {
  method <- match.arg(method)
  n <- length(x)
  K <- kernel_matrix(kernel, x)
  K <- (K + t(K)) / 2
  if (!all(is.finite(K))) stop("non-finite covariance matrix")
  Z <- matrix(stats::rnorm(n * n_samples), n, n_samples)

  if (method == "cholesky") {
    scale <- max(mean(abs(diag(K))), .Machine$double.eps)
    for (jitter in c(0, 1e-10, 1e-8, 1e-6, 1e-4) * scale) {
      R <- tryCatch(chol(if (jitter > 0) K + diag(jitter, n) else K),
                    error = function(e) NULL)
      if (!is.null(R)) return(crossprod(R, Z))
    }
  }
  e <- eigen(K, symmetric = TRUE)
  e$vectors %*% (sqrt(pmax(e$values, 0)) * Z)
}
