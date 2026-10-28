# Internal helpers -------------------------------------------------------------

# Equally spaced grid on [0, 1] (numpy.linspace(0, 1, n))
.linspace <- function(n) {
  if (n <= 1L) return(0)
  seq(0, 1, length.out = n)
}

# Log-uniform draws
.rlogunif <- function(n, lo, hi) exp(stats::runif(n, log(lo), log(hi)))

# Safe single-element sampling (avoids the sample(5) pitfall)
.pick <- function(x, size = 1L, replace = FALSE, prob = NULL) {
  x[sample.int(length(x), size = size, replace = replace, prob = prob)]
}

.fmt <- function(x) formatC(x, digits = 4, format = "g")

.safe_sd <- function(v) {
  s <- stats::sd(v, na.rm = TRUE)
  if (!is.finite(s) || s < .Machine$double.eps) 1 else s
}

# Standardize
.std <- function(v) (v - mean(v, na.rm = TRUE)) / .safe_sd(v)

# Rescale `new` to have the mean/sd of `old`
.rescale_like <- function(new, old) {
  mean(old, na.rm = TRUE) + .safe_sd(old) * .std(new)
}

.softplus <- function(z) ifelse(z > 30, z, log1p(exp(z)))

# Apply a numeric -> numeric function to a vector or univariate `ts`,
# preserving time attributes and extra attributes (e.g. "kernel").
.with_ts <- function(y, f) {
  if (stats::is.ts(y)) {
    if (is.matrix(y)) stop("expects a univariate series; use augment_series() for 'mts'")
    tsp_ <- stats::tsp(y)
    keep <- attributes(y)
    keep[c("tsp", "class", "dim", "dimnames", "names")] <- NULL
    out <- stats::ts(f(as.numeric(y)), start = tsp_[1], frequency = tsp_[3])
    for (nm in names(keep)) attr(out, nm) <- keep[[nm]]
    out
  } else {
    f(as.numeric(y))
  }
}

# Copy non-time attributes from one object to another
.copy_attrs <- function(to, from) {
  keep <- attributes(from)
  keep[c("tsp", "class", "dim", "dimnames", "names")] <- NULL
  for (nm in names(keep)) attr(to, nm) <- keep[[nm]]
  to
}

.default_start <- function(frequency) {
  if (frequency %in% c(7, 24)) c(1, 1) else if (frequency == 1) 2000 else c(2000, 1)
}
