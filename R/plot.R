#' Plot a collection of series
#'
#' Quick grid of plots for a list of `ts` objects (e.g. from
#' [kernel_synth()]), titled with their frequency and augmentations.
#'
#' @param series A `ts` or a list of `ts`.
#' @param n Maximum number of series to plot.
#' @param ncol Number of columns of the grid.
#' @param ... Passed to [graphics::plot()].
#'
#' @return `series`, invisibly.
#' @export
#' @examples
#' ys <- kernel_synth(n = 6, length = 200, seed = 1)
#' plot_series(ys)
plot_series <- function(series, n = 9L, ncol = 3L, ...) {
  if (stats::is.ts(series)) series <- list(series)
  m <- min(n, length(series))
  op <- graphics::par(mfrow = c(ceiling(m / ncol), min(ncol, m)),
                      mar = c(2, 2, 2, 0.5), cex.main = 0.8)
  on.exit(graphics::par(op))
  for (i in seq_len(m)) {
    y <- series[[i]]
    aug <- unlist(attr(y, "augmentations"))
    ttl <- sprintf("#%d f=%s%s", i, stats::frequency(y),
                   if (length(aug)) paste0(" [", paste(unique(aug), collapse = ","), "]") else "")
    graphics::plot(y, main = ttl, xlab = "", ylab = "", plot.type = "single",
                   col = if (NCOL(y) > 1) seq_len(NCOL(y)) else 1, ...)
  }
  invisible(series)
}
