#' kernelsynth: diverse synthetic time series from random GP kernel compositions
#'
#' An R port and extension of KernelSynth (Chronos, Ansari et al., 2024),
#' aimed at building diverse synthetic corpora for pretraining and
#' meta-learning of forecasting models. All generators return `ts` objects.
#'
#' Main entry points:
#' * [kernel_synth()] / [KernelSynth] (S3 / R6): enriched generator.
#' * [generate_time_series()]: faithful port of the original algorithm.
#' * [kernels], [changepoint()], [warp_kernel()], [simulate.ks_kernel()]:
#'   kernel algebra and GP sampling.
#' * [augment_series()] and [augmentations]: post-hoc augmentations.
#' * [ts_mixup()], [make_tasks()], [meta_features()]: pretraining and
#'   meta-learning utilities.
#'
#' KernelSynth was created by the Chronos authors: Abdul Fatir Ansari, Lorenzo Stella, Caner Turkmen, Xiyuan Zhang, Pedro Mercado, Huibin Shen, Oleksandr Shchur, Syama Sundar Rangapuram, Sebastian Pineda Arango, Shubham Kapoor, Jasper Zschiegner, Danielle C. Maddix, Hao Wang, Michael W. Mahoney, Kari Torkkola, Andrew Gordon Wilson, Michael Bohlke-Schneider, Yuyang Wang.
#'
#' @references Ansari, A. F., Stella, L., Turkmen, C., et al. (2024).
#'   Chronos: Learning the Language of Time Series.
#'   Transactions on Machine Learning Research. \doi{10.48550/arXiv.2403.07815}
#'
#' @importFrom R6 R6Class
#' @keywords internal
"_PACKAGE"
