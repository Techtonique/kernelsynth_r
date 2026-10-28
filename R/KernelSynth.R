#' KernelSynth R6 generator
#'
#' @description
#' R6 interface to [kernel_synth()]. The generator stores its configuration;
#' successive calls to `$generate()` produce new, reproducible batches when a
#' `seed` is given (the seed is incremented at each call).
#'
#' @export
#' @examples
#' gen <- KernelSynth$new(length = 120, frequency = 12, seed = 123)
#' gen
#' ys <- gen$generate(5)
#' tasks <- gen$tasks(5, h = 12)
#' k <- gen$sample_kernel()
KernelSynth <- R6::R6Class(
  "KernelSynth",
  public = list(
    #' @field length Series length or range `c(min, max)`.
    length = NULL,
    #' @field max_kernels Maximum number of base kernels.
    max_kernels = NULL,
    #' @field bank `"extended"` or `"chronos"`.
    bank = NULL,
    #' @field frequency Frequency (or candidates); `NULL` for random.
    frequency = NULL,
    #' @field randomize Randomize kernel hyperparameters.
    randomize = NULL,
    #' @field augment Apply augmentations.
    augment = NULL,
    #' @field control Augmentation probabilities.
    control = NULL,
    #' @field n_vars Number of variables per series.
    n_vars = NULL,
    #' @field p_changepoint Probability of changepoint combinations.
    p_changepoint = NULL,
    #' @field weights Kernel group weights.
    weights = NULL,
    #' @field seed Master seed (or `NULL`).
    seed = NULL,

    #' @description Create a new generator. Arguments as in [kernel_synth()].
    #' @param length,max_kernels,bank,frequency,randomize,augment,control,n_vars,p_changepoint,weights,seed
    #'   See [kernel_synth()].
    initialize = function(length = 1024L, max_kernels = 5L,
                          bank = c("extended", "chronos"), frequency = NULL,
                          randomize = TRUE, augment = TRUE,
                          control = augmentation_control(), n_vars = 1L,
                          p_changepoint = 0.1, weights = group_weights(),
                          seed = NULL) {
      self$length <- length
      self$max_kernels <- max_kernels
      self$bank <- match.arg(bank)
      self$frequency <- frequency
      self$randomize <- randomize
      self$augment <- augment
      self$control <- control
      self$n_vars <- n_vars
      self$p_changepoint <- p_changepoint
      self$weights <- weights
      self$seed <- seed
      invisible(self)
    },

    #' @description Generate series.
    #' @param n Number of series.
    #' @param n_cores Number of cores.
    #' @param as_mts Bind into a single `mts` when possible.
    #' @return A `ts` or a list of `ts`.
    generate = function(n = 1L, n_cores = 1L, as_mts = FALSE) {
      kernel_synth(n = n, length = self$length, max_kernels = self$max_kernels,
                   bank = self$bank, frequency = self$frequency,
                   randomize = self$randomize, augment = self$augment,
                   control = self$control, n_vars = self$n_vars,
                   p_changepoint = self$p_changepoint, weights = self$weights,
                   seed = private$next_seed(), n_cores = n_cores, as_mts = as_mts)
    },

    #' @description Generate series and split them into forecasting tasks.
    #' @param n Number of series.
    #' @param h Horizon (or function of the frequency).
    #' @param ... Passed to `$generate()`.
    #' @return A list of tasks, see [make_tasks()].
    tasks = function(n = 1L, h, ...) {
      ys <- self$generate(n, ...)
      make_tasks(if (n == 1L) list(ys) else ys, h)
    },

    #' @description Draw a random composite kernel from the configured bank.
    #' @param frequency Frequency used to build the extended bank.
    #' @param length Series length used to build the bank.
    #' @return A `ks_kernel`.
    sample_kernel = function(frequency = NULL, length = NULL) {
      len <- if (is.null(length)) max(self$length) else length
      freq <- if (!is.null(frequency)) frequency
              else if (!is.null(self$frequency)) .pick(self$frequency) else NULL
      b <- if (self$bank == "chronos") kernel_bank_chronos(len)
           else kernel_bank_extended(len, freq, self$randomize)
      sample_composite_kernel(b, self$max_kernels,
                              p_changepoint = if (self$bank == "chronos") 0 else self$p_changepoint,
                              weights = self$weights)
    },

    #' @description Print the configuration.
    #' @param ... Unused.
    print = function(...) {
      cat("<KernelSynth generator>\n")
      cat("  bank        :", self$bank, "\n")
      cat("  length      :", paste(self$length, collapse = "-"), "\n")
      cat("  frequency   :", if (is.null(self$frequency)) "random (1, 4, 12, 52, 7, 24)"
                             else paste(self$frequency, collapse = ", "), "\n")
      cat("  max_kernels :", self$max_kernels, "\n")
      cat("  n_vars      :", self$n_vars, "\n")
      cat("  augment     :", self$augment, "\n")
      cat("  seed        :", if (is.null(self$seed)) "none" else self$seed, "\n")
      invisible(self)
    }
  ),
  private = list(
    calls = 0L,
    next_seed = function() {
      private$calls <- private$calls + 1L
      if (is.null(self$seed)) NULL else self$seed + private$calls - 1L
    }
  )
)
