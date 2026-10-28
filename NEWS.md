# kernelsynth 0.1.0

* Initial version.
* Faithful port of Chronos' KernelSynth (`generate_time_series()`,
  `kernel_bank_chronos()`, `random_binary_map()`, `sample_gp_prior()`).
* Enriched generator `kernel_synth()` (S3) and `KernelSynth` (R6):
  frequency-aware extended kernel bank, randomized hyperparameters,
  changepoint combinations, input warping, variable lengths, multivariate
  series (linear model of coregionalization), augmentations, parallel and
  reproducible generation. Returns `ts` objects.
* Kernel algebra (`+`, `*`), `simulate()` method for kernels.
* `ts_mixup()`, `make_tasks()` and `meta_features()` for pretraining and
  meta-learning.
