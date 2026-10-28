# kernelsynth

Diverse synthetic time series from random Gaussian process kernel compositions, for pretraining and meta-learning of forecasting models. An R port and extension of [KernelSynth](https://github.com/amazon-science/chronos-forecasting/blob/main/scripts/kernel-synth.py) (Chronos, Ansari et al., 2024) with S3 and R6 interfaces. **All generators return `ts` objects.**

## 1. Install

```R
install.packages("remotes")
remotes::install_github("Techtonique/kernelsynth_r")
```

No Python needed: kernels and GP sampling are implemented in base R (only `R6` as a dependency).

## 2. Use it

```R
library(kernelsynth)

# Original KernelSynth (faithful port of kernel-synth.py)
set.seed(123)
y <- generate_time_series(max_kernels = 5, length = 1024)
attr(y, "kernel")
plot(y)

# S3 interface: enriched generator
ys <- kernel_synth(n = 100, length = c(64, 512), seed = 42)  # list of ts
plot_series(ys)
ks_info(ys[[1]])                                              # kernel, augmentations, ...

# monthly series only, as a single mts
m <- kernel_synth(n = 10, length = 120, frequency = 12, as_mts = TRUE, seed = 1)

# multivariate (3 variables, hourly)
mv <- kernel_synth(length = 168, frequency = 24, n_vars = 3, seed = 2)

# R6 interface
gen <- KernelSynth$new(length = 144, frequency = 12, seed = 2026)
batch <- gen$generate(1000, n_cores = 4)
tasks <- gen$tasks(100, h = 12)          # x / xx / h, M4metalearning-style

# Kernel algebra + simulate()
k <- exp_sine_squared(periodicity = 12 / 143) * rbf(0.5) + dot_product(0)
y <- simulate(k, length = 144, frequency = 12, seed = 3)

# Pretraining & meta-learning utilities
mixed <- ts_mixup(ys, n = 50, k_max = 3)  # TSMixup
mf <- meta_features(ys)                   # trend/seasonal strength, entropy, ...
```

## 3. What is added to KernelSynth

| | KernelSynth (Chronos) | kernelsynth |
|---|---|---|
| Kernel bank | 33 fixed kernels | + Matern (1/2, 3/2, 5/2), Brownian, polynomial, cosine, spectral mixture, locally periodic; grouped and weighted |
| Hyperparameters | fixed | random (log-uniform priors) |
| Periods | fixed, frequency-agnostic | aligned with the `ts` frequency (1, 4, 12, 52, 7, 24, or any) |
| Combinations | `+`, `*` | `+`, `*`, changepoints, input warping |
| Post-processing | none | trends, level shifts, positivity / multiplicative, intermittency, counts, heteroscedasticity, noise, outliers, scaling, missing values |
| Shape | length 1024, univariate | variable lengths, multivariate (coregionalization) |
| Extras | | TSMixup, task splits, meta-features, reproducible parallel generation |

See the vignette: `vignette("kernelsynth")`.

## Acknowledgments

KernelSynth was created by the Chronos authors: Abdul Fatir Ansari, Lorenzo Stella, Caner Turkmen, Xiyuan Zhang, Pedro Mercado, Huibin Shen, Oleksandr Shchur, Syama Sundar Rangapuram, Sebastian Pineda Arango, Shubham Kapoor, Jasper Zschiegner, Danielle C. Maddix, Hao Wang, Michael W. Mahoney, Kari Torkkola, Andrew Gordon Wilson, Michael Bohlke-Schneider, Yuyang Wang.

This software contains code derived from [Chronos](https://github.com/amazon-science/chronos-forecasting) (`scripts/kernel-synth.py`), licensed under the Apache License 2.0; kernels follow [scikit-learn](https://scikit-learn.org/)'s parameterization. See the [NOTICE](NOTICE) file for details.

Ansari, A. F., Stella, L., Turkmen, C., et al. (2024). Chronos: Learning the Language of Time Series. *Transactions on Machine Learning Research*.
