# spikedGLHT

An R package for testing general linear hypotheses under high-dimensional
multivariate regression with spiked noise covariance, implementing
[Li, Aue, Paul and Peng (2024), JASA](https://doi.org/10.1080/01621459.2023.2278825).

The package computes regularized likelihood ratio (LR), Lawley–Hotelling (LH),
and Bartlett–Nanda–Pillai (BNP) tests with minimax weights. It supports normal
and parametric bootstrap calibration, KN and PLY spike selection, and nonzero
null hypotheses. It works when the response dimension exceeds the sample size.

## Installation

From a checkout of this package directory:

```sh
Rscript -e 'install.packages(c("quadprog", "RMTstat"), repos="https://cloud.r-project.org")'
R CMD INSTALL .
```

Install from GitHub:

```r
install.packages("remotes")
remotes::install_github("HaoranLi/spikedGLHT")
```
RMTstat is optional when specifying `K` or using `method = "PLY"`.

## Two-sample example

```r
library(spikedGLHT)
set.seed(123)
p <- 50
N <- 100
X <- rbind(c(rep(1, 40), rep(0, 60)),
           c(rep(0, 40), rep(1, 60)))
C <- matrix(c(1, -1), ncol = 1)
Y <- matrix(rnorm(p * N), p, N)
Y[1:3, ] <- sweep(Y[1:3, ], 1, sqrt(c(20, 12, 8)), "*")

fit <- spiked_glht(Y, X, C, K = 3, rho = 0.5)
fit
fit$p.value
fit$weights

# Parametric bootstrap, recommended for skewed null distributions.
set.seed(124)
boot <- spiked_glht(Y, X, C, K = 3, rho = 0.5,
                    calibration = "bootstrap", B = 1999)
boot

# Estimate the number of spikes, searching up to 5.
auto <- spiked_glht(Y, X, C, method = "KN", K_max = 5)
```

Columns are subjects: `Y` is p × N, `X` is m × N, and `C` is m × q.
Include any intercept in `X`. The null is `B %*% C = 0`; supply a p × q
matrix `D` to test `B %*% C = D`. Multiple columns of `C` form a joint test.
These functions expect numeric matrices, not formulas or data frames.

`rho` is the paper's varrho: the largest allowed signal fraction in the
nonspiked subspace. Choose it based on the alternative of interest; the
simulation study considers 0.2, 0.5, and 0.8. It is not a significance level.
`alpha` controls the rejection threshold.

`statistic` contains standardized statistics, `raw.statistic` contains the
unstandardized statistics, and `p.value` uses the requested calibration.
`reject` compares the raw statistic with its critical value. Bootstrap p-values
use a plus-one correction; rejection uses the paper's empirical quantile rule,
so they can differ near a Monte Carlo boundary. The three reported tests are
not adjusted for multiplicity. Use `set.seed()` for bootstrap reproducibility.

The theory assumes independent Gaussian errors, a common spiked covariance,
fixed m, q and K, and distinct supercritical spikes. Calibration is approximate.
K=0 uses a documented isotropic extension; rho is inapplicable in that case.
The bootstrap representation requires p ≥ 2K.

## Development

```sh
R CMD build .
R CMD check --no-manual spikedGLHT_0.1.0.tar.gz
```

Tests use base R and run through `R CMD check`. GitHub Actions runs the package
checks on pushes and pull requests. A small empirical-size simulation is in
`inst/examples/empirical-size.R`. The original paper and research scripts are
kept outside this package; no HCP data or manuscript files are distributed.

## Citation and license

Run `citation("spikedGLHT")` for the full JASA citation.
The package is distributed under GPL-3. Author and maintainer metadata use the
published article's author list and corresponding-author email; verify the
contact details before release.
