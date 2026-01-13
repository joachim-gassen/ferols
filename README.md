# ferols: A Fixed-Effects Robust M Estimator with Huber Loss

**ferols** provides an experimental implementation of fixed-effects linear regression
with **Huber M-estimation** using **iteratively reweighted least squares (IRLS)**. 
It is build on and designed to integrate tightly with the
[`fixest`](https://lrberge.github.io/fixest/) ecosystem.

The package focuses on settings with high-dimensional fixed effects (e.g. panel
data with unit and time effects) and offers robust inference via
Huber ψ/φ sandwich variance estimators, including clustered standard errors.

This package is under active development.  
The API, implementation details, and numerical behavior may change.
Use with care and **do not rely on it to generate reproducible empirical results yet**.

Standard fixed-effects estimators are sensitive to outliers and heavy-tailed
error distributions. While `fixest` provides fast and reliable estimation
for large linear models with fixed effects, it does not currently offer
robust estimators.

`ferols` builds on our recent work

> Joachim Gassen and David Veenman (2026): Estimation Precision and Robust Inference in Archival Research, SSRN Working Paper, http://dx.doi.org/10.2139/ssrn.4975569.

and combines robust M-estimation (specifically Huber loss) with:

- high-dimensional fixed-effect absorption,
- fast estimation via `fixest::feols`,
- and fixest-style variance–covariance interfaces (`vcov()`, clustering, etc.).

The design goal is for `ferols` to **feel like a natural extension of fixest**, not a
stand-alone modeling framework.

---

## Features

- Fixed-effects linear regression with Huber M-estimation
- Uses OLS to estimate scale in the first step
- IRWLS estimation using `fixest::feols()` as a workhorse
- One-way clustered Huber sandwich standard errors
- Fixest-compatible `vcov()`, `summary()`, and `print()` methods

---

## Not yet implemented (sorted by priority)

- Alternative starting points for scale estimation
- Multi-way clustering
- Alternative loss functions

---

## Installation (development version)

`ferols` is not on CRAN. To install the current development version from GitHub:

```r
install.packages("remotes")
remotes::install_github("joachim-gassen/ferols")
```

## Typical use

```r
library(ferols)
mod <- ferols(y ~ x1 + x2 | id + year, data = my_data, vcov = ~ id)
summary(mod)
```
