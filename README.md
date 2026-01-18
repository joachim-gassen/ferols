# ferols: A Fixed-Effects Robust M Estimator with Huber Loss

**ferols** provides an experimental implementation for a fixed-effects linear 
regression estimator using **Huber M-estimation** with 
**iteratively reweighted least squares (IRLS)**. It is inspired by 
the [package `robtwfe` by David Veenman](https://github.com/dveenman/robtwfe). 
It is build on and designed to integrate tightly with the
[`fixest`](https://lrberge.github.io/fixest/) ecosystem.

It builds on our recent work

> Joachim Gassen and David Veenman (2026): Estimation Precision and 
> Robust Inference in Archival Research, SSRN Working Paper,
> http://dx.doi.org/10.2139/ssrn.4975569.

and combines robust M-estimation (Huber loss) with:

- high-dimensional fixed-effect absorption,
- fast estimation via `fixest::feols`,
- several algorithms to estimate scale in the first step,
- and a `fixest` style variance–covariance interface allowing for
  one-way clustered Huber sandwich standard errors.


## Why?

Standard fixed-effects estimators are sensitive to outliers and heavy-tailed
error distributions. While `fixest` provides fast and reliable estimation
for large linear models with fixed effects, it does not currently offer
robust estimators.

The `ferols` package is a first step to fill this gap. It is inteded to be 
useful for setting with high-dimensional fixedeffects (e.g. panel data with 
unit and time effects) and offers robustinference via Huber ψ/φ sandwich 
variance estimators, including clustered standard errors.


## Disclaimer 

`ferols` is not meant to replace established R packages for robust regression
as these provide much more flexible, tested, and rigorous implementations. 
Rather it aims to fill a narrow gap as these packages currently do not allow
to absorb fixed effects during estimation and might thus run in 
performance/convergence problems in setting with many fixed effects.

This package is under active development.  
The API, implementation details, and numerical behavior may and probably will
change. Use with care and 
**do not rely on it to generate reproducible empirical results yet**.


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

## Not yet implemented (sorted by priority)

- Multi-way clustering
- Alternative loss functions


