# ffarfima

Memory estimation, order selection and forecasting for long-memory functional
time series, modelled as functional ARFIMA (FARFIMA) processes.

## Installation

```r
install.packages(c("FCVAR", "far", "ftsa", "longmemo", "rainbow", "sde", "vars", "withr", "wwntests"))
install.packages("ffarfima_0.0.0.9000.tar.gz", repos = NULL, type = "source")
library(ffarfima)
vignette("ffarfima")
```

## Data

Curves go in a matrix with **grid points in rows and time in columns**, for
example ages by years. `fts` objects from rainbow and ftsa also work.

## Quick start

```r
library(ffarfima)
library(demography)

x <- log10(smooth.demogdata(extract.years(fr.mort, 1901:2006))$rate$female[1:101, ])
blocks <- split_curves(x, train = 0.6, valid = 0.2)

fit <- ffarfima(blocks$train, blocks$valid, h = 10, level = c(0.8, 0.95), trim = 20, seed = 1)
fit

fc <- predict(fit, test = blocks$test)
fc
plot(fc, h = 1, band = "AC", xlab = "Age", ylab = "log10 death rate")
plot(fc, h = 3, target_test_curve = 11, band = "AC", xlab = "Age", ylab = "log10 death rate")
```

`ffarfima()` uses the training block to estimate the model and the validation
block to calibrate the bands for horizons 1 to `h`. Without `test`, `predict()`
forecasts beyond the last validation curve. With `test`, it also reports RMSPE,
MAPE, ECP, CPD, interval score and WIS on the test block, and
`plot(target_test_curve = k)` shows the forecast of the k-th test curve made
`h` steps earlier, against the observed curve.

## Choices

| Step | Argument | Options |
|---|---|---|
| Memory parameter | `d` | `"bc_lpwn2"` (default), `"bc_lpwn1"`, `"lpwn2"`, `"lpwn1"`, `"lw"`, `"dfa"`, or a number |
| Orders (p, q) | `order` | `"ic_full"` (default), `"ic_auto"`, `"bj"`, `"cv"`, or `c(p, q)` |
| Search range | `p`, `q` | e.g. `p = 0:5, q = 0:5` |
| Tuning | `memory_args`, `order_args` | e.g. `list(m_const = 0.75, B = 100)`, `list(n_lags = 20, alpha = 0.05)` |
| Horizon and levels | `h`, `level`, `window` | e.g. `h = 10, level = c(0.8, 0.95), window = "expanding"` |
| Bands | `band` in `ffarfima()` | `"SR"` scaled residual, `"NS"` split conformal, `"AC"` adaptive conformal |

A FAR(1) benchmark on d = 1 is fitted with the same interface:

```r
bench <- ffarfima(blocks$train, blocks$valid, d = 1, method = "far", trim = 20)
predict(bench, test = blocks$test)
```
