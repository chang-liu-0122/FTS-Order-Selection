#' ffarfima: forecasting long-memory functional time series
#'
#' Curves are stored as a matrix with grid points (for example ages) in rows
#' and time (for example years) in columns. `fts` objects from the rainbow and
#' ftsa packages are also accepted.
#'
#' @section Workflow:
#' 1. [split_curves()] splits the series into training, validation and
#'    optional test blocks.
#' 2. [ffarfima()] estimates the memory parameter `d` on the training block
#'    ([estimate_memory()]), selects the orders `(p, q)` ([select_order()]),
#'    fits the operators, and calibrates the prediction bands on the
#'    validation block for horizons 1 to `h`.
#' 3. [predict.ffarfima()] forecasts beyond the data and, when a test block is
#'    given, evaluates the forecasts on it.
#'
#' @examples
#' sim <- sim_ffarfima(n = 300, d = 0.3, p = 1, q = 1, seed = 1)
#' blocks <- split_curves(sim$X, train = 0.6, valid = 0.2)
#' fit <- ffarfima(blocks$train, blocks$valid, d = "lpwn2", order = "ic_full",
#'                 p = 0:2, q = 0:2, h = 5, trim = 50)
#' fc <- predict(fit, test = blocks$test, band = "AC")
#' fc
#' plot(fc, h = 1)
#' plot(fc, h = 3, target_test_curve = 10)
#' @keywords internal
#' @importFrom stats ar cov lm.fit lsfit optim optimize predict sd var median coef fitted residuals
"_PACKAGE"
