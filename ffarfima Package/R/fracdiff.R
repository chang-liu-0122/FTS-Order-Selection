#' Fractional differencing of a functional time series
#'
#' Applies the fractional difference operator \eqn{(1 - L)^d} to each grid point.
#' A negative `d` integrates.
#'
#' @param x Numeric matrix of curves, with grid points in rows and time in columns.
#' @param d Memory parameter.
#' @return A matrix with the same dimensions as `x`.
#' @examples
#' x <- sim_ffarfima(n = 100, d = 0.3, p = 1, q = 0, seed = 1)$X
#' y <- frac_diff(x, 0.3)
#' max(abs(frac_diff(y, -0.3) - x))
#' @export
frac_diff <- function(x, d) {
  x <- as.matrix(x)
  t(FCVAR::FracDiff(t(x), d = d))
}
