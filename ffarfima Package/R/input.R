as_curves <- function(x, name, n_grid = NULL) {
  if (inherits(x, "fts")) {
    y <- as.matrix(x$y)
    rownames(y) <- x$x
    x <- y
  }

  if (is.data.frame(x)) {
    x <- as.matrix(x)
  }

  if (!is.matrix(x) || !is.numeric(x)) {
    stop(sprintf("`%s` must be a numeric matrix with grid points in rows and curves (time) in columns, or an `fts` object.",
                 name), call. = FALSE)
  }

  if (!is.null(n_grid) && nrow(x) != n_grid) {
    stop(sprintf("`%s` has %d grid points (rows) but the training curves have %d.", name, nrow(x), n_grid),
         call. = FALSE)
  }

  bad <- sum(!is.finite(x))

  if (bad > 0) {
    stop(sprintf(paste("`%s` has %d missing or infinite values (for example log rates of zero).",
                       "Impute or smooth them before fitting."), name, bad), call. = FALSE)
  }

  x
}

check_trim <- function(n_curves, trim, d) {
  if (d != 0 && n_curves - trim < 10) {
    stop(sprintf(paste("`trim` = %d curves are removed after fractional differencing, leaving %d of the %d",
                       "training curves (at least 10 are needed). Use a smaller `trim`, for example",
                       "trim = %d. If the data are time by grid rather than grid by time, transpose them."),
                 trim, max(0L, n_curves - trim), n_curves, max(1L, floor(n_curves / 3))), call. = FALSE)
  }
}

check_follow <- function(prev, nxt, prev_name, next_name) {
  a <- suppressWarnings(as.numeric(colnames(prev)))
  b <- suppressWarnings(as.numeric(colnames(nxt)))
  if (length(a) < 2 || length(b) < 1 || any(!is.finite(a)) || any(!is.finite(b))) {
    return(invisible(NULL))
  }
  step <- a[length(a)] - a[length(a) - 1]
  if (abs(b[1] - (a[length(a)] + step)) > 1e-8) {
    warning(sprintf("`%s` starts at %s, but the curve after the end of %s (%s) would be %s. The blocks must follow each other in time.",
                    next_name, colnames(nxt)[1], prev_name, colnames(prev)[ncol(prev)], format(a[length(a)] + step)),
            call. = FALSE)
  }
  invisible(NULL)
}

check_args <- function(args, fun, name, exclude) {
  allowed <- setdiff(names(formals(fun)), exclude)
  unknown <- setdiff(names(args), allowed)
  if (length(unknown) > 0 || (length(args) > 0 && is.null(names(args)))) {
    stop(sprintf("Unknown entries in `%s`: %s. Allowed: %s.", name,
                 paste(if (length(unknown)) unknown else "<unnamed>", collapse = ", "),
                 paste(allowed, collapse = ", ")), call. = FALSE)
  }
}

progress_start <- function(verbose, n, text) {
  if (!verbose) {
    return(NULL)
  }
  if (length(text) > 0 && nzchar(text)) {
    message(text)
  }
  if (interactive() && n > 1) utils::txtProgressBar(max = n, style = 3) else NULL
}

progress_step <- function(bar, i) {
  if (!is.null(bar)) {
    utils::setTxtProgressBar(bar, i)
  }
}

progress_end <- function(bar) {
  if (!is.null(bar)) {
    close(bar)
  }
}

#' Split curves into training, validation and test blocks
#'
#' @param x Curves, with grid points in rows and time in columns, or an `fts`
#'   object.
#' @param train,valid Proportions of the curves in the training and validation
#'   blocks; the remaining curves form the test block. Each block length is
#'   rounded down, as in `floor(T * train)`.
#' @return A list with matrices `train`, `valid` and `test` (`test` is `NULL`
#'   when `train + valid = 1`).
#' @examples
#' x <- sim_ffarfima(n = 100, d = 0.2, p = 1, q = 0, no_grid = 21, seed = 1)$X
#' blocks <- split_curves(x, train = 0.6, valid = 0.2)
#' sapply(blocks, ncol)
#' @export
split_curves <- function(x, train = 0.6, valid = 0.2) {
  x <- as_curves(x, "x")

  if (train <= 0 || valid <= 0 || train + valid > 1) {
    stop("`train` and `valid` must be positive proportions with train + valid <= 1.", call. = FALSE)
  }

  train_end <- floor(ncol(x) * train)
  valid_end <- train_end + floor(ncol(x) * valid)

  list(
    train = x[, seq_len(train_end), drop = FALSE],
    valid = x[, (train_end + 1):valid_end, drop = FALSE],
    test = if (valid_end < ncol(x)) x[, (valid_end + 1):ncol(x), drop = FALSE] else NULL
  )
}
