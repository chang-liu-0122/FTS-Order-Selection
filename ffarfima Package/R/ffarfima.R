#' Fit a functional ARFIMA model and calibrate its prediction bands
#'
#' Uses the training curves to estimate the memory parameter, select the FARMA
#' orders of the fractionally differenced curves and fit the operators. When
#' validation curves are supplied, rolling forecasts are then made from every
#' validation origin up to `h` steps ahead, and their errors calibrate the
#' prediction bands for each horizon and level. Forecasts and evaluation are
#' produced by [predict.ffarfima()].
#'
#' @param train Training curves: a matrix with grid points in rows and curves
#'   (time) in columns, or an `fts` object.
#' @param valid Validation curves that follow `train`. Needed to calibrate the
#'   bands and for `order = "cv"`. Without it, only point forecasts are made.
#' @param d Memory estimator passed to [estimate_memory()] (`"bc_lpwn2"`,
#'   `"bc_lpwn1"`, `"lpwn2"`, `"lpwn1"`, `"lw"`, `"dfa"`), or a known value.
#' @param order Order selection method passed to [select_order()]
#'   (`"ic_full"`, `"ic_auto"`, `"bj"`, `"cv"`), or a known `c(p, q)`.
#' @param p,q Candidate orders searched by the selection method.
#' @param method `"farma"` for the FARMA operator estimator, or `"far"` for a
#'   FAR(1) benchmark fitted by [far::far()] with the number of components
#'   chosen by `far_pve`. The benchmark ignores `order`.
#' @param h Maximum forecast horizon. The bands are calibrated for horizons 1
#'   to `h`, and `order = "cv"` scores the same horizons. The validation block
#'   needs at least `h + 1` curves.
#' @param level Nominal coverage levels of the bands, for example `c(0.8, 0.95)`.
#' @param band Prediction bands to calibrate: any of `"SR"` (scaled residual),
#'   `"NS"` (split conformal) and `"AC"` (adaptive conformal). This is the place
#'   to choose the conformal method; [predict.ffarfima()] uses these bands.
#' @param window `"expanding"` or `"fixed"` estimation window at each rolling
#'   origin, used for calibration, cross-validation and the test evaluation.
#'   The fixed window has the length of the training block.
#' @param valid_eval Optional curves that replace `valid` as the observed values
#'   when calibrating the bands and scoring cross-validation, for example raw
#'   curves when the model is fitted to smoothed ones.
#' @param Npc Number of principal components used by the FARMA estimator.
#' @param trim Number of initial fractionally differenced curves discarded
#'   before estimation, to remove the start-up effect of the truncated filter.
#' @param far_pve Proportion of variance explained that sets the number of
#'   components for `method = "far"`.
#' @param memory_args List of further arguments for [estimate_memory()], such
#'   as `m_const`, `cov_type` or `B`.
#' @param order_args List of further arguments for [select_order()], such as
#'   `n_lags`, `alpha` or `cv_score`.
#' @param eta_grid Candidate step sizes for the adaptive conformal band.
#' @param zeta_grid Candidate multipliers for the scaled residual band.
#' @param seed Optional seed for the bootstrap in the memory estimator and the
#'   white-noise test of the Box-Jenkins method.
#' @param verbose Whether to report progress.
#' @return An object of class `ffarfima` with
#' * `d` and `order`: the memory parameter and the orders `c(p = , q = )`;
#' * `phi`, `theta`: the AR and MA operators, grid by grid by lag;
#' * `fitted`, `fitted_differenced`, `residuals`: the in-sample one-step-ahead
#'   fit on the scale of the data and on the differenced scale, and the
#'   residuals (see [fitted.ffarfima()]);
#' * `center`: the training mean curve, subtracted before differencing and
#'   added back to the forecasts;
#' * `memory`, `order_selection`, `calibration`: the details of each step;
#' * `settings` and `data`: the settings and curves used.
#' @examples
#' sim <- sim_ffarfima(n = 300, d = 0.3, p = 1, q = 0, seed = 1)
#' blocks <- split_curves(sim$X)
#' fit <- ffarfima(blocks$train, blocks$valid, d = "lpwn2", order = "ic_full",
#'                 p = 0:2, q = 0:2, h = 5, trim = 50)
#' fit
#' @export
ffarfima <- function(train,
                     valid = NULL,
                     d = c("bc_lpwn2", "bc_lpwn1", "lpwn2", "lpwn1", "lw", "dfa"),
                     order = c("ic_full", "ic_auto", "bj", "cv"),
                     p = 0:5,
                     q = 0:5,
                     method = c("farma", "far"),
                     h = 10,
                     level = c(0.8, 0.95),
                     band = c("SR", "NS", "AC"),
                     window = c("expanding", "fixed"),
                     valid_eval = NULL,
                     Npc = 1,
                     trim = 100,
                     far_pve = 0.95,
                     memory_args = list(),
                     order_args = list(),
                     eta_grid = c(0, 0.005, 0.01, 0.02, 0.05, 0.1),
                     zeta_grid = seq(0.1, 50, length.out = 2496),
                     seed = NULL,
                     verbose = TRUE) {
  method <- match.arg(method)
  window <- match.arg(window)
  band <- match.arg(band, several.ok = TRUE)
  train <- as_curves(train, "train")
  h <- as.integer(h)

  if (length(h) != 1L || is.na(h) || h < 1L) {
    stop("`h` must be a single positive integer.", call. = FALSE)
  }

  if (!is.numeric(level) || any(level <= 0 | level >= 1)) {
    stop("`level` must lie strictly between 0 and 1, for example level = c(0.8, 0.95).", call. = FALSE)
  }

  if (!is.null(valid)) {
    valid <- as_curves(valid, "valid", nrow(train))
    check_follow(train, valid, "the training curves", "valid")
    if (ncol(valid) - h + 1L < 2L) {
      stop(sprintf(paste("With h = %d the validation block needs at least %d curves, so that two forecast",
                         "origins calibrate the longest horizon; `valid` has %d. Use h <= %d or a longer",
                         "validation block."), h, h + 1L, ncol(valid), ncol(valid) - 1L), call. = FALSE)
    }
  }

  if (!is.null(valid_eval)) {
    if (is.null(valid)) {
      stop("`valid_eval` needs `valid`.", call. = FALSE)
    }
    valid_eval <- as_curves(valid_eval, "valid_eval", nrow(train))
    if (ncol(valid_eval) != ncol(valid)) {
      stop("`valid_eval` must have the same number of curves as `valid`.", call. = FALSE)
    }
  }

  if (method == "far" && !missing(order)) {
    warning("`order` is ignored for the FAR(1) benchmark.", call. = FALSE)
  }

  check_args(memory_args, estimate_memory, "memory_args", c("x", "method", "seed", "verbose"))
  check_args(order_args, select_order, "order_args",
             c("x", "d", "method", "p", "q", "Npc", "trim", "valid", "valid_eval", "h", "window", "seed", "verbose"))
  check_trim(ncol(train), trim, if (is.numeric(d)) d else 1)

  if (!is.null(seed)) {
    withr::local_seed(seed)
  }

  memory <- NULL

  if (is.numeric(d)) {
    d_hat <- as.numeric(d)
  } else {
    memory <- do.call(estimate_memory, c(list(x = train, method = match.arg(d), verbose = verbose),
                                         memory_args))
    d_hat <- memory$d

    if (!isTRUE(attr(memory, "diagnostics")$admissible)) {
      stop(sprintf(paste("The bias-corrected memory estimate is below zero (d = %.4f).",
                         "Choose another estimator with `d =` or supply d directly."), d_hat),
           call. = FALSE)
    }
  }

  if (length(d_hat) != 1L || !is.finite(d_hat)) {
    stop("`d` must be an estimator name or a single finite value.", call. = FALSE)
  }

  order_selection <- NULL

  if (method == "far") {
    p_hat <- 1L
    q_hat <- 0L
  } else if (is.numeric(order)) {
    if (length(order) != 2L || any(order < 0) || all(order == 0)) {
      stop("A known `order` must be c(p, q) with p + q > 0.", call. = FALSE)
    }
    p_hat <- as.integer(order[1])
    q_hat <- as.integer(order[2])
  } else {
    order_selection <- do.call(
      select_order,
      c(list(x = train, d = d_hat, method = match.arg(order), p = p, q = q, Npc = Npc,
             trim = trim, valid = valid, valid_eval = valid_eval, h = h, window = window,
             verbose = verbose), order_args)
    )
    p_hat <- as.integer(order_selection$p)
    q_hat <- as.integer(order_selection$q)
  }

  history <- origin_history(train, ncol(train), d_hat, trim, "expanding", ncol(train))

  est <- if (method == "farma") {
    farma_estimation(history$Y_train, p = p_hat, q = q_hat, m = NULL, Npc = Npc)
  } else {
    forecast_far_path(history$Y_train, h = 1, far_pve = far_pve)
  }

  fitted_curves <- NULL
  fitted_differenced <- NULL
  residual_curves <- NULL

  if (method == "farma") {
    tryCatch(validate_farma_fit(est, history$Y_train),
             farma_unstable_fit = function(e) warning(conditionMessage(e), call. = FALSE))
    cols <- (if (d_hat == 0) 0L else trim) + seq_len(ncol(history$Y_train))
    labels <- list(rownames(train), colnames(train)[cols])
    residual_curves <- est$residuals
    fitted_differenced <- est$fitted_original
    dimnames(residual_curves) <- labels
    dimnames(fitted_differenced) <- labels
    fitted_curves <- train[, cols, drop = FALSE] - residual_curves
  }

  calibration <- NULL

  if (!is.null(valid)) {
    X <- cbind(train, valid)
    X_eval <- cbind(train, if (is.null(valid_eval)) valid else valid_eval)
    roller <- rolling_forecaster(X, d_hat, p_hat, q_hat, method, Npc, far_pve, trim, window, ncol(train))
    calibration <- calibrate_bands(roller$forecast, X_eval, ncol(train), ncol(X), h,
                                   round(1 - level, 10), eta_grid, zeta_grid, verbose)
    failures <- roller$log()
    if (!is.null(failures) && nrow(failures) > 0) {
      warning(sprintf("%d validation forecasts failed and were skipped. The first error was: %s",
                      nrow(failures), failures$message[1]), call. = FALSE)
    }
  }

  center <- history$mean_curve
  names(center) <- rownames(train)

  structure(
    list(
      d = d_hat,
      order = c(p = p_hat, q = q_hat),
      phi = if (method == "farma") est$phi else NULL,
      theta = if (method == "farma") est$theta else NULL,
      fitted = fitted_curves,
      fitted_differenced = fitted_differenced,
      residuals = residual_curves,
      center = center,
      memory = memory,
      order_selection = order_selection,
      calibration = calibration,
      settings = list(method = method, h = h, level = level, band = band, window = window,
                      trim = trim, Npc = Npc, far_pve = far_pve),
      data = list(train = train, valid = valid, valid_eval = valid_eval)
    ),
    class = "ffarfima",
    estimator = est
  )
}

memory_label <- function(m) {
  labels <- c(bc_lpwn2 = "BC-LPWN2", bc_lpwn1 = "BC-LPWN1", lpwn2 = "LPWN2",
              lpwn1 = "LPWN1", lw = "LW", dfa = "DFA")
  ops <- c(sample = "empirical covariance", lrc_RS = "long-run covariance (RS)",
           lrc_LRS = "long-run covariance (LRS)")
  paste0(labels[[m$method]], " on the leading score of the ", ops[[m$cov_type]],
         ", bandwidth T^", m$m_const)
}

order_label <- function(o) {
  labels <- c(ic_full = "functional IC, full search", ic_auto = "functional IC, automatic search",
              bj = "functional Box-Jenkins", cv = "cross-validation")
  extra <- switch(o$method,
                  bj = paste0(", ", o$n_lags, " lags"),
                  cv = paste0(", h = 1..", o$h, ", ", o$window, " window"),
                  "")
  paste0(labels[[o$method]], extra, ", ", o$n_candidates, " candidates")
}

#' @export
print.ffarfima <- function(x, ...) {
  set <- x$settings
  dat <- x$data
  cat("FARFIMA(p = ", x$order[["p"]], ", d = ", format(x$d, digits = 4), ", q = ", x$order[["q"]], ")\n", sep = "")
  cat("  Data:      ", nrow(dat$train), " grid points x ", ncol(dat$train), " training curves",
      if (!is.null(dat$valid)) paste0(", ", ncol(dat$valid), " validation curves"), "\n", sep = "")
  cat("  Memory:    ", if (is.null(x$memory)) "supplied" else memory_label(x$memory), "\n", sep = "")
  if (!is.null(x$memory) && is.finite(x$memory$bias)) {
    cat("             uncorrected d = ", format(x$memory$d_initial, digits = 4),
        ", bootstrap bias = ", format(x$memory$bias, digits = 3), "\n", sep = "")
  }
  if (!is.null(x$memory)) {
    cat(memory_notes(x$memory, "             "), sep = "")
  }
  cat("  Orders:    ", if (set$method == "far") "fixed at (1, 0) for the FAR benchmark"
      else if (is.null(x$order_selection)) "supplied" else order_label(x$order_selection), "\n", sep = "")
  cat("  Operators: ", if (set$method == "far") {
    paste0("FAR(1) via far::far with ", attr(x, "estimator")$kn, " components (", 100 * set$far_pve, "% variance)")
  } else {
    paste0("FARMA estimator with ", set$Npc, " principal component", if (set$Npc > 1) "s")
  }, ", first ", set$trim, " differenced curves trimmed\n", sep = "")
  cat("  Bands:     ", if (is.null(x$calibration)) {
    "none (no validation block), point forecasts only"
  } else {
    paste0(paste(set$band, collapse = ", "), " at ", paste(level_label(set$level), collapse = ", "),
           " for h = 1..", set$h, ", calibrated on ", ncol(dat$valid), " validation curves (", set$window, " window)")
  }, "\n", sep = "")
  invisible(x)
}

#' @export
summary.ffarfima <- function(object, ...) {
  print(object)
  if (!is.null(object$order_selection)) {
    cat("\nCandidates evaluated:\n")
    print(object$order_selection$table, row.names = FALSE, digits = 4)
  }
  invisible(object)
}

#' @export
coef.ffarfima <- function(object, ...) {
  list(d = object$d, phi = object$phi, theta = object$theta)
}

#' Fitted curves and residuals of a FARFIMA fit
#'
#' The in-sample one-step-ahead fit, for the training curves used in the
#' estimation (after the first `trim` fractionally differenced curves are
#' dropped). The first `max(p, q)` of these are `NA`, because the model needs
#' that many earlier curves. Both scales are uncentred, so the curves are
#' directly comparable with the data.
#'
#' * `fitted(fit)` is the fit on the scale of the data: each curve is the
#'   one-step-ahead forecast from the earlier curves, and
#'   `fitted(fit) + residuals(fit)` equals the training curves.
#' * `fitted(fit, scale = "differenced")` is the fit of the short-memory FARMA
#'   part, on the fractionally differenced scale.
#'
#' The residuals are the same on both scales, because the fractional
#' difference filter starts with coefficient one.
#'
#' @param object An `ffarfima` fit with `method = "farma"`.
#' @param scale `"original"` for the scale of the data, or `"differenced"` for
#'   the fractionally differenced scale.
#' @param ... Not used.
#' @return A matrix with grid points in rows and curves in columns.
#' @name fitted.ffarfima
#' @export
fitted.ffarfima <- function(object, scale = c("original", "differenced"), ...) {
  scale <- match.arg(scale)
  if (object$settings$method == "far") {
    stop("Fitted values are only available for the FARMA method.", call. = FALSE)
  }
  if (scale == "original") object$fitted else object$fitted_differenced
}

#' @rdname fitted.ffarfima
#' @export
residuals.ffarfima <- function(object, ...) {
  if (object$settings$method == "far") {
    stop("Residuals are only available for the FARMA method.", call. = FALSE)
  }
  object$residuals
}
