#' Forecast a fitted FARFIMA model, and evaluate it on a test block
#'
#' Forecasts the curves `h` steps beyond the last curve supplied, with the
#' prediction bands calibrated in [ffarfima()]. The band types are chosen with
#' the `band` argument of [ffarfima()]:
#' * `"SR"`, scaled residual: the forecast plus or minus the smallest multiple
#'   of the pointwise validation error standard deviation that reaches the
#'   nominal coverage on the validation block;
#' * `"NS"`, split conformal: the same scale with a width set by the conformal
#'   quantile of normalised non-conformity scores;
#' * `"AC"`, adaptive conformal: the `"NS"` band whose miscoverage level is
#'   updated after each observed error, with a step size chosen on the second
#'   half of the validation block.
#'
#' When `test` is supplied, rolling forecasts are made from every test origin,
#' with `d` and `(p, q)` held fixed and the operators re-estimated at each
#' origin, and evaluated against the observed test curves. The final forecast
#' then starts after the last test curve; otherwise after the last validation
#' curve.
#'
#' @param object An `ffarfima` fit.
#' @param test Optional test curves that follow the validation curves.
#' @param h Maximum forecast horizon, at most the `h` used in [ffarfima()].
#' @param band Band types to return, by default those chosen in [ffarfima()].
#'   Only bands calibrated by the fit can be requested.
#' @param test_eval Optional curves that replace `test` as the observed values
#'   in the evaluation, for example raw curves when the model is fitted to
#'   smoothed ones.
#' @param verbose Whether to report progress.
#' @param ... Not used.
#' @return An object of class `ffarfima_forecast` containing
#' * `point_forecast`: forecast curves beyond the data, grid by `h`;
#' * `bands`: for each band type and level (named like `"80%"`), the `lower`
#'   and `upper` matrices, grid by `h`;
#' * `evaluation` (when `test` is supplied): three tables, `point` (RMSPE and
#'   MAPE by horizon), `interval` (ECP, CPD and interval score by band, level
#'   and horizon) and `wis` (weighted interval score by band and horizon);
#' * `origin`: the time of the last curve, after which the forecast starts;
#' * `settings`: the horizon, levels, bands and window used;
#' * `model`: the fit.
#'
#' Use `plot()` to draw one horizon, either beyond the data or for a curve of
#' the test block (`target_test_curve`), and `as.data.frame()` for long
#' tables, including every forecast made during the test evaluation.
#' @examples
#' sim <- sim_ffarfima(n = 300, d = 0.3, p = 1, q = 0, seed = 1)
#' blocks <- split_curves(sim$X)
#' fit <- ffarfima(blocks$train, blocks$valid, d = 0.3, order = c(1, 0), h = 5, trim = 50)
#' fc <- predict(fit, test = blocks$test)
#' fc
#' head(as.data.frame(fc))
#' @export
predict.ffarfima <- function(object,
                             test = NULL,
                             h = object$settings$h,
                             band = object$settings$band,
                             test_eval = NULL,
                             verbose = TRUE,
                             ...) {
  moved <- intersect(names(list(...)), c("valid", "level", "window", "valid_eval", "eta_grid", "zeta_grid"))
  if (length(moved) > 0) {
    stop(sprintf("`%s` is set in ffarfima(), not predict().", moved[1]), call. = FALSE)
  }

  h_max <- as.integer(h)
  calibrated <- !is.null(object$calibration)
  fit_band <- object$settings$band
  band <- match.arg(band, c("SR", "NS", "AC"), several.ok = TRUE)
  if (any(!band %in% fit_band)) {
    stop(sprintf("Band %s was not chosen in ffarfima(), which calibrated %s. Refit with band = c(...) to add it.",
                 paste(setdiff(band, fit_band), collapse = ", "), paste(fit_band, collapse = ", ")), call. = FALSE)
  }

  if (length(h_max) != 1L || is.na(h_max) || h_max < 1L) {
    stop("`h` must be a single positive integer.", call. = FALSE)
  }

  if (calibrated && h_max > object$settings$h) {
    stop(sprintf("The bands were calibrated for h <= %d in ffarfima(). Use h <= %d, or refit with a larger h.",
                 object$settings$h, object$settings$h), call. = FALSE)
  }

  train <- object$data$train
  valid <- object$data$valid

  if (!is.null(test)) {
    test <- as_curves(test, "test", nrow(train))
    if (is.null(valid)) {
      check_follow(train, test, "the training curves", "test")
    } else {
      check_follow(valid, test, "the validation curves", "test")
    }
  }

  if (!is.null(test_eval)) {
    if (is.null(test)) {
      stop("`test_eval` needs `test`.", call. = FALSE)
    }
    test_eval <- as_curves(test_eval, "test_eval", nrow(train))
    if (ncol(test_eval) != ncol(test)) {
      stop("`test_eval` must have the same number of curves as `test`.", call. = FALSE)
    }
  }

  X <- cbind(train, valid, test)
  X_eval <- cbind(
    train,
    if (is.null(object$data$valid_eval)) valid else object$data$valid_eval,
    if (is.null(test_eval)) test else test_eval
  )

  set <- object$settings
  alpha_set <- round(1 - set$level, 10)
  T_all <- ncol(X)
  p_grid <- nrow(X)
  train_end <- ncol(train)
  valid_end <- train_end + if (is.null(valid)) 0L else ncol(valid)
  has_time <- !is.null(colnames(X)) && all(nzchar(colnames(X)))
  labels <- if (has_time) colnames(X) else as.character(seq_len(T_all))
  roller <- rolling_forecaster(X, object$d, object$order[["p"]], object$order[["q"]], set$method, set$Npc,
                               set$far_pve, set$trim, set$window, train_end)
  calibration <- object$calibration
  aci_err <- if (calibrated) lapply(calibration, function(cal) cal$aci_err) else NULL
  evaluation <- NULL
  test_forecasts <- NULL

  if (!is.null(test)) {
    test_steps <- seq.int(valid_end, T_all - 1L)
    sse <- numeric(h_max)
    sae <- numeric(h_max)
    n_pt <- numeric(h_max)
    test_forecasts <- list()
    acc <- list()
    state <- list()

    for (alpha in alpha_set) {
      a_name <- as.character(alpha)
      acc[[a_name]] <- lapply(seq_len(h_max), function(hh) list(SR = c(0, 0, 0), NS = c(0, 0, 0), AC = c(0, 0, 0)))
      state[[a_name]] <- lapply(seq_len(h_max), function(hh) {
        list(a_t = alpha, origin = integer(0), err = numeric(0), used = logical(0))
      })
    }

    if (verbose) {
      message(sprintf("Evaluating on %d test origins", length(test_steps)))
    }

    bar <- progress_start(verbose && interactive(), length(test_steps), character(0))

    for (t0 in test_steps) {
      progress_step(bar, t0 - valid_end + 1L)
      h_here <- min(h_max, T_all - t0)
      pred <- roller$forecast(t0, h_here, "test")

      if (is.null(pred)) {
        next
      }

      truth <- X_eval[, (t0 + 1):(t0 + h_here), drop = FALSE]
      err <- pred - truth
      ok <- is.finite(err)
      err[!ok] <- 0
      sse[seq_len(h_here)] <- sse[seq_len(h_here)] + colSums(err^2)
      sae[seq_len(h_here)] <- sae[seq_len(h_here)] + colSums(abs(err))
      n_pt[seq_len(h_here)] <- n_pt[seq_len(h_here)] + colSums(ok)
      record <- list(origin = labels[t0], target = labels[(t0 + 1):(t0 + h_here)],
                     target_test_curve = (t0 + 1):(t0 + h_here) - valid_end,
                     point_forecast = pred, observed = truth, bands = list())

      for (alpha in alpha_set) {
        if (!calibrated) {
          next
        }

        a_name <- as.character(alpha)
        cal <- calibration[[a_name]]
        lo <- list(SR = pred * NA, NS = pred * NA, AC = pred * NA)
        up <- lo

        for (hh in seq_len(h_here)) {
          gamma_h <- as.numeric(cal$gamma[[hh]])
          f_h <- pred[, hh]
          y_h <- truth[, hh]

          if (length(gamma_h) != length(f_h)) {
            next
          }

          a <- acc[[a_name]][[hh]]

          if (is.finite(cal$zeta[hh])) {
            hw <- cal$zeta[hh] * gamma_h
            e_h <- f_h - y_h
            keep <- is.finite(e_h) & is.finite(hw)
            if (any(keep)) {
              a$SR <- a$SR + c(sum(abs(e_h[keep]) <= hw[keep]),
                               sum(2 * hw[keep] + 2 / alpha * pmax(abs(e_h[keep]) - hw[keep], 0)),
                               sum(keep))
            }
            lo$SR[, hh] <- f_h - hw
            up$SR[, hh] <- f_h + hw
          }

          band_ns <- conformal_interval(f_h, gamma_h, cal$Q_hat[hh])
          a$NS <- a$NS + band_tally(y_h, band_ns, alpha)
          lo$NS[, hh] <- band_ns$lower
          up$NS[, hh] <- band_ns$upper

          st <- state[[a_name]][[hh]]
          for (j in which(!st$used & st$origin + hh <= t0)) {
            st$a_t <- st$a_t + cal$eta * (alpha - st$err[j])
            st$used[j] <- TRUE
          }

          band_ac <- conformal_interval(f_h, gamma_h, aci_quantile(cal$scores[[hh]], st$a_t))
          a$AC <- a$AC + band_tally(y_h, band_ac, alpha)
          lo$AC[, hh] <- band_ac$lower
          up$AC[, hh] <- band_ac$upper

          keep <- is.finite(y_h) & is.finite(band_ac$lower) & is.finite(band_ac$upper)
          if (any(keep)) {
            st$origin <- c(st$origin, t0)
            st$err <- c(st$err, mean(y_h[keep] < band_ac$lower[keep] | y_h[keep] > band_ac$upper[keep]))
            st$used <- c(st$used, FALSE)
          }

          state[[a_name]][[hh]] <- st
          acc[[a_name]][[hh]] <- a
        }

        for (m in band) {
          record$bands[[m]][[level_label(1 - alpha)]] <- list(lower = lo[[m]], upper = up[[m]])
        }
      }

      test_forecasts[[length(test_forecasts) + 1L]] <- record
    }

    progress_end(bar)

    if (calibrated) {
      aci_err <- lapply(state, function(s) lapply(s, function(st) st$err))
    }

    point <- data.frame(h = seq_len(h_max), RMSPE = sqrt(sse / n_pt), MAPE = sae / n_pt)
    interval <- NULL
    wis <- NULL

    if (calibrated) {
      interval <- do.call(rbind, lapply(alpha_set, function(alpha) {
        do.call(rbind, lapply(band, function(m) {
          do.call(rbind, lapply(seq_len(h_max), function(hh) {
            v <- acc[[as.character(alpha)]][[hh]][[m]]
            ecp <- if (v[3] > 0) v[1] / v[3] else NA_real_
            data.frame(band = m, level = 1 - alpha, h = hh, ECP = ecp,
                       CPD = abs(ecp - (1 - alpha)), IS = if (v[3] > 0) v[2] / v[3] else NA_real_)
          }))
        }))
      }))

      wis <- do.call(rbind, lapply(band, function(m) {
        do.call(rbind, lapply(seq_len(h_max), function(hh) {
          rows <- interval[interval$band == m & interval$h == hh, ]
          data.frame(band = m, h = hh,
                     WIS = (0.5 * point$MAPE[hh] + sum((1 - rows$level) / 2 * rows$IS)) /
                       (nrow(rows) + 0.5))
        }))
      }))
    }

    evaluation <- list(point = point, interval = interval, wis = wis)
  }

  if (verbose && !calibrated) {
    message("No validation block in the fit: point forecasts only.")
  }

  mean_future <- roller$forecast(T_all, h_max, "forecast")
  if (is.null(mean_future)) {
    stop("The final forecast failed; see the fit for details.", call. = FALSE)
  }
  dimnames(mean_future) <- list(rownames(train), future_labels(colnames(X), h_max))
  bands <- list()

  if (calibrated) {
    for (alpha in alpha_set) {
      a_name <- as.character(alpha)
      cal <- calibration[[a_name]]
      lev <- level_label(1 - alpha)
      lower <- upper <- list(SR = mean_future * NA, NS = mean_future * NA, AC = mean_future * NA)

      for (hh in seq_len(h_max)) {
        gamma_h <- as.numeric(cal$gamma[[hh]])
        if (length(gamma_h) != p_grid) {
          next
        }
        f_h <- mean_future[, hh]
        lower$SR[, hh] <- f_h - cal$zeta[hh] * gamma_h
        upper$SR[, hh] <- f_h + cal$zeta[hh] * gamma_h
        b <- conformal_interval(f_h, gamma_h, cal$Q_hat[hh])
        lower$NS[, hh] <- b$lower
        upper$NS[, hh] <- b$upper
        err_h <- aci_err[[a_name]][[hh]]
        a_t <- alpha + cal$eta * sum(alpha - err_h[is.finite(err_h)])
        b <- conformal_interval(f_h, gamma_h, aci_quantile(cal$scores[[hh]], a_t))
        lower$AC[, hh] <- b$lower
        upper$AC[, hh] <- b$upper
      }

      for (m in band) {
        bands[[m]][[lev]] <- list(lower = lower[[m]], upper = upper[[m]])
      }
    }
  }

  failures <- roller$log()
  if (!is.null(failures) && nrow(failures) > 0) {
    warning(sprintf("%d forecast origins failed and were skipped. The first error was: %s",
                    nrow(failures), failures$message[1]), call. = FALSE)
  }

  structure(
    list(
      point_forecast = mean_future,
      bands = bands,
      evaluation = evaluation,
      origin = labels[T_all],
      settings = list(h = h_max, level = set$level, band = if (calibrated) band else character(0),
                      window = set$window),
      model = object
    ),
    class = "ffarfima_forecast",
    test_forecasts = test_forecasts,
    n_test = if (is.null(test)) 0L else ncol(test),
    last_curve = X[, T_all],
    has_time = has_time,
    failures = failures
  )
}

level_label <- function(level) {
  paste0(format(100 * level, trim = TRUE), "%")
}

future_labels <- function(time, h) {
  t_num <- suppressWarnings(as.numeric(time))
  steps <- diff(t_num)
  if (length(t_num) > 1 && all(is.finite(t_num)) && all(abs(steps - steps[1]) < 1e-8)) {
    format(t_num[length(t_num)] + steps[1] * seq_len(h), trim = TRUE)
  } else {
    paste0("h", seq_len(h))
  }
}

band_tally <- function(y, band, alpha) {
  keep <- is.finite(y) & is.finite(band$lower) & is.finite(band$upper)
  if (!any(keep)) {
    return(c(0, 0, 0))
  }
  c(sum(y[keep] >= band$lower[keep] & y[keep] <= band$upper[keep]),
    sum(interval_score_cells(y[keep], band$lower[keep], band$upper[keep], alpha)),
    sum(keep))
}

band_names <- c(SR = "scaled residual", NS = "split conformal", AC = "adaptive conformal")

#' @export
print.ffarfima_forecast <- function(x, ...) {
  m <- x$model
  set <- x$settings
  cat("Forecasts from FARFIMA(p = ", m$order[["p"]], ", d = ", format(m$d, digits = 4), ", q = ", m$order[["q"]],
      ")\n", sep = "")
  cat("  Horizons:  h = 1..", set$h, " after curve ", x$origin, " (", set$window, " window)\n", sep = "")
  if (length(set$band) > 0) {
    cat("  Bands:     ", paste0(set$band, " (", band_names[set$band], ")", collapse = ", "), "\n", sep = "")
    cat("  Levels:    ", paste(level_label(set$level), collapse = ", "), "\n", sep = "")
  } else {
    cat("  Bands:     none (the fit has no validation block)\n")
  }

  if (!is.null(x$evaluation)) {
    ev <- x$evaluation
    cat("\nTest evaluation over ", attr(x, "n_test"), " test curves, averaged over horizons:\n", sep = "")
    cat("  RMSPE = ", format(mean(ev$point$RMSPE), digits = 4),
        ", MAPE = ", format(mean(ev$point$MAPE), digits = 4), "\n", sep = "")
    if (!is.null(ev$interval)) {
      tab <- stats::aggregate(ev$interval[, c("ECP", "CPD", "IS")],
                              by = ev$interval[, c("band", "level")], FUN = mean)
      tab <- tab[order(tab$level, match(tab$band, c("SR", "NS", "AC"))), ]
      tab$level <- level_label(tab$level)
      print(tab, row.names = FALSE, digits = 4)
      w <- tapply(ev$wis$WIS, ev$wis$band, mean)[set$band]
      cat("  WIS: ", paste0(names(w), " = ", format(w, digits = 4), collapse = ", "), "\n", sep = "")
    }
  }

  invisible(x)
}

#' @export
summary.ffarfima_forecast <- function(object, ...) {
  print(object)
  if (!is.null(object$evaluation)) {
    cat("\nPoint accuracy by horizon:\n")
    print(object$evaluation$point, row.names = FALSE, digits = 4)
  }
  invisible(object)
}

#' Long tables of forecasts and bands
#'
#' @param x An `ffarfima_forecast` object.
#' @param row.names,optional Not used.
#' @param test `FALSE` for the forecast beyond the data; `TRUE` for the
#'   forecasts made during the test evaluation, with the observed curves.
#' @param ... Not used.
#' @return A data frame with one row per grid point, horizon, band and level.
#'   Columns are `band`, `level`, `h`, `target` (the time being forecast),
#'   `grid`, `point_forecast`, `lower` and `upper`. With `test = TRUE` there are
#'   also `origin` (the time the forecast was made from), `target_test_curve`
#'   (the position in the test block of the curve being forecast) and
#'   `observed`.
#' @export
as.data.frame.ffarfima_forecast <- function(x, row.names = NULL, optional = FALSE, test = FALSE, ...) {
  grid <- if (is.null(rownames(x$point_forecast))) seq_len(nrow(x$point_forecast)) else rownames(x$point_forecast)
  long <- function(pf, bands, target) {
    base <- data.frame(h = rep(seq_along(target), each = length(grid)), target = rep(target, each = length(grid)),
                       grid = rep(grid, length(target)), point_forecast = as.vector(pf))
    if (length(bands) == 0) {
      return(base)
    }
    do.call(rbind, lapply(names(bands), function(m) {
      do.call(rbind, lapply(names(bands[[m]]), function(lev) {
        b <- bands[[m]][[lev]]
        data.frame(band = m, level = as.numeric(sub("%", "", lev)) / 100, base,
                   lower = as.vector(b$lower), upper = as.vector(b$upper))
      }))
    }))
  }
  if (!test) {
    return(long(x$point_forecast, x$bands, colnames(x$point_forecast)))
  }
  records <- attr(x, "test_forecasts")
  if (is.null(records)) {
    stop("This forecast has no test evaluation; call predict() with `test`.", call. = FALSE)
  }
  do.call(rbind, lapply(records, function(r) {
    out <- data.frame(origin = r$origin, long(r$point_forecast, r$bands, r$target))
    out$target_test_curve <- r$target_test_curve[out$h]
    out$observed <- rep(as.vector(r$observed), length.out = nrow(out))
    out
  }))
}

#' Plot a forecast curve with its prediction band
#'
#' Without `target_test_curve`, draws the forecast `h` steps beyond the last
#' curve. With `target_test_curve = k`, draws the forecast of the `k`-th curve
#' of the test block that was made `h` steps earlier, together with that
#' observed curve.
#'
#' @param x An `ffarfima_forecast` object.
#' @param h Horizon to display.
#' @param band Band type, one of those returned by [predict.ffarfima()].
#' @param level Coverage level of the band.
#' @param target_test_curve Optional position in the test block of the curve
#'   being forecast, counting from 1, between `h` and the number of test
#'   curves. The forecast was made from the curve `h` positions earlier.
#' @param argvals Grid values for the horizontal axis. Defaults to the row
#'   names of the curves when they are numeric (for example ages).
#' @param xlab,ylab,main Axis labels and title.
#' @param legend Position of the legend, for example `"topleft"`, or `NULL`
#'   for no legend.
#' @param ... Further graphical parameters passed to [graphics::plot()].
#' @return `x`, invisibly.
#' @export
plot.ffarfima_forecast <- function(x, h = 1, band = x$settings$band[1], level = x$settings$level[1],
                                   target_test_curve = NULL, argvals = NULL, xlab = "Grid", ylab = "Value",
                                   main = NULL, legend = "bottomright", ...) {
  set <- x$settings
  has_band <- length(set$band) > 0
  has_time <- isTRUE(attr(x, "has_time"))
  if (has_band && !band %in% set$band) {
    stop(sprintf("`band` must be one of the bands returned by predict(): %s.",
                 paste(set$band, collapse = ", ")), call. = FALSE)
  }
  if (h < 1 || h > set$h) {
    stop(sprintf("`h` must be between 1 and %d.", set$h), call. = FALSE)
  }

  if (is.null(target_test_curve)) {
    f <- x$point_forecast[, h]
    b <- if (has_band) x$bands[[band]][[level_label(level)]] else NULL
    lo <- if (is.null(b)) NULL else b$lower[, h]
    hi <- if (is.null(b)) NULL else b$upper[, h]
    ref <- attr(x, "last_curve")
    ref_label <- if (has_time) paste("last observed:", x$origin) else "last observed curve"
    ref_lty <- 2
    when <- if (has_time) {
      paste0(colnames(x$point_forecast)[h], " forecast from ", x$origin, " (h = ", h, ")")
    } else {
      paste0("h = ", h, " after the last curve")
    }
  } else {
    records <- attr(x, "test_forecasts")
    if (is.null(records)) {
      stop("`target_test_curve` needs a forecast evaluated on a test block; call predict() with `test`.",
           call. = FALSE)
    }
    n_test <- attr(x, "n_test")
    k <- target_test_curve
    if (length(k) != 1 || !is.numeric(k) || k < h || k > n_test) {
      stop(sprintf("A %d-step forecast exists for target test curves %d to %d.", h, h, n_test), call. = FALSE)
    }
    hit <- Filter(function(r) length(r$target_test_curve) >= h && r$target_test_curve[h] == k, records)
    if (length(hit) == 0) {
      stop(sprintf("The %d-step forecast of test curve %d failed.", h, k), call. = FALSE)
    }
    r <- hit[[1]]
    f <- r$point_forecast[, h]
    b <- if (has_band) r$bands[[band]][[level_label(level)]] else NULL
    lo <- if (is.null(b)) NULL else b$lower[, h]
    hi <- if (is.null(b)) NULL else b$upper[, h]
    ref <- r$observed[, h]
    ref_label <- if (has_time) paste("observed:", r$target[h]) else paste("observed test curve", k)
    ref_lty <- 1
    when <- if (has_time) {
      paste0(r$target[h], " forecast from ", r$origin, " (test curve ", k, ")")
    } else {
      paste0("Test curve ", k, ", forecast ", h, " step", if (h > 1) "s", " ahead")
    }
  }

  if (has_band && is.null(b)) {
    stop(sprintf("`level` must be one of %s.", paste(set$level, collapse = ", ")), call. = FALSE)
  }

  if (is.null(argvals)) {
    argvals <- suppressWarnings(as.numeric(rownames(x$point_forecast)))
    if (length(argvals) == 0 || any(!is.finite(argvals))) {
      argvals <- seq_len(nrow(x$point_forecast))
    }
  }

  if (is.null(main)) {
    main <- if (has_band) paste0(when, ", ", level_label(level), " ", band) else when
  }

  yr <- range(c(lo, hi, f, ref), finite = TRUE)
  graphics::plot(argvals, f, type = "n", ylim = yr, xlab = xlab, ylab = ylab, main = main, ...)
  if (has_band) {
    ok <- is.finite(lo) & is.finite(hi)
    graphics::polygon(c(argvals[ok], rev(argvals[ok])), c(lo[ok], rev(hi[ok])),
                      col = grDevices::adjustcolor("steelblue", 0.3), border = NA)
  }
  graphics::lines(argvals, ref, col = if (ref_lty == 1) "black" else "grey50", lty = ref_lty)
  graphics::lines(argvals, f, col = "steelblue4", lwd = 2)
  if (!is.null(legend)) {
    labs <- c(ref_label, "forecast", if (has_band) paste(level_label(level), band, "band"))
    graphics::legend(legend, labs,
                     col = c(if (ref_lty == 1) "black" else "grey50", "steelblue4",
                             if (has_band) grDevices::adjustcolor("steelblue", 0.3)),
                     lty = c(ref_lty, 1, if (has_band) NA), lwd = c(1, 2, if (has_band) NA),
                     pch = c(NA, NA, if (has_band) 15), pt.cex = 2, bty = "n")
  }
  invisible(x)
}
