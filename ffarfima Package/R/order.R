#' Select the short-memory orders of a FARFIMA model
#'
#' The training curves are centred and fractionally differenced with `d`, and
#' FARMA(p, q) candidates are compared on the differenced series.
#'
#' Methods:
#' * `"ic_full"`: adaptive functional Hannan-Quinn criterion over every candidate.
#' * `"ic_auto"`: the same criterion with a local search from a few starting
#'   models; the noise level is taken from the largest candidate.
#' * `"bj"`: functional Box-Jenkins. Candidates are tried in increasing `p + q`
#'   (larger `p` first) and the first whose residuals pass the functional
#'   portmanteau white-noise test is returned.
#' * `"cv"`: forecast-based cross-validation. Forecasts are made from every
#'   validation origin for the horizons whose target lies in the validation
#'   block, and candidates are ranked by the mean over horizons of the h-step
#'   error.
#'
#' @param x Training curves, with grid points in rows and time in columns.
#' @param d Memory parameter, or an `ffarfima_memory` object.
#' @param method `"ic_full"`, `"ic_auto"`, `"bj"` or `"cv"`.
#' @param p,q Candidate autoregressive and moving-average orders. The pair
#'   (0, 0) is always excluded.
#' @param Npc Number of principal components used by the FARMA estimator.
#' @param trim Number of initial fractionally differenced curves discarded.
#' @param n_lags Number of residual autocovariance lags. For the criteria the
#'   default is `floor(log(T))` of each candidate's residual series; for
#'   `"bj"` it is 20.
#' @param alpha Significance level of the white-noise test (`"bj"`).
#' @param iid Whether the white-noise test assumes iid rather than weak white
#'   noise under the null (`"bj"`).
#' @param starts Data frame of starting models with columns `p` and `q`
#'   (`"ic_auto"`). Defaults to (2, 2), (1, 0) and (0, 1).
#' @param valid Validation curves following `x` (`"cv"`).
#' @param valid_eval Optional curves that replace `valid` as the observed values
#'   when scoring the cross-validation forecasts (`"cv"`).
#' @param h Maximum forecast horizon (`"cv"`).
#' @param window `"expanding"` or `"fixed"` estimation window (`"cv"`).
#' @param cv_score `"rmse"` or `"mae"` (`"cv"`).
#' @param seed Optional seed. The white-noise test of `"bj"` uses Monte Carlo
#'   integration, so its p-values vary between runs without one.
#' @param verbose Whether to report progress.
#' @return An object of class `ffarfima_order` with the selected `p` and `q`
#'   and a `table` of all candidates evaluated.
#' @examples
#' x <- sim_ffarfima(n = 300, d = 0, p = 1, q = 0, seed = 1)$X
#' select_order(x, d = 0, method = "ic_full", p = 0:2, q = 0:2)
#' @export
select_order <- function(x,
                         d,
                         method = c("ic_full", "ic_auto", "bj", "cv"),
                         p = 0:5,
                         q = 0:5,
                         Npc = 1,
                         trim = 100,
                         n_lags = NULL,
                         alpha = 0.05,
                         iid = FALSE,
                         starts = NULL,
                         valid = NULL,
                         valid_eval = NULL,
                         h = 10,
                         window = c("expanding", "fixed"),
                         cv_score = c("rmse", "mae"),
                         seed = NULL,
                         verbose = FALSE) {
  method <- match.arg(method)
  window <- match.arg(window)
  cv_score <- match.arg(cv_score)
  x <- as_curves(x, "x")

  if (inherits(d, "ffarfima_memory")) {
    d <- d$d
  }

  if (length(d) != 1L || !is.finite(d)) {
    stop("`d` must be a single finite value.", call. = FALSE)
  }

  check_trim(ncol(x), trim, d)

  if (!is.null(seed)) {
    withr::local_seed(seed)
  }

  p <- sort(unique(as.integer(p)))
  q <- sort(unique(as.integer(q)))
  candidates <- expand.grid(q = q, p = p)[, c("p", "q")]
  candidates <- candidates[!(candidates$p == 0 & candidates$q == 0), ]
  rownames(candidates) <- NULL

  if (nrow(candidates) == 0L || any(candidates < 0)) {
    stop("`p` and `q` must be non-negative and give at least one candidate.", call. = FALSE)
  }

  out <- if (method == "cv") {
    if (is.null(valid)) {
      stop("`valid` is required for cross-validation.", call. = FALSE)
    }
    select_order_cv(x, as_curves(valid, "valid", nrow(x)), d, candidates, Npc, trim, h, window,
                    cv_score, verbose,
                    if (is.null(valid_eval)) NULL else as_curves(valid_eval, "valid_eval", nrow(x)))
  } else {
    Y <- origin_history(x, ncol(x), d, trim, "expanding", ncol(x))$Y_train

    if (verbose) {
      message(sprintf("Selecting (p, q) by %s over %d candidates",
                      c(ic_full = "functional IC", ic_auto = "functional IC (automatic search)",
                        bj = "functional Box-Jenkins")[[method]], nrow(candidates)))
    }

    switch(
      method,
      ic_full = select_order_ic_full(Y, candidates, Npc, n_lags),
      ic_auto = select_order_ic_auto(Y, candidates, Npc, n_lags, starts),
      bj = select_order_bj(Y, candidates, Npc, if (is.null(n_lags)) 20L else n_lags, alpha, iid)
    )
  }

  out$p <- as.integer(out$p)
  out$q <- as.integer(out$q)
  structure(c(list(method = method, d = d, n_candidates = nrow(candidates)), out),
            class = "ffarfima_order")
}

ic_fit_term <- function(Y, p, q, Npc, n_lags) {
  fit <- tryCatch(farma_estimation(Y, p = p, q = q, m = NULL, Npc = Npc),
                  error = function(e) NULL)

  if (is.null(fit)) {
    return(list(F_T = Inf, converged = FALSE))
  }

  E <- fit$residuals[, (max(p, q) + 1):ncol(Y), drop = FALSE]
  E <- t(scale(t(E), center = TRUE, scale = FALSE))

  if (any(!is.finite(E)) || ncol(E) < 3) {
    return(list(F_T = Inf, converged = FALSE))
  }

  H <- if (is.null(n_lags)) floor(log(ncol(E))) else n_lags

  list(F_T = temporal_dependence(E, H), converged = TRUE)
}

select_order_ic_full <- function(Y, candidates, Npc, n_lags) {
  tab <- do.call(rbind, lapply(seq_len(nrow(candidates)), function(i) {
    res <- ic_fit_term(Y, candidates$p[i], candidates$q[i], Npc, n_lags)
    data.frame(p = candidates$p[i], q = candidates$q[i], F_T = res$F_T,
               converged = res$converged)
  }))

  ok <- is.finite(tab$F_T) & tab$converged

  if (!any(ok)) {
    stop(sprintf(paste("No candidate FARMA model could be fitted to the %d x %d differenced training curves.",
                       "Check that the data are grid points by time."), nrow(Y), ncol(Y)), call. = FALSE)
  }

  delta_hat <- min(tab$F_T[ok])
  tab$penalty <- log(log(ncol(Y))) * delta_hat * (tab$p + tab$q)
  tab$IC <- tab$F_T + tab$penalty
  best <- which.min(ifelse(is.finite(tab$IC), tab$IC, Inf))

  list(p = tab$p[best], q = tab$q[best], delta_hat = delta_hat, table = tab)
}

select_order_ic_auto <- function(Y, candidates, Npc, n_lags, starts) {
  in_set <- function(pp, qq) {
    paste(pp, qq) %in% paste(candidates$p, candidates$q)
  }

  big <- ic_fit_term(Y, max(candidates$p), max(candidates$q), Npc, n_lags)

  if (!isTRUE(big$converged) || !is.finite(big$F_T)) {
    stop("The largest candidate model could not be fitted, so the noise level is unavailable.",
         call. = FALSE)
  }

  delta_hat <- big$F_T
  big_fit <- tryCatch(farma_estimation(Y, p = max(candidates$p), q = max(candidates$q), m = NULL, Npc = Npc),
                      error = function(e) NULL)
  if (!is.null(big_fit)) {
    dg <- farma_fit_diagnostics(big_fit, Y)
    if (isTRUE(dg$ar_radius >= 1) || isTRUE(dg$ma_radius >= 1) || !isTRUE(dg$residual_ratio <= 50)) {
      warning(sprintf(paste("The largest candidate FARMA(%d, %d) is unstable, so its residual dependence is",
                            "not a usable noise level for \"ic_auto\". Use smaller p and q ranges, or",
                            "order = \"ic_full\"."), max(candidates$p), max(candidates$q)), call. = FALSE)
    }
  }
  penalty_rate <- log(log(ncol(Y))) * delta_hat
  tried <- list()

  get_ic <- function(pp, qq) {
    res <- ic_fit_term(Y, pp, qq, Npc, n_lags)
    ic <- if (isTRUE(res$converged) && is.finite(res$F_T)) {
      res$F_T + penalty_rate * (pp + qq)
    } else {
      Inf
    }
    row <- data.frame(p = pp, q = qq, F_T = res$F_T, converged = res$converged,
                      penalty = penalty_rate * (pp + qq), IC = ic)
    tried[[length(tried) + 1L]] <<- row
    row
  }

  if (is.null(starts)) {
    starts <- data.frame(p = c(2, 1, 0), q = c(2, 0, 1))
  }

  starts <- data.frame(p = as.integer(starts$p), q = as.integer(starts$q))
  starts <- starts[in_set(starts$p, starts$q), , drop = FALSE]

  if (nrow(starts) == 0L) {
    stop("No starting model lies in the candidate set.", call. = FALSE)
  }

  start_tab <- do.call(rbind, lapply(seq_len(nrow(starts)), function(i) {
    get_ic(starts$p[i], starts$q[i])
  }))
  start_tab <- start_tab[is.finite(start_tab$IC) & start_tab$converged, , drop = FALSE]

  if (nrow(start_tab) == 0L) {
    stop("None of the starting models could be fitted.", call. = FALSE)
  }

  current <- start_tab[which.min(start_tab$IC), ]
  path <- data.frame(step = 0, p = current$p, q = current$q, IC = current$IC)
  step <- 0

  repeat {
    moves <- expand.grid(dp = c(-1, 0, 1), dq = c(-1, 0, 1))
    moves <- moves[!(moves$dp == 0 & moves$dq == 0), ]
    nb <- data.frame(p = as.integer(current$p + moves$dp), q = as.integer(current$q + moves$dq))
    nb <- nb[in_set(nb$p, nb$q), , drop = FALSE]

    if (nrow(nb) == 0L) {
      break
    }

    nb_tab <- do.call(rbind, lapply(seq_len(nrow(nb)), function(i) get_ic(nb$p[i], nb$q[i])))
    nb_tab <- nb_tab[is.finite(nb_tab$IC) & nb_tab$converged, , drop = FALSE]

    if (nrow(nb_tab) == 0L || min(nb_tab$IC) >= current$IC) {
      break
    }

    current <- nb_tab[which.min(nb_tab$IC), ]
    step <- step + 1
    path <- rbind(path, data.frame(step = step, p = current$p, q = current$q, IC = current$IC))
  }

  tab <- do.call(rbind, tried)
  tab <- tab[!duplicated(tab[, c("p", "q")]), ]
  tab <- tab[order(tab$p, tab$q), ]
  rownames(tab) <- NULL

  list(p = current$p, q = current$q, delta_hat = delta_hat, table = tab, search_path = path)
}

residual_wn_pvalue <- function(E, lag, iid) {
  wwntests::fport_test(f_data = E, test = "multi-lag", lag = lag, iid = iid,
                       suppress_print_output = TRUE)$p_value
}

select_order_bj <- function(Y, candidates, Npc, n_lags, alpha, iid) {
  candidates <- candidates[order(candidates$p + candidates$q, -candidates$p), ]
  lag_test <- max(1, min(n_lags, ncol(Y) - 1))
  rows <- list()

  for (i in seq_len(nrow(candidates))) {
    pp <- candidates$p[i]
    qq <- candidates$q[i]

    p_value <- tryCatch({
      fit <- farma_estimation(Y, p = pp, q = qq, m = NULL, Npc = Npc)
      E <- fit$residuals[, (max(pp, qq) + 1):ncol(Y), drop = FALSE]
      E <- t(scale(t(E), center = TRUE, scale = FALSE))
      if (any(!is.finite(E))) NA_real_ else residual_wn_pvalue(E, lag_test, iid)
    }, error = function(e) NA_real_)

    passed <- is.finite(p_value) && p_value > alpha
    rows[[i]] <- data.frame(p = pp, q = qq, p_value = p_value, passed = passed)

    if (passed) {
      return(list(p = pp, q = qq, n_lags = lag_test, table = do.call(rbind, rows)))
    }
  }

  warning("No candidate passed the white-noise test; the largest candidate is returned.",
          call. = FALSE)

  list(p = pp, q = qq, n_lags = lag_test, table = do.call(rbind, rows))
}

select_order_cv <- function(train, valid, d, candidates, Npc, trim, h, window, cv_score, verbose,
                            valid_eval = NULL) {
  X <- cbind(train, valid)
  X_eval <- cbind(train, if (is.null(valid_eval)) valid else valid_eval)
  train_end <- ncol(train)
  valid_end <- ncol(X)
  h <- as.integer(h)

  if (train_end <= trim) {
    stop("The training sample is too short after removing `trim` curves.", call. = FALSE)
  }

  if (valid_end - h < train_end) {
    stop(sprintf("Cross-validation with h = %d needs at least %d validation curves; `valid` has %d.",
                 h, h, ncol(valid)), call. = FALSE)
  }

  origins <- seq(train_end, valid_end - 1)
  expected <- vapply(seq_len(h), function(hh) sum(origins + hh <= valid_end), numeric(1))
  rows <- list()
  bar <- progress_start(verbose, nrow(candidates), sprintf(
    "Selecting (p, q) by cross-validation: %d candidates x %d origins", nrow(candidates), length(origins)))

  for (i in seq_len(nrow(candidates))) {
    pp <- candidates$p[i]
    qq <- candidates$q[i]
    sse <- numeric(h)
    sae <- numeric(h)
    n <- numeric(h)
    n_ok <- numeric(h)

    for (t0 in origins) {
      h_here <- min(h, valid_end - t0)
      err <- tryCatch({
        history <- origin_history(X, t0, d, trim, window, train_end)
        if (any(!is.finite(history$Y_history)) || any(!is.finite(history$Y_train)) ||
            ncol(history$Y_train) <= max(pp, qq) + 5) {
          NULL
        } else {
          out <- forecast_original_scale(
            Y_train = history$Y_train, d_hat = d, h = h, method = "farma",
            p_est = pp, q_est = qq, Npc = Npc,
            Y_history_for_integration = history$Y_history
          )
          pred <- out$farma$X_future[, seq_len(h_here), drop = FALSE] +
            matrix(history$mean_curve, nrow = nrow(X), ncol = h_here)
          e <- X_eval[, (t0 + 1):(t0 + h_here), drop = FALSE] - pred
          if (any(!is.finite(e))) NULL else e
        }
      }, error = function(e) NULL)

      if (!is.null(err)) {
        k <- seq_len(h_here)
        sse[k] <- sse[k] + colSums(err^2)
        sae[k] <- sae[k] + colSums(abs(err))
        n[k] <- n[k] + nrow(err)
        n_ok[k] <- n_ok[k] + 1
      }
    }

    rmse <- sqrt(sse / n)
    mae <- sae / n
    complete <- all(n_ok == expected) && all(is.finite(rmse))

    rows[[i]] <- data.frame(
      p = pp,
      q = qq,
      RMSPE = if (complete) mean(rmse) else Inf,
      MAPE = if (complete) mean(mae) else Inf,
      n_origin = n_ok[1]
    )
    progress_step(bar, i)
  }

  progress_end(bar)
  tab <- do.call(rbind, rows)
  score <- if (cv_score == "rmse") tab$RMSPE else tab$MAPE
  tab <- tab[order(score, tab$p + tab$q, tab$p, tab$q), ]
  rownames(tab) <- NULL
  score <- if (cv_score == "rmse") tab$RMSPE else tab$MAPE

  if (!is.finite(score[1])) {
    stop("No candidate produced forecasts at every validation origin.", call. = FALSE)
  }

  list(p = tab$p[1], q = tab$q[1], score = score[1], h = h, window = window,
       cv_score = cv_score, table = tab)
}

#' @export
print.ffarfima_order <- function(x, ...) {
  cat("Selected (p, q) = (", x$p, ", ", x$q, ") at d = ", format(x$d, digits = 4), "\n", sep = "")
  cat("  ", order_label(x), "\n", sep = "")
  invisible(x)
}
