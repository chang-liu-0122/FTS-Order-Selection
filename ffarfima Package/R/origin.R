origin_history <- function(X, t0, d, trim, window, window_size) {
  X_history <- X[, 1:t0, drop = FALSE]

  X_mean_sample <- if (window == "expanding") {
    X_history
  } else {
    X[, max(1, t0 - window_size + 1):t0, drop = FALSE]
  }

  mean_curve <- rowMeans(X_mean_sample)
  X_history_centered <- sweep(X_history, 1, mean_curve, "-")

  if (d == 0) {
    Y_history <- X_history_centered
    Y_estimation <- Y_history
  } else {
    Y_history <- frac_difference_history(X_history_centered, d)

    if (ncol(Y_history) <= trim) {
      stop("The history is too short after removing `trim` fractionally differenced curves.",
           call. = FALSE)
    }

    Y_estimation <- Y_history[, (trim + 1):ncol(Y_history), drop = FALSE]
  }

  Y_train <- if (window == "expanding") {
    Y_estimation
  } else {
    Y_estimation[, max(1, ncol(Y_estimation) - window_size + 1):ncol(Y_estimation), drop = FALSE]
  }

  list(Y_history = Y_history, Y_train = Y_train, mean_curve = mean_curve)
}

forecast_at_origin <- function(X, t0, h, d, p, q, engine, Npc, far_pve, trim,
                               window, window_size, residual_ratio_limit = 50) {
  history <- origin_history(X, t0, d, trim, window, window_size)

  out <- forecast_original_scale(
    Y_train = history$Y_train,
    d_hat = d,
    h = h,
    method = engine,
    p_est = p,
    q_est = q,
    Npc = Npc,
    far_pve = far_pve,
    Y_history_for_integration = history$Y_history,
    farma_residual_ratio_limit = residual_ratio_limit
  )

  out[[engine]]$X_future + matrix(history$mean_curve, nrow = nrow(X), ncol = h)
}
