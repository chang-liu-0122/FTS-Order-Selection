select_order_cv_forecast_full <- function(
    X,
    h_max = 10,
    train_ratio = 0.6,
    valid_ratio = 0.2,
    p_max = 5,
    q_max = 5,
    Npc = 1,
    m_const = 0.75,
    R_noise = 0,
    B = 100,
    burn_in = 200,
    fracdiff_trim = 100,
    seed = 1,
    window_type = c("expanding", "fixed"),
    memory_cov_type = "sample",
    memory_cov_band_const = 1
) {

  window_type <- match.arg(window_type)
  memory_cov_type <- match.arg(
    memory_cov_type,
    c("sample", "lrc_RS", "lrc_LRS")
  )

  X <- as.matrix(X)

  if (any(!is.finite(X))) {
    stop("X contains NA, NaN or Inf values.")
  }

  TT <- ncol(X)
  p_grid <- nrow(X)

  h_max <- as.integer(h_max[1])

  train_end <- floor(TT * train_ratio)
  valid_end <- train_end + floor(TT * valid_ratio)
  window_size <- train_end

  if (train_end <= fracdiff_trim) {
    stop("The initial training sample is too short after applying fracdiff_trim.")
  }

  if (valid_end - h_max < train_end) {
    stop("Validation period is too short for h_max.")
  }

  X_initial_train <- X[, 1:train_end, drop = FALSE]

  mean_initial <- rowMeans(X_initial_train)

  X_initial_train_centered <- sweep(
    X_initial_train,
    MARGIN = 1,
    STATS = mean_initial,
    FUN = "-"
  )

  set.seed(seed)

  d_result <- estimate_d_BC_LPWN2(
    X_initial_train_centered,
    m_const = m_const,
    R_noise = R_noise,
    B = B,
    burn_in = burn_in,
    memory_cov_type = memory_cov_type,
    memory_cov_band_const = memory_cov_band_const
  )

  if (isTRUE(d_result$d_BC_LPWN2_rejected)) {
    stop(sprintf(
      paste(
        "Rejected inadmissible BC-LPWN2 estimate:",
        "unrestricted d = %.6f is below zero."
      ),
      d_result$d_BC_LPWN2_unrestricted
    ))
  }

  d_hat <- d_result$d_BC_LPWN2

  d_hat <- as.numeric(d_hat)

  if (length(d_hat) != 1 || !is.finite(d_hat)) {
    stop("Estimated d is not a single finite value.")
  }

  valid_steps <- seq(train_end, valid_end - h_max)

  get_centered_history <- function(t0) {

    X_history <- X[, 1:t0, drop = FALSE]

    if (window_type == "expanding") {
      X_mean_sample <- X_history
    } else {
      raw_start_idx <- max(1, t0 - window_size + 1)
      X_mean_sample <- X[, raw_start_idx:t0, drop = FALSE]
    }

    mean_curve <- rowMeans(X_mean_sample)

    X_history_centered <- sweep(
      X_history,
      MARGIN = 1,
      STATS = mean_curve,
      FUN = "-"
    )

    list(
      X_history_centered = X_history_centered,
      mean_curve = mean_curve
    )
  }

  make_Y_train <- function(X_history_centered) {

    if (d_hat == 0) {
      Y_history <- X_history_centered
    } else {
      Y_history <- frac_difference_history(
        X_history = X_history_centered,
        d_hat = d_hat
      )
    }

    if (d_hat == 0) {
      Y_estimation_history <- Y_history
    } else {
      Y_estimation_history <- Y_history[
        ,
        (fracdiff_trim + 1):ncol(Y_history),
        drop = FALSE
      ]
    }

    if (window_type == "expanding") {
      Y_train <- Y_estimation_history
    } else {
      start_idx <- max(1, ncol(Y_estimation_history) - window_size + 1)

      Y_train <- Y_estimation_history[
        ,
        start_idx:ncol(Y_estimation_history),
        drop = FALSE
      ]
    }

    list(
      Y_history = Y_history,
      Y_train = Y_train
    )
  }

  candidate_grid <- expand.grid(
    p = 0:p_max,
    q = 0:q_max
  ) |>
    dplyr::filter(!(p == 0 & q == 0))

  cv_detail_list <- list()
  error_list <- list()

  counter <- 1
  err_counter <- 1

  for (ii in seq_len(nrow(candidate_grid))) {

    p_now <- candidate_grid$p[ii]
    q_now <- candidate_grid$q[ii]

    for (k in seq_along(valid_steps)) {

      t0 <- valid_steps[k]

      centered_obj <- get_centered_history(t0 = t0)

      train_obj <- tryCatch(
        make_Y_train(centered_obj$X_history_centered),
        error = function(e) e
      )

      if (inherits(train_obj, "error")) {
        error_list[[err_counter]] <- data.frame(
          p = p_now,
          q = q_now,
          t0 = t0,
          stage = "make_Y_train",
          error_message = train_obj$message
        )
        err_counter <- err_counter + 1
        next
      }

      if (any(!is.finite(train_obj$Y_history)) ||
          any(!is.finite(train_obj$Y_train))) {
        error_list[[err_counter]] <- data.frame(
          p = p_now,
          q = q_now,
          t0 = t0,
          stage = "nonfinite_Y",
          error_message = "Y_history or Y_train contains non-finite values."
        )
        err_counter <- err_counter + 1
        next
      }

      if (ncol(train_obj$Y_train) <= max(p_now, q_now) + 5) {
        error_list[[err_counter]] <- data.frame(
          p = p_now,
          q = q_now,
          t0 = t0,
          stage = "short_Y_train",
          error_message = "Y_train is too short for this candidate order."
        )
        err_counter <- err_counter + 1
        next
      }

      out <- tryCatch(
        forecast_original_scale(
          Y_train = train_obj$Y_train,
          d_hat = d_hat,
          h = h_max,
          method = "farma",
          p_est = p_now,
          q_est = q_now,
          Npc = Npc,
          Y_history_for_integration = train_obj$Y_history,
          return_fit = FALSE
        ),
        error = function(e) e
      )

      if (inherits(out, "error")) {
        error_list[[err_counter]] <- data.frame(
          p = p_now,
          q = q_now,
          t0 = t0,
          stage = "forecast_original_scale",
          error_message = out$message
        )
        err_counter <- err_counter + 1
        next
      }

      if (is.null(out$farma) || is.null(out$farma$X_future)) {
        error_list[[err_counter]] <- data.frame(
          p = p_now,
          q = q_now,
          t0 = t0,
          stage = "forecast_original_scale",
          error_message = "forecast_original_scale did not return farma$X_future."
        )
        err_counter <- err_counter + 1
        next
      }

      pred_centered_mat <- as.matrix(out$farma$X_future)

      if (nrow(pred_centered_mat) != p_grid ||
          ncol(pred_centered_mat) != h_max) {
        error_list[[err_counter]] <- data.frame(
          p = p_now,
          q = q_now,
          t0 = t0,
          stage = "forecast_dimension",
          error_message = "Forecast dimension does not match X and h_max."
        )
        err_counter <- err_counter + 1
        next
      }

      pred_future_mat <- pred_centered_mat +
        matrix(
          centered_obj$mean_curve,
          nrow = p_grid,
          ncol = h_max,
          byrow = FALSE
        )

      true_future_mat <- X[, (t0 + 1):(t0 + h_max), drop = FALSE]

      error_mat <- true_future_mat - pred_future_mat

      if (any(!is.finite(error_mat))) {
        error_list[[err_counter]] <- data.frame(
          p = p_now,
          q = q_now,
          t0 = t0,
          stage = "nonfinite_error",
          error_message = "Forecast error contains non-finite values."
        )
        err_counter <- err_counter + 1
        next
      }

      horizon_df <- tibble::tibble(
        p = p_now,
        q = q_now,
        t0 = t0,
        h = seq_len(h_max),
        SSE = apply(error_mat, 2, function(x) sum(x^2, na.rm = TRUE)),
        SAE = apply(error_mat, 2, function(x) sum(abs(x), na.rm = TRUE)),
        n_forecast = apply(error_mat, 2, function(x) sum(is.finite(x)))
      ) |>
        dplyr::mutate(
          RMSE = sqrt(SSE / n_forecast),
          MAE = SAE / n_forecast
        )

      cv_detail_list[[counter]] <- horizon_df
      counter <- counter + 1
    }
  }

  cv_detail <- dplyr::bind_rows(cv_detail_list)
  error_log <- dplyr::bind_rows(error_list)

  if (nrow(cv_detail) == 0) {
    return(
      list(
        model_summary = tibble::tibble(
          T = TT,
          train_end = train_end,
          valid_end = valid_end,
          h_max = h_max,
          window_type = window_type,
          d_hat = d_hat,
          p_hat = NA_integer_,
          q_hat = NA_integer_,
          selection_score = Inf,
          fracdiff_trim = fracdiff_trim,
          m_const = m_const,
          R_noise = R_noise,
          B = B,
          burn_in = burn_in,
          memory_cov_type = memory_cov_type,
          memory_cov_band_const = memory_cov_band_const
        ),
        cv_overall = tibble::tibble(),
        cv_horizon = tibble::tibble(),
        cv_detail = cv_detail,
        cv_matrix = matrix(Inf, nrow = p_max + 1, ncol = q_max + 1),
        error_log = error_log
      )
    )
  }

  cv_horizon <- cv_detail |>
    dplyr::group_by(p, q, h) |>
    dplyr::summarise(
      SSE = sum(SSE, na.rm = TRUE),
      SAE = sum(SAE, na.rm = TRUE),
      n_forecast = sum(n_forecast, na.rm = TRUE),
      n_origin = dplyr::n_distinct(t0),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      RMSE = sqrt(SSE / n_forecast),
      MAE = SAE / n_forecast
    )

  cv_overall <- cv_horizon |>
    dplyr::group_by(p, q) |>
    dplyr::summarise(
      selection_score = mean(RMSE, na.rm = TRUE),
      mean_MAE = mean(MAE, na.rm = TRUE),
      n_valid_horizons = sum(is.finite(RMSE)),
      min_n_origin = min(n_origin, na.rm = TRUE),
      min_RMSE = min(RMSE, na.rm = TRUE),
      max_RMSE = max(RMSE, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      n_origin_expected = length(valid_steps),
      selection_score = dplyr::if_else(
        n_valid_horizons == h_max &
          min_n_origin == n_origin_expected,
        selection_score,
        Inf
      )
    ) |>
    dplyr::right_join(candidate_grid, by = c("p", "q")) |>
    dplyr::mutate(
      selection_score = dplyr::if_else(
        is.na(selection_score),
        Inf,
        selection_score
      ),
      mean_MAE = dplyr::if_else(
        is.na(mean_MAE),
        Inf,
        mean_MAE
      ),
      n_valid_horizons = dplyr::if_else(
        is.na(n_valid_horizons),
        0L,
        as.integer(n_valid_horizons)
      ),
      min_n_origin = dplyr::if_else(
        is.na(min_n_origin),
        0L,
        as.integer(min_n_origin)
      ),
      n_origin_expected = dplyr::if_else(
        is.na(n_origin_expected),
        length(valid_steps),
        as.integer(n_origin_expected)
      )
    ) |>
    dplyr::arrange(selection_score, p + q, p, q)

  best_row <- cv_overall |>
    dplyr::filter(is.finite(selection_score)) |>
    dplyr::slice(1)

  if (nrow(best_row) == 0) {
    p_hat <- NA_integer_
    q_hat <- NA_integer_
    best_score <- Inf
  } else {
    p_hat <- as.integer(best_row$p)
    q_hat <- as.integer(best_row$q)
    best_score <- as.numeric(best_row$selection_score)
  }

  cv_matrix <- matrix(
    Inf,
    nrow = p_max + 1,
    ncol = q_max + 1,
    dimnames = list(
      paste0("p=", 0:p_max),
      paste0("q=", 0:q_max)
    )
  )

  for (i in seq_len(nrow(cv_overall))) {
    p_i <- cv_overall$p[i]
    q_i <- cv_overall$q[i]
    cv_matrix[p_i + 1, q_i + 1] <- cv_overall$selection_score[i]
  }

  list(
    model_summary = tibble::tibble(
      T = TT,
      train_end = train_end,
      valid_end = valid_end,
      h_max = h_max,
      window_type = window_type,
      d_hat = d_hat,
      p_hat = p_hat,
      q_hat = q_hat,
      selection_score = best_score,
      fracdiff_trim = fracdiff_trim,
      m_const = m_const,
      R_noise = R_noise,
      B = B,
      burn_in = burn_in,
      memory_cov_type = memory_cov_type,
      memory_cov_band_const = memory_cov_band_const
    ),
    cv_overall = cv_overall,
    cv_horizon = cv_horizon,
    cv_detail = cv_detail,
    cv_matrix = cv_matrix,
    error_log = error_log
  )
}
