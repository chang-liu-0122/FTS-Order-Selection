rolling_split_conformal_original_ecp_compact <- function(X,
                                                         d_hat,
                                                         d_hat_far = NULL,
                                                         seed = 1,
                                                         p_est = 1,
                                                         q_est = 1,
                                                         h_max = 10,
                                                         train_ratio = 0.6,
                                                         valid_ratio = 0.2,
                                                         alpha_set = c(0.2, 0.05),
                                                         method = c("farma", "far", "both"),
                                                         window_type = c("expanding", "fixed"),
                                                         Npc = 1,
                                                         far_pve = 0.95,
                                                         zeta_grid = seq(0.1, 10, length.out = 500),
                                                         fracdiff_trim = 100)
{
  method <- match.arg(method)
  window_type <- match.arg(window_type)

  X <- as.matrix(X)

  if (any(!is.finite(X))) {
    stop("X contains NA, NaN or Inf values.")
  }

  if (is.null(d_hat_far)) {
    d_hat_far <- d_hat
  }

  d_hat <- as.numeric(d_hat)
  d_hat_far <- as.numeric(d_hat_far)
  h_max <- as.integer(h_max[1])

  T <- ncol(X)
  p_grid <- nrow(X)

  train_end <- floor(T * train_ratio)
  valid_end <- train_end + floor(T * valid_ratio)
  window_size <- train_end

  if (train_end < 1L ||
      valid_end <= train_end ||
      valid_end >= T) {
    stop(
      paste(
        "train_ratio and valid_ratio must produce",
        "non-empty training, validation and test blocks."
      )
    )
  }

  if (h_max < 1L) {
    stop("h_max must be at least 1.")
  }

  methods_used <- if (method == "both") {
    c("farma", "far")
  } else {
    method
  }

  d_method <- list(farma = d_hat, far = d_hat_far)

  # All potential origins are retained.
  #
  # At each origin, h_here determines how many horizons remain
  # inside the validation or test block.
  valid_steps <- seq.int(from = train_end, to = valid_end - 1L)

  test_steps <- seq.int(from = valid_end, to = T - 1L)

  get_centered_history <- function(t0)
  {
    X_history <- X[, 1:t0, drop = FALSE]

    if (window_type == "expanding") {
      X_mean_sample <- X_history
    } else {
      raw_start_idx <- max(1, t0 - window_size + 1)

      X_mean_sample <-
        X[, raw_start_idx:t0, drop = FALSE]
    }

    mean_curve <- rowMeans(X_mean_sample)

    X_history_centered <- sweep(X_history,
                                MARGIN = 1,
                                STATS = mean_curve,
                                FUN = "-")

    list(X_history_centered = X_history_centered, mean_curve = mean_curve)
  }

  make_Y_train <- function(X_history_centered, d_current)
  {
    if (d_current == 0) {
      Y_history <- X_history_centered
    } else {
      Y_history <- frac_difference_history(X_history = X_history_centered, d_hat = d_current)
    }

    if (d_current == 0) {
      Y_estimation_history <- Y_history
    } else {
      if (ncol(Y_history) <= fracdiff_trim) {
        stop(paste(
          "The available history is too short",
          "after applying fracdiff_trim."
        ))
      }

      Y_estimation_history <-
        Y_history[, (fracdiff_trim + 1):ncol(Y_history), drop = FALSE]
    }

    if (window_type == "expanding") {
      Y_train <- Y_estimation_history
    } else {
      start_idx <- max(1, ncol(Y_estimation_history) - window_size + 1)

      Y_train <-
        Y_estimation_history[, start_idx:ncol(Y_estimation_history), drop = FALSE]
    }

    list(Y_history = Y_history, Y_train = Y_train)
  }

  error_list <- list()

  # E_valid[[method]][[h]] stores one column per available
  # validation forecast origin for horizon h.
  E_valid <- setNames(lapply(methods_used, function(m)
  {
    lapply(seq_len(h_max), function(h)
    {
      matrix(numeric(0), nrow = p_grid, ncol = 0)
    })
  }), methods_used)

  # ----------------------------------------------------------
  # Validation forecasts
  # ----------------------------------------------------------

  for (k in seq_along(valid_steps))
  {
    t0 <- valid_steps[k]

    # The final target must remain inside the validation block.
    #
    # Consequently, h = 1 uses every validation origin, whereas
    # longer horizons use progressively fewer origins.
    h_here <- min(h_max, valid_end - t0)

    if (h_here < 1L) {
      next
    }

    centered_obj <- get_centered_history(t0 = t0)

    true_future_mat <-
      X[, (t0 + 1):(t0 + h_here), drop = FALSE]

    for (m in methods_used)
    {
      d_current <- d_method[[m]]

      train_obj <- make_Y_train(X_history_centered =
                                  centered_obj$X_history_centered,
                                d_current = d_current)

      out <- tryCatch(
        forecast_original_scale(
          Y_train = train_obj$Y_train,
          d_hat = d_current,
          h = h_here,
          method = m,
          p_est = p_est,
          q_est = q_est,
          Npc = Npc,
          far_pve = far_pve,
          Y_history_for_integration =
            train_obj$Y_history,
          return_fit = FALSE
        ),
        error = function(e)
          e
      )

      if (inherits(out, "error"))
      {
        error_list[[length(error_list) + 1]] <-
          data.frame(
            seed = seed,
            window_type = window_type,
            method = m,
            stage = "validation_forecast",
            t0 = t0,
            error_message = out$message
          )

        next
      }

      if (is.null(out[[m]])) {
        next
      }

      pred_centered_mat <- out[[m]]$X_future

      pred_future_mat <-
        pred_centered_mat +
        matrix(
          centered_obj$mean_curve,
          nrow = p_grid,
          ncol = h_here,
          byrow = FALSE
        )

      for (h in seq_len(h_here))
      {
        # Validation convention: observed minus forecast.
        e_h <-
          true_future_mat[, h] -
          pred_future_mat[, h]

        E_valid[[m]][[h]] <-
          cbind(E_valid[[m]][[h]], e_h)
      }
    }
  }

  # ----------------------------------------------------------
  # Pointwise ECP calibration
  # ----------------------------------------------------------

  calibration <- list()
  calibration_list <- list()

  for (m in methods_used)
  {
    calibration[[m]] <- list()

    for (alpha in alpha_set)
    {
      alpha_name <- as.character(alpha)

      calibration[[m]][[alpha_name]] <-
        list(gamma = vector("list", h_max),
             zeta = rep(NA_real_, h_max))

      for (h in seq_len(h_max))
      {
        E_h <- E_valid[[m]][[h]]

        if (ncol(E_h) == 0) {
          next
        }

        # calibrate_zeta_ecp retains your pointwise definition:
        #
        # mean(abs(E_h) <= zeta * gamma_h)
        cal <- calibrate_zeta_ecp(E_valid = E_h,
                                  alpha = alpha,
                                  zeta_grid = zeta_grid)

        calibration[[m]][[alpha_name]]$gamma[[h]] <-
          cal$gamma

        calibration[[m]][[alpha_name]]$zeta[h] <-
          cal$zeta
      }

      calibration_list[[length(calibration_list) + 1]] <-
        tibble::tibble(
          seed = seed,
          window_type = window_type,
          method = m,
          alpha = alpha,
          d_used = d_method[[m]],
          gamma = list(calibration[[m]][[alpha_name]]$gamma),
          zeta = list(calibration[[m]][[alpha_name]]$zeta)
        )
    }
  }

  conformal_accum <- setNames(lapply(methods_used, function(m)
  {
    setNames(lapply(as.character(alpha_set), function(alpha_name)
    {
      lapply(seq_len(h_max), function(h)
      {
        list(
          cover_count = 0,
          score_sum = 0,
          n_interval = 0
        )
      })
    }),
    as.character(alpha_set))
  }), methods_used)

  error_curve_list <- list()

  # ----------------------------------------------------------
  # Test forecasts
  # ----------------------------------------------------------

  for (k in seq_along(test_steps))
  {
    t0 <- test_steps[k]

    # The final target must remain inside the observed test block.
    h_here <- min(h_max, T - t0)

    if (h_here < 1L) {
      next
    }

    centered_obj <- get_centered_history(t0 = t0)

    true_future_mat <-
      X[, (t0 + 1):(t0 + h_here), drop = FALSE]

    for (m in methods_used)
    {
      d_current <- d_method[[m]]

      train_obj <- make_Y_train(X_history_centered =
                                  centered_obj$X_history_centered,
                                d_current = d_current)

      out <- tryCatch(
        forecast_original_scale(
          Y_train = train_obj$Y_train,
          d_hat = d_current,
          h = h_here,
          method = m,
          p_est = p_est,
          q_est = q_est,
          Npc = Npc,
          far_pve = far_pve,
          Y_history_for_integration =
            train_obj$Y_history,
          return_fit = FALSE
        ),
        error = function(e)
          e
      )

      if (inherits(out, "error"))
      {
        error_list[[length(error_list) + 1]] <-
          data.frame(
            seed = seed,
            window_type = window_type,
            method = m,
            stage = "testing_forecast",
            t0 = t0,
            error_message = out$message
          )

        next
      }

      if (is.null(out[[m]])) {
        next
      }

      pred_centered_mat <- out[[m]]$X_future

      pred_future_mat <-
        pred_centered_mat +
        matrix(
          centered_obj$mean_curve,
          nrow = p_grid,
          ncol = h_here,
          byrow = FALSE
        )

      # Test convention retained from your original function.
      # The sign does not affect pointwise coverage or interval score
      # because both use abs(error).
      error_future_mat <-
        pred_future_mat -
        true_future_mat

      error_curve_list[[length(error_curve_list) + 1]] <-
        tibble::tibble(
          seed = seed,
          method = m,
          window_type = window_type,
          step = k,
          t0 = t0,
          d_used = d_current,
          error = list(error_future_mat)
        )

      for (alpha in alpha_set)
      {
        alpha_name <- as.character(alpha)

        for (h in seq_len(h_here))
        {
          error_h <-
            as.numeric(error_future_mat[, h])

          gamma_h <-
            as.numeric(calibration[[m]][[alpha_name]]$gamma[[h]])

          zeta_h <-
            calibration[[m]][[alpha_name]]$zeta[h]

          if (length(gamma_h) != length(error_h) ||
              !is.finite(zeta_h)) {
            next
          }

          half_width_h <-
            zeta_h * gamma_h

          ok <-
            is.finite(error_h) &
            is.finite(half_width_h)

          error_h <- error_h[ok]
          half_width_h <- half_width_h[ok]

          if (length(error_h) == 0) {
            next
          }

          # Standard interval score for the symmetric interval
          # [forecast - half_width, forecast + half_width].
          score_h <-
            2 * half_width_h +
            2 / alpha *
            pmax(abs(error_h) - half_width_h, 0)

          accum <-
            conformal_accum[[m]][[alpha_name]][[h]]

          accum$cover_count <-
            accum$cover_count +
            sum(abs(error_h) <= half_width_h)

          accum$score_sum <-
            accum$score_sum +
            sum(score_h)

          accum$n_interval <-
            accum$n_interval +
            length(error_h)

          conformal_accum[[m]][[alpha_name]][[h]] <- accum
        }
      }
    }
  }

  # ----------------------------------------------------------
  # Pointwise test summaries
  # ----------------------------------------------------------

  conformal_summary_list <- list()

  for (m in methods_used)
  {
    for (alpha in alpha_set)
    {
      alpha_name <- as.character(alpha)

      for (h in seq_len(h_max))
      {
        accum <-
          conformal_accum[[m]][[alpha_name]][[h]]

        if (accum$n_interval == 0) {
          next
        }

        ECP <-
          accum$cover_count /
          accum$n_interval

        CPD <-
          abs(ECP - (1 - alpha))

        Score <-
          accum$score_sum /
          accum$n_interval

        conformal_summary_list[[length(conformal_summary_list) + 1]] <-
          tibble::tibble(
            seed = seed,
            method = m,
            window_type = window_type,
            alpha = alpha,
            h = h,
            d_used = d_method[[m]],
            ECP = ECP,
            CPD = CPD,
            Score = Score
          )
      }
    }
  }

  list(
    error_detail =
      dplyr::bind_rows(error_curve_list),

    calibration_detail =
      dplyr::bind_rows(calibration_list),

    conformal_summary_detail =
      dplyr::bind_rows(conformal_summary_list),

    error_log =
      dplyr::bind_rows(error_list)
  )
}



# rolling_split_conformal_original_ecp_compact <- function(X,
#                                                          d_hat,
#                                                          d_hat_far = NULL,
#                                                          seed = 1,
#                                                          p_est = 1,
#                                                          q_est = 1,
#                                                          h_max = 10,
#                                                          train_ratio = 0.6,
#                                                          valid_ratio = 0.2,
#                                                          alpha_set = c(0.2, 0.05),
#                                                          method = c("farma", "far", "both"),
#                                                          window_type = c("expanding", "fixed"),
#                                                          Npc = 1,
#                                                          far_pve = 0.95,
#                                                          zeta_grid = seq(0.1, 10, length.out = 500),
#                                                          fracdiff_trim = 100)
# {
#   method <- match.arg(method)
#   window_type <- match.arg(window_type)
#
#   X <- as.matrix(X)
#
#   if (any(!is.finite(X)))
#   {
#     stop("X contains NA, NaN or Inf values.")
#   }
#
#   if (is.null(d_hat_far)) {
#     d_hat_far <- d_hat
#   }
#
#   d_hat <- as.numeric(d_hat)
#   d_hat_far <- as.numeric(d_hat_far)
#
#   h_max <- as.integer(h_max[1])
#
#   T <- ncol(X)
#   p_grid <- nrow(X)
#
#   train_end <- floor(T * train_ratio)
#   valid_end <- train_end + floor(T * valid_ratio)
#   window_size <- train_end
#
#   methods_used <- if (method == "both") {
#     c("farma", "far")
#   } else {
#     method
#   }
#
#   d_method <- list(farma = d_hat, far = d_hat_far)
#
#   valid_steps <- seq(train_end, valid_end - h_max)
#   test_steps <- seq(valid_end + 1, T - h_max)
#
#   get_centered_history <- function(t0)
#   {
#     X_history <- X[, 1:t0, drop = FALSE]
#
#     if (window_type == "expanding")
#     {
#       X_mean_sample <- X_history
#     } else {
#       raw_start_idx <- max(1, t0 - window_size + 1)
#
#       X_mean_sample <- X[, raw_start_idx:t0, drop = FALSE]
#     }
#
#     mean_curve <- rowMeans(X_mean_sample)
#
#     X_history_centered <- sweep(X_history,
#                                 MARGIN = 1,
#                                 STATS = mean_curve,
#                                 FUN = "-")
#
#     list(X_history_centered = X_history_centered, mean_curve = mean_curve)
#   }
#
#   make_Y_train <- function(X_history_centered, d_current)
#   {
#     if (d_current == 0) {
#       Y_history <- X_history_centered
#     } else {
#       Y_history <- frac_difference_history(X_history = X_history_centered, d_hat = d_current)
#     }
#
#     if (d_current == 0) {
#       Y_estimation_history <- Y_history
#     } else {
#       Y_estimation_history <-
#         Y_history[, (fracdiff_trim + 1):ncol(Y_history), drop = FALSE]
#     }
#
#     if (window_type == "expanding") {
#       Y_train <- Y_estimation_history
#     } else {
#       start_idx <- max(1, ncol(Y_estimation_history) - window_size + 1)
#
#       Y_train <- Y_estimation_history[, start_idx:ncol(Y_estimation_history), drop = FALSE]
#     }
#
#     list(Y_history = Y_history, Y_train = Y_train)
#   }
#
#   error_list <- list()
#
#   E_valid <- setNames(lapply(methods_used, function(m)
#   {
#     lapply(seq_len(h_max), function(h)
#     {
#       matrix(numeric(0), nrow = p_grid, ncol = 0)
#     })
#   }), methods_used)
#
#   for (k in seq_along(valid_steps))
#   {
#     t0 <- valid_steps[k]
#
#     centered_obj <- get_centered_history(t0 = t0)
#
#     true_future_mat <- X[, (t0 + 1):(t0 + h_max), drop = FALSE]
#
#     for (m in methods_used)
#     {
#       d_current <- d_method[[m]]
#
#       train_obj <- make_Y_train(X_history_centered = centered_obj$X_history_centered,
#                                 d_current = d_current)
#
#       out <- tryCatch(
#         forecast_original_scale(
#           Y_train = train_obj$Y_train,
#           d_hat = d_current,
#           h = h_max,
#           method = m,
#           p_est = p_est,
#           q_est = q_est,
#           Npc = Npc,
#           far_pve = far_pve,
#           Y_history_for_integration = train_obj$Y_history,
#           return_fit = FALSE
#         ),
#         error = function(e)
#           e
#       )
#
#       if (inherits(out, "error"))
#       {
#         error_list[[length(error_list) + 1]] <-
#           data.frame(
#             seed = seed,
#             window_type = window_type,
#             method = m,
#             stage = "validation_forecast",
#             t0 = t0,
#             error_message = out$message
#           )
#
#         next
#       }
#
#       if (is.null(out[[m]])) {
#         next
#       }
#
#       pred_centered_mat <- out[[m]]$X_future
#
#       pred_future_mat <- pred_centered_mat +
#         matrix(
#           centered_obj$mean_curve,
#           nrow = p_grid,
#           ncol = h_max,
#           byrow = FALSE
#         )
#
#       for (h in seq_len(h_max))
#       {
#         e_h <- true_future_mat[, h] - pred_future_mat[, h]
#
#         E_valid[[m]][[h]] <-
#           cbind(E_valid[[m]][[h]], e_h)
#       }
#     }
#   }
#
#   calibration <- list()
#   calibration_list <- list()
#
#   for (m in methods_used)
#   {
#     calibration[[m]] <- list()
#
#     for (alpha in alpha_set)
#     {
#       alpha_name <- as.character(alpha)
#
#       calibration[[m]][[alpha_name]] <-
#         list(gamma = vector("list", h_max),
#              zeta = rep(NA_real_, h_max))
#
#       for (h in seq_len(h_max))
#       {
#         E_h <- E_valid[[m]][[h]]
#
#         if (ncol(E_h) == 0) {
#           next
#         }
#
#         cal <- calibrate_zeta_ecp(E_valid = E_h,
#                                   alpha = alpha,
#                                   zeta_grid = zeta_grid)
#
#         calibration[[m]][[alpha_name]]$gamma[[h]] <- cal$gamma
#         calibration[[m]][[alpha_name]]$zeta[h] <- cal$zeta
#       }
#
#       calibration_list[[length(calibration_list) + 1]] <-
#         tibble::tibble(
#           seed = seed,
#           window_type = window_type,
#           method = m,
#           alpha = alpha,
#           d_used = d_method[[m]],
#           gamma = list(calibration[[m]][[alpha_name]]$gamma),
#           zeta = list(calibration[[m]][[alpha_name]]$zeta)
#         )
#     }
#   }
#
#   conformal_accum <- setNames(lapply(methods_used, function(m)
#   {
#     setNames(lapply(as.character(alpha_set), function(alpha_name)
#     {
#       lapply(seq_len(h_max), function(h)
#       {
#         list(
#           cover_count = 0,
#           score_sum = 0,
#           n_interval = 0
#         )
#       })
#     }),
#     as.character(alpha_set))
#   }), methods_used)
#
#   error_curve_list <- list()
#
#   for (k in seq_along(test_steps))
#   {
#     t0 <- test_steps[k]
#
#     centered_obj <- get_centered_history(t0 = t0)
#
#     true_future_mat <- X[, (t0 + 1):(t0 + h_max), drop = FALSE]
#
#     for (m in methods_used)
#     {
#       d_current <- d_method[[m]]
#
#       train_obj <- make_Y_train(X_history_centered = centered_obj$X_history_centered,
#                                 d_current = d_current)
#
#       out <- tryCatch(
#         forecast_original_scale(
#           Y_train = train_obj$Y_train,
#           d_hat = d_current,
#           h = h_max,
#           method = m,
#           p_est = p_est,
#           q_est = q_est,
#           Npc = Npc,
#           far_pve = far_pve,
#           Y_history_for_integration = train_obj$Y_history,
#           return_fit = FALSE
#         ),
#         error = function(e)
#           e
#       )
#
#       if (inherits(out, "error"))
#       {
#         error_list[[length(error_list) + 1]] <-
#           data.frame(
#             seed = seed,
#             window_type = window_type,
#             method = m,
#             stage = "testing_forecast",
#             t0 = t0,
#             error_message = out$message
#           )
#
#         next
#       }
#
#       if (is.null(out[[m]])) {
#         next
#       }
#
#       pred_centered_mat <- out[[m]]$X_future
#
#       pred_future_mat <- pred_centered_mat +
#         matrix(
#           centered_obj$mean_curve,
#           nrow = p_grid,
#           ncol = h_max,
#           byrow = FALSE
#         )
#
#       error_future_mat <- pred_future_mat - true_future_mat
#
#       error_curve_list[[length(error_curve_list) + 1]] <-
#         tibble::tibble(
#           seed = seed,
#           method = m,
#           window_type = window_type,
#           step = k,
#           t0 = t0,
#           d_used = d_current,
#           error = list(error_future_mat)
#         )
#
#       for (alpha in alpha_set)
#       {
#         alpha_name <- as.character(alpha)
#
#         for (h in seq_len(h_max))
#         {
#           error_h <- as.numeric(error_future_mat[, h])
#
#           gamma_h <- as.numeric(calibration[[m]][[alpha_name]]$gamma[[h]])
#
#           zeta_h <- calibration[[m]][[alpha_name]]$zeta[h]
#
#           if (length(gamma_h) != length(error_h) ||
#               !is.finite(zeta_h)) {
#             next
#           }
#
#           half_width_h <- zeta_h * gamma_h
#
#           ok <-
#             is.finite(error_h) &
#             is.finite(half_width_h)
#
#           error_h <- error_h[ok]
#           half_width_h <- half_width_h[ok]
#
#           if (length(error_h) == 0) {
#             next
#           }
#
#           score_h <-
#             2 * half_width_h +
#             2 / alpha *
#             pmax(abs(error_h) - half_width_h, 0)
#
#           accum <- conformal_accum[[m]][[alpha_name]][[h]]
#
#           accum$cover_count <-
#             accum$cover_count +
#             sum(abs(error_h) <= half_width_h)
#
#           accum$score_sum <-
#             accum$score_sum +
#             sum(score_h)
#
#           accum$n_interval <-
#             accum$n_interval +
#             length(error_h)
#
#           conformal_accum[[m]][[alpha_name]][[h]] <- accum
#         }
#       }
#     }
#   }
#
#   conformal_summary_list <- list()
#
#   for (m in methods_used)
#   {
#     for (alpha in alpha_set)
#     {
#       alpha_name <- as.character(alpha)
#
#       for (h in seq_len(h_max))
#       {
#         accum <- conformal_accum[[m]][[alpha_name]][[h]]
#
#         if (accum$n_interval == 0) {
#           next
#         }
#
#         ECP <- accum$cover_count / accum$n_interval
#         CPD <- abs(ECP - (1 - alpha))
#         Score <- accum$score_sum / accum$n_interval
#
#         conformal_summary_list[[length(conformal_summary_list) + 1]] <-
#           tibble::tibble(
#             seed = seed,
#             method = m,
#             window_type = window_type,
#             alpha = alpha,
#             h = h,
#             d_used = d_method[[m]],
#             ECP = ECP,
#             CPD = CPD,
#             Score = Score
#           )
#       }
#     }
#   }
#
#   list(
#     error_detail =
#       dplyr::bind_rows(error_curve_list),
#
#     calibration_detail =
#       dplyr::bind_rows(calibration_list),
#
#     conformal_summary_detail =
#       dplyr::bind_rows(conformal_summary_list),
#
#     error_log =
#       dplyr::bind_rows(error_list)
#   )
# }




parallel_sim_forecasting_FARMA_compare_compact <- function(n_rep = 100,
                                                           grid,
                                                           h_max = 10,
                                                           train_ratio = 0.6,
                                                           valid_ratio = 0.2,
                                                           alpha_set = c(0.2, 0.05),
                                                           window_types = c("expanding", "fixed"),
                                                           Npc = 1,
                                                           no_grid = 101,
                                                           far_pve = 0.95,
                                                           zeta_grid = seq(0.1, 10, length.out = 500),
                                                           fracdiff_trim = 100)
{
  if (!"scenario_id" %in% names(grid))
  {
    grid <- grid %>%
      dplyr::mutate(scenario_id = dplyr::row_number(),
                    .before = 1)
  }

  progressr::with_progress({
    seeds <- seq_len(n_rep)

    p_prog <- progressr::progressor(steps =
                                      length(seeds) *
                                      nrow(grid) *
                                      length(window_types))

    results <- future.apply::future_lapply(seeds, function(s)
    {
      error_curve_list <- list()
      calibration_list <- list()
      conformal_summary_list <- list()
      error_list <- list()

      for (i in seq_len(nrow(grid)))
      {
        scenario_i <- grid$scenario_id[i]
        n_i <- grid$T[i]
        d_i <- grid$d[i]
        ar_i <- grid$ar_strength[i]
        ma_i <- grid$ma_strength[i]
        p_i <- grid$p[i]
        q_i <- grid$q[i]

        sim_seed <- s * 100000 + i

        RNGkind(
          kind = "Mersenne-Twister",
          normal.kind = "Inversion",
          sample.kind = "Rejection"
        )

        sim <- tryCatch(
          sim_FARMA(
            n = n_i,
            d = d_i,
            no_grid = no_grid,
            seed_number = sim_seed,
            p = p_i,
            q = q_i,
            ar_strength = ar_i,
            ma_strength = ma_i,
            lag_strength = "equal"
          ),
          error = function(e)
            e
        )

        if (inherits(sim, "error"))
        {
          error_list[[length(error_list) + 1]] <-
            data.frame(
              scenario_id = scenario_i,
              n = n_i,
              d_true = d_i,
              ar_strength = ar_i,
              ma_strength = ma_i,
              p_true = p_i,
              q_true = q_i,
              seed = s,
              sim_seed = sim_seed,
              window_type = NA_character_,
              stage = "simulation",
              t0 = NA_integer_,
              error_message = sim$message
            )

          for (wtype in window_types)
          {
            p_prog()
          }

          next
        }

        X <- sim$X

        for (wtype in window_types)
        {
          out <- tryCatch(
            rolling_split_conformal_original_ecp_compact(
              X = X,
              d_hat = d_i,
              d_hat_far = d_i,
              seed = s,
              p_est = p_i,
              q_est = q_i,
              h_max = h_max,
              train_ratio = train_ratio,
              valid_ratio = valid_ratio,
              alpha_set = alpha_set,
              method = "both",
              window_type = wtype,
              Npc = Npc,
              far_pve = far_pve,
              zeta_grid = zeta_grid,
              fracdiff_trim = fracdiff_trim
            ),
            error = function(e)
              e
          )

          if (inherits(out, "error"))
          {
            error_list[[length(error_list) + 1]] <-
              data.frame(
                scenario_id = scenario_i,
                n = n_i,
                d_true = d_i,
                ar_strength = ar_i,
                ma_strength = ma_i,
                p_true = p_i,
                q_true = q_i,
                seed = s,
                sim_seed = sim_seed,
                window_type = wtype,
                stage = "rolling_forecast",
                t0 = NA_integer_,
                error_message = out$message
              )

            p_prog()

            next
          }

          metadata <- list(
            scenario_id = scenario_i,
            n = n_i,
            d_true = d_i,
            ar_strength = ar_i,
            ma_strength = ma_i,
            p_true = p_i,
            q_true = q_i,
            sim_seed = sim_seed
          )

          if (nrow(out$error_detail) > 0)
          {
            error_curve_list[[length(error_curve_list) + 1]] <-
              out$error_detail %>%
              dplyr::mutate(!!!metadata, .before = 1)
          }

          if (nrow(out$calibration_detail) > 0)
          {
            calibration_list[[length(calibration_list) + 1]] <-
              out$calibration_detail %>%
              dplyr::mutate(!!!metadata, .before = 1)
          }

          if (nrow(out$conformal_summary_detail) > 0)
          {
            conformal_summary_list[[length(conformal_summary_list) + 1]] <-
              out$conformal_summary_detail %>%
              dplyr::mutate(!!!metadata, .before = 1)
          }

          if (nrow(out$error_log) > 0)
          {
            error_list[[length(error_list) + 1]] <-
              out$error_log %>%
              dplyr::mutate(!!!metadata, .before = 1)
          }

          rm(out)
          gc()

          p_prog()
        }

        rm(X, sim)
        gc()
      }

      list(
        error_detail =
          dplyr::bind_rows(error_curve_list),

        calibration_detail =
          dplyr::bind_rows(calibration_list),

        conformal_summary_detail =
          dplyr::bind_rows(conformal_summary_list),

        error_log =
          dplyr::bind_rows(error_list)
      )
    }, future.seed = TRUE)

    list(
      grid = grid,

      error_detail =
        dplyr::bind_rows(lapply(results, function(x)
          x$error_detail)),

      calibration_detail =
        dplyr::bind_rows(lapply(results, function(x)
          x$calibration_detail)),

      conformal_summary_detail =
        dplyr::bind_rows(lapply(results, function(x)
          x$conformal_summary_detail)),

      error_log =
        dplyr::bind_rows(lapply(results, function(x)
          x$error_log))
    )
  })
}


extract_required_metrics_from_compact_order_method <- function(compact_out, rmse_ratio_threshold = 10)
{
  error_detail <- compact_out$error_detail
  conformal_summary_detail <- compact_out$conformal_summary_detail
  estimate_detail <- compact_out$estimate_detail

  if (nrow(error_detail) == 0) {
    stop("error_detail is empty.")
  }

  if (nrow(conformal_summary_detail) == 0) {
    stop("conformal_summary_detail is empty.")
  }

  if (!"order_selection_method" %in% names(error_detail)) {
    error_detail$order_selection_method <- "adaptive_ic"
  }

  if (!"order_selection_method" %in% names(conformal_summary_detail)) {
    conformal_summary_detail$order_selection_method <- "adaptive_ic"
  }

  if (!is.null(estimate_detail) && nrow(estimate_detail) > 0) {
    if (!"order_selection_method" %in% names(estimate_detail)) {
      estimate_detail$order_selection_method <- "adaptive_ic"
    }
  }

  point_by_rep_raw <- error_detail %>%
    dplyr::mutate(horizon_stats =
                    purrr::map(error, function(error_mat)
                    {
                      tibble::tibble(
                        h = seq_len(ncol(error_mat)),
                        SSE = apply(error_mat, 2, function(x)
                        {
                          x <- x[is.finite(x)]
                          sum(x^2)
                        }),
                        SAE = apply(error_mat, 2, function(x)
                        {
                          x <- x[is.finite(x)]
                          sum(abs(x))
                        }),
                        n_forecast = apply(error_mat, 2, function(x)
                        {
                          sum(is.finite(x))
                        })
                      )
                    })) %>%
    dplyr::select(
      scenario_id,
      n,
      d_true,
      ar_strength,
      ma_strength,
      p_true,
      q_true,
      seed,
      sim_seed,
      order_selection_method,
      window_type,
      method,
      horizon_stats
    ) %>%
    tidyr::unnest(horizon_stats) %>%
    dplyr::group_by(
      scenario_id,
      n,
      d_true,
      ar_strength,
      ma_strength,
      p_true,
      q_true,
      seed,
      sim_seed,
      order_selection_method,
      window_type,
      method,
      h
    ) %>%
    dplyr::summarise(
      SSE = sum(SSE, na.rm = TRUE),
      SAE = sum(SAE, na.rm = TRUE),
      n_forecast = sum(n_forecast, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(RMSE = sqrt(SSE / n_forecast),
                  MAE = SAE / n_forecast)

  far_benchmark <- point_by_rep_raw %>%
    dplyr::filter(method == "far") %>%
    dplyr::group_by(
      scenario_id,
      n,
      d_true,
      ar_strength,
      ma_strength,
      p_true,
      q_true,
      seed,
      sim_seed,
      window_type,
      h
    ) %>%
    dplyr::summarise(
      RMSE_far = dplyr::first(RMSE),
      MAE_far = dplyr::first(MAE),
      .groups = "drop"
    )

  farma_by_rep <- point_by_rep_raw %>%
    dplyr::filter(method == "farma") %>%
    dplyr::select(
      scenario_id,
      n,
      d_true,
      ar_strength,
      ma_strength,
      p_true,
      q_true,
      seed,
      sim_seed,
      order_selection_method,
      window_type,
      h,
      RMSE_farma = RMSE,
      MAE_farma = MAE
    )

  ratio_by_rep <- farma_by_rep %>%
    dplyr::left_join(
      far_benchmark,
      by = c(
        "scenario_id",
        "n",
        "d_true",
        "ar_strength",
        "ma_strength",
        "p_true",
        "q_true",
        "seed",
        "sim_seed",
        "window_type",
        "h"
      )
    ) %>%
    dplyr::mutate(RMSE_ratio = RMSE_farma / RMSE_far,
                  MAE_ratio = MAE_farma / MAE_far) %>%
    dplyr::arrange(scenario_id,
                   seed,
                   sim_seed,
                   order_selection_method,
                   window_type,
                   h)

  removed_ratio_rows <- ratio_by_rep %>%
    dplyr::filter(!is.finite(RMSE_ratio) |
                    RMSE_ratio > rmse_ratio_threshold)

  ratio_by_rep_clean <- ratio_by_rep %>%
    dplyr::filter(is.finite(RMSE_ratio), RMSE_ratio <= rmse_ratio_threshold)

  ratio_summary <- ratio_by_rep_clean %>%
    dplyr::group_by(
      scenario_id,
      n,
      d_true,
      ar_strength,
      ma_strength,
      p_true,
      q_true,
      order_selection_method,
      window_type,
      h
    ) %>%
    dplyr::summarise(
      RMSE_far = mean(RMSE_far, na.rm = TRUE),
      RMSE_farma = mean(RMSE_farma, na.rm = TRUE),
      MAE_far = mean(MAE_far, na.rm = TRUE),
      MAE_farma = mean(MAE_farma, na.rm = TRUE),
      RMSE_ratio = RMSE_farma / RMSE_far,
      MAE_ratio = MAE_farma / MAE_far,
      n_rep_used = dplyr::n(),
      .groups = "drop"
    ) %>%
    dplyr::arrange(scenario_id, order_selection_method, window_type, h)

  bad_keys <- removed_ratio_rows %>%
    dplyr::distinct(scenario_id,
                    seed,
                    sim_seed,
                    order_selection_method,
                    window_type,
                    h)

  conformal_summary_detail_clean <- conformal_summary_detail

  if (nrow(bad_keys) > 0 &&
      "h" %in% names(conformal_summary_detail_clean)) {
    conformal_summary_detail_clean <- conformal_summary_detail_clean %>%
      dplyr::anti_join(
        bad_keys,
        by = c(
          "scenario_id",
          "seed",
          "sim_seed",
          "order_selection_method",
          "window_type",
          "h"
        )
      )
  }

  interval_summary <- conformal_summary_detail_clean %>%
    dplyr::group_by(
      scenario_id,
      n,
      d_true,
      ar_strength,
      ma_strength,
      p_true,
      q_true,
      order_selection_method,
      window_type,
      method,
      alpha,
      h
    ) %>%
    dplyr::summarise(
      ECP = mean(ECP, na.rm = TRUE),
      CPD = mean(CPD, na.rm = TRUE),
      Score = mean(Score, na.rm = TRUE),
      n_rep_used = dplyr::n(),
      .groups = "drop"
    ) %>%
    dplyr::arrange(scenario_id,
                   order_selection_method,
                   window_type,
                   alpha,
                   h,
                   method)

  if (!is.null(estimate_detail) && nrow(estimate_detail) > 0) {
    d_summary <- estimate_detail %>%
      dplyr::group_by(
        scenario_id,
        n,
        d_true,
        ar_strength,
        ma_strength,
        p_true,
        q_true,
        order_selection_method
      ) %>%
      dplyr::summarise(
        mean_d_hat = mean(d_hat, na.rm = TRUE),
        sd_d_hat = sd(d_hat, na.rm = TRUE),
        mean_d_bias = mean(d_hat - d_true, na.rm = TRUE),
        mean_abs_d_bias = mean(abs(d_hat - d_true), na.rm = TRUE),
        RMSE_d = sqrt(mean((d_hat - d_true)^2, na.rm = TRUE)),
        n_rep = dplyr::n(),
        .groups = "drop"
      ) %>%
      dplyr::arrange(scenario_id, order_selection_method)

    order_count <- estimate_detail %>%
      dplyr::group_by(
        scenario_id,
        n,
        d_true,
        ar_strength,
        ma_strength,
        p_true,
        q_true,
        order_selection_method,
        p_est,
        q_est
      ) %>%
      dplyr::summarise(count = dplyr::n(), .groups = "drop")

    order_total <- estimate_detail %>%
      dplyr::group_by(
        scenario_id,
        n,
        d_true,
        ar_strength,
        ma_strength,
        p_true,
        q_true,
        order_selection_method
      ) %>%
      dplyr::summarise(total = dplyr::n(), .groups = "drop")

    order_summary <- order_count %>%
      dplyr::left_join(
        order_total,
        by = c(
          "scenario_id",
          "n",
          "d_true",
          "ar_strength",
          "ma_strength",
          "p_true",
          "q_true",
          "order_selection_method"
        )
      ) %>%
      dplyr::mutate(proportion = count / total) %>%
      dplyr::arrange(scenario_id,
                     order_selection_method,
                     dplyr::desc(proportion),
                     p_est,
                     q_est)

    order_accuracy <- estimate_detail %>%
      dplyr::group_by(
        scenario_id,
        n,
        d_true,
        ar_strength,
        ma_strength,
        p_true,
        q_true,
        order_selection_method
      ) %>%
      dplyr::summarise(
        order_accuracy = mean(p_est == p_true &
                                q_est == q_true, na.rm = TRUE),
        p_accuracy = mean(p_est == p_true, na.rm = TRUE),
        q_accuracy = mean(q_est == q_true, na.rm = TRUE),
        mean_p_est = mean(p_est, na.rm = TRUE),
        mean_q_est = mean(q_est, na.rm = TRUE),
        n_rep = dplyr::n(),
        .groups = "drop"
      ) %>%
      dplyr::arrange(scenario_id, order_selection_method)

  } else {
    d_summary <- tibble::tibble()
    order_summary <- tibble::tibble()
    order_accuracy <- tibble::tibble()

  }

  list(
    ratio_by_rep = ratio_by_rep,
    removed_ratio_rows = removed_ratio_rows,
    ratio_by_rep_clean = ratio_by_rep_clean,
    ratio_summary = ratio_summary,
    interval_summary = interval_summary,
    d_summary = d_summary,
    order_summary = order_summary,
    order_accuracy = order_accuracy
  )
}
