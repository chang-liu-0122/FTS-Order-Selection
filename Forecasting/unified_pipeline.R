# Unified pipeline for one series: estimate d, select (p, q) by IC and CV, then forecast
unified_pipeline <- function(X,
                             h_max = 10,
                             train_ratio = 0.6,
                             valid_ratio = 0.2,
                             alpha_set = c(0.2, 0.05),
                             window_types = c("expanding", "fixed"),
                             p_max = 5,
                             q_max = 5,
                             Npc = 1,
                             far_pve = 0.95,
                             zeta_grid = seq(0.1, 50, length.out = 2496),
                             m_const = 0.75,
                             R_noise = 0,
                             B = 100,
                             burn_in = 200,
                             fracdiff_trim = 100,
                             seed = 1,
                             far_benchmark_d = NULL,
                             order_selection_method = c("adaptive_ic", "cv", "both"),
                             cv_order_window_type = c("expanding", "fixed"),
                             memory_cov_type = "sample",
                             memory_cov_band_const = 1,
                             X_eval = NULL) {
  order_selection_method <- match.arg(order_selection_method)
  cv_order_window_type <- match.arg(cv_order_window_type, several.ok = TRUE)
  order_selection_method <- match.arg(order_selection_method)
  memory_cov_type <- match.arg(memory_cov_type, c("sample", "lrc_RS", "lrc_LRS"))

  order_methods <- if (order_selection_method == "both") {
    c("adaptive_ic", "cv")
  } else {
    order_selection_method
  }

  window_types <- match.arg(window_types, c("expanding", "fixed"), several.ok = TRUE)

  X <- as.matrix(X)

  if (any(!is.finite(X))) {
    stop("X contains NA, NaN or Inf values.")
  }

  X_eval <- if (is.null(X_eval))
    X
  else
    as.matrix(X_eval)

  if (any(!is.finite(X_eval))) {
    stop("X_eval contains NA, NaN or Inf values.")
  }

  if (!identical(dim(X_eval), dim(X))) {
    stop("X_eval must have the same dimensions as X.")
  }

  TT <- ncol(X)
  train_end <- floor(TT * train_ratio)
  valid_end <- train_end + floor(TT * valid_ratio)

  if (train_end <= fracdiff_trim) {
    stop("The initial training sample is too short after applying fracdiff_trim.")
  }

  # Label results with the order-selection method and the selected (p, q)
  add_order_tag <- function(df,
                            os_method,
                            os_window,
                            os_score,
                            p_est,
                            q_est) {
    if (is.null(df) || nrow(df) == 0) {
      return(tibble::tibble())
    }

    df |>
      dplyr::mutate(
        order_selection_method = os_method,
        order_selection_window_type = os_window,
        order_selection_score = os_score,
        p_est = p_est,
        q_est = q_est,
        .before = 1
      )
  }

  # Rolling forecasts in each window type for a selected d and (p, q)
  run_rolling_given_order <- function(d_hat,
                                      d_hat_far,
                                      p_est,
                                      q_est,
                                      os_method,
                                      os_window,
                                      os_score) {
    error_curve_list <- list()
    calibration_list <- list()
    conformal_summary_list <- list()
    error_list <- list()

    for (wtype in window_types) {
      out_w <- tryCatch(
        rolling_forecast(
          X = X,
          d_hat = d_hat,
          d_hat_far = d_hat_far,
          seed = seed,
          p_est = p_est,
          q_est = q_est,
          h_max = h_max,
          train_ratio = train_ratio,
          valid_ratio = valid_ratio,
          alpha_set = alpha_set,
          method = "both",
          window_type = wtype,
          Npc = Npc,
          far_pve = far_pve,
          zeta_grid = zeta_grid,
          fracdiff_trim = fracdiff_trim,
          X_eval = X_eval
        ),
        error = function(e)
          e
      )

      if (inherits(out_w, "error")) {
        error_list[[length(error_list) + 1]] <-
          tibble::tibble(
            order_selection_method = os_method,
            order_selection_window_type = os_window,
            order_selection_score = os_score,
            p_est = p_est,
            q_est = q_est,
            window_type = wtype,
            stage = "rolling_forecast",
            t0 = NA_integer_,
            error_message = out_w$message
          )
        next
      }

      if (nrow(out_w$error_detail) > 0) {
        error_curve_list[[length(error_curve_list) + 1]] <-
          add_order_tag(out_w$error_detail,
                        os_method,
                        os_window,
                        os_score,
                        p_est,
                        q_est)
      }

      if (nrow(out_w$calibration_detail) > 0) {
        calibration_list[[length(calibration_list) + 1]] <-
          add_order_tag(out_w$calibration_detail,
                        os_method,
                        os_window,
                        os_score,
                        p_est,
                        q_est)
      }

      if (nrow(out_w$conformal_summary_detail) > 0) {
        conformal_summary_list[[length(conformal_summary_list) + 1]] <-
          add_order_tag(
            out_w$conformal_summary_detail,
            os_method,
            os_window,
            os_score,
            p_est,
            q_est
          )
      }

      if (nrow(out_w$error_log) > 0) {
        error_list[[length(error_list) + 1]] <-
          add_order_tag(out_w$error_log,
                        os_method,
                        os_window,
                        os_score,
                        p_est,
                        q_est)
      }
    }

    list(
      error_detail = dplyr::bind_rows(error_curve_list),
      calibration_detail = dplyr::bind_rows(calibration_list),
      conformal_summary_detail = dplyr::bind_rows(conformal_summary_list),
      error_log = dplyr::bind_rows(error_list)
    )
  }

  # BC-LPWN estimate of d on the centred initial training curves
  estimate_d_once <- function() {
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
      X = X_initial_train_centered,
      m_const = m_const,
      R_noise = R_noise,
      B = B,
      burn_in = burn_in,
      memory_cov_type = memory_cov_type,
      memory_cov_band_const = memory_cov_band_const
    )

    d_hat <- as.numeric(d_result$d_BC_LPWN2)

    if (isTRUE(d_result$d_BC_LPWN2_rejected)) {
      stop(structure(
        class = c("lpwn_inadmissible_d", "error", "condition"),
        list(
          message = sprintf(
            paste(
              "Rejected inadmissible BC-LPWN2 estimate:",
              "unrestricted d = %.6f is below zero."
            ),
            d_result$d_BC_LPWN2_unrestricted
          ),
          call = NULL
        )
      ))
    }

    if (length(d_hat) != 1L || !is.finite(d_hat)) {
      stop("Estimated d is not a single finite value.")
    }

    list(
      d_hat = d_hat,
      X_initial_train_centered = X_initial_train_centered,
      memory_cov_type = memory_cov_type,
      memory_cov_band_const = memory_cov_band_const
    )
  }

  # Select (p, q) by the adaptive functional IC, then forecast
  run_adaptive_ic <- function() {
    d_obj <- estimate_d_once()
    d_hat <- d_obj$d_hat
    X_initial_train_centered <- d_obj$X_initial_train_centered

    if (is.null(far_benchmark_d)) {
      d_hat_far <- d_hat
    } else {
      d_hat_far <- as.numeric(far_benchmark_d)
    }

    if (length(d_hat_far) != 1 || !is.finite(d_hat_far)) {
      stop("far_benchmark_d must be NULL or a single finite value.")
    }

    if (d_hat == 0) {
      X_diff_initial_train <- X_initial_train_centered
    } else {
      X_diff_initial_train <- FracDiff(t(X_initial_train_centered), d_hat) |>
        t()
    }

    if (any(!is.finite(X_diff_initial_train))) {
      stop("Fractionally differenced training curves contain NA or Inf.")
    }

    if (d_hat == 0) {
      X_diff_for_order_selection <- X_diff_initial_train
    } else {
      X_diff_for_order_selection <-
        X_diff_initial_train[, (fracdiff_trim + 1):ncol(X_diff_initial_train), drop = FALSE]
    }

    order_out <- select_order_far_full(
      X = X_diff_for_order_selection,
      p_max = p_max,
      q_max = q_max,
      Npc = Npc
    )

    p_est <- as.integer(order_out$p_hat)
    q_est <- as.integer(order_out$q_hat)

    if (length(p_est) != 1 ||
        length(q_est) != 1 ||
        !is.finite(p_est) ||
        !is.finite(q_est)) {
      stop("Adaptive IC order selection did not return valid p_est and q_est.")
    }

    rolling_out <- run_rolling_given_order(
      d_hat = d_hat,
      d_hat_far = d_hat_far,
      p_est = p_est,
      q_est = q_est,
      os_method = "adaptive_ic",
      os_window = NA_character_,
      os_score = NA_real_
    )

    list(
      model_summary = tibble::tibble(
        order_selection_method = "adaptive_ic",
        order_selection_window_type = NA_character_,
        order_selection_score = NA_real_,
        T = TT,
        train_end = train_end,
        valid_end = valid_end,
        test_start = valid_end + 1,
        centering = "expanding_or_fixed_training_window_mean",
        d_hat = d_hat,
        d_hat_far = d_hat_far,
        p_est = p_est,
        q_est = q_est,
        fracdiff_trim = fracdiff_trim
      ),
      order_selection_detail = order_out,
      error_detail = rolling_out$error_detail,
      calibration_detail = rolling_out$calibration_detail,
      conformal_summary_detail = rolling_out$conformal_summary_detail,
      error_log = rolling_out$error_log
    )
  }

  # Select d and (p, q) by CV for one window type, then forecast
  run_cv_one_window <- function(cv_wtype) {
    cv_order <- select_order_cv_forecast_full(
      X = X,
      X_eval = X_eval,
      h_max = h_max,
      train_ratio = train_ratio,
      valid_ratio = valid_ratio,
      p_max = p_max,
      q_max = q_max,
      Npc = Npc,
      m_const = m_const,
      R_noise = R_noise,
      B = B,
      burn_in = burn_in,
      fracdiff_trim = fracdiff_trim,
      seed = seed,
      window_type = cv_wtype,
      memory_cov_type = memory_cov_type,
      memory_cov_band_const = memory_cov_band_const
    )

    d_hat <- as.numeric(cv_order$model_summary$d_hat[1])
    p_est <- as.integer(cv_order$model_summary$p_hat[1])
    q_est <- as.integer(cv_order$model_summary$q_hat[1])
    order_selection_score <- as.numeric(cv_order$model_summary$selection_score[1])

    if (length(d_hat) != 1 || !is.finite(d_hat)) {
      stop("CV order selection did not return a valid d_hat.")
    }

    if (length(p_est) != 1 ||
        length(q_est) != 1 ||
        !is.finite(p_est) ||
        !is.finite(q_est)) {
      stop("CV order selection did not return valid p_est and q_est.")
    }

    if (is.null(far_benchmark_d)) {
      d_hat_far <- d_hat
    } else {
      d_hat_far <- as.numeric(far_benchmark_d)
    }

    if (length(d_hat_far) != 1 || !is.finite(d_hat_far)) {
      stop("far_benchmark_d must be NULL or a single finite value.")
    }

    os_method <- paste0("cv_", cv_wtype)

    rolling_out <- run_rolling_given_order(
      d_hat = d_hat,
      d_hat_far = d_hat_far,
      p_est = p_est,
      q_est = q_est,
      os_method = os_method,
      os_window = cv_wtype,
      os_score = order_selection_score
    )

    list(
      model_summary = tibble::tibble(
        order_selection_method = os_method,
        order_selection_window_type = cv_wtype,
        order_selection_score = order_selection_score,
        T = TT,
        train_end = train_end,
        valid_end = valid_end,
        test_start = valid_end + 1,
        centering = "expanding_or_fixed_training_window_mean",
        d_hat = d_hat,
        d_hat_far = d_hat_far,
        p_est = p_est,
        q_est = q_est,
        fracdiff_trim = fracdiff_trim
      ),
      order_selection_detail = cv_order,
      error_detail = rolling_out$error_detail,
      calibration_detail = rolling_out$calibration_detail,
      conformal_summary_detail = rolling_out$conformal_summary_detail,
      error_log = rolling_out$error_log
    )
  }

  out_list <- list()

  if ("adaptive_ic" %in% order_methods) {
    out_list[["adaptive_ic"]] <- run_adaptive_ic()
  }

  if ("cv" %in% order_methods) {
    for (cv_wtype in cv_order_window_type) {
      out_list[[paste0("cv_", cv_wtype)]] <- run_cv_one_window(cv_wtype)
    }
  }

  list(
    model_summary = dplyr::bind_rows(lapply(out_list, function(z)
      z$model_summary)),
    order_selection_detail = lapply(out_list, function(z)
      z$order_selection_detail),
    error_detail = dplyr::bind_rows(lapply(out_list, function(z)
      z$error_detail)),
    calibration_detail = dplyr::bind_rows(lapply(out_list, function(z)
      z$calibration_detail)),
    conformal_summary_detail = dplyr::bind_rows(lapply(out_list, function(z)
      z$conformal_summary_detail)),
    error_log = dplyr::bind_rows(lapply(out_list, function(z)
      z$error_log))
  )
}





# Unified simulation: simulate each scenario and run the unified pipeline in parallel
simulate_unified <- function(n_rep = 100,
                             grid,
                             h_max = 10,
                             train_ratio = 0.6,
                             valid_ratio = 0.2,
                             alpha_set = c(0.2, 0.05),
                             window_types = c("expanding", "fixed"),
                             p_max = 5,
                             q_max = 5,
                             Npc = 1,
                             no_grid = 101,
                             far_pve = 0.95,
                             zeta_grid = seq(0.1, 50, length.out = 2496),
                             m_const = 0.75,
                             R_noise = 0,
                             B = 100,
                             burn_in = 200,
                             lag_strength = "equal",
                             fracdiff_trim = 150,
                             far_benchmark_d = 1,
                             order_selection_method = c("adaptive_ic", "cv", "both"),
                             cv_order_window_type = c("expanding", "fixed")) {
  order_selection_method <- match.arg(order_selection_method)
  cv_order_window_type <- match.arg(cv_order_window_type, several.ok = TRUE)
  window_types <- match.arg(window_types, c("expanding", "fixed"), several.ok = TRUE)

  if (!"scenario_id" %in% names(grid)) {
    grid <- grid |>
      dplyr::mutate(scenario_id = dplyr::row_number(),
                    .before = 1)
  }

  progressr::with_progress({
    seeds <- seq_len(n_rep)

    p_prog <- progressr::progressor(steps = length(seeds) * nrow(grid))

    results <- future.apply::future_lapply(seeds, function(s) {
      estimate_list <- list()
      error_curve_list <- list()
      conformal_summary_list <- list()
      error_list <- list()

      for (i in seq_len(nrow(grid))) {
        scenario_i <- grid$scenario_id[i]
        n_i <- grid$T[i]
        d_i <- grid$d[i]
        ar_i <- grid$ar_strength[i]
        ma_i <- grid$ma_strength[i]
        p_i <- grid$p[i]
        q_i <- grid$q[i]

        # Keyed on the scenario identity rather than the row position, so a
        # subset of the grid reproduces the same data as the full grid.
        sim_seed <- s * 100000 + scenario_i

        sim_out <- tryCatch({
          sim <- sim_FARMA(
            n = n_i,
            d = d_i,
            no_grid = no_grid,
            seed_number = sim_seed,
            p = p_i,
            q = q_i,
            ar_strength = ar_i,
            ma_strength = ma_i,
            lag_strength = lag_strength
          )

          sim$X
        }, error = function(e)
          e)

        if (inherits(sim_out, "error")) {
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
              order_selection_method = order_selection_method,
              window_type = NA_character_,
              stage = "simulation",
              t0 = NA_integer_,
              error_message = sim_out$message
            )

          p_prog()
          next
        }

        X <- sim_out

        out <- tryCatch(
          unified_pipeline(
            X = X,
            h_max = h_max,
            train_ratio = train_ratio,
            valid_ratio = valid_ratio,
            alpha_set = alpha_set,
            window_types = window_types,
            p_max = p_max,
            q_max = q_max,
            Npc = Npc,
            far_pve = far_pve,
            zeta_grid = zeta_grid,
            m_const = m_const,
            R_noise = R_noise,
            B = B,
            burn_in = burn_in,
            fracdiff_trim = fracdiff_trim,
            seed = 500000 + sim_seed,
            far_benchmark_d = far_benchmark_d,
            order_selection_method = order_selection_method,
            cv_order_window_type = cv_order_window_type
          ),
          error = function(e)
            e
        )

        if (inherits(out, "error")) {
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
              order_selection_method = order_selection_method,
              window_type = NA_character_,
              stage = "centered_unified_pipeline",
              t0 = NA_integer_,
              error_message = out$message
            )

          p_prog()
          next
        }

        estimate_list[[length(estimate_list) + 1]] <-
          out$model_summary |>
          dplyr::mutate(
            scenario_id = scenario_i,
            n = n_i,
            d_true = d_i,
            d_bias = d_hat - d_i,
            ar_strength = ar_i,
            ma_strength = ma_i,
            p_true = p_i,
            q_true = q_i,
            order_correct = p_est == p_i & q_est == q_i,
            seed = s,
            sim_seed = sim_seed,
            .before = 1
          )

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

        if (nrow(out$error_detail) > 0) {
          error_curve_list[[length(error_curve_list) + 1]] <-
            compact_forecast_errors(out$error_detail) |>
            dplyr::mutate(!!!metadata, .before = 1)
        }

        if (nrow(out$conformal_summary_detail) > 0) {
          conformal_summary_list[[length(conformal_summary_list) + 1]] <-
            out$conformal_summary_detail |>
            dplyr::mutate(!!!metadata, .before = 1)
        }

        if (nrow(out$error_log) > 0) {
          error_list[[length(error_list) + 1]] <-
            out$error_log |>
            dplyr::mutate(!!!metadata, .before = 1)
        }

        rm(X, out, sim_out)
        gc()

        p_prog()
      }

      list(
        estimate_detail = dplyr::bind_rows(estimate_list),
        error_detail = dplyr::bind_rows(error_curve_list),
        conformal_summary_detail = dplyr::bind_rows(conformal_summary_list),
        error_log = dplyr::bind_rows(error_list)
      )
    }, future.seed = TRUE)

    list(
      grid = grid,
      estimate_detail = dplyr::bind_rows(lapply(results, function(z)
        z$estimate_detail)),
      error_detail = dplyr::bind_rows(lapply(results, function(z)
        z$error_detail)),
      conformal_summary_detail = dplyr::bind_rows(lapply(results, function(z)
        z$conformal_summary_detail)),
      error_log = dplyr::bind_rows(lapply(results, function(z)
        z$error_log))
    )
  })
}
