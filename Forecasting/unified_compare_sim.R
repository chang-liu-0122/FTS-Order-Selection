rolling_forecasting_pipeline_compact_order_method <- function(X,
                                                              h_max = 10,
                                                              train_ratio = 0.6,
                                                              valid_ratio = 0.2,
                                                              alpha_set = c(0.2, 0.05),
                                                              window_types = c("expanding", "fixed"),
                                                              p_max = 5,
                                                              q_max = 5,
                                                              Npc = 1,
                                                              far_pve = 0.95,
                                                              zeta_grid = seq(0.1, 10, length.out = 500),
                                                              m_const = 0.65,
                                                              R_noise = 0,
                                                              B = 100,
                                                              burn_in = 200,
                                                              fracdiff_trim = 100,
                                                              seed = 1,
                                                              far_benchmark_d = NULL,
                                                              order_selection_method = c("adaptive_ic", "cv", "both"),
                                                              cv_order_window_type = c("expanding", "fixed")) {
  order_selection_method <- match.arg(order_selection_method)
  cv_order_window_type <- match.arg(cv_order_window_type, several.ok = TRUE)

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

  TT <- ncol(X)
  train_end <- floor(TT * train_ratio)
  valid_end <- train_end + floor(TT * valid_ratio)

  if (train_end <= fracdiff_trim) {
    stop("The initial training sample is too short after applying fracdiff_trim.")
  }

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
        rolling_split_conformal_original_ecp_compact(
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
          fracdiff_trim = fracdiff_trim
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

    d_hat <- estimate_d_BC_LPWN2(
      X_initial_train_centered,
      m_const = m_const,
      R_noise = R_noise,
      B = B,
      burn_in = burn_in
    )$d_BC_LPWN2

    d_hat <- as.numeric(d_hat)

    if (length(d_hat) != 1 || !is.finite(d_hat)) {
      stop("Estimated d is not a single finite value.")
    }

    list(d_hat = d_hat, X_initial_train_centered = X_initial_train_centered)
  }

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

  run_cv_one_window <- function(cv_wtype) {
    cv_order <- select_order_cv_forecast_full(
      X = X,
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
      window_type = cv_wtype
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

parallel_repeat_full_forecasting_pipeline_compact_order_method <- function(n_rep = 100,
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
                                                                           zeta_grid = seq(0.1, 10, length.out = 500),
                                                                           m_const = 0.65,
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
      calibration_list <- list()
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

        sim_seed <- s * 100000 + i

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
          rolling_forecasting_pipeline_compact_order_method(
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
            out$error_detail |>
            dplyr::mutate(!!!metadata, .before = 1)
        }

        if (nrow(out$calibration_detail) > 0) {
          calibration_list[[length(calibration_list) + 1]] <-
            out$calibration_detail |>
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
        calibration_detail = dplyr::bind_rows(calibration_list),
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
      calibration_detail = dplyr::bind_rows(lapply(results, function(z)
        z$calibration_detail)),
      conformal_summary_detail = dplyr::bind_rows(lapply(results, function(z)
        z$conformal_summary_detail)),
      error_log = dplyr::bind_rows(lapply(results, function(z)
        z$error_log))
    )
  })
}




rolling_forecasting_pipeline_empirical <- function(X,
                                                   h_max = 10,
                                                   train_ratio = 0.6,
                                                   valid_ratio = 0.2,
                                                   alpha_set = c(0.2, 0.05),
                                                   window_types = c("expanding", "fixed"),
                                                   p_max = 5,
                                                   q_max = 5,
                                                   Npc = 1,
                                                   far_pve = 0.95,
                                                   zeta_grid = seq(0.1, 10, length.out = 500),
                                                   m_const = 0.65,
                                                   R_noise = 0,
                                                   B = 100,
                                                   burn_in = 200,
                                                   fracdiff_trim = 100,
                                                   seed = 1,
                                                   far_benchmark_d = NULL)
{
  X <- as.matrix(X)

  if (any(!is.finite(X))) {
    stop("X contains NA, NaN or Inf values.")
  }

  T <- ncol(X)
  p_grid <- nrow(X)

  h_max <- as.integer(h_max[1])

  train_end <- floor(T * train_ratio)

  valid_end <-
    train_end +
    floor(T * valid_ratio)

  if (train_end <= fracdiff_trim) {
    stop(paste(
      "The initial training sample is too short",
      "after applying fracdiff_trim."
    ))
  }

  if (valid_end <= train_end || valid_end >= T) {
    stop(
      paste(
        "train_ratio and valid_ratio must produce",
        "non-empty validation and test blocks."
      )
    )
  }

  if (h_max < 1L) {
    stop("h_max must be at least 1.")
  }

  # ==========================================================
  # Initial training sample
  # ==========================================================

  X_initial_train <-
    X[, 1:train_end, drop = FALSE]

  mean_initial <-
    rowMeans(X_initial_train)

  X_initial_train_centered <-
    sweep(X_initial_train,
          MARGIN = 1,
          STATS = mean_initial,
          FUN = "-")

  # ==========================================================
  # Estimate the memory parameter
  # ==========================================================

  set.seed(seed)

  d_hat <-
    estimate_d_BC_LPWN2(
      X_initial_train_centered,
      m_const = m_const,
      R_noise = R_noise,
      B = B,
      burn_in = burn_in
    )$d_BC_LPWN2

  d_hat <- as.numeric(d_hat)

  if (length(d_hat) != 1 ||
      !is.finite(d_hat)) {
    stop("Estimated d is not a single finite value.")
  }

  if (is.null(far_benchmark_d)) {
    d_hat_far <- d_hat
  } else {
    d_hat_far <-
      as.numeric(far_benchmark_d)
  }

  if (length(d_hat_far) != 1 ||
      !is.finite(d_hat_far)) {
    stop(paste("far_benchmark_d must be NULL", "or a single finite value."))
  }

  # ==========================================================
  # Fractionally difference the initial training sample
  # ==========================================================

  if (d_hat == 0) {
    X_diff_initial_train <-
      X_initial_train_centered
  } else {
    X_diff_initial_train <-
      FracDiff(t(X_initial_train_centered), d_hat) %>%
      t()
  }

  if (any(!is.finite(X_diff_initial_train))) {
    stop(paste(
      "Fractionally differenced training curves",
      "contain NA or Inf."
    ))
  }

  if (d_hat == 0) {
    X_diff_for_order_selection <-
      X_diff_initial_train
  } else {
    X_diff_for_order_selection <-
      X_diff_initial_train[, (fracdiff_trim + 1):ncol(X_diff_initial_train), drop = FALSE]
  }

  # ==========================================================
  # Select the FARMA order using the initial training sample
  # ==========================================================

  order_out <-
    select_order_far_full(
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
    stop(paste(
      "Order selection did not return",
      "valid values for p_est and q_est."
    ))
  }

  # ==========================================================
  # Rolling forecasting for one window type
  # ==========================================================

  rolling_split_centered <- function(window_type)
  {
    window_type <-
      match.arg(window_type, c("expanding", "fixed"))

    window_size <- train_end

    methods_used <-
      c("farma", "far")

    d_method <-
      list(farma = d_hat, far = d_hat_far)

    # --------------------------------------------------------
    # Horizon-specific forecast origins
    #
    # At each origin, h_here determines how many forecast
    # horizons remain inside the validation or test block.
    # --------------------------------------------------------

    valid_steps <-
      seq.int(from = train_end, to = valid_end - 1L)

    test_steps <-
      seq.int(from = valid_end, to = T - 1L)

    get_centered_history <- function(t0)
    {
      X_history <-
        X[, 1:t0, drop = FALSE]

      if (window_type == "expanding") {
        X_mean_sample <- X_history
      } else {
        raw_start_idx <-
          max(1, t0 - window_size + 1)

        X_mean_sample <-
          X[, raw_start_idx:t0, drop = FALSE]
      }

      mean_curve <-
        rowMeans(X_mean_sample)

      X_history_centered <-
        sweep(X_history,
              MARGIN = 1,
              STATS = mean_curve,
              FUN = "-")

      list(
        X_history = X_history,
        X_history_centered =
          X_history_centered,
        mean_curve = mean_curve
      )
    }

    make_Y_train <- function(X_history_centered, d_current)
    {
      if (d_current == 0) {
        Y_history <-
          X_history_centered
      } else {
        Y_history <-
          frac_difference_history(X_history =
                                    X_history_centered, d_hat = d_current)
      }

      if (d_current == 0) {
        Y_estimation_history <-
          Y_history
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
        Y_train <-
          Y_estimation_history
      } else {
        start_idx <-
          max(1, ncol(Y_estimation_history) -
                window_size + 1)

        Y_train <-
          Y_estimation_history[, start_idx:ncol(Y_estimation_history), drop = FALSE]
      }

      list(Y_history = Y_history, Y_train = Y_train)
    }

    error_list <- list()

    # E_valid[[method]][[h]] contains one column per
    # available validation origin for horizon h.
    E_valid <-
      setNames(lapply(methods_used, function(m)
      {
        lapply(seq_len(h_max), function(h)
        {
          matrix(numeric(0), nrow = p_grid, ncol = 0)
        })
      }), methods_used)

    # ========================================================
    # Validation forecasts
    # ========================================================

    for (k in seq_along(valid_steps))
    {
      t0 <- valid_steps[k]

      # Keep all targets within the validation block.
      h_here <-
        min(h_max, valid_end - t0)

      if (h_here < 1L) {
        next
      }

      centered_obj <-
        get_centered_history(t0 = t0)

      true_future_mat <-
        X[, (t0 + 1):(t0 + h_here), drop = FALSE]

      for (m in methods_used)
      {
        d_current <-
          d_method[[m]]

        train_obj <-
          tryCatch(
            make_Y_train(
              X_history_centered =
                centered_obj$X_history_centered,
              d_current =
                d_current
            ),
            error = function(e)
              e
          )

        if (inherits(train_obj, "error"))
        {
          error_list[[length(error_list) + 1]] <-
            data.frame(
              seed = seed,
              window_type = window_type,
              method = m,
              stage =
                "validation_training",
              t0 = t0,
              error_message =
                train_obj$message
            )

          next
        }

        out <-
          tryCatch(
            forecast_original_scale(
              Y_train =
                train_obj$Y_train,
              d_hat =
                d_current,
              h =
                h_here,
              method =
                m,
              p_est =
                p_est,
              q_est =
                q_est,
              Npc =
                Npc,
              far_pve =
                far_pve,
              Y_history_for_integration =
                train_obj$Y_history,
              return_fit =
                FALSE
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
              stage =
                "validation_forecast",
              t0 = t0,
              error_message =
                out$message
            )

          next
        }

        if (is.null(out[[m]])) {
          next
        }

        pred_centered_mat <-
          out[[m]]$X_future

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
          e_h <-
            true_future_mat[, h] -
            pred_future_mat[, h]

          E_valid[[m]][[h]] <-
            cbind(E_valid[[m]][[h]], e_h)
        }
      }
    }

    # ========================================================
    # Pointwise ECP calibration
    # ========================================================

    calibration <- list()
    calibration_list <- list()

    for (m in methods_used)
    {
      calibration[[m]] <- list()

      for (alpha in alpha_set)
      {
        alpha_name <-
          as.character(alpha)

        calibration[[m]][[alpha_name]] <-
          list(gamma =
                 vector("list", h_max),
               zeta =
                 rep(NA_real_, h_max))

        for (h in seq_len(h_max))
        {
          E_h <-
            E_valid[[m]][[h]]

          if (ncol(E_h) == 0) {
            next
          }

          # This retains the pointwise ECP definition.
          cal <-
            calibrate_zeta_ecp(E_valid = E_h,
                               alpha = alpha,
                               zeta_grid =
                                 zeta_grid)

          calibration[[m]][[alpha_name]]$gamma[[h]] <-
            cal$gamma

          calibration[[m]][[alpha_name]]$zeta[h] <-
            cal$zeta
        }

        calibration_list[[length(calibration_list) + 1]] <-
          tibble::tibble(
            seed = seed,
            window_type =
              window_type,
            method = m,
            alpha = alpha,
            d_used =
              d_method[[m]],
            gamma =
              list(calibration[[m]][[alpha_name]]$gamma),
            zeta =
              list(calibration[[m]][[alpha_name]]$zeta)
          )
      }
    }

    # ========================================================
    # Test accumulators
    # ========================================================

    conformal_accum <-
      setNames(lapply(methods_used, function(m)
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

    # ========================================================
    # Test forecasts
    # ========================================================

    for (k in seq_along(test_steps))
    {
      t0 <- test_steps[k]

      # Keep all targets inside the observed test block.
      h_here <-
        min(h_max, T - t0)

      if (h_here < 1L) {
        next
      }

      centered_obj <-
        get_centered_history(t0 = t0)

      true_future_mat <-
        X[, (t0 + 1):(t0 + h_here), drop = FALSE]

      for (m in methods_used)
      {
        d_current <-
          d_method[[m]]

        train_obj <-
          tryCatch(
            make_Y_train(
              X_history_centered =
                centered_obj$X_history_centered,
              d_current =
                d_current
            ),
            error = function(e)
              e
          )

        if (inherits(train_obj, "error"))
        {
          error_list[[length(error_list) + 1]] <-
            data.frame(
              seed = seed,
              window_type = window_type,
              method = m,
              stage =
                "testing_training",
              t0 = t0,
              error_message =
                train_obj$message
            )

          next
        }

        out <-
          tryCatch(
            forecast_original_scale(
              Y_train =
                train_obj$Y_train,
              d_hat =
                d_current,
              h =
                h_here,
              method =
                m,
              p_est =
                p_est,
              q_est =
                q_est,
              Npc =
                Npc,
              far_pve =
                far_pve,
              Y_history_for_integration =
                train_obj$Y_history,
              return_fit =
                FALSE
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
              stage =
                "testing_forecast",
              t0 = t0,
              error_message =
                out$message
            )

          next
        }

        if (is.null(out[[m]])) {
          next
        }

        pred_centered_mat <-
          out[[m]]$X_future

        pred_future_mat <-
          pred_centered_mat +
          matrix(
            centered_obj$mean_curve,
            nrow = p_grid,
            ncol = h_here,
            byrow = FALSE
          )

        error_future_mat <-
          pred_future_mat -
          true_future_mat

        error_curve_list[[length(error_curve_list) + 1]] <-
          tibble::tibble(
            seed = seed,
            method = m,
            window_type =
              window_type,
            step = k,
            t0 = t0,
            d_used =
              d_current,
            error =
              list(error_future_mat)
          )

        for (alpha in alpha_set)
        {
          alpha_name <-
            as.character(alpha)

          for (h in seq_len(h_here))
          {
            error_h <-
              as.numeric(error_future_mat[, h])

            gamma_h <-
              as.numeric(calibration[[m]][[alpha_name]]$gamma[[h]])

            zeta_h <-
              calibration[[m]][[alpha_name]]$zeta[h]

            if (length(gamma_h) !=
                length(error_h) ||
                !is.finite(zeta_h)) {
              next
            }

            half_width_h <-
              zeta_h * gamma_h

            ok <-
              is.finite(error_h) &
              is.finite(half_width_h)

            error_h <-
              error_h[ok]

            half_width_h <-
              half_width_h[ok]

            if (length(error_h) == 0) {
              next
            }

            score_h <-
              2 * half_width_h +
              2 / alpha *
              pmax(abs(error_h) -
                     half_width_h, 0)

            accum <-
              conformal_accum[[m]][[alpha_name]][[h]]

            accum$cover_count <-
              accum$cover_count +
              sum(abs(error_h) <=
                    half_width_h)

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

    # ========================================================
    # Pointwise test summaries
    # ========================================================

    conformal_summary_list <- list()

    for (m in methods_used)
    {
      for (alpha in alpha_set)
      {
        alpha_name <-
          as.character(alpha)

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
            abs(ECP -
                  (1 - alpha))

          Score <-
            accum$score_sum /
            accum$n_interval

          conformal_summary_list[[length(conformal_summary_list) + 1]] <-
            tibble::tibble(
              seed = seed,
              method = m,
              window_type =
                window_type,
              alpha = alpha,
              h = h,
              d_used =
                d_method[[m]],
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

  # ==========================================================
  # Run expanding and/or fixed windows
  # ==========================================================

  error_curve_list <- list()
  calibration_list <- list()
  conformal_summary_list <- list()
  error_list <- list()

  for (wtype in window_types)
  {
    out <-
      tryCatch(
        rolling_split_centered(window_type = wtype),
        error = function(e)
          e
      )

    if (inherits(out, "error"))
    {
      error_list[[length(error_list) + 1]] <-
        tibble::tibble(
          window_type = wtype,
          stage =
            "rolling_forecast",
          t0 = NA_integer_,
          error_message =
            out$message
        )

      next
    }

    if (nrow(out$error_detail) > 0) {
      error_curve_list[[length(error_curve_list) + 1]] <-
        out$error_detail
    }

    if (nrow(out$calibration_detail) > 0) {
      calibration_list[[length(calibration_list) + 1]] <-
        out$calibration_detail
    }

    if (nrow(out$conformal_summary_detail) > 0) {
      conformal_summary_list[[length(conformal_summary_list) + 1]] <-
        out$conformal_summary_detail
    }

    if (nrow(out$error_log) > 0) {
      error_list[[length(error_list) + 1]] <-
        out$error_log
    }
  }

  # ==========================================================
  # Return
  # ==========================================================

  list(
    model_summary =
      tibble::tibble(
        T = T,
        train_end = train_end,
        valid_end = valid_end,

        # First test observation is valid_end + 1.
        # The first forecast origin is valid_end.
        test_start = valid_end + 1,

        centering =
          paste("expanding_or_fixed", "training_window_mean", sep = "_"),

        d_hat = d_hat,
        d_hat_far = d_hat_far,
        p_est = p_est,
        q_est = q_est,
        fracdiff_trim =
          fracdiff_trim
      ),

    order_selection_detail =
      order_out,

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



# rolling_forecasting_pipeline_empirical <- function(X,
#                                                    h_max = 10,
#                                                    train_ratio = 0.6,
#                                                    valid_ratio = 0.2,
#                                                    alpha_set = c(0.2, 0.05),
#                                                    window_types = c("expanding", "fixed"),
#                                                    p_max = 5,
#                                                    q_max = 5,
#                                                    Npc = 1,
#                                                    far_pve = 0.95,
#                                                    zeta_grid = seq(0.1, 10, length.out = 500),
#                                                    m_const = 0.65,
#                                                    R_noise = 0,
#                                                    B = 100,
#                                                    burn_in = 200,
#                                                    fracdiff_trim = 100,
#                                                    seed = 1,
#                                                    far_benchmark_d = NULL)
# {
#   X <- as.matrix(X)
#
#   if (any(!is.finite(X)))
#   {
#     stop("X contains NA, NaN or Inf values.")
#   }
#
#   T <- ncol(X)
#   p_grid <- nrow(X)
#
#   train_end <- floor(T * train_ratio)
#
#   valid_end <- train_end +
#     floor(T * valid_ratio)
#
#   if (train_end <= fracdiff_trim)
#   {
#     stop("The initial training sample is too short after applying fracdiff_trim.")
#   }
#
#   X_initial_train <- X[, 1:train_end, drop = FALSE]
#
#   mean_initial <- rowMeans(X_initial_train)
#
#   X_initial_train_centered <- sweep(X_initial_train,
#                                     MARGIN = 1,
#                                     STATS = mean_initial,
#                                     FUN = "-")
#
#   set.seed(seed)
#
#   d_hat <- estimate_d_BC_LPWN2(
#     X_initial_train_centered,
#     m_const = m_const,
#     R_noise = R_noise,
#     B = B,
#     burn_in = burn_in
#   )$d_BC_LPWN2
#
#   d_hat <- as.numeric(d_hat)
#
#   if (length(d_hat) != 1 ||
#       !is.finite(d_hat))
#   {
#     stop("Estimated d is not a single finite value.")
#   }
#
#   if (is.null(far_benchmark_d)) {
#     d_hat_far <- d_hat
#   } else {
#     d_hat_far <- as.numeric(far_benchmark_d)
#   }
#
#   if (length(d_hat_far) != 1 ||
#       !is.finite(d_hat_far))
#   {
#     stop("far_benchmark_d must be NULL or a single finite value.")
#   }
#
#   if (d_hat == 0)
#   {
#     X_diff_initial_train <- X_initial_train_centered
#   } else {
#     X_diff_initial_train <- FracDiff(t(X_initial_train_centered), d_hat) %>%
#       t()
#   }
#
#   if (any(!is.finite(X_diff_initial_train)))
#   {
#     stop("Fractionally differenced training curves contain NA or Inf.")
#   }
#
#   if (d_hat == 0)
#   {
#     X_diff_for_order_selection <- X_diff_initial_train
#   } else {
#     X_diff_for_order_selection <-
#       X_diff_initial_train[, (fracdiff_trim + 1):ncol(X_diff_initial_train), drop = FALSE]
#   }
#
#   order_out <- select_order_far_full(
#     X = X_diff_for_order_selection,
#     p_max = p_max,
#     q_max = q_max,
#     Npc = Npc
#   )
#
#   p_est <- as.integer(order_out$p_hat)
#   q_est <- as.integer(order_out$q_hat)
#
#   if (length(p_est) != 1 ||
#       length(q_est) != 1 ||
#       !is.finite(p_est) ||
#       !is.finite(q_est))
#   {
#     stop("Order selection did not return valid values for p_est and q_est.")
#   }
#
#   rolling_split_centered <- function(window_type)
#   {
#     window_type <- match.arg(window_type, c("expanding", "fixed"))
#
#     h_max <- as.integer(h_max[1])
#
#     window_size <- train_end
#
#     methods_used <- c("farma", "far")
#
#     d_method <- list(farma = d_hat, far = d_hat_far)
#
#     valid_steps <- seq(train_end, valid_end - h_max)
#
#     test_steps <- seq(valid_end + 1, T - h_max)
#
#     get_centered_history <- function(t0)
#     {
#       X_history <- X[, 1:t0, drop = FALSE]
#
#       if (window_type == "expanding")
#       {
#         X_mean_sample <- X_history
#       } else {
#         raw_start_idx <- max(1, t0 - window_size + 1)
#
#         X_mean_sample <- X[, raw_start_idx:t0, drop = FALSE]
#       }
#
#       mean_curve <- rowMeans(X_mean_sample)
#
#       X_history_centered <- sweep(X_history,
#                                   MARGIN = 1,
#                                   STATS = mean_curve,
#                                   FUN = "-")
#
#       list(
#         X_history = X_history,
#         X_history_centered = X_history_centered,
#         mean_curve = mean_curve
#       )
#     }
#
#     make_Y_train <- function(X_history_centered, d_current)
#     {
#       if (d_current == 0) {
#         Y_history <- X_history_centered
#       } else {
#         Y_history <- frac_difference_history(X_history = X_history_centered, d_hat = d_current)
#       }
#
#       if (d_current == 0) {
#         Y_estimation_history <- Y_history
#       } else {
#         Y_estimation_history <-
#           Y_history[, (fracdiff_trim + 1):ncol(Y_history), drop = FALSE]
#       }
#
#       if (window_type == "expanding") {
#         Y_train <- Y_estimation_history
#       } else {
#         start_idx <- max(1, ncol(Y_estimation_history) - window_size + 1)
#
#         Y_train <- Y_estimation_history[, start_idx:ncol(Y_estimation_history), drop = FALSE]
#       }
#
#       list(Y_history = Y_history, Y_train = Y_train)
#     }
#
#     error_list <- list()
#
#     E_valid <- setNames(lapply(methods_used, function(m)
#     {
#       lapply(seq_len(h_max), function(h)
#       {
#         matrix(numeric(0), nrow = p_grid, ncol = 0)
#       })
#     }), methods_used)
#
#     for (k in seq_along(valid_steps))
#     {
#       t0 <- valid_steps[k]
#
#       centered_obj <- get_centered_history(t0 = t0)
#
#       true_future_mat <- X[, (t0 + 1):(t0 + h_max), drop = FALSE]
#
#       for (m in methods_used)
#       {
#         d_current <- d_method[[m]]
#
#         train_obj <- make_Y_train(X_history_centered = centered_obj$X_history_centered,
#                                   d_current = d_current)
#
#         out <- tryCatch(
#           forecast_original_scale(
#             Y_train = train_obj$Y_train,
#             d_hat = d_current,
#             h = h_max,
#             method = m,
#             p_est = p_est,
#             q_est = q_est,
#             Npc = Npc,
#             far_pve = far_pve,
#             Y_history_for_integration = train_obj$Y_history,
#             return_fit = FALSE
#           ),
#           error = function(e)
#             e
#         )
#
#         if (inherits(out, "error"))
#         {
#           error_list[[length(error_list) + 1]] <-
#             data.frame(
#               seed = seed,
#               window_type = window_type,
#               method = m,
#               stage = "validation_forecast",
#               t0 = t0,
#               error_message = out$message
#             )
#
#           next
#         }
#
#         if (is.null(out[[m]])) {
#           next
#         }
#
#         pred_centered_mat <- out[[m]]$X_future
#
#         pred_future_mat <- pred_centered_mat +
#           matrix(
#             centered_obj$mean_curve,
#             nrow = p_grid,
#             ncol = h_max,
#             byrow = FALSE
#           )
#
#         for (h in seq_len(h_max))
#         {
#           e_h <- true_future_mat[, h] -
#             pred_future_mat[, h]
#
#           E_valid[[m]][[h]] <-
#             cbind(E_valid[[m]][[h]], e_h)
#         }
#       }
#     }
#
#     calibration <- list()
#     calibration_list <- list()
#
#     for (m in methods_used)
#     {
#       calibration[[m]] <- list()
#
#       for (alpha in alpha_set)
#       {
#         alpha_name <- as.character(alpha)
#
#         calibration[[m]][[alpha_name]] <-
#           list(gamma = vector("list", h_max),
#                zeta = rep(NA_real_, h_max))
#
#         for (h in seq_len(h_max))
#         {
#           E_h <- E_valid[[m]][[h]]
#
#           if (ncol(E_h) == 0) {
#             next
#           }
#
#           cal <- calibrate_zeta_ecp(E_valid = E_h,
#                                     alpha = alpha,
#                                     zeta_grid = zeta_grid)
#
#           calibration[[m]][[alpha_name]]$gamma[[h]] <-
#             cal$gamma
#
#           calibration[[m]][[alpha_name]]$zeta[h] <-
#             cal$zeta
#         }
#
#         calibration_list[[length(calibration_list) + 1]] <-
#           tibble::tibble(
#             seed = seed,
#             window_type = window_type,
#             method = m,
#             alpha = alpha,
#             d_used = d_method[[m]],
#             gamma = list(calibration[[m]][[alpha_name]]$gamma),
#             zeta = list(calibration[[m]][[alpha_name]]$zeta)
#           )
#       }
#     }
#
#     conformal_accum <- setNames(lapply(methods_used, function(m)
#     {
#       setNames(lapply(as.character(alpha_set), function(alpha_name)
#       {
#         lapply(seq_len(h_max), function(h)
#         {
#           list(
#             cover_count = 0,
#             score_sum = 0,
#             n_interval = 0
#           )
#         })
#       }),
#       as.character(alpha_set))
#     }), methods_used)
#
#     error_curve_list <- list()
#
#     for (k in seq_along(test_steps))
#     {
#       t0 <- test_steps[k]
#
#       centered_obj <- get_centered_history(t0 = t0)
#
#       true_future_mat <- X[, (t0 + 1):(t0 + h_max), drop = FALSE]
#
#       for (m in methods_used)
#       {
#         d_current <- d_method[[m]]
#
#         train_obj <- make_Y_train(X_history_centered = centered_obj$X_history_centered,
#                                   d_current = d_current)
#
#         out <- tryCatch(
#           forecast_original_scale(
#             Y_train = train_obj$Y_train,
#             d_hat = d_current,
#             h = h_max,
#             method = m,
#             p_est = p_est,
#             q_est = q_est,
#             Npc = Npc,
#             far_pve = far_pve,
#             Y_history_for_integration = train_obj$Y_history,
#             return_fit = FALSE
#           ),
#           error = function(e)
#             e
#         )
#
#         if (inherits(out, "error"))
#         {
#           error_list[[length(error_list) + 1]] <-
#             data.frame(
#               seed = seed,
#               window_type = window_type,
#               method = m,
#               stage = "testing_forecast",
#               t0 = t0,
#               error_message = out$message
#             )
#
#           next
#         }
#
#         if (is.null(out[[m]])) {
#           next
#         }
#
#         pred_centered_mat <- out[[m]]$X_future
#
#         pred_future_mat <- pred_centered_mat +
#           matrix(
#             centered_obj$mean_curve,
#             nrow = p_grid,
#             ncol = h_max,
#             byrow = FALSE
#           )
#
#         error_future_mat <-
#           pred_future_mat -
#           true_future_mat
#
#         error_curve_list[[length(error_curve_list) + 1]] <-
#           tibble::tibble(
#             seed = seed,
#             method = m,
#             window_type = window_type,
#             step = k,
#             t0 = t0,
#             d_used = d_current,
#             error = list(error_future_mat)
#           )
#
#         for (alpha in alpha_set)
#         {
#           alpha_name <- as.character(alpha)
#
#           for (h in seq_len(h_max))
#           {
#             error_h <- as.numeric(error_future_mat[, h])
#
#             gamma_h <- as.numeric(calibration[[m]][[alpha_name]]$gamma[[h]])
#
#             zeta_h <-
#               calibration[[m]][[alpha_name]]$zeta[h]
#
#             if (length(gamma_h) != length(error_h) ||
#                 !is.finite(zeta_h)) {
#               next
#             }
#
#             half_width_h <- zeta_h * gamma_h
#
#             ok <-
#               is.finite(error_h) &
#               is.finite(half_width_h)
#
#             error_h <- error_h[ok]
#             half_width_h <- half_width_h[ok]
#
#             if (length(error_h) == 0) {
#               next
#             }
#
#             score_h <-
#               2 * half_width_h +
#               2 / alpha *
#               pmax(abs(error_h) - half_width_h, 0)
#
#             accum <-
#               conformal_accum[[m]][[alpha_name]][[h]]
#
#             accum$cover_count <-
#               accum$cover_count +
#               sum(abs(error_h) <= half_width_h)
#
#             accum$score_sum <-
#               accum$score_sum +
#               sum(score_h)
#
#             accum$n_interval <-
#               accum$n_interval +
#               length(error_h)
#
#             conformal_accum[[m]][[alpha_name]][[h]] <-
#               accum
#           }
#         }
#       }
#     }
#
#     conformal_summary_list <- list()
#
#     for (m in methods_used)
#     {
#       for (alpha in alpha_set)
#       {
#         alpha_name <- as.character(alpha)
#
#         for (h in seq_len(h_max))
#         {
#           accum <- conformal_accum[[m]][[alpha_name]][[h]]
#
#           if (accum$n_interval == 0) {
#             next
#           }
#
#           ECP <- accum$cover_count /
#             accum$n_interval
#
#           CPD <- abs(ECP -
#                        (1 - alpha))
#
#           Score <- accum$score_sum /
#             accum$n_interval
#
#           conformal_summary_list[[length(conformal_summary_list) + 1]] <-
#             tibble::tibble(
#               seed = seed,
#               method = m,
#               window_type = window_type,
#               alpha = alpha,
#               h = h,
#               d_used = d_method[[m]],
#               ECP = ECP,
#               CPD = CPD,
#               Score = Score
#             )
#         }
#       }
#     }
#
#     list(
#       error_detail =
#         dplyr::bind_rows(error_curve_list),
#
#       calibration_detail =
#         dplyr::bind_rows(calibration_list),
#
#       conformal_summary_detail =
#         dplyr::bind_rows(conformal_summary_list),
#
#       error_log =
#         dplyr::bind_rows(error_list)
#     )
#   }
#
#   error_curve_list <- list()
#   calibration_list <- list()
#   conformal_summary_list <- list()
#   error_list <- list()
#
#   for (wtype in window_types)
#   {
#     out <- tryCatch(
#       rolling_split_centered(window_type = wtype),
#       error = function(e)
#         e
#     )
#
#     if (inherits(out, "error"))
#     {
#       error_list[[length(error_list) + 1]] <-
#         tibble::tibble(
#           window_type = wtype,
#           stage = "rolling_forecast",
#           t0 = NA_integer_,
#           error_message = out$message
#         )
#
#       next
#     }
#
#     if (nrow(out$error_detail) > 0)
#     {
#       error_curve_list[[length(error_curve_list) + 1]] <-
#         out$error_detail
#     }
#
#     if (nrow(out$calibration_detail) > 0)
#     {
#       calibration_list[[length(calibration_list) + 1]] <-
#         out$calibration_detail
#     }
#
#     if (nrow(out$conformal_summary_detail) > 0)
#     {
#       conformal_summary_list[[length(conformal_summary_list) + 1]] <-
#         out$conformal_summary_detail
#     }
#
#     if (nrow(out$error_log) > 0)
#     {
#       error_list[[length(error_list) + 1]] <-
#         out$error_log
#     }
#   }
#
#   list(
#     model_summary =
#       tibble::tibble(
#         T = T,
#         train_end = train_end,
#         valid_end = valid_end,
#         test_start = valid_end + 1,
#         centering = "expanding_or_fixed_training_window_mean",
#         d_hat = d_hat,
#         d_hat_far = d_hat_far,
#         p_est = p_est,
#         q_est = q_est,
#         fracdiff_trim = fracdiff_trim
#       ),
#
#     order_selection_detail =
#       order_out,
#
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
