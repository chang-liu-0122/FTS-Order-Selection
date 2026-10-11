forecast_farma_path <- function(fit, Y_train, h)
{
  Yc <- t(scale(t(Y_train), center = TRUE, scale = FALSE))
  d_grid <- nrow(Yc)

  p <- if (is.null(fit$phi))
    0
  else
    dim(fit$phi)[3]
  q <- if (is.null(fit$theta))
    0
  else
    dim(fit$theta)[3]

  Y_hist <- Yc
  eps_hist <- fit$residuals
  eps_hist[is.na(eps_hist)] <- 0

  Y_future <- matrix(NA_real_, nrow = d_grid, ncol = h)

  for (step in 1:h)
  {
    ar_part <- rep(0, d_grid)

    if (p > 0)
    {
      for (j in 1:p)
      {
        ar_part <- ar_part +
          fit$phi[, , j] %*% Y_hist[, ncol(Y_hist) - j + 1]
      }
    }

    ma_part <- rep(0, d_grid)

    if (q > 0)
    {
      for (k in 1:q)
      {
        ma_part <- ma_part +
          fit$theta[, , k] %*% eps_hist[, ncol(eps_hist) - k + 1]
      }
    }

    y_new_centered <- ar_part + ma_part
    y_new <- as.numeric(y_new_centered + fit$mean)

    Y_future[, step] <- y_new

    Y_hist <- cbind(Y_hist, y_new_centered)

    eps_hist <- cbind(eps_hist, matrix(0, nrow = d_grid, ncol = 1))
  }

  Y_future
}

forecast_far_path <- function(Y_train, h, far_pve = 0.95)
{
  d_grid <- nrow(Y_train)

  mean_far <- rowMeans(Y_train, na.rm = TRUE)

  Yc <- sweep(Y_train, 1, mean_far, "-")

  eig_far <- eigen(stats::cov(t(Yc)))

  pve_far <- cumsum(eig_far$values) / sum(eig_far$values)

  kn_far <- which(pve_far >= far_pve)[1]

  fit_far <- far::far(far::as.fdata(Yc), kn = kn_far)

  y_current <- Yc[, ncol(Yc)]

  Y_future <- matrix(NA_real_, nrow = d_grid, ncol = h)

  for (step in 1:h)
  {
    y_next_centered <- as.numeric(stats::predict(fit_far, newdata = far::as.fdata(matrix(
      y_current, ncol = 1
    )))[[1]])

    y_next <- y_next_centered + mean_far

    Y_future[, step] <- y_next

    y_current <- y_next_centered
  }

  list(
    kn = kn_far,
    fit = fit_far,
    mean = mean_far,
    Y_future = Y_future
  )
}

frac_difference_history <- function(X_history, d_hat)
{
  t(FCVAR::FracDiff(x = t(X_history), d = d_hat))
}

frac_integrate_forecast_with_history <- function(Y_history, Y_future, d_hat)
{
  Y_all <- cbind(Y_history, Y_future)

  X_all <- t(FCVAR::FracDiff(x = t(Y_all), d = -d_hat))

  h <- ncol(Y_future)

  X_all[, (ncol(X_all) - h + 1):ncol(X_all), drop = FALSE]
}

farma_companion_radius <- function(operator_array, leading_sign = 1)
{
  if (is.null(operator_array)) {
    return(NA_real_)
  }

  operator_array <- as.array(operator_array)
  operator_dim <- dim(operator_array)

  if (length(operator_dim) != 3L || operator_dim[1] != operator_dim[2]) {
    stop("FARMA operators must be a square matrix array.")
  }

  if (any(!is.finite(operator_array))) {
    return(Inf)
  }

  state_dim <- operator_dim[1]
  n_lags <- operator_dim[3]
  companion <- matrix(0, nrow = state_dim * n_lags,
                      ncol = state_dim * n_lags)

  first_row <- do.call(
    cbind,
    lapply(seq_len(n_lags), function(lag_i) {
      leading_sign * operator_array[, , lag_i]
    })
  )
  companion[seq_len(state_dim), ] <- first_row

  if (n_lags > 1L) {
    lower_rows <- state_dim + seq_len(state_dim * (n_lags - 1L))
    previous_states <- seq_len(state_dim * (n_lags - 1L))
    companion[lower_rows, previous_states] <-
      diag(state_dim * (n_lags - 1L))
  }

  max(Mod(eigen(companion, only.values = TRUE)$values))
}

farma_fit_diagnostics <- function(fit, Y_train)
{
  phi_for_check <- if (!is.null(fit$phi_reduced)) {
    fit$phi_reduced
  } else {
    fit$phi
  }

  theta_for_check <- if (!is.null(fit$theta_reduced)) {
    fit$theta_reduced
  } else {
    fit$theta
  }

  residual_values <- as.numeric(fit$residuals)
  residual_values <- residual_values[is.finite(residual_values)]
  training_sd <- stats::sd(as.numeric(Y_train), na.rm = TRUE)

  residual_ratio <- if (length(residual_values) == 0L ||
                        !is.finite(training_sd) || training_sd <= 0) {
    Inf
  } else {
    max(abs(residual_values)) / training_sd
  }

  list(
    ar_radius = farma_companion_radius(phi_for_check, leading_sign = 1),
    ma_radius = farma_companion_radius(theta_for_check, leading_sign = -1),
    residual_ratio = residual_ratio
  )
}

validate_farma_fit <- function(fit,
                               Y_train,
                               ar_radius_limit = 1,
                               ma_radius_limit = 1,
                               residual_ratio_limit = 50)
{
  diagnostics <- farma_fit_diagnostics(fit = fit, Y_train = Y_train)
  rejection_reasons <- character(0)

  if (!is.null(fit$phi) &&
      (!is.finite(diagnostics$ar_radius) ||
       diagnostics$ar_radius >= ar_radius_limit)) {
    rejection_reasons <- c(
      rejection_reasons,
      sprintf("AR companion radius %.4f >= %.4f",
              diagnostics$ar_radius, ar_radius_limit)
    )
  }

  if (!is.null(fit$theta) &&
      (!is.finite(diagnostics$ma_radius) ||
       diagnostics$ma_radius >= ma_radius_limit)) {
    rejection_reasons <- c(
      rejection_reasons,
      sprintf("MA companion radius %.4f >= %.4f",
              diagnostics$ma_radius, ma_radius_limit)
    )
  }

  if (!is.finite(diagnostics$residual_ratio) ||
      diagnostics$residual_ratio > residual_ratio_limit) {
    rejection_reasons <- c(
      rejection_reasons,
      sprintf("maximum residual ratio %.4f > %.4f",
              diagnostics$residual_ratio, residual_ratio_limit)
    )
  }

  if (length(rejection_reasons) > 0L) {

    stop(
      structure(
        class = c("farma_unstable_fit", "error", "condition"),
        list(
          message = paste("Rejected unstable FARMA fit:",
                          paste(rejection_reasons, collapse = "; ")),
          call = NULL
        )
      )
    )
  }

  invisible(diagnostics)
}

forecast_original_scale <- function(Y_train,
                                    d_hat,
                                    h,
                                    method = c("farma", "far", "both"),
                                    p_est = NULL,
                                    q_est = NULL,
                                    Npc = 1,
                                    far_pve = 0.95,
                                    Y_history_for_integration = NULL,
                                    return_fit = FALSE,
                                    farma_ar_radius_limit = 1,
                                    farma_ma_radius_limit = 1,
                                    farma_residual_ratio_limit = 50)
{
  method <- match.arg(method)

  if (is.null(Y_history_for_integration)) {
    Y_history_for_integration <- Y_train
  }

  out <- list()

  if (method %in% c("farma", "both"))
  {
    if (is.null(p_est) || is.null(q_est)) {
      stop("For method = 'farma' or 'both', p_est and q_est must be provided.")
    }

    fit_farma <- farma_estimation(
      Y_train,
      p = p_est,
      q = q_est,
      m = NULL,
      Npc = Npc
    )

    fit_diagnostics <- validate_farma_fit(
      fit = fit_farma,
      Y_train = Y_train,
      ar_radius_limit = farma_ar_radius_limit,
      ma_radius_limit = farma_ma_radius_limit,
      residual_ratio_limit = farma_residual_ratio_limit
    )

    Y_future_farma <- forecast_farma_path(fit = fit_farma,
                                          Y_train = Y_train,
                                          h = h)

    X_future_farma <- frac_integrate_forecast_with_history(Y_history = Y_history_for_integration,
                                                           Y_future = Y_future_farma,
                                                           d_hat = d_hat)

    out$farma <- list(Y_future = Y_future_farma, X_future = X_future_farma)

    if (return_fit) {
      out$farma$fit <- fit_farma
      out$farma$diagnostics <- fit_diagnostics
    }
  }

  if (method %in% c("far", "both"))
  {
    far_forecast <- forecast_far_path(Y_train = Y_train,
                                      h = h,
                                      far_pve = far_pve)

    Y_future_far <- far_forecast$Y_future

    X_future_far <- frac_integrate_forecast_with_history(Y_history = Y_history_for_integration,
                                                         Y_future = Y_future_far,
                                                         d_hat = d_hat)

    out$far <- list(kn = far_forecast$kn,
                    Y_future = Y_future_far,
                    X_future = X_future_far)

    if (return_fit) {
      out$far$fit <- far_forecast$fit
      out$far$mean <- far_forecast$mean
    }
  }

  out
}

