interval_score <- function(holdout, lb, ub, alpha) {
  lb_ind <- ifelse(holdout < lb, 1, 0)
  ub_ind <- ifelse(holdout > ub, 1, 0)
  score  <- (ub - lb) +
    2 / alpha * ((lb - holdout) * lb_ind + (holdout - ub) * ub_ind)
  cover <- 1 - (length(which(lb_ind == 1)) + length(which(ub_ind == 1))) /
    length(holdout)
  cpd <- abs(cover - (1 - alpha))
  return(c(cover, cpd, mean(score)))
}



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

  eig_far <- eigen(cov(t(Yc)))

  pve_far <- cumsum(eig_far$values) / sum(eig_far$values)

  kn_far <- which(pve_far >= far_pve)[1]

  fit_far <- far(as.fdata(Yc), kn = kn_far)

  y_current <- Yc[, ncol(Yc)]

  Y_future <- matrix(NA_real_, nrow = d_grid, ncol = h)

  for (step in 1:h)
  {
    y_next_centered <- as.numeric(predict(fit_far, newdata = as.fdata(matrix(
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
  FracDiff(x = t(X_history), d = d_hat) %>%
    t()
}




frac_integrate_forecast_with_history <- function(Y_history, Y_future, d_hat)
{
  Y_all <- cbind(Y_history, Y_future)

  X_all <- FracDiff(x = t(Y_all), d = -d_hat) %>%
    t()

  h <- ncol(Y_future)

  X_all[, (ncol(X_all) - h + 1):ncol(X_all), drop = FALSE]
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
                                    return_fit = FALSE)
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

    Y_future_farma <- forecast_farma_path(fit = fit_farma,
                                          Y_train = Y_train,
                                          h = h)

    X_future_farma <- frac_integrate_forecast_with_history(Y_history = Y_history_for_integration,
                                                           Y_future = Y_future_farma,
                                                           d_hat = d_hat)

    out$farma <- list(Y_future = Y_future_farma, X_future = X_future_farma)

    if (return_fit) {
      out$farma$fit <- fit_farma
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



calibrate_zeta_ecp <- function(E_valid,
                               alpha,
                               zeta_grid = seq(0.1, 10, length.out = 500))
{
  E_valid <- as.matrix(E_valid)

  E_centered <- t(scale(t(E_valid), center = TRUE, scale = FALSE))

  gamma <- apply(E_centered, 1, sd, na.rm = TRUE)

  gamma[!is.finite(gamma) | gamma <= 0] <- NA_real_

  ecp_grid <- sapply(zeta_grid, function(z)
  {
    threshold <- matrix(
      z * gamma,
      nrow = nrow(E_valid),
      ncol = ncol(E_valid),
      byrow = FALSE
    )

    mean(abs(E_valid) <= threshold, na.rm = TRUE)
  })

  cand <- zeta_grid[ecp_grid >= (1 - alpha)]

  if (length(cand) == 0) {
    zeta <- NA_real_
  } else {
    zeta <- min(cand)
  }

  list(
    gamma = gamma,
    zeta = zeta,
    ecp_grid = ecp_grid,
    zeta_grid = zeta_grid,
    target = 1 - alpha,
    n_valid = ncol(E_valid),
    n_grid = nrow(E_valid)
  )
}
