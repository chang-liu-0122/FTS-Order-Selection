calibrate_zeta_ecp <- function(E_valid,
                               alpha,
                               zeta_grid = seq(0.1, 50, length.out = 2496))
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

conformal_base_interval <- function(forecast, gamma) {
  lower <- forecast - gamma
  upper <- forecast + gamma

  list(lower = lower, upper = upper, width = upper - lower)
}

conformal_scores <- function(truth, forecast, gamma) {
  base <- conformal_base_interval(forecast, gamma)

  s <- pmax(base$lower - truth, truth - base$upper) / base$width

  sort(s[is.finite(s)])
}

conformal_quantile <- function(scores_sorted, alpha) {
  n <- length(scores_sorted)

  if (n == 0L) {
    return(NA_real_)
  }

  k <- ceiling((n + 1) * (1 - alpha))

  if (k > n) {
    return(Inf)
  }

  scores_sorted[k]
}

conformal_interval <- function(forecast, gamma, Q) {
  base <- conformal_base_interval(forecast, gamma)

  list(lower = base$lower - Q * base$width,
       upper = base$upper + Q * base$width)
}

interval_score_cells <- function(truth, lower, upper, alpha) {
  (upper - lower) +
    (2 / alpha) * (lower - truth) * (truth < lower) +
    (2 / alpha) * (truth - upper) * (truth > upper)
}

aci_gamma <- function(E) {
  E <- as.matrix(E)

  gamma <- apply(t(scale(t(E), center = TRUE, scale = FALSE)), 1, sd, na.rm = TRUE)

  gamma[!is.finite(gamma) | gamma <= 0] <- NA_real_

  gamma
}

aci_quantile <- function(scores_sorted, a_t) {
  n <- length(scores_sorted)

  if (n == 0L || !is.finite(a_t)) {
    return(NA_real_)
  }

  k <- ceiling((n + 1) * (1 - a_t))

  if (k < 1L) {
    return(-0.5)
  }

  scores_sorted[min(k, n)]
}

aci_run_sequence <- function(scores_sorted, gamma, forecast, truth, origins, h,
                             alpha, eta) {
  n_seq <- length(origins)
  err <- rep(NA_real_, n_seq)
  used <- rep(FALSE, n_seq)
  a_t <- alpha

  covered <- 0
  n_cells <- 0
  score_sum <- 0

  for (i in seq_len(n_seq)) {
    for (j in which(!used & origins + h <= origins[i])) {
      if (is.finite(err[j])) {
        a_t <- a_t + eta * (alpha - err[j])
      }

      used[j] <- TRUE
    }

    band <- conformal_interval(forecast[[i]], gamma, aci_quantile(scores_sorted, a_t))
    y <- truth[[i]]
    keep <- is.finite(y) & is.finite(band$lower) & is.finite(band$upper)

    if (!any(keep)) {
      next
    }

    inside <- y[keep] >= band$lower[keep] & y[keep] <= band$upper[keep]
    err[i] <- mean(!inside)

    covered <- covered + sum(inside)
    n_cells <- n_cells + sum(keep)
    score_sum <- score_sum +
      sum(interval_score_cells(y[keep], band$lower[keep], band$upper[keep], alpha))
  }

  c(covered = covered, n = n_cells, score_sum = score_sum)
}

aci_select_eta <- function(F_by_h, Y_by_h, origins_by_h, split_origin, alpha,
                           eta_grid = c(0, 0.005, 0.01, 0.02, 0.05, 0.1)) {
  eta_grid <- sort(unique(eta_grid))

  prepared <- list()

  for (h in seq_along(F_by_h)) {
    Fh <- F_by_h[[h]]
    Yh <- Y_by_h[[h]]
    o <- origins_by_h[[h]]

    if (is.null(Fh) || length(o) == 0L || ncol(Fh) != length(o)) {
      next
    }

    cal <- o < split_origin
    tun <- o >= split_origin

    if (sum(cal) < 2L || sum(tun) < 1L) {
      next
    }

    gamma <- aci_gamma(Yh[, cal, drop = FALSE] - Fh[, cal, drop = FALSE])
    scores <- conformal_scores(Yh[, cal, drop = FALSE], Fh[, cal, drop = FALSE], gamma)

    if (length(scores) == 0L) {
      next
    }

    prepared[[length(prepared) + 1L]] <- list(
      h = h,
      gamma = gamma,
      scores = scores,
      forecast = lapply(which(tun), function(c) Fh[, c]),
      truth = lapply(which(tun), function(c) Yh[, c]),
      origins = o[tun]
    )
  }

  criterion <- vapply(eta_grid, function(eta) {
    is_h <- vapply(prepared, function(x) {
      r <- aci_run_sequence(x$scores, x$gamma, x$forecast, x$truth,
                            x$origins, x$h, alpha, eta)
      if (r[["n"]] > 0) r[["score_sum"]] / r[["n"]] else NA_real_
    }, numeric(1))

    if (any(is.finite(is_h))) mean(is_h, na.rm = TRUE) else NA_real_
  }, numeric(1))

  if (!any(is.finite(criterion))) {
    return(list(eta = 0, criterion = criterion, eta_grid = eta_grid))
  }

  rounded <- round(criterion, 10)
  best <- which(rounded == min(rounded, na.rm = TRUE))[1]

  list(eta = eta_grid[best], criterion = criterion, eta_grid = eta_grid)
}

