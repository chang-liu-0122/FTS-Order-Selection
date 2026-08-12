lagMat = function(Z, k)
{
  if (k > 0) {
    do.call(what = cbind,
            args = c(rep(NA, times = k),
                     list(Z[,-(ncol(Z)-(k-1):0), drop=FALSE])))
  } else if (k < 0) {
    do.call(what = cbind,
            args = c(list(Z[,-(1:abs(k)), drop=FALSE]),
                     rep(NA, times = abs(k))))
  } else {
    Z
  }
}


lagMat_AR <- function(Z, k)
{
  if (k > 0) {
    cbind(matrix(NA, nrow(Z), k),
          Z[, 1:(ncol(Z)-k), drop=FALSE])
  } else {
    Z
  }
}



farma_estimation <- function(X, p = 1, q = 0, m = NULL, Npc = 1)
{
  n <- ncol(X)
  d <- nrow(X)

  if (p < 0 || q < 0) stop("p and q must be >= 0")
  if (p == 0 && q == 0) stop("At least one of p or q must be > 0")

  if (is.null(m)) m <- min(n - 2, 15 + p + q)


  # CENTER DATA
  Xm <- rowMeans(X)


  Xc <- t(scale(t(X), center = TRUE, scale = FALSE))

  # ---- FPCA ----
  Xe <- eigen(cov(t(Xc)))$vectors
  Npc <- min(Npc, d)

  # CASE 1: FAR(p)
  if (q == 0)
  {
    Xs <- crossprod(Xe[, 1:Npc, drop = FALSE], Xc)
    Ws <- do.call(rbind,
                  lapply(1:p, lagMat_AR,
                         Z = Xs[1:Npc, , drop = FALSE]))

    Y  <- Xc[, (p+1):n, drop = FALSE]
    Ws <- Ws[, (p+1):n, drop = FALSE]

    D <- Y %*% t(Ws) / (n - p)
    G <- Ws %*% t(Ws) / (n - p)

    A_hat <- D %*% solve(G)

    dim(A_hat) <- c(d, Npc, p)

    A_hat <- apply(A_hat, 3,
                   function(M) M %*% t(Xe[, 1:Npc, drop = FALSE]))

    dim(A_hat) <- c(d, d, p)

    # ---- recursion ----
    fitted <- matrix(NA_real_, d, n)
    resid  <- matrix(NA_real_, d, n)

    for (t in (p+1):n)
    {
      pred <- 0
      for (j in 1:p) {
        pred <- pred + A_hat[,,j] %*% Xc[, t-j]
      }

      fitted[, t] <- pred
      resid[, t]  <- Xc[, t] - pred
    }

    fitted_original <- fitted + Xm

    return(list(
      phi = A_hat,
      theta = NULL,
      fitted = fitted,
      fitted_original = fitted_original,
      residuals = resid,
      mean = Xm
    ))
  }

  # CASE 2: FARMA(p,q)

  # ---- Step 1: long AR approx ----
  Xs <- crossprod(Xe, Xc)

  Ws <- do.call(rbind,
                lapply(1:m, lagMat, Z = Xs[1:Npc, , drop=FALSE]))

  D <- cov(t(Xs), t(Ws), use = "pairwise.complete.obs")
  opshat <- D %*% solve(cov(t(Ws), use = "pairwise.complete.obs"))

  epshats <- Xs - opshat %*% Ws
  epshats[is.na(epshats)] <- 0

  # =======================
  # Step 2: estimate FARMA(p,q) operators (Phi, Theta)
  # =======================
  if (p == 0) {
    Ws <- do.call(rbind,
                  lapply(1:q, lagMat, Z = epshats[1:Npc, , drop=FALSE]))
  } else {
    Ws <- rbind(
      do.call(rbind, lapply(1:q, lagMat, Z = epshats[1:Npc, , drop=FALSE])),
      do.call(rbind, lapply(1:p, lagMat, Z = Xs[1:Npc, , drop=FALSE]))
    )
  }

  D <- cov(t(Xc), t(Ws), use = "pairwise.complete.obs")
  opshat <- D %*% solve(cov(t(Ws), use = "pairwise.complete.obs"))

  # reshape and lift back to function space
  dim(opshat) <- c(d, Npc, ncol(opshat) / Npc)
  opshat <- apply(opshat, 3, "%*%", t(Xe[, 1:Npc, drop=FALSE]))
  dim(opshat) <- c(d, d, p + q)

  thetahat <- opshat[, , 1:q, drop=FALSE]
  if (p > 0) {
    phihat <- opshat[, , (q + 1):(q + p), drop=FALSE]
  } else {
    phihat <- NULL
  }

  # =======================
  # Step 3: compute innovations/residuals via recursion
  # Model: X_t = sum Phi_j X_{t-j} + eps_t + sum Theta_k eps_{t-k}
  # =======================
  maxLag <- max(p, q)
  eps_hat <- matrix(0, nrow = d, ncol = n)
  fitted  <- matrix(NA_real_, nrow = d, ncol = n)

  # start recursion at t = maxLag + 1
  for (t in (maxLag + 1):n) {
    ar_part <- 0
    if (p > 0) {
      for (j in 1:p) {
        ar_part <- ar_part + phihat[, , j] %*% Xc[, t - j]
      }
    }

    ma_part <- 0
    if (q > 0) {
      for (k in 1:q) {
        ma_part <- ma_part + thetahat[, , k] %*% eps_hat[, t - k]
      }
    }

    # fitted mean (excluding current innovation)
    fitted[, t] <- ar_part + ma_part

    # innovation
    eps_hat[, t] <- Xc[, t] - fitted[, t]
  }

  eps_hat[, 1:maxLag] <- NA
  fitted[, 1:maxLag]  <- NA

  fitted_original <- fitted + Xm

  return(list(
    phi = phihat,
    theta = thetahat,
    fitted = fitted,
    fitted_original = fitted_original,
    residuals = eps_hat,
    mean = Xm
  ))
}
