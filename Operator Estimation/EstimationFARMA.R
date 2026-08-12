################################################################################
#
# Estimation of functional ARMA models
# - Implementation and comparison of methods by Kuenzer and Kuehnert
# - All functional time series are p x n matrices
# - Without any guarantee that there are no errors
#
# Author:  Thomas Kuenzer
# Version: 2024-12-17
#
############################################################################## #

# Method: Kuenzer ---------------------------------------------------------

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





# Method by Kuenzer
hannanrissanen = function (X0, p = 1, q = 1, m = 15 + p + q, Npc = 1)
{
  if (q == 0) stop("No pure AR process allowed.")

  n = ncol(X0)
  d = nrow(X0)

  Npc = min(Npc, d)

  m = min(m, n-2)

  Xm = rowMeans(X0)

  X0 = X0 - Xm

  Xe = eigen(cov(t(X0)))$vectors
  Xs = crossprod(Xe, X0)

  # Step 1

  Ws = do.call(what = rbind,
               args = lapply(1:m, lagMat, Z = Xs[1:Npc, , drop=FALSE]))

  D = cov(t(Xs), t(Ws), use = "pairwise.complete.obs")

  opshat = D %*% solve(cov(t(Ws), use = "pairwise.complete.obs"))

  epshats = Xs - opshat %*% Ws
  epshats[is.na(epshats)] = 0

  # Step 2

  if (p == 0) {
    Ws = do.call(what = rbind,
                 args = lapply(1:q, lagMat, Z = epshats[1:Npc, , drop=FALSE]))
  } else {
    Ws = rbind(do.call(what = rbind,
                       args = lapply(1:q, lagMat, Z = epshats[1:Npc, , drop=FALSE])),
               do.call(what = rbind,
                       args = lapply(1:p, lagMat, Z = Xs[1:Npc, , drop=FALSE])))
  }

  D = cov(t(X0), t(Ws), use = "pairwise.complete.obs")

  opshat = D %*% solve(cov(t(Ws), use = "pairwise.complete.obs"))

  dim(opshat) = c(d, Npc, ncol(opshat)/Npc)
  opshat = apply(opshat, 3, "%*%", t(Xe[, 1:Npc , drop=FALSE]))
  dim(opshat) = c(d, d, p+q)

  thetahat = opshat[,,1:q, drop=FALSE]
  if (p > 0) {
    phihat = opshat[,,q + 1:p, drop=FALSE]
  } else {
    phihat = NULL
  }

  list(phi = phihat, theta = thetahat)
}




# Method: Kuehnert --------------------------------------------------------

blockmatrix = function (A, inds)
{
  M = apply(inds, 1, function(i) A[,,i])
  dim(M) = c(nrow(A), ncol(A)*ncol(inds), nrow(inds))
  M = aperm(M, c(1,3,2))
  dim(M) = c(nrow(M)*ncol(M), dim(M)[3])
  M
}

pseudoinverse = function (A, Npc, delta = NULL, delta1 = NULL, deltaK = 1)
{
  eigs = eigen(tcrossprod(A))
  if (Npc < 1) {
    Npc = min(which(cumsum(eigs$values) >= Npc * sum(eigs$values)),
              length(eigs$values))
  }
  if (is.null(delta)) {
    if (!is.null(delta1)) {
      delta = eigs$values[1] * delta1
    } else {
      delta = eigs$values[Npc] * deltaK
    }
  }
  evs = eigs$vectors[,1:Npc, drop=FALSE]
  B = crossprod(evs, A) # t(evs) %*% A
  # t(A) %*% evs %*% solve(B%*%t(B) + diag(delta, Npc)) %*% t(evs)
  crossprod(A, evs) %*% tcrossprod(solve(tcrossprod(B) + diag(delta, Npc)), evs)
}




kuehnert = function (X0, p = 1, q = 1,
                     L = 15+p+2*q,     # max lag for estimation of Pi
                     K = NULL,   #Npc1 # Npc for inversion of Gamma
                     vn = NULL,        # -"-
                     M = 1,      #Npc2 # Npc for inversion of Pi
                     thetan = NULL)
{
  if (q == 0) stop("No pure AR process allowed.")

  n = ncol(X0)
  d = nrow(X0)

  s = p+q

  L = max(L, s+q-1)

  if (is.null(K)) {
    # K = ceiling(10* (n/L)^(1/3))
    K = ceiling(L*M)
  }

  K = min(K, d*L)

  M = min(M, q*d)

  Xm = rowMeans(X0)

  X0 = X0 - Xm

  Ws = do.call(what = rbind,
               args = lapply(1:L, lagMat, Z = X0))

  SL = cov(t(Ws), use = "pairwise.complete.obs")

  # scalef = norm(SL[1:d,1:d], "2")
  # if (is.null(vn)) {
  #   vn = scalef^2 / sqrt(n)
  # }
  # if (is.null(thetan)) {
  #   thetan = scalef^2 / sqrt(n)
  # }

  SLi = pseudoinverse(A = SL, Npc = K, delta = vn, deltaK = 1/K)

  SLK = cov(t(X0), t(Ws), use = "pairwise.complete.obs")
  PiLK = SLK %*% SLi

  dim(PiLK) = c(nrow(PiLK), nrow(PiLK), L)

  pisq = PiLK[,,(s+q-1):s, drop=FALSE]
  dim(pisq) = c(nrow(pisq), ncol(pisq)*dim(pisq)[3])

  if (p > 0) {
    pinds = outer(q:1, q:1, "+") + p - 2
    Pisq = blockmatrix(A = PiLK, inds = pinds)
    Pisqi = pseudoinverse(A = Pisq, Npc = M, delta = thetan, delta1 = 1/sqrt(K))

    beta = -pisq %*% Pisqi
    dim(beta) = c(d,d, q)

    alpha = PiLK[,,1:p, drop=FALSE]
    for (i in 1:p) {
      if (i <= q)  {
        alpha[,,i] = alpha[,,i] - beta[,,i]
      }
      if (i > 1) {
        for (j in 1:min(i-1,q)) {
          alpha[,,i] = alpha[,,i] + beta[,,j] %*% PiLK[,,i-j]
        }
      }
    }
  } else {
    if (q == 1) {
      beta = pisq
      dim(beta) = c(d,d, q)
      alpha = NULL
    } else {
      stop("p = 0 and q > 1 not implemented")
    }
  }

  list(phi = alpha, theta = beta, pi = PiLK)
}




farma_estimation_general <- function(
    X0, p = 1, q = 1,
    L = 15 + p + 2*q,
    K = NULL,
    vn = NULL,
    M = 1,
    thetan = NULL,
    p_large = 10
) {

  n <- ncol(X0)
  d <- nrow(X0)

  lagMat <- function(Z, k) {
    Z[, (k+1):n, drop = FALSE]
  }

  Xm <- rowMeans(X0)
  X0c <- X0 - Xm

  ############################################################
  # CASE 1: FAR(p)
  ############################################################
  if (q == 0) {

    Ws <- do.call(rbind, lapply(1:p, function(k) lagMat(X0c, k)))
    Y  <- X0c[, (p+1):n, drop = FALSE]
    Ws <- Ws[, 1:(n-p), drop = FALSE]

    Gamma <- cov(t(Ws))
    gamma <- cov(t(Y), t(Ws))

    alpha <- gamma %*% solve(Gamma)
    dim(alpha) <- c(d, d, p)

    # residuals
    E <- matrix(0, d, n)
    for (t in (p+1):n) {
      pred <- 0
      for (i in 1:p) {
        pred <- pred + alpha[,,i] %*% X0c[,t-i]
      }
      E[,t] <- X0c[,t] - pred
    }

    return(list(phi = alpha, theta = NULL, residuals = E))
  }

  ############################################################
  # CASE 2: MA(q)
  ############################################################
  if (p == 0 && q > 1) {

    # Use large AR approximation
    Ws <- do.call(rbind, lapply(1:p_large, function(k) lagMat(X0c, k)))
    Y  <- X0c[, (p_large+1):n, drop = FALSE]
    Ws <- Ws[, 1:(n-p_large), drop = FALSE]

    Gamma <- cov(t(Ws))
    gamma <- cov(t(Y), t(Ws))

    Pi <- gamma %*% solve(Gamma)
    dim(Pi) <- c(d, d, p_large)

    beta <- array(0, dim = c(d, d, q))

    for (i in 1:q) {
      beta[,,i] <- -Pi[,,i]
      if (i > 1) {
        for (j in 1:(i-1)) {
          beta[,,i] <- beta[,,i] - beta[,,j] %*% Pi[,,i-j]
        }
      }
    }

    # residuals
    E <- matrix(0, d, n)
    for (t in (q+1):n) {
      pred <- 0
      for (j in 1:q) {
        pred <- pred + beta[,,j] %*% E[,t-j]
      }
      E[,t] <- X0c[,t] - pred
    }

    return(list(phi = NULL, theta = beta, residuals = E))
  }

  ############################################################
  # CASE 3: FARMA(p,q)
  ############################################################
  fit <- kuehnert(
    X0 = X0,
    p = p,
    q = q,
    L = L,
    K = K,
    vn = vn,
    M = M,
    thetan = thetan
  )


  return(list(phi = alpha, theta = beta))
}
