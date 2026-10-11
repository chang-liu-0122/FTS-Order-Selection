pengFit <- function(x,
                    levels = 50,
                    minnpts = 3,
                    cut.off = 10^c(0.7, 2.5),
                    method = c("mean", "median")) {

  x <- as.vector(x)
  n <- length(x)

  increment <- (log10(n / minnpts)) / levels
  M <- floor(10^((1:levels) * increment))
  M <- M[M > 2]

  stats <- match.fun(method[1])

  PENG <- numeric(length(M))

  for (k in seq_along(M)) {

    m <- M[k]
    nCols <- n %/% m

    X <- matrix(x[1:(m * nCols)], ncol = nCols)

    Y <- apply(X, 2, cumsum)

    V <- numeric(nCols)
    t_mat <- cbind(1, 1:m)

    for (i in 1:nCols) {
      res <- lm.fit(t_mat, Y[, i])$residuals
      V[i] <- var(res)
    }

    PENG[k] <- stats(V)
  }

  wt <- trunc((sign((M - cut.off[1]) * (cut.off[2] - M)) + 1) / 2)

  fit <- lsfit(log10(M), log10(PENG), wt)

  beta <- fit$coef[[2]]
  H <- beta / 2

  return(list(
    H = H,
    beta = beta,
    M = M,
    PENG = PENG
  ))
}

peng_general <- function(
    score,
    levels = 50,
    minnpts = 3,
    cut.off = 10^c(0.7, 2.5),
    method = "mean",
    tuning_const = 1,
    m_const = 0.65
) {

  get_d <- function(x) {
    fit <- pengFit(
      x,
      levels = levels,
      minnpts = minnpts,
      cut.off = cut.off,
      method = method
    )
    d <- fit$H - 0.5
    return(d)
  }

  tuning_parameter = floor(1 + length(score)^(m_const))
  n <- length(score)
  d_est <- get_d(score)

  gamma_n <- tuning_const / log(tuning_parameter)

  if (d_est < (1 - gamma_n)) {

    d_final <- d_est
    ijk <- 1

  } else {

    ijk <- 1

    while (d_est >= (1 - gamma_n)) {

      score <- diff(score, differences = 1)
      ijk <- ijk + 1

      d_est <- get_d(score)
    }

    d_final <- d_est + (ijk - 1)
  }

  return(d_final)
}

