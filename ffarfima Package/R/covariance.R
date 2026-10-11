cov_LRS = function(fun_dat, band_const = 1){

  T = ncol(fun_dat)
  p_dim = nrow(fun_dat)
  m = min(T - 1, p_dim) * band_const

  X_bar = rowMeans(fun_dat, na.rm = TRUE)
  center_dat = sweep(fun_dat, 1, X_bar)

  gamma_l <- function(lag, T)
  {
    gamma_lag_sum = 0
    if(lag >= 0)
    {
      for(ij in 1:(T-lag))
      {
        gamma_lag_sum = gamma_lag_sum + t(as.matrix(center_dat[,ij]) %*%         t(as.matrix(center_dat[,(ij+lag)])))
      }
    }
    else
    {
      for(ij in 1:(T - abs(lag)))
      {
        gamma_lag_sum = gamma_lag_sum + t(as.matrix(center_dat[,ij + abs(lag)]) %*% t(as.matrix(center_dat[,ij])))
      }
    }
    return(gamma_lag_sum/T)
  }

  C = 0
  index = seq(-m, m, by = 1)
  for(k in 1:length(index))
  {
    C = C + (m - abs(index[k])) * gamma_l(lag = index[k], T = T)
  }

  return(C)
}

get_covariance <- function(X, type, band_const = 1) {

  Tn <- ncol(X)
  Xc <- t(scale(t(X), center = TRUE, scale = FALSE))

  if (type == "sample") {
    return((Xc %*% t(Xc)) / Tn)
  }

  if (type == "lrc_RS") {
    return(ftsa::long_run_covariance_estimation(X))
  }

  if (type == "lrc_LRS") {
    return(cov_LRS(X, band_const = band_const))
  }

  if (type == "plrc_RS") {
    sample <- (Xc %*% t(Xc)) / Tn
    lrc <- ftsa::long_run_covariance_estimation(X)
    return(lrc - sample)
  }

  if (type == "plrc_LRS") {
    sample <- (Xc %*% t(Xc)) / Tn
    lrc <- cov_LRS(X, band_const = band_const)
    return(lrc - sample)
  }

  stop("Unknown covariance type")
}

leading_pc <- function(X, C) {
  Xc <- t(scale(t(X), center = TRUE, scale = FALSE))
  eig <- eigen(C)
  psi1 <- eig$vectors[, 1]

  xi <- as.numeric(t(psi1) %*% Xc)

  return(list(
    psi1 = psi1,
    scores = xi
  ))
}

