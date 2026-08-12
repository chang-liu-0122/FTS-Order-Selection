library(ftsa)
library(progressr)



cov_LRS = function(fun_dat, band_const = 1){

  T = ncol(fun_dat)
  p_dim = nrow(fun_dat)
  m = min(T - 1, p_dim) * band_const

  X_bar = rowMeans(fun_dat, na.rm = TRUE)
  center_dat = sweep(fun_dat, 1, X_bar)

  # calculating long-run covariance for a given lag

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

  # calculating the sum of long-run covariance for all lags (m equals the total number of grid       points)

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

  # Step 2: scores
  xi <- as.numeric(t(psi1) %*% Xc)

  return(list(
    psi1 = psi1,
    scores = xi
  ))
}


estimate_d_with_cov <- function(X, cov_type, method = c("whittle", "peng")) {

  method <- match.arg(method)

  # compute only needed covariance
  C <- get_covariance(X, cov_type)

  pc <- leading_pc(X, C)
  scores <- pc$scores

  if (method == "whittle") {
    d_hat <- local_whittle_general(scores, m_const = 0.65)
  }

  if (method == "peng") {
    d_hat <- peng_general(scores, method = "mean")
  }

  return(as.numeric(d_hat))
}


parallel_sim_d_est_FAR <- function(n_rep = 100,
                                   cov_types = c("sample", "lrc_RS", "lrc_LRS")) {

  with_progress({

    seeds <- 1:n_rep
    p <- progressor(steps = length(seeds) * nrow(grid))

    results <- future_lapply(seeds, function(s) {

      set.seed(s)
      library(fracdiff)

      res_list <- lapply(1:nrow(grid), function(i) {

        T_i  <- grid$T[i]
        ar_i <- grid$ar_strength[i]
        d_i  <- grid$d[i]

        X <- sim_FARMA(
          n = T_i,
          d = d_i,
          no_grid = 101,
          seed_number = s,
          p = 1,
          q = 0,
          ar_strength = ar_i
        )$X

        # base output
        out <- data.frame(
          T = T_i,
          ar_strength = ar_i,
          d_true = d_i,
          seed = s
        )

        # ---- sample ----
        if ("sample" %in% cov_types) {
          out$sample_whittle <- estimate_d_with_cov(X, "sample", "whittle")
          out$sample_peng    <- estimate_d_with_cov(X, "sample", "peng")
        }

        # ---- LRC RS ----
        if ("lrc_RS" %in% cov_types) {
          out$lrc_RS_whittle <- estimate_d_with_cov(X, "lrc_RS", "whittle")
          out$lrc_RS_peng    <- estimate_d_with_cov(X, "lrc_RS", "peng")
        }

        # ---- LRC LRS ----
        if ("lrc_LRS" %in% cov_types) {
          out$lrc_LRS_whittle <- estimate_d_with_cov(X, "lrc_LRS", "whittle")
          out$lrc_LRS_peng    <- estimate_d_with_cov(X, "lrc_LRS", "peng")
        }

        # ---- punctured ----
        if ("plrc_RS" %in% cov_types) {
          out$plrc_RS_whittle <- estimate_d_with_cov(X, "plrc_RS", "whittle")
          out$plrc_RS_peng    <- estimate_d_with_cov(X, "plrc_RS", "peng")
        }

        if ("plrc_LRS" %in% cov_types) {
          out$plrc_LRS_whittle <- estimate_d_with_cov(X, "plrc_LRS", "whittle")
          out$plrc_LRS_peng    <- estimate_d_with_cov(X, "plrc_LRS", "peng")
        }

        p()
        out
      })

      do.call(rbind, res_list)
    })

    do.call(rbind, results)
  })
}




parallel_sim_d_est_FARMA <- function(n_rep = 100,
                                     cov_types = c("sample", "lrc_RS", "lrc_LRS")) {

  with_progress({

    seeds <- 1:n_rep
    p <- progressor(steps = length(seeds) * nrow(grid))

    results <- future_lapply(seeds, function(s) {

      set.seed(s)
      library(fracdiff)

      res_list <- lapply(1:nrow(grid), function(i) {

        T_i <- grid$T[i]
        d_i <- grid$d[i]
        strength_ar_i <- grid$strength_ar[i]
        strength_ma_i <- grid$strength_ma[i]

        X <- sim_FARMA(
          n = T_i,
          d = d_i,
          no_grid = 101,
          seed_number = s,
          p = 1,
          q = 1,
          ar_strength = strength_ar_i,
          ma_strength = strength_ma_i
        )$X

        # base output
        out <- data.frame(
          T = T_i,
          d_true = d_i,
          strength_ar = strength_ar_i,
          strength_ma = strength_ma_i,
          seed = s
        )

        # ---- sample ----
        if ("sample" %in% cov_types) {
          out$sample_whittle <- estimate_d_with_cov(X, "sample", "whittle")
          out$sample_peng    <- estimate_d_with_cov(X, "sample", "peng")
        }

        # ---- LRC RS ----
        if ("lrc_RS" %in% cov_types) {
          out$lrc_RS_whittle <- estimate_d_with_cov(X, "lrc_RS", "whittle")
          out$lrc_RS_peng    <- estimate_d_with_cov(X, "lrc_RS", "peng")
        }

        # ---- LRC LRS ----
        if ("lrc_LRS" %in% cov_types) {
          out$lrc_LRS_whittle <- estimate_d_with_cov(X, "lrc_LRS", "whittle")
          out$lrc_LRS_peng    <- estimate_d_with_cov(X, "lrc_LRS", "peng")
        }

        # ---- optional punctured ----
        if ("plrc_RS" %in% cov_types) {
          out$plrc_RS_whittle <- estimate_d_with_cov(X, "plrc_RS", "whittle")
          out$plrc_RS_peng    <- estimate_d_with_cov(X, "plrc_RS", "peng")
        }

        if ("plrc_LRS" %in% cov_types) {
          out$plrc_LRS_whittle <- estimate_d_with_cov(X, "plrc_LRS", "whittle")
          out$plrc_LRS_peng    <- estimate_d_with_cov(X, "plrc_LRS", "peng")
        }

        p()
        out
      })

      do.call(rbind, res_list)
    })

    do.call(rbind, results)
  })
}
