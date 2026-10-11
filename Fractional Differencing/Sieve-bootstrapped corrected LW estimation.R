compute_bias_LPWN_multi <- function(
    X,
    d_current,
    m_const = 0.75,
    R_short_vec = c(1, 2),
    R_noise = 0,
    B = 100,
    burn_in = 200,
    memory_cov_type = "sample",
    memory_cov_band_const = 1
) {
  memory_cov_type <- match.arg(
    memory_cov_type,
    c("sample", "lrc_RS", "lrc_LRS")
  )

  K <- length(R_short_vec)

  X_diff <- t(FracDiff(t(X), d = d_current))

  X_bootstrap <- sieve_bootstrap(
    fun_dat = t(X_diff),
    ncomp_porder_selection = "CPV_AICC",
    CPV_percent = 0.95,
    VAR_type = "none",
    B = B,
    burn_in = burn_in
  )

  d_boot_mat <- matrix(NA_real_, nrow = B, ncol = K)

  for (b in seq_len(B)) {
    Xb_diff <- X_bootstrap$X_boot[, , b]

    Xb <- t(FracDiff(t(Xb_diff), d = -d_current))

    cov_mat <- get_covariance(
      Xb,
      type = memory_cov_type,
      band_const = memory_cov_band_const
    )

    lead <- leading_pc(Xb, cov_mat)
    scores <- lead$scores

    for (k in seq_len(K)) {
      d_boot_mat[b, k] <- tryCatch(
        LPWN_general(
          scores,
          m_const = m_const,
          R_short = R_short_vec[k],
          R_noise = R_noise
        ),
        error = function(e) NA_real_
      )
    }
  }

  bias_vec <- colMeans(d_boot_mat, na.rm = TRUE) - d_current
  sd_vec <- apply(d_boot_mat - d_current, 2, sd, na.rm = TRUE)

  list(
    bias = bias_vec,
    sd = sd_vec,
    boot = d_boot_mat
  )
}


estimate_d_BC_LPWN2 <- function(
    X,
    m_const = 0.75,
    R_noise = 0,
    B = 100,
    burn_in = 200,
    memory_cov_type = "sample",
    memory_cov_band_const = 1
) {
  memory_cov_type <- match.arg(
    memory_cov_type,
    c("sample", "lrc_RS", "lrc_LRS")
  )

  cov_mat <- get_covariance(
    X,
    type = memory_cov_type,
    band_const = memory_cov_band_const
  )

  lead <- leading_pc(X, cov_mat)
  scores <- lead$scores

  d_LPWN2 <- LPWN_general(
    scores,
    m_const = m_const,
    R_short = 2,
    R_noise = R_noise
  )

  out <- compute_bias_LPWN_multi(
    X = X,
    d_current = d_LPWN2,
    m_const = m_const,
    R_short_vec = 2,
    R_noise = R_noise,
    B = B,
    burn_in = burn_in,
    memory_cov_type = memory_cov_type,
    memory_cov_band_const = memory_cov_band_const
  )

  d_BC_LPWN2_unrestricted <- d_LPWN2 - out$bias[1]

  # The bias correction is not box constrained. Reject estimates below zero,
  # which fall outside the admissible range of the LPWN construction, while
  # retaining the unrestricted value for seed-level diagnostics.
  d_BC_LPWN2_rejected <-
    !is.finite(d_BC_LPWN2_unrestricted) || d_BC_LPWN2_unrestricted < 0

  d_BC_LPWN2 <- if (d_BC_LPWN2_rejected) {
    NA_real_
  } else {
    as.numeric(d_BC_LPWN2_unrestricted)
  }

  list(
    d_LPWN2 = d_LPWN2,
    d_BC_LPWN2 = d_BC_LPWN2,
    d_BC_LPWN2_unrestricted = as.numeric(d_BC_LPWN2_unrestricted),
    d_BC_LPWN2_rejected = d_BC_LPWN2_rejected,
    bias_LPWN2 = out$bias[1],
    sd_LPWN2 = out$sd[1],
    bootstrap_estimates = out$boot[, 1],
    memory_cov_type = memory_cov_type,
    memory_cov_band_const = memory_cov_band_const,
    leading_eigenfunction = lead$psi1,
    leading_scores = lead$scores
  )
}




estimate_d_LPWN_hybrid <- function(
     X,
     m_const = 0.75,
     R_noise = 0,
     B = 100,
     burn_in = 200
 ) {

   # FPCA
   Xc <- t(scale(t(X), center = TRUE, scale = FALSE))
   cov_mat <- get_covariance(X, type = "sample", band_const = 1)
   lead <- leading_pc(X, cov_mat)

   scores <- lead$scores

   # LW estimates
   d_LW <- local_whittle_general(scores, m_const)

   # Analytical estimates
   d_R1 <- LPWN_general(scores, m_const, R_short = 1, R_noise)
   d_R2 <- LPWN_general(scores, m_const, R_short = 2, R_noise)

   # ONE bootstrap for both
   out <- compute_bias_LPWN_multi(
     X,
     d_current = d_R1,
     m_const = m_const,
     R_short_vec = c(1, 2),
     R_noise = R_noise,
     B = B,
     burn_in = burn_in
   )

   # bias correction
   d_R1_bc <- d_R1 - out$bias[1]
   d_R2_bc <- d_R2 - out$bias[2]

   return(list(
     d_LW = d_LW,
     d_R1 = d_R1,
     d_R2 = d_R2,
     d_R1_bc = d_R1_bc,
     d_R2_bc = d_R2_bc
   ))
 }



 parallel_sim_LPWN_FAR <- function(n_rep = 100) {

   with_progress({

     seeds <- 1:n_rep
     p <- progressor(steps = length(seeds) * nrow(grid))

     results <- future_lapply(seeds, function(s) {

       set.seed(s)

       res_list <- lapply(1:nrow(grid), function(i) {

         out <- tryCatch({

           T_i  <- grid$T[i]
           ar_i <- grid$ar_strength[i]
           d_i  <- grid$d[i]

           sim <- sim_FARMA(
             n = T_i,
             d = d_i,
             no_grid = 101,
             seed_number = s * 1000 + i,
             p = 1,
             q = 0,
             ar_strength = ar_i,
             ma_strength = 0.5
           )

           X <- sim$X
           res_r <- estimate_d_LPWN_hybrid(X, B = 399)

           data.frame(
             T = T_i,
             ar_strength = ar_i,
             d_true = d_i,
             seed = s,
             d_LW = res_r$d_LW,
             d_R1 = res_r$d_R1,
             d_R2 = res_r$d_R2,
             d_R1_bc = res_r$d_R1_bc,
             d_R2_bc = res_r$d_R2_bc
           )

         }, error = function(e) {

         #   just return NA row if error
           data.frame(
             T = grid$T[i],
             ar_strength = grid$ar_strength[i],
             d_true = grid$d[i],
             seed = s,
             d_LW = NA,
             d_R1 = NA,
             d_R2 = NA,
             d_R1_bc = NA,
             d_R2_bc = NA
           )
         })

         p()
         out
       })

       do.call(rbind, res_list)

     })

     do.call(rbind, results)
   })
 }



parallel_sim_LPWN_FARMA <- function(n_rep = 100) {

  with_progress({

    seeds <- 1:n_rep
    p <- progressor(steps = length(seeds) * nrow(grid))

    results <- future_lapply(seeds, function(s) {

      set.seed(s)

      res_list <- lapply(1:nrow(grid), function(i) {

        out <- tryCatch({

          T_i  <- grid$T[i]
          ar_i <- grid$ar_strength[i]
          ma_i <- grid$ma_strength[i]
          d_i  <- grid$d[i]

          sim <- sim_FARMA(
            n = T_i,
            d = d_i,
            no_grid = 101,
            seed_number = s * 1000 + i,
            p = 1,
            q = 1,
            ar_strength = ar_i,
            ma_strength = ma_i
          )

          X <- sim$X
          res_r <- estimate_d_LPWN_hybrid(X, B = 399)

          data.frame(
            T = T_i,
            ar_strength = ar_i,
            ma_strength = ma_i,
            d_true = d_i,
            seed = s,
            d_LW = res_r$d_LW,
            d_R1 = res_r$d_R1,
            d_R2 = res_r$d_R2,
            d_R1_bc = res_r$d_R1_bc,
            d_R2_bc = res_r$d_R2_bc
          )

        }, error = function(e) {

          data.frame(
            T = grid$T[i],
            ar_strength = grid$ar_strength[i],
            ma_strength = grid$ma_strength[i],
            d_true = grid$d[i],
            seed = s,
            d_LW = NA,
            d_R1 = NA,
            d_R2 = NA,
            d_R1_bc = NA,
            d_R2_bc = NA
          )
        })

        p()
        out
      })

      do.call(rbind, res_list)

    })

    do.call(rbind, results)
  })
}



