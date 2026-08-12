##################
# load R packages
##################

require(MASS)
require(fda)
require(ftsa)
require(orthopolynom)
require(doMC)
require(longmemo)

sim_ex2 <- function(n, d, seed_number, A_structure = c("diag", "band"), 
                    basis_structure = c("polynomials", "fourier"), basis_permute = c("TRUE", "FALSE"), 
                    cutt = 3, lbnumber = 21, no_discrete_points = 101, burnin_points = 100)
{
    basis_permute = match.arg(basis_permute)
    basis_structure = match.arg(basis_structure)
    set.seed(122 + seed_number)
    grid = seq(0, 1, length.out = no_discrete_points)
    
    # simulate non-stationary part

    Sigma_mat = matrix(c(1, 0.5, 0.5, 0.5, 1, 0.5, 0.5, 0.5, 1), 3, 3)
    step_1 = t(mvrnorm(n + burnin_points, mu = rep(0, 3), Sigma = Sigma_mat))

    # calculating the beta values

    beta_val = vector("numeric", n + burnin_points)
    for(i in 1:(n + burnin_points))
    {
        beta_val[i] = exp(lgamma(i + d) - lgamma(i + 1) - lgamma(d))
    }
    rm(i)

    eta_star = matrix(NA, cutt, n + burnin_points)
    for(t_val in 1:(n + burnin_points))
    {
        inner_sum = matrix(NA, cutt, t_val)
        index = t_val:1
        for(ik in 1:t_val)
        {
            inner_sum[,ik] = (beta_val[ik] * step_1[,index[ik]])
        }
        eta_star[,t_val] = rowSums(inner_sum)
        rm(index); rm(inner_sum)
    }

    # simulate stationary part

    no_cutt = lbnumber - cutt
    step_2 = t(mvrnorm(n + burnin_points, mu = rep(0, no_cutt), Sigma = diag(no_cutt)))

    if(A_structure == "band")
    {
        A = matrix(NA, no_cutt, no_cutt)
        for(ij in 1:no_cutt)
        {
            for(ik in 1:no_cutt)
            {
                if(abs(ij - ik) <= 3)
                {
                    A[ij,ik] = runif(1, min = -0.3, max = 0.3)
                }
                else
                {
                    A[ij,ik] = 0
                }
            }
        }
    }
    if(A_structure == "diag")
    {
        A = diag(runif(no_cutt, -0.5, 0.5))
    }
        
    eta_star_stationary = matrix(NA, no_cutt, n + burnin_points)
    eta_star_stationary[,1] = step_2[,1]
    for(iw in 1:(n + (burnin_points - 1)))
    {
        eta_star_stationary[,iw + 1] = A %*% eta_star_stationary[,iw] + step_2[,iw + 1]
    }

    # combine coefficients (long-memory + short_memory)

    coef_fn = rbind(eta_star, eta_star_stationary)
    
    if(basis_structure == "polynomials")
    {
        normalized.p.list <- legendre.polynomials(n = (lbnumber - 1), normalized = TRUE)
        basis_fn <- as.matrix(as.data.frame(polynomial.values(polynomials = normalized.p.list, x = grid)))
        if(basis_permute == "TRUE")
        {
            rpert = append(sample(1:cutt, size = length(1:cutt)), sample((cutt+1):lbnumber, size = length((cutt+1):lbnumber)))
        }
        if(basis_permute == "FALSE")
        {
            rpert = 1:lbnumber
        }
        eval_Xfd = basis_fn[,rpert] %*% coef_fn
    }
    if(basis_structure == "fourier")
    {
        basis_fn <- create.fourier.basis(rangeval = c(0, 1), nbasis = 21)
        Xfd <- fd(coef_fn, basis_fn)
        eval_Xfd <- eval.fd(grid, Xfd)
    }
    return(eval_Xfd[,(burnin_points + 1):(n + burnin_points)])
}


##########
# n = 500
##########

# A_structure = "diag"

registerDoMC(36)
ex2_first_d_0.75_n500_bandconst_1_const_tuning_0.4_diag = foreach(iwk = 1:1000) %dopar% H_estimate(fun_dat = sim_ex2(n = 500, d = 0.75, seed_number = iwk, A_structure = "diag", basis_permute = "TRUE", no_discrete_points = 101),
                                                                                             band_const = 1, tuning_const = 0.4, delta_val = 1e-2, 
                                                                                             choice_score = "first", taper_choice = "none")

ex2_first_d_1.0_n500_bandconst_1_const_tuning_0.4_diag = foreach(iwk = 1:1000) %dopar% H_estimate(fun_dat = sim_ex2(n = 500, d = 1.0, seed_number = iwk, A_structure = "diag", basis_permute = "TRUE", no_discrete_points = 101),
                                                                                             band_const = 1, tuning_const = 0.4, delta_val = 1e-2, 
                                                                                             choice_score = "first", taper_choice = "none")

ex2_first_d_1.25_n500_bandconst_1_const_tuning_0.4_diag = foreach(iwk = 1:1000) %dopar% H_estimate(fun_dat = sim_ex2(n = 500, d = 1.25, seed_number = iwk, A_structure = "diag", basis_permute = "TRUE", no_discrete_points = 101),
                                                                                             band_const = 1, tuning_const = 0.4, delta_val = 1e-2, 
                                                                                             choice_score = "first", taper_choice = "none")


ex2_first_d_1.25_n500_bandconst_1_const_tuning_0.4_diag_Velasco = foreach(iwk = 1:1000) %dopar% H_estimate(fun_dat = sim_ex2(n = 500, d = 1.25, seed_number = iwk, A_structure = "diag", basis_permute = "TRUE", no_discrete_points = 101),
                                                                                                   band_const = 1, tuning_const = 0.4, delta_val = 1e-2, 
                                                                                                   choice_score = "first", taper_choice = "Velasco")

ex2_first_d_1.25_n500_bandconst_1_const_tuning_0.4_diag_burnin500 = foreach(iwk = 1:1000) %dopar% H_estimate(fun_dat = sim_ex2(n = 500, d = 1.25, seed_number = iwk, A_structure = "diag", basis_permute = "TRUE", no_discrete_points = 101, burnin_points = 500),
                                                                                                   band_const = 1, tuning_const = 0.4, delta_val = 1e-2, 
                                                                                                   choice_score = "first", taper_choice = "none")



ex2_mat_first_d_0.75_n500_bandconst_1_const_tuning_0.4_diag = ex2_mat_first_d_1.0_n500_bandconst_1_const_tuning_0.4_diag = ex2_mat_first_d_1.25_n500_bandconst_1_const_tuning_0.4_diag = 
ex2_mat_first_d_1.25_n500_bandconst_1_const_tuning_0.4_diag_Velasco = ex2_mat_first_d_1.25_n500_bandconst_1_const_tuning_0.4_diag_burnin500 = matrix(NA, 1000, 3)
for(ik in 1:1000)
{
    dum = ex2_first_d_0.75_n500_bandconst_1_const_tuning_0.4_diag[[ik]]
    ex2_mat_first_d_0.75_n500_bandconst_1_const_tuning_0.4_diag[ik,] = c(dum$d_est, dum$ncomp, dum$diff_order)
    rm(dum) 
    
    dum = ex2_first_d_1.0_n500_bandconst_1_const_tuning_0.4_diag[[ik]]
    ex2_mat_first_d_1.0_n500_bandconst_1_const_tuning_0.4_diag[ik,] = c(dum$d_est, dum$ncomp, dum$diff_order)
    rm(dum)
    
    dum = ex2_first_d_1.25_n500_bandconst_1_const_tuning_0.4_diag[[ik]]
    ex2_mat_first_d_1.25_n500_bandconst_1_const_tuning_0.4_diag[ik,] = c(dum$d_est, dum$ncomp, dum$diff_order)
    rm(dum)
    
    dum = ex2_first_d_1.25_n500_bandconst_1_const_tuning_0.4_diag_Velasco[[ik]]
    ex2_mat_first_d_1.25_n500_bandconst_1_const_tuning_0.4_diag_Velasco[ik,] = c(dum$d_est, dum$ncomp, dum$diff_order)
    rm(dum)
    
    dum = ex2_first_d_1.25_n500_bandconst_1_const_tuning_0.4_diag_burnin500[[ik]]
    ex2_mat_first_d_1.25_n500_bandconst_1_const_tuning_0.4_diag_burnin500[ik,] = c(dum$d_est, dum$ncomp, dum$diff_order)
    rm(dum); print(ik); rm(ik)
}

polynomial_d_0.75_n500_diag_summary = round(c(mean(ex2_mat_first_d_0.75_n500_bandconst_1_const_tuning_0.4_diag[,1]), median(ex2_mat_first_d_0.75_n500_bandconst_1_const_tuning_0.4_diag[,1])), 4) # 0.7907 0.7965
polynomial_d_0.75_n500_diag_bias = round(mean(ex2_mat_first_d_0.75_n500_bandconst_1_const_tuning_0.4_diag[,1] - 0.75), 4)     # 0.0407
polynomial_d_0.75_n500_diag_variance = round(var(ex2_mat_first_d_0.75_n500_bandconst_1_const_tuning_0.4_diag[,1]), 4)             # 0.0049
polynomial_d_0.75_n500_diag_MSE = round(mean((ex2_mat_first_d_0.75_n500_bandconst_1_const_tuning_0.4_diag[,1] - 0.75)^2), 4) # 0.0066

polynomial_d_1.0_n500_diag_summary = round(c(mean(ex2_mat_first_d_1.0_n500_bandconst_1_const_tuning_0.4_diag[,1]), median(ex2_mat_first_d_1.0_n500_bandconst_1_const_tuning_0.4_diag[,1])), 4)   # 0.9584 0.9584
polynomial_d_1.0_n500_diag_bias = round(mean(ex2_mat_first_d_1.0_n500_bandconst_1_const_tuning_0.4_diag[,1] - 1.0), 4)     # -0.0416
polynomial_d_1.0_n500_diag_variance = round(var(ex2_mat_first_d_1.0_n500_bandconst_1_const_tuning_0.4_diag[,1]), 4)            # 0.0055
polynomial_d_1.0_n500_diag_MSE = round(mean((ex2_mat_first_d_1.0_n500_bandconst_1_const_tuning_0.4_diag[,1] - 1.0)^2), 4) # 0.0073

polynomial_d_1.25_n500_diag_summary = round(c(mean(ex2_mat_first_d_1.25_n500_bandconst_1_const_tuning_0.4_diag[,1]), median(ex2_mat_first_d_1.25_n500_bandconst_1_const_tuning_0.4_diag[,1])), 4) # 1.1695 1.1702
polynomial_d_1.25_n500_diag_bias = round(mean(ex2_mat_first_d_1.25_n500_bandconst_1_const_tuning_0.4_diag[,1] - 1.25), 4)     # -0.0805
polynomial_d_1.25_n500_diag_variance = round(var(ex2_mat_first_d_1.25_n500_bandconst_1_const_tuning_0.4_diag[,1]), 4)             # 0.0058
polynomial_d_1.25_n500_diag_MSE = round(mean((ex2_mat_first_d_1.25_n500_bandconst_1_const_tuning_0.4_diag[,1] - 1.25)^2), 4) # 0.0123

# round(c(mean(ex2_mat_first_d_1.25_n500_bandconst_1_const_tuning_0.4_diag_Velasco[,1]), median(ex2_mat_first_d_1.25_n500_bandconst_1_const_tuning_0.4_diag_Velasco[,1])), 4) # 1.1796 1.1840
# round(c(mean(ex2_mat_first_d_1.25_n500_bandconst_1_const_tuning_0.4_diag_burnin500[,1]), median(ex2_mat_first_d_1.25_n500_bandconst_1_const_tuning_0.4_diag_burnin500[,1])), 4) # 1.1655 1.1706

polynomial_d_0.75_n500_diag_table = table(ex2_mat_first_d_0.75_n500_bandconst_1_const_tuning_0.4_diag[,2]) # 451 547 1 1 
polynomial_d_1.0_n500_diag_table = table(ex2_mat_first_d_1.0_n500_bandconst_1_const_tuning_0.4_diag[,2])  # 237 760 3
polynomial_d_1.25_n500_diag_table = table(ex2_mat_first_d_1.25_n500_bandconst_1_const_tuning_0.4_diag[,2]) # 286 703 11

# table(ex2_mat_first_d_1.25_n500_bandconst_1_const_tuning_0.4_diag_burnin500[,2]) # 300 694 6


# A_structure = "band"

registerDoMC(36)
ex2_first_d_0.75_n500_bandconst_1_const_tuning_0.4 = foreach(iwk = 1:1000) %dopar% H_estimate(fun_dat = sim_ex2(n = 500, d = 0.75, seed_number = iwk, A_structure = "band", basis_permute = "TRUE"),
                                                                                             band_const = 1, tuning_const = 0.4, delta_val = 1e-2, 
                                                                                             choice_score = "first", taper_choice = "none")

ex2_first_d_1.0_n500_bandconst_1_const_tuning_0.4 = foreach(iwk = 1:1000) %dopar% H_estimate(fun_dat = sim_ex2(n = 500, d = 1.0, seed_number = iwk, A_structure = "band", basis_permute = "TRUE"),
                                                                                             band_const = 1, tuning_const = 0.4, delta_val = 1e-2, 
                                                                                             choice_score = "first", taper_choice = "none")

ex2_first_d_1.25_n500_bandconst_1_const_tuning_0.4 = foreach(iwk = 1:1000) %dopar% H_estimate(fun_dat = sim_ex2(n = 500, d = 1.25, seed_number = iwk, A_structure = "band", basis_permute = "TRUE"),
                                                                                             band_const = 1, tuning_const = 0.4, delta_val = 1e-2, 
                                                                                             choice_score = "first", taper_choice = "none")


ex2_mat_first_d_0.75_n500_bandconst_1_const_tuning_0.4 = ex2_mat_first_d_1.0_n500_bandconst_1_const_tuning_0.4 = ex2_mat_first_d_1.25_n500_bandconst_1_const_tuning_0.4 = matrix(NA, 1000, 3)
for(ik in 1:1000)
{
    dum = ex2_first_d_0.75_n500_bandconst_1_const_tuning_0.4[[ik]]
    ex2_mat_first_d_0.75_n500_bandconst_1_const_tuning_0.4[ik,] = c(dum$d_est, dum$ncomp, dum$diff_order)
    rm(dum) 
    
    dum = ex2_first_d_1.0_n500_bandconst_1_const_tuning_0.4[[ik]]
    ex2_mat_first_d_1.0_n500_bandconst_1_const_tuning_0.4[ik,] = c(dum$d_est, dum$ncomp, dum$diff_order)
    rm(dum)
    
    dum = ex2_first_d_1.25_n500_bandconst_1_const_tuning_0.4[[ik]]
    ex2_mat_first_d_1.25_n500_bandconst_1_const_tuning_0.4[ik,] = c(dum$d_est, dum$ncomp, dum$diff_order)
    rm(dum); print(ik); rm(ik)
}

polynomial_d_0.75_n500_band_summary = round(c(mean(ex2_mat_first_d_0.75_n500_bandconst_1_const_tuning_0.4[,1]), median(ex2_mat_first_d_0.75_n500_bandconst_1_const_tuning_0.4[,1])), 4) # 0.7886 0.7917
polynomial_d_0.75_n500_band_bias = round(mean(ex2_mat_first_d_0.75_n500_bandconst_1_const_tuning_0.4[,1] - 0.75), 4) # 0.0386
polynomial_d_0.75_n500_band_variance = round(var(ex2_mat_first_d_0.75_n500_bandconst_1_const_tuning_0.4[,1]), 4) # 0.0046
polynomial_d_0.75_n500_band_MSE = round(mean((ex2_mat_first_d_0.75_n500_bandconst_1_const_tuning_0.4[,1] - 0.75)^2), 4) # 0.0061

polynomial_d_1.0_n500_band_summary = round(c(mean(ex2_mat_first_d_1.0_n500_bandconst_1_const_tuning_0.4[,1]), median(ex2_mat_first_d_1.0_n500_bandconst_1_const_tuning_0.4[,1])), 4)   # 0.9555 0.9538
polynomial_d_1.0_n500_band_bias = round(mean(ex2_mat_first_d_1.0_n500_bandconst_1_const_tuning_0.4[,1] - 1.0), 4) # -0.0445
polynomial_d_1.0_n500_band_variance = round(var(ex2_mat_first_d_1.0_n500_bandconst_1_const_tuning_0.4[,1]), 4) # 0.0055
polynomial_d_1.0_n500_band_MSE = round(mean((ex2_mat_first_d_1.0_n500_bandconst_1_const_tuning_0.4[,1] - 1.0)^2), 4) # 0.0075

polynomial_d_1.25_n500_band_summary = round(c(mean(ex2_mat_first_d_1.25_n500_bandconst_1_const_tuning_0.4[,1]), median(ex2_mat_first_d_1.25_n500_bandconst_1_const_tuning_0.4[,1])), 4) # 1.1671 1.1688
polynomial_d_1.25_n500_band_bias = round(mean(ex2_mat_first_d_1.25_n500_bandconst_1_const_tuning_0.4[,1] - 1.25), 4) # -0.0829
polynomial_d_1.25_n500_band_variance = round(var(ex2_mat_first_d_1.25_n500_bandconst_1_const_tuning_0.4[,1]), 4) # 0.0057
polynomial_d_1.25_n500_band_MSE = round(mean((ex2_mat_first_d_1.25_n500_bandconst_1_const_tuning_0.4[,1] - 1.25)^2), 4) # 0.0126

polynomial_d_0.75_n500_band_table = table(ex2_mat_first_d_0.75_n500_bandconst_1_const_tuning_0.4[,2]) # 456 544
polynomial_d_1.0_n500_band_table = table(ex2_mat_first_d_1.0_n500_bandconst_1_const_tuning_0.4[,2])  # 239 759 2
polynomial_d_1.25_n500_band_table = table(ex2_mat_first_d_1.25_n500_bandconst_1_const_tuning_0.4[,2]) # 270 719 11

###########
# n = 1000
###########

# A_structure = "diag"

registerDoMC(24)
ex2_first_d_0.75_n1000_bandconst_1_const_tuning_0.4_diag = foreach(iwk = 1:1000) %dopar% H_estimate(fun_dat = sim_ex2(n = 1000, d = 0.75, seed_number = iwk, A_structure = "diag", basis_permute = "TRUE"),
                                                                                              band_const = 1, tuning_const = 0.4, delta_val = 1e-2,
                                                                                              choice_score = "first", taper_choice = "none")

ex2_first_d_1.0_n1000_bandconst_1_const_tuning_0.4_diag = foreach(iwk = 1:1000) %dopar% H_estimate(fun_dat = sim_ex2(n = 1000, d = 1.0, seed_number = iwk, A_structure = "diag", basis_permute = "TRUE"),
                                                                                              band_const = 1, tuning_const = 0.4, delta_val = 1e-2, 
                                                                                              choice_score = "first", taper_choice = "none")

ex2_first_d_1.25_n1000_bandconst_1_const_tuning_0.4_diag = foreach(iwk = 1:1000) %dopar% H_estimate(fun_dat = sim_ex2(n = 1000, d = 1.25, seed_number = iwk, A_structure = "diag", basis_permute = "TRUE"),
                                                                                              band_const = 1, tuning_const = 0.4, delta_val = 1e-2, 
                                                                                              choice_score = "first", taper_choice = "none")


ex2_mat_first_d_0.75_n1000_bandconst_1_const_tuning_0.4_diag = ex2_mat_first_d_1.0_n1000_bandconst_1_const_tuning_0.4_diag = ex2_mat_first_d_1.25_n1000_bandconst_1_const_tuning_0.4_diag = matrix(NA, 1000, 3)
for(ik in 1:1000)
{
    dum = ex2_first_d_0.75_n1000_bandconst_1_const_tuning_0.4_diag[[ik]] 
    ex2_mat_first_d_0.75_n1000_bandconst_1_const_tuning_0.4_diag[ik,] = c(dum$d_est, dum$ncomp, dum$diff_order)
    rm(dum)
    
    dum = ex2_first_d_1.0_n1000_bandconst_1_const_tuning_0.4_diag[[ik]] 
    ex2_mat_first_d_1.0_n1000_bandconst_1_const_tuning_0.4_diag[ik,] = c(dum$d_est, dum$ncomp, dum$diff_order)
    rm(dum)
    
    dum = ex2_first_d_1.25_n1000_bandconst_1_const_tuning_0.4_diag[[ik]] 
    ex2_mat_first_d_1.25_n1000_bandconst_1_const_tuning_0.4_diag[ik,] = c(dum$d_est, dum$ncomp, dum$diff_order)
    rm(dum); print(ik); rm(ik)
}

polynomial_d_0.75_n1000_diag_summary = round(c(mean(ex2_mat_first_d_0.75_n1000_bandconst_1_const_tuning_0.4_diag[,1]), median(ex2_mat_first_d_0.75_n1000_bandconst_1_const_tuning_0.4_diag[,1])), 4) # 0.7948 0.7957
polynomial_d_0.75_n1000_diag_bias = round(mean(ex2_mat_first_d_0.75_n1000_bandconst_1_const_tuning_0.4_diag[,1] - 0.75), 4) # 0.0448
polynomial_d_0.75_n1000_diag_variance = round(var(ex2_mat_first_d_0.75_n1000_bandconst_1_const_tuning_0.4_diag[,1]), 4) # 0.0032
polynomial_d_0.75_n1000_diag_MSE = round(mean((ex2_mat_first_d_0.75_n1000_bandconst_1_const_tuning_0.4_diag[,1] - 0.75)^2), 4) # 0.0052

polynomial_d_1.0_n1000_diag_summary = round(c(mean(ex2_mat_first_d_1.0_n1000_bandconst_1_const_tuning_0.4_diag[,1]),  median(ex2_mat_first_d_1.0_n1000_bandconst_1_const_tuning_0.4_diag[,1])), 4)  # 0.9703 0.9726
polynomial_d_1.0_n1000_diag_bias = round(mean(ex2_mat_first_d_1.0_n1000_bandconst_1_const_tuning_0.4_diag[,1] - 1.0), 4) # -0.0297
polynomial_d_1.0_n1000_diag_variance = round(var(ex2_mat_first_d_1.0_n1000_bandconst_1_const_tuning_0.4_diag[,1]), 4) # 0.0036
polynomial_d_1.0_n1000_diag_MSE = round(mean((ex2_mat_first_d_1.0_n1000_bandconst_1_const_tuning_0.4_diag[,1] - 1.0)^2), 4) # 0.0044

polynomial_d_1.25_n1000_diag_summary = round(c(mean(ex2_mat_first_d_1.25_n1000_bandconst_1_const_tuning_0.4_diag[,1]), median(ex2_mat_first_d_1.25_n1000_bandconst_1_const_tuning_0.4_diag[,1])), 4) # 1.1902 1.1920
polynomial_d_1.25_n1000_diag_bias = round(mean(ex2_mat_first_d_1.25_n1000_bandconst_1_const_tuning_0.4_diag[,1] - 1.25), 4) # -0.0598
polynomial_d_1.25_n1000_diag_variance = round(var(ex2_mat_first_d_1.25_n1000_bandconst_1_const_tuning_0.4_diag[,1]), 4) # 0.0035
polynomial_d_1.25_n1000_diag_MSE = round(mean((ex2_mat_first_d_1.25_n1000_bandconst_1_const_tuning_0.4_diag[,1] - 1.25)^2), 4) # 0.0071

polynomial_d_0.75_n1000_diag_table = table(ex2_mat_first_d_0.75_n1000_bandconst_1_const_tuning_0.4_diag[,2]) # 264 736
polynomial_d_1.0_n1000_diag_table = table(ex2_mat_first_d_1.0_n1000_bandconst_1_const_tuning_0.4_diag[,2])  # 152 846 2
polynomial_d_1.25_n1000_diag_table = table(ex2_mat_first_d_1.25_n1000_bandconst_1_const_tuning_0.4_diag[,2]) # 246 739 15

# A_structure = "band"

ex2_first_d_0.75_n1000_bandconst_1_const_tuning_0.4 = foreach(iwk = 1:1000) %dopar% H_estimate(fun_dat = sim_ex2(n = 1000, d = 0.75, seed_number = iwk, A_structure = "band", basis_permute = "TRUE"),
                                                                                             band_const = 1, tuning_const = 0.4, delta_val = 1e-2,
                                                                                             choice_score = "first", taper_choice = "none")

ex2_first_d_1.0_n1000_bandconst_1_const_tuning_0.4 = foreach(iwk = 1:1000) %dopar% H_estimate(fun_dat = sim_ex2(n = 1000, d = 1.0, seed_number = iwk, A_structure = "band", basis_permute = "TRUE"),
                                                                                             band_const = 1, tuning_const = 0.4, delta_val = 1e-2, 
                                                                                             choice_score = "first", taper_choice = "none")

ex2_first_d_1.25_n1000_bandconst_1_const_tuning_0.4 = foreach(iwk = 1:1000) %dopar% H_estimate(fun_dat = sim_ex2(n = 1000, d = 1.25, seed_number = iwk, A_structure = "band", basis_permute = "TRUE"),
                                                                                             band_const = 1, tuning_const = 0.4, delta_val = 1e-2, 
                                                                                             choice_score = "first", taper_choice = "none")

ex2_mat_first_d_0.75_n1000_bandconst_1_const_tuning_0.4 = ex2_mat_first_d_1.0_n1000_bandconst_1_const_tuning_0.4 = ex2_mat_first_d_1.25_n1000_bandconst_1_const_tuning_0.4 = matrix(NA, 1000, 3)
for(ik in 1:1000)
{
    dum = ex2_first_d_0.75_n1000_bandconst_1_const_tuning_0.4[[ik]]
    ex2_mat_first_d_0.75_n1000_bandconst_1_const_tuning_0.4[ik,] = c(dum$d_est, dum$ncomp, dum$diff_order)
    rm(dum) 
    
    dum = ex2_first_d_1.0_n1000_bandconst_1_const_tuning_0.4[[ik]]
    ex2_mat_first_d_1.0_n1000_bandconst_1_const_tuning_0.4[ik,] = c(dum$d_est, dum$ncomp, dum$diff_order)
    rm(dum) 
    
    dum = ex2_first_d_1.25_n1000_bandconst_1_const_tuning_0.4[[ik]]
    ex2_mat_first_d_1.25_n1000_bandconst_1_const_tuning_0.4[ik,] = c(dum$d_est, dum$ncomp, dum$diff_order)
    rm(dum); print(ik); rm(ik)
}

polynomial_d_0.75_n1000_band_summary = round(c(mean(ex2_mat_first_d_0.75_n1000_bandconst_1_const_tuning_0.4[,1]), median(ex2_mat_first_d_0.75_n1000_bandconst_1_const_tuning_0.4[,1])), 4) # 0.7937 0.7962
polynomial_d_0.75_n1000_band_bias = round(mean(ex2_mat_first_d_0.75_n1000_bandconst_1_const_tuning_0.4[,1] - 0.75), 4) # 0.0437
polynomial_d_0.75_n1000_band_variance = round(var(ex2_mat_first_d_0.75_n1000_bandconst_1_const_tuning_0.4[,1]), 4) # 0.0033
polynomial_d_0.75_n1000_band_MSE = round(mean((ex2_mat_first_d_0.75_n1000_bandconst_1_const_tuning_0.4[,1] - 0.75)^2), 4) # 0.0052

polynomial_d_1.0_n1000_band_summary = round(c(mean(ex2_mat_first_d_1.0_n1000_bandconst_1_const_tuning_0.4[,1]), median(ex2_mat_first_d_1.0_n1000_bandconst_1_const_tuning_0.4[,1])), 4)   # 0.9697 0.9721
polynomial_d_1.0_n1000_band_bias = round(mean(ex2_mat_first_d_1.0_n1000_bandconst_1_const_tuning_0.4[,1] - 1.0), 4) # -0.0303
polynomial_d_1.0_n1000_band_variance = round(var(ex2_mat_first_d_1.0_n1000_bandconst_1_const_tuning_0.4[,1]), 4) # 0.0035
polynomial_d_1.0_n1000_band_MSE = round(mean((ex2_mat_first_d_1.0_n1000_bandconst_1_const_tuning_0.4[,1] - 1.0)^2), 4) # 0.0044

polynomial_d_1.25_n1000_band_summary = round(c(mean(ex2_mat_first_d_1.25_n1000_bandconst_1_const_tuning_0.4[,1]), median(ex2_mat_first_d_1.25_n1000_bandconst_1_const_tuning_0.4[,1])), 4) # 1.1892 1.1896
polynomial_d_1.25_n1000_band_bias = round(mean(ex2_mat_first_d_1.25_n1000_bandconst_1_const_tuning_0.4[,1] - 1.25), 4) # -0.0608
polynomial_d_1.25_n1000_band_variance = round(var(ex2_mat_first_d_1.25_n1000_bandconst_1_const_tuning_0.4[,1]), 4) # 0.0035
polynomial_d_1.25_n1000_band_MSE = round(mean((ex2_mat_first_d_1.25_n1000_bandconst_1_const_tuning_0.4[,1] - 1.25)^2), 4) # 0.0072

polynomial_d_0.75_n1000_band_table = table(ex2_mat_first_d_0.75_n1000_bandconst_1_const_tuning_0.4[,2]) # 280 720
polynomial_d_1.0_n1000_band_table = table(ex2_mat_first_d_1.0_n1000_bandconst_1_const_tuning_0.4[,2])  # 158 839 3
polynomial_d_1.25_n1000_band_table = table(ex2_mat_first_d_1.25_n1000_bandconst_1_const_tuning_0.4[,2]) # 244 750 6


###########
# n = 2000
###########

# A_structure = "diag"

registerDoMC(36)
ex2_first_d_0.75_n2000_bandconst_1_const_tuning_0.4_diag = foreach(iwk = 1:1000) %dopar% H_estimate(fun_dat = sim_ex2(n = 2000, d = 0.75, seed_number = iwk, A_structure = "diag", basis_permute = "TRUE"),
                                                                                              band_const = 1, tuning_const = 0.4, delta_val = 1e-2, 
                                                                                              choice_score = "first", taper_choice = "none")

ex2_first_d_1.0_n2000_bandconst_1_const_tuning_0.4_diag = foreach(iwk = 1:1000) %dopar% H_estimate(fun_dat = sim_ex2(n = 2000, d = 1.0, seed_number = iwk, A_structure = "diag", basis_permute = "TRUE"),
                                                                                              band_const = 1, tuning_const = 0.4, delta_val = 1e-2, 
                                                                                              choice_score = "first", taper_choice = "none")

ex2_first_d_1.25_n2000_bandconst_1_const_tuning_0.4_diag = foreach(iwk = 1:1000) %dopar% H_estimate(fun_dat = sim_ex2(n = 2000, d = 1.25, seed_number = iwk, A_structure = "diag", basis_permute = "TRUE"),
                                                                                              band_const = 1, tuning_const = 0.4, delta_val = 1e-2, 
                                                                                              choice_score = "first", taper_choice = "none")

ex2_mat_first_d_0.75_n2000_bandconst_1_const_tuning_0.4_diag = ex2_mat_first_d_1.0_n2000_bandconst_1_const_tuning_0.4_diag = ex2_mat_first_d_1.25_n2000_bandconst_1_const_tuning_0.4_diag = matrix(NA, 1000, 3)
for(ik in 1:1000)
{
    dum = ex2_first_d_0.75_n2000_bandconst_1_const_tuning_0.4_diag[[ik]]
    ex2_mat_first_d_0.75_n2000_bandconst_1_const_tuning_0.4_diag[ik,] = c(dum$d_est, dum$ncomp, dum$diff_order)
    rm(dum) 
    
    dum = ex2_first_d_1.0_n2000_bandconst_1_const_tuning_0.4_diag[[ik]]
    ex2_mat_first_d_1.0_n2000_bandconst_1_const_tuning_0.4_diag[ik,] = c(dum$d_est, dum$ncomp, dum$diff_order)
    rm(dum)
    
    dum = ex2_first_d_1.25_n2000_bandconst_1_const_tuning_0.4_diag[[ik]]
    ex2_mat_first_d_1.25_n2000_bandconst_1_const_tuning_0.4_diag[ik,] = c(dum$d_est, dum$ncomp, dum$diff_order)
    rm(dum); print(ik); rm(ik)
}

polynomial_d_0.75_n2000_diag_summary = round(c(mean(ex2_mat_first_d_0.75_n2000_bandconst_1_const_tuning_0.4_diag[,1]), median(ex2_mat_first_d_0.75_n2000_bandconst_1_const_tuning_0.4_diag[,1])), 4) # 0.7924 0.7931
polynomial_d_0.75_n2000_diag_bias = round(mean(ex2_mat_first_d_0.75_n2000_bandconst_1_const_tuning_0.4_diag[,1] - 0.75), 4) # 0.0424
polynomial_d_0.75_n2000_diag_variance = round(var(ex2_mat_first_d_0.75_n2000_bandconst_1_const_tuning_0.4_diag[,1]), 4) # 0.0022
polynomial_d_0.75_n2000_diag_MSE = round(mean((ex2_mat_first_d_0.75_n2000_bandconst_1_const_tuning_0.4_diag[,1] - 0.75)^2), 4) # 0.0040

polynomial_d_1.0_n2000_diag_summary = round(c(mean(ex2_mat_first_d_1.0_n2000_bandconst_1_const_tuning_0.4_diag[,1]), median(ex2_mat_first_d_1.0_n2000_bandconst_1_const_tuning_0.4_diag[,1])), 4)   # 0.9816 0.9838
polynomial_d_1.0_n2000_diag_bias = round(mean(ex2_mat_first_d_1.0_n2000_bandconst_1_const_tuning_0.4_diag[,1] - 1.0), 4) # -0.0184
polynomial_d_1.0_n2000_diag_variance = round(var(ex2_mat_first_d_1.0_n2000_bandconst_1_const_tuning_0.4_diag[,1]), 4) # 0.0020
polynomial_d_1.0_n2000_diag_MSE = round(mean((ex2_mat_first_d_1.0_n2000_bandconst_1_const_tuning_0.4_diag[,1] - 1.0)^2), 4) # 0.0023

polynomial_d_1.25_n2000_diag_summary = round(c(mean(ex2_mat_first_d_1.25_n2000_bandconst_1_const_tuning_0.4_diag[,1]), median(ex2_mat_first_d_1.25_n2000_bandconst_1_const_tuning_0.4_diag[,1])), 4) # 1.2079 1.2095
polynomial_d_1.25_n2000_diag_bias = round(mean(ex2_mat_first_d_1.25_n2000_bandconst_1_const_tuning_0.4_diag[,1] - 1.25), 4) # -0.0421
polynomial_d_1.25_n2000_diag_variance = round(var(ex2_mat_first_d_1.25_n2000_bandconst_1_const_tuning_0.4_diag[,1]), 4) # 0.0019
polynomial_d_1.25_n2000_diag_MSE = round(mean((ex2_mat_first_d_1.25_n2000_bandconst_1_const_tuning_0.4_diag[,1] - 1.25)^2), 4) # 0.0037

polynomial_d_0.75_n2000_diag_table = table(ex2_mat_first_d_0.75_n2000_bandconst_1_const_tuning_0.4_diag[,2]) # 114 886
polynomial_d_1.0_n2000_diag_table = table(ex2_mat_first_d_1.0_n2000_bandconst_1_const_tuning_0.4_diag[,2])  # 111 886 3
polynomial_d_1.25_n2000_diag_table = table(ex2_mat_first_d_1.25_n2000_bandconst_1_const_tuning_0.4_diag[,2]) # 220 771 9


# A_structure = "band"

ex2_first_d_0.75_n2000_bandconst_1_const_tuning_0.4 = foreach(iwk = 1:1000) %dopar% H_estimate(fun_dat = sim_ex2(n = 2000, d = 0.75, seed_number = iwk, A_structure = "band", basis_permute = "TRUE"),
                                                                                              band_const = 1, tuning_const = 0.4, delta_val = 1e-2, 
                                                                                              choice_score = "first", taper_choice = "none")

ex2_first_d_1.0_n2000_bandconst_1_const_tuning_0.4 = foreach(iwk = 1:1000) %dopar% H_estimate(fun_dat = sim_ex2(n = 2000, d = 1.0, seed_number = iwk, A_structure = "band", basis_permute = "TRUE"),
                                                                                              band_const = 1, tuning_const = 0.4, delta_val = 1e-2, 
                                                                                              choice_score = "first", taper_choice = "none")

ex2_first_d_1.25_n2000_bandconst_1_const_tuning_0.4 = foreach(iwk = 1:1000) %dopar% H_estimate(fun_dat = sim_ex2(n = 2000, d = 1.25, seed_number = iwk, A_structure = "band", basis_permute = "TRUE"),
                                                                                              band_const = 1, tuning_const = 0.4, delta_val = 1e-2, 
                                                                                              choice_score = "first", taper_choice = "none")

ex2_mat_first_d_0.75_n2000_bandconst_1_const_tuning_0.4 = ex2_mat_first_d_1.0_n2000_bandconst_1_const_tuning_0.4 = ex2_mat_first_d_1.25_n2000_bandconst_1_const_tuning_0.4 = matrix(NA, 1000, 3)
for(ik in 1:1000)
{
    dum = ex2_first_d_0.75_n2000_bandconst_1_const_tuning_0.4[[ik]]
    ex2_mat_first_d_0.75_n2000_bandconst_1_const_tuning_0.4[ik,] = c(dum$d_est, dum$ncomp, dum$diff_order)
    rm(dum) 
    
    dum = ex2_first_d_1.0_n2000_bandconst_1_const_tuning_0.4[[ik]]
    ex2_mat_first_d_1.0_n2000_bandconst_1_const_tuning_0.4[ik,] = c(dum$d_est, dum$ncomp, dum$diff_order)
    rm(dum)
    
    dum = ex2_first_d_1.25_n2000_bandconst_1_const_tuning_0.4[[ik]]
    ex2_mat_first_d_1.25_n2000_bandconst_1_const_tuning_0.4[ik,] = c(dum$d_est, dum$ncomp, dum$diff_order)
    rm(dum); print(ik); rm(ik)
}

polynomial_d_0.75_n2000_band_summary = round(c(mean(ex2_mat_first_d_0.75_n2000_bandconst_1_const_tuning_0.4[,1]), median(ex2_mat_first_d_0.75_n2000_bandconst_1_const_tuning_0.4[,1])), 4) # 0.7924 0.7929
polynomial_d_0.75_n2000_band_bias = round(mean(ex2_mat_first_d_0.75_n2000_bandconst_1_const_tuning_0.4[,1] - 0.75), 4) # 0.0424
polynomial_d_0.75_n2000_band_variance = round(var(ex2_mat_first_d_0.75_n2000_bandconst_1_const_tuning_0.4[,1]), 4) # 0.0023
polynomial_d_0.75_n2000_band_MSE = round(mean((ex2_mat_first_d_0.75_n2000_bandconst_1_const_tuning_0.4[,1] - 0.75)^2), 4) # 0.0040

polynomial_d_1.0_n2000_band_summary = round(c(mean(ex2_mat_first_d_1.0_n2000_bandconst_1_const_tuning_0.4[,1]), median(ex2_mat_first_d_1.0_n2000_bandconst_1_const_tuning_0.4[,1])), 4)   # 0.9817 0.9849
polynomial_d_1.0_n2000_band_bias = round(mean(ex2_mat_first_d_1.0_n2000_bandconst_1_const_tuning_0.4[,1] - 1.0), 4) # -0.0183
polynomial_d_1.0_n2000_band_variance = round(var(ex2_mat_first_d_1.0_n2000_bandconst_1_const_tuning_0.4[,1]), 4) # 0.0021
polynomial_d_1.0_n2000_band_MSE = round(mean((ex2_mat_first_d_1.0_n2000_bandconst_1_const_tuning_0.4[,1] - 1.0)^2), 4) # 0.0024

polynomial_d_1.25_n2000_band_summary = round(c(mean(ex2_mat_first_d_1.25_n2000_bandconst_1_const_tuning_0.4[,1]), median(ex2_mat_first_d_1.25_n2000_bandconst_1_const_tuning_0.4[,1])), 4) # 1.2081 1.2093
polynomial_d_1.25_n2000_band_bias = round(mean(ex2_mat_first_d_1.25_n2000_bandconst_1_const_tuning_0.4[,1] - 1.25), 4) # -0.0419
polynomial_d_1.25_n2000_band_variance = round(var(ex2_mat_first_d_1.25_n2000_bandconst_1_const_tuning_0.4[,1]), 4) # 0.0019
polynomial_d_1.25_n2000_band_MSE = round(mean((ex2_mat_first_d_1.25_n2000_bandconst_1_const_tuning_0.4[,1] - 1.25)^2), 4) # 0.0036

polynomial_d_0.75_n2000_band_table = table(ex2_mat_first_d_0.75_n2000_bandconst_1_const_tuning_0.4[,2]) # 128 872
polynomial_d_1.0_n2000_band_table = table(ex2_mat_first_d_1.0_n2000_bandconst_1_const_tuning_0.4[,2])  # 111 882 7
polynomial_d_1.25_n2000_band_table = table(ex2_mat_first_d_1.25_n2000_bandconst_1_const_tuning_0.4[,2]) # 224 759 17

##########
# summary
##########

require(xtable)
xtable(rbind(c(polynomial_d_0.75_n500_diag_summary, polynomial_d_0.75_n500_diag_bias, polynomial_d_0.75_n500_diag_variance, polynomial_d_0.75_n500_diag_MSE),
c(polynomial_d_0.75_n1000_diag_summary, polynomial_d_0.75_n1000_diag_bias, polynomial_d_0.75_n1000_diag_variance, polynomial_d_0.75_n1000_diag_MSE),
c(polynomial_d_0.75_n2000_diag_summary, polynomial_d_0.75_n2000_diag_bias, polynomial_d_0.75_n2000_diag_variance, polynomial_d_0.75_n2000_diag_MSE)), digits = 4)

xtable(rbind(c(polynomial_d_1.0_n500_diag_summary,  polynomial_d_1.0_n500_diag_bias,  polynomial_d_1.0_n500_diag_variance,  polynomial_d_1.0_n500_diag_MSE),
c(polynomial_d_1.0_n1000_diag_summary,  polynomial_d_1.0_n1000_diag_bias,  polynomial_d_1.0_n1000_diag_variance,  polynomial_d_1.0_n1000_diag_MSE),
c(polynomial_d_1.0_n2000_diag_summary,  polynomial_d_1.0_n2000_diag_bias,  polynomial_d_1.0_n2000_diag_variance,  polynomial_d_1.0_n2000_diag_MSE)), digits = 4)

xtable(rbind(c(polynomial_d_1.25_n500_diag_summary, polynomial_d_1.25_n500_diag_bias, polynomial_d_1.25_n500_diag_variance, polynomial_d_1.25_n500_diag_MSE),
c(polynomial_d_1.25_n1000_diag_summary, polynomial_d_1.25_n1000_diag_bias, polynomial_d_1.25_n1000_diag_variance, polynomial_d_1.25_n1000_diag_MSE),
c(polynomial_d_1.25_n2000_diag_summary, polynomial_d_1.25_n2000_diag_bias, polynomial_d_1.25_n2000_diag_variance, polynomial_d_1.25_n2000_diag_MSE)), digits = 4)

xtable(rbind(c(polynomial_d_0.75_n500_band_summary, polynomial_d_0.75_n500_band_bias, polynomial_d_0.75_n500_band_variance, polynomial_d_0.75_n500_band_MSE),
c(polynomial_d_0.75_n1000_band_summary, polynomial_d_0.75_n1000_band_bias, polynomial_d_0.75_n1000_band_variance, polynomial_d_0.75_n1000_band_MSE),
c(polynomial_d_0.75_n2000_band_summary, polynomial_d_0.75_n2000_band_bias, polynomial_d_0.75_n2000_band_variance, polynomial_d_0.75_n2000_band_MSE)), digits = 4)

xtable(rbind(c(polynomial_d_1.0_n500_band_summary,  polynomial_d_1.0_n500_band_bias,  polynomial_d_1.0_n500_band_variance,  polynomial_d_1.0_n500_band_MSE),
c(polynomial_d_1.0_n1000_band_summary,  polynomial_d_1.0_n1000_band_bias,  polynomial_d_1.0_n1000_band_variance,  polynomial_d_1.0_n1000_band_MSE),
c(polynomial_d_1.0_n2000_band_summary,  polynomial_d_1.0_n2000_band_bias,  polynomial_d_1.0_n2000_band_variance,  polynomial_d_1.0_n2000_band_MSE)), digits = 4)

xtable(rbind(c(polynomial_d_1.25_n500_band_summary, polynomial_d_1.25_n500_band_bias, polynomial_d_1.25_n500_band_variance, polynomial_d_1.25_n500_band_MSE),
c(polynomial_d_1.25_n1000_band_summary, polynomial_d_1.25_n1000_band_bias, polynomial_d_1.25_n1000_band_variance, polynomial_d_1.25_n1000_band_MSE),
c(polynomial_d_1.25_n2000_band_summary, polynomial_d_1.25_n2000_band_bias, polynomial_d_1.25_n2000_band_variance, polynomial_d_1.25_n2000_band_MSE)), digits = 4)

# number of retained principal components

polynomial_d_0.75_n500_diag_table
polynomial_d_0.75_n1000_diag_table
polynomial_d_0.75_n2000_diag_table

polynomial_d_1.0_n500_diag_table
polynomial_d_1.0_n1000_diag_table
polynomial_d_1.0_n2000_diag_table

polynomial_d_1.25_n500_diag_table
polynomial_d_1.25_n1000_diag_table
polynomial_d_1.25_n2000_diag_table


polynomial_d_0.75_n500_band_table
polynomial_d_0.75_n1000_band_table
polynomial_d_0.75_n2000_band_table

polynomial_d_1.0_n500_band_table
polynomial_d_1.0_n1000_band_table
polynomial_d_1.0_n2000_band_table

polynomial_d_1.25_n500_band_table
polynomial_d_1.25_n1000_band_table
polynomial_d_1.25_n2000_band_table

