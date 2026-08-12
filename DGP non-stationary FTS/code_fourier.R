#############################
# basis function = "fourier"
#############################

# n = 500

registerDoMC(36)
ex2_first_d_0.75_n500_bandconst_1_const_tuning_0.4_diag_fourier = foreach(iwk = 1:1000) %dopar% H_estimate(fun_dat = sim_ex2(n = 500, d = 0.75, seed_number = iwk, A_structure = "diag", basis_structure = "fourier"),
                                                                                                          band_const = 1, tuning_const = 0.4, delta_val = 1e-2, 
                                                                                                          choice_score = "first", taper_choice = "none")

ex2_first_d_1.0_n500_bandconst_1_const_tuning_0.4_diag_fourier = foreach(iwk = 1:1000) %dopar% H_estimate(fun_dat = sim_ex2(n = 500, d = 1.0, seed_number = iwk, A_structure = "diag", basis_structure = "fourier"),
                                                                                                          band_const = 1, tuning_const = 0.4, delta_val = 1e-2, 
                                                                                                          choice_score = "first", taper_choice = "none")

ex2_first_d_1.25_n500_bandconst_1_const_tuning_0.4_diag_fourier = foreach(iwk = 1:1000) %dopar% H_estimate(fun_dat = sim_ex2(n = 500, d = 1.25, seed_number = iwk, A_structure = "diag", basis_structure = "fourier"),
                                                                                                          band_const = 1, tuning_const = 0.4, delta_val = 1e-2, 
                                                                                                          choice_score = "first", taper_choice = "none")


# n = 1000

registerDoMC(36)
ex2_first_d_0.75_n1000_bandconst_1_const_tuning_0.4_diag_fourier = foreach(iwk = 1:1000) %dopar% H_estimate(fun_dat = sim_ex2(n = 1000, d = 0.75, seed_number = iwk, A_structure = "diag", basis_structure = "fourier"),
                                                                                                           band_const = 1, tuning_const = 0.4, delta_val = 1e-2, 
                                                                                                           choice_score = "first", taper_choice = "none")

ex2_first_d_1.0_n1000_bandconst_1_const_tuning_0.4_diag_fourier = foreach(iwk = 1:1000) %dopar% H_estimate(fun_dat = sim_ex2(n = 1000, d = 1.0, seed_number = iwk, A_structure = "diag", basis_structure = "fourier"),
                                                                                                           band_const = 1, tuning_const = 0.4, delta_val = 1e-2, 
                                                                                                           choice_score = "first", taper_choice = "none")

ex2_first_d_1.25_n1000_bandconst_1_const_tuning_0.4_diag_fourier = foreach(iwk = 1:1000) %dopar% H_estimate(fun_dat = sim_ex2(n = 1000, d = 1.25, seed_number = iwk, A_structure = "diag", basis_structure = "fourier"),
                                                                                                           band_const = 1, tuning_const = 0.4, delta_val = 1e-2, 
                                                                                                           choice_score = "first", taper_choice = "none")

# n = 2000

registerDoMC(36)
ex2_first_d_0.75_n2000_bandconst_1_const_tuning_0.4_diag_fourier = foreach(iwk = 1:1000) %dopar% H_estimate(fun_dat = sim_ex2(n = 2000, d = 0.75, seed_number = iwk, A_structure = "diag", basis_structure = "fourier"),
                                                                                                           band_const = 1, tuning_const = 0.4, delta_val = 1e-2, 
                                                                                                           choice_score = "first", taper_choice = "none")

ex2_first_d_1.0_n2000_bandconst_1_const_tuning_0.4_diag_fourier = foreach(iwk = 1:1000) %dopar% H_estimate(fun_dat = sim_ex2(n = 2000, d = 1.0, seed_number = iwk, A_structure = "diag", basis_structure = "fourier"),
                                                                                                           band_const = 1, tuning_const = 0.4, delta_val = 1e-2, 
                                                                                                           choice_score = "first", taper_choice = "none")

ex2_first_d_1.25_n2000_bandconst_1_const_tuning_0.4_diag_fourier = foreach(iwk = 1:1000) %dopar% H_estimate(fun_dat = sim_ex2(n = 2000, d = 1.25, seed_number = iwk, A_structure = "diag", basis_structure = "fourier"),
                                                                                                           band_const = 1, tuning_const = 0.4, delta_val = 1e-2, 
                                                                                                           choice_score = "first", taper_choice = "none")

# n = 500

ex2_mat_first_d_0.75_n500_bandconst_1_const_tuning_0.4_diag_no_fourier = ex2_mat_first_d_1.0_n500_bandconst_1_const_tuning_0.4_diag_no_fourier = ex2_mat_first_d_1.25_n500_bandconst_1_const_tuning_0.4_diag_no_fourier = matrix(NA, 1000, 3)
for(ik in 1:1000)
{
  dum = ex2_first_d_0.75_n500_bandconst_1_const_tuning_0.4_diag_fourier[[ik]]
  ex2_mat_first_d_0.75_n500_bandconst_1_const_tuning_0.4_diag_no_fourier[ik,] = c(dum$d_est, dum$ncomp, dum$diff_order)
  rm(dum) 
  
  dum = ex2_first_d_1.0_n500_bandconst_1_const_tuning_0.4_diag_fourier[[ik]]
  ex2_mat_first_d_1.0_n500_bandconst_1_const_tuning_0.4_diag_no_fourier[ik,] = c(dum$d_est, dum$ncomp, dum$diff_order)
  rm(dum)
  
  dum = ex2_first_d_1.25_n500_bandconst_1_const_tuning_0.4_diag_fourier[[ik]]
  ex2_mat_first_d_1.25_n500_bandconst_1_const_tuning_0.4_diag_no_fourier[ik,] = c(dum$d_est, dum$ncomp, dum$diff_order)
  rm(dum); print(ik); rm(ik)
}

round(c(mean(ex2_mat_first_d_0.75_n500_bandconst_1_const_tuning_0.4_diag_no_fourier[,1]), median(ex2_mat_first_d_0.75_n500_bandconst_1_const_tuning_0.4_diag_no_fourier[,1])), 4) # 
round(c(mean(ex2_mat_first_d_1.0_n500_bandconst_1_const_tuning_0.4_diag_no_fourier[,1]), median(ex2_mat_first_d_1.0_n500_bandconst_1_const_tuning_0.4_diag_no_fourier[,1])), 4)   # 
round(c(mean(ex2_mat_first_d_1.25_n500_bandconst_1_const_tuning_0.4_diag_no_fourier[,1]), median(ex2_mat_first_d_1.25_n500_bandconst_1_const_tuning_0.4_diag_no_fourier[,1])), 4) # 

table(ex2_mat_first_d_0.75_n500_bandconst_1_const_tuning_0.4_diag_no_fourier[,2]) # 248 251 501
table(ex2_mat_first_d_1.0_n500_bandconst_1_const_tuning_0.4_diag_no_fourier[,2])  # 155 131 714 
table(ex2_mat_first_d_1.25_n500_bandconst_1_const_tuning_0.4_diag_no_fourier[,2]) # 326 188 486 

# n = 1000

ex2_mat_first_d_0.75_n1000_bandconst_1_const_tuning_0.4_diag_no_fourier = ex2_mat_first_d_1.0_n1000_bandconst_1_const_tuning_0.4_diag_no_fourier = ex2_mat_first_d_1.25_n1000_bandconst_1_const_tuning_0.4_diag_no_fourier = matrix(NA, 1000, 3)
for(ik in 1:1000)
{
  dum = ex2_first_d_0.75_n1000_bandconst_1_const_tuning_0.4_diag_fourier[[ik]]
  ex2_mat_first_d_0.75_n1000_bandconst_1_const_tuning_0.4_diag_no_fourier[ik,] = c(dum$d_est, dum$ncomp, dum$diff_order)
  rm(dum) 
  
  dum = ex2_first_d_1.0_n1000_bandconst_1_const_tuning_0.4_diag_fourier[[ik]]
  ex2_mat_first_d_1.0_n1000_bandconst_1_const_tuning_0.4_diag_no_fourier[ik,] = c(dum$d_est, dum$ncomp, dum$diff_order)
  rm(dum)
  
  dum = ex2_first_d_1.25_n1000_bandconst_1_const_tuning_0.4_diag_fourier[[ik]]
  ex2_mat_first_d_1.25_n1000_bandconst_1_const_tuning_0.4_diag_no_fourier[ik,] = c(dum$d_est, dum$ncomp, dum$diff_order)
  rm(dum); print(ik); rm(ik)
}

round(c(mean(ex2_mat_first_d_0.75_n1000_bandconst_1_const_tuning_0.4_diag_no_fourier[,1]), median(ex2_mat_first_d_0.75_n1000_bandconst_1_const_tuning_0.4_diag_no_fourier[,1])), 4) # 
round(c(mean(ex2_mat_first_d_1.0_n1000_bandconst_1_const_tuning_0.4_diag_no_fourier[,1]), median(ex2_mat_first_d_1.0_n1000_bandconst_1_const_tuning_0.4_diag_no_fourier[,1])), 4)   # 
round(c(mean(ex2_mat_first_d_1.25_n1000_bandconst_1_const_tuning_0.4_diag_no_fourier[,1]), median(ex2_mat_first_d_1.25_n1000_bandconst_1_const_tuning_0.4_diag_no_fourier[,1])), 4) # 

table(ex2_mat_first_d_0.75_n1000_bandconst_1_const_tuning_0.4_diag_no_fourier[,2]) # 79 38 883
table(ex2_mat_first_d_1.0_n1000_bandconst_1_const_tuning_0.4_diag_no_fourier[,2])  # 85 26 889  
table(ex2_mat_first_d_1.25_n1000_bandconst_1_const_tuning_0.4_diag_no_fourier[,2]) # 243 99 658  

# n = 2000

ex2_mat_first_d_0.75_n2000_bandconst_1_const_tuning_0.4_diag_no_fourier = ex2_mat_first_d_1.0_n2000_bandconst_1_const_tuning_0.4_diag_no_fourier = ex2_mat_first_d_1.25_n2000_bandconst_1_const_tuning_0.4_diag_no_fourier = matrix(NA, 1000, 3)
for(ik in 1:1000)
{
  dum = ex2_first_d_0.75_n2000_bandconst_1_const_tuning_0.4_diag_fourier[[ik]]
  ex2_mat_first_d_0.75_n2000_bandconst_1_const_tuning_0.4_diag_no_fourier[ik,] = c(dum$d_est, dum$ncomp, dum$diff_order)
  rm(dum) 
  
  dum = ex2_first_d_1.0_n2000_bandconst_1_const_tuning_0.4_diag_fourier[[ik]]
  ex2_mat_first_d_1.0_n2000_bandconst_1_const_tuning_0.4_diag_no_fourier[ik,] = c(dum$d_est, dum$ncomp, dum$diff_order)
  rm(dum)
  
  dum = ex2_first_d_1.25_n2000_bandconst_1_const_tuning_0.4_diag_fourier[[ik]]
  ex2_mat_first_d_1.25_n2000_bandconst_1_const_tuning_0.4_diag_no_fourier[ik,] = c(dum$d_est, dum$ncomp, dum$diff_order)
  rm(dum); print(ik); rm(ik)
}

round(c(mean(ex2_mat_first_d_0.75_n2000_bandconst_1_const_tuning_0.4_diag_no_fourier[,1]), median(ex2_mat_first_d_0.75_n2000_bandconst_1_const_tuning_0.4_diag_no_fourier[,1])), 4) # 
round(c(mean(ex2_mat_first_d_1.0_n2000_bandconst_1_const_tuning_0.4_diag_no_fourier[,1]), median(ex2_mat_first_d_1.0_n2000_bandconst_1_const_tuning_0.4_diag_no_fourier[,1])), 4)   # 
round(c(mean(ex2_mat_first_d_1.25_n2000_bandconst_1_const_tuning_0.4_diag_no_fourier[,1]), median(ex2_mat_first_d_1.25_n2000_bandconst_1_const_tuning_0.4_diag_no_fourier[,1])), 4) # 

table(ex2_mat_first_d_0.75_n2000_bandconst_1_const_tuning_0.4_diag_no_fourier[,2]) # 7 1 992
table(ex2_mat_first_d_1.0_n2000_bandconst_1_const_tuning_0.4_diag_no_fourier[,2])  # 45 7 948
table(ex2_mat_first_d_1.25_n2000_bandconst_1_const_tuning_0.4_diag_no_fourier[,2]) # 217 65 718

