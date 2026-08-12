##########
# n = 500
##########

# A_structure = "diag"

registerDoMC(24)
ex2_first_d_0.75_n500_bandconst_1_const_tuning_0.4_diag_no_basis_permute = foreach(iwk = 1:1000) %dopar% H_estimate(fun_dat = sim_ex2(n = 500, d = 0.75, seed_number = iwk, A_structure = "diag", basis_permute = "FALSE"),
                                                                                             band_const = 1, tuning_const = 0.4, delta_val = 1e-2, 
                                                                                             choice_score = "first", taper_choice = "none")

ex2_first_d_1.0_n500_bandconst_1_const_tuning_0.4_diag_no_basis_permute = foreach(iwk = 1:1000) %dopar% H_estimate(fun_dat = sim_ex2(n = 500, d = 1.0, seed_number = iwk, A_structure = "diag", basis_permute = "FALSE"),
                                                                                             band_const = 1, tuning_const = 0.4, delta_val = 1e-2, 
                                                                                             choice_score = "first", taper_choice = "none")

ex2_first_d_1.25_n500_bandconst_1_const_tuning_0.4_diag_no_basis_permute = foreach(iwk = 1:1000) %dopar% H_estimate(fun_dat = sim_ex2(n = 500, d = 1.255, seed_number = iwk, A_structure = "diag", basis_permute = "FALSE"),
                                                                                             band_const = 1, tuning_const = 0.4, delta_val = 1e-2, 
                                                                                             choice_score = "first", taper_choice = "none")


ex2_mat_first_d_0.75_n500_bandconst_1_const_tuning_0.4_diag_no_basis_permute = ex2_mat_first_d_1.0_n500_bandconst_1_const_tuning_0.4_diag_no_basis_permute = ex2_mat_first_d_1.25_n500_bandconst_1_const_tuning_0.4_diag_no_basis_permute = matrix(NA, 1000, 3)
for(ik in 1:1000)
{
    dum = ex2_first_d_0.75_n500_bandconst_1_const_tuning_0.4_diag_no_basis_permute[[ik]]
    ex2_mat_first_d_0.75_n500_bandconst_1_const_tuning_0.4_diag_no_basis_permute[ik,] = c(dum$d_est, dum$ncomp, dum$diff_order)
    rm(dum) 
    
    dum = ex2_first_d_1.0_n500_bandconst_1_const_tuning_0.4_diag_no_basis_permute[[ik]]
    ex2_mat_first_d_1.0_n500_bandconst_1_const_tuning_0.4_diag_no_basis_permute[ik,] = c(dum$d_est, dum$ncomp, dum$diff_order)
    rm(dum)
    
    dum = ex2_first_d_1.25_n500_bandconst_1_const_tuning_0.4_diag_no_basis_permute[[ik]]
    ex2_mat_first_d_1.25_n500_bandconst_1_const_tuning_0.4_diag_no_basis_permute[ik,] = c(dum$d_est, dum$ncomp, dum$diff_order)
    rm(dum); print(ik); rm(ik)
}

round(c(mean(ex2_mat_first_d_0.75_n500_bandconst_1_const_tuning_0.4_diag_no_basis_permute[,1]), median(ex2_mat_first_d_0.75_n500_bandconst_1_const_tuning_0.4_diag_no_basis_permute[,1])), 4) # 0.7898 0.7952
round(c(mean(ex2_mat_first_d_1.0_n500_bandconst_1_const_tuning_0.4_diag_no_basis_permute[,1]), median(ex2_mat_first_d_1.0_n500_bandconst_1_const_tuning_0.4_diag_no_basis_permute[,1])), 4)   # 0.9575 0.9590
round(c(mean(ex2_mat_first_d_1.25_n500_bandconst_1_const_tuning_0.4_diag_no_basis_permute[,1]), median(ex2_mat_first_d_1.25_n500_bandconst_1_const_tuning_0.4_diag_no_basis_permute[,1])), 4) # 1.1719 1.1734

table(ex2_mat_first_d_0.75_n500_bandconst_1_const_tuning_0.4_diag_no_basis_permute[,2]) # 446 553 1
table(ex2_mat_first_d_1.0_n500_bandconst_1_const_tuning_0.4_diag_no_basis_permute[,2])  # 220 777 3
table(ex2_mat_first_d_1.25_n500_bandconst_1_const_tuning_0.4_diag_no_basis_permute[,2]) # 263 728 9


# A_structure = "band"

registerDoMC(24)
ex2_first_d_0.75_n500_bandconst_1_const_tuning_0.4_no_basis_permute = foreach(iwk = 1:1000) %dopar% H_estimate(fun_dat = sim_ex2(n = 500, d = 0.75, seed_number = iwk, A_structure = "band", basis_permute = "FALSE"),
                                                                                             band_const = 1, tuning_const = 0.4, delta_val = 1e-2, 
                                                                                             choice_score = "first", taper_choice = "none")

ex2_first_d_1.0_n500_bandconst_1_const_tuning_0.4_no_basis_permute = foreach(iwk = 1:1000) %dopar% H_estimate(fun_dat = sim_ex2(n = 500, d = 1.0, seed_number = iwk, A_structure = "band", basis_permute = "FALSE"),
                                                                                             band_const = 1, tuning_const = 0.4, delta_val = 1e-2, 
                                                                                             choice_score = "first", taper_choice = "none")

ex2_first_d_1.25_n500_bandconst_1_const_tuning_0.4_no_basis_permute = foreach(iwk = 1:1000) %dopar% H_estimate(fun_dat = sim_ex2(n = 500, d = 1.25, seed_number = iwk, A_structure = "band", basis_permute = "FALSE"),
                                                                                             band_const = 1, tuning_const = 0.4, delta_val = 1e-2, 
                                                                                             choice_score = "first", taper_choice = "none")


ex2_mat_first_d_0.75_n500_bandconst_1_const_tuning_0.4_no_basis_permute = ex2_mat_first_d_1.0_n500_bandconst_1_const_tuning_0.4_no_basis_permute = ex2_mat_first_d_1.25_n500_bandconst_1_const_tuning_0.4_no_basis_permute = matrix(NA, 1000, 3)
for(ik in 1:1000)
{
    dum = ex2_first_d_0.75_n500_bandconst_1_const_tuning_0.4_no_basis_permute[[ik]]
    ex2_mat_first_d_0.75_n500_bandconst_1_const_tuning_0.4_no_basis_permute[ik,] = c(dum$d_est, dum$ncomp, dum$diff_order)
    rm(dum) 
    
    dum = ex2_first_d_1.0_n500_bandconst_1_const_tuning_0.4_no_basis_permute[[ik]]
    ex2_mat_first_d_1.0_n500_bandconst_1_const_tuning_0.4_no_basis_permute[ik,] = c(dum$d_est, dum$ncomp, dum$diff_order)
    rm(dum)
    
    dum = ex2_first_d_1.25_n500_bandconst_1_const_tuning_0.4_no_basis_permute[[ik]]
    ex2_mat_first_d_1.25_n500_bandconst_1_const_tuning_0.4_no_basis_permute[ik,] = c(dum$d_est, dum$ncomp, dum$diff_order)
    rm(dum); print(ik); rm(ik)
}

round(c(mean(ex2_mat_first_d_0.75_n500_bandconst_1_const_tuning_0.4_no_basis_permute[,1]), median(ex2_mat_first_d_0.75_n500_bandconst_1_const_tuning_0.4_no_basis_permute[,1])), 4) # 0.7893 0.7947
round(c(mean(ex2_mat_first_d_1.0_n500_bandconst_1_const_tuning_0.4_no_basis_permute[,1]), median(ex2_mat_first_d_1.0_n500_bandconst_1_const_tuning_0.4_no_basis_permute[,1])), 4)   # 0.9572 0.9600
round(c(mean(ex2_mat_first_d_1.25_n500_bandconst_1_const_tuning_0.4_no_basis_permute[,1]), median(ex2_mat_first_d_1.25_n500_bandconst_1_const_tuning_0.4_no_basis_permute[,1])), 4) # 1.1670 1.1689

table(ex2_mat_first_d_0.75_n500_bandconst_1_const_tuning_0.4_no_basis_permute[,2]) # 455 542 3
table(ex2_mat_first_d_1.0_n500_bandconst_1_const_tuning_0.4_no_basis_permute[,2])  # 230 767 3 
table(ex2_mat_first_d_1.25_n500_bandconst_1_const_tuning_0.4_no_basis_permute[,2]) # 261 730 9

###########
# n = 1000
###########

# A_structure = "diag"

registerDoMC(24)
ex2_first_d_0.75_n1000_bandconst_1_const_tuning_0.4_diag_no_basis_permute = foreach(iwk = 1:1000) %dopar% H_estimate(fun_dat = sim_ex2(n = 1000, d = 0.75, seed_number = iwk, A_structure = "diag", basis_permute = "FALSE"),
                                                                                              band_const = 1, tuning_const = 0.4, delta_val = 1e-2,
                                                                                              choice_score = "first", taper_choice = "none")

ex2_first_d_1.0_n1000_bandconst_1_const_tuning_0.4_diag_no_basis_permute = foreach(iwk = 1:1000) %dopar% H_estimate(fun_dat = sim_ex2(n = 1000, d = 1.0, seed_number = iwk, A_structure = "diag", basis_permute = "FALSE"),
                                                                                              band_const = 1, tuning_const = 0.4, delta_val = 1e-2, 
                                                                                              choice_score = "first", taper_choice = "none")

ex2_first_d_1.25_n1000_bandconst_1_const_tuning_0.4_diag_no_basis_permute = foreach(iwk = 1:1000) %dopar% H_estimate(fun_dat = sim_ex2(n = 1000, d = 1.25, seed_number = iwk, A_structure = "diag", basis_permute = "FALSE"),
                                                                                              band_const = 1, tuning_const = 0.4, delta_val = 1e-2, 
                                                                                              choice_score = "first", taper_choice = "none")


ex2_mat_first_d_0.75_n1000_bandconst_1_const_tuning_0.4_diag_no_basis_permute = ex2_mat_first_d_1.0_n1000_bandconst_1_const_tuning_0.4_diag_no_basis_permute = ex2_mat_first_d_1.25_n1000_bandconst_1_const_tuning_0.4_diag_no_basis_permute = matrix(NA, 1000, 3)
for(ik in 1:1000)
{
    dum = ex2_first_d_0.75_n1000_bandconst_1_const_tuning_0.4_diag_no_basis_permute[[ik]] 
    ex2_mat_first_d_0.75_n1000_bandconst_1_const_tuning_0.4_diag_no_basis_permute[ik,] = c(dum$d_est, dum$ncomp, dum$diff_order)
    rm(dum)
    
    dum = ex2_first_d_1.0_n1000_bandconst_1_const_tuning_0.4_diag_no_basis_permute[[ik]] 
    ex2_mat_first_d_1.0_n1000_bandconst_1_const_tuning_0.4_diag_no_basis_permute[ik,] = c(dum$d_est, dum$ncomp, dum$diff_order)
    rm(dum)
    
    dum = ex2_first_d_1.25_n1000_bandconst_1_const_tuning_0.4_diag_no_basis_permute[[ik]] 
    ex2_mat_first_d_1.25_n1000_bandconst_1_const_tuning_0.4_diag_no_basis_permute[ik,] = c(dum$d_est, dum$ncomp, dum$diff_order)
    rm(dum); print(ik); rm(ik)
}

round(c(mean(ex2_mat_first_d_0.75_n1000_bandconst_1_const_tuning_0.4_diag_no_basis_permute[,1]), median(ex2_mat_first_d_0.75_n1000_bandconst_1_const_tuning_0.4_diag_no_basis_permute[,1])), 4) # 0.7948 0.7983
round(c(mean(ex2_mat_first_d_1.0_n1000_bandconst_1_const_tuning_0.4_diag_no_basis_permute[,1]),  median(ex2_mat_first_d_1.0_n1000_bandconst_1_const_tuning_0.4_diag_no_basis_permute[,1])), 4)  # 0.9702 0.9726
round(c(mean(ex2_mat_first_d_1.25_n1000_bandconst_1_const_tuning_0.4_diag_no_basis_permute[,1]), median(ex2_mat_first_d_1.25_n1000_bandconst_1_const_tuning_0.4_diag_no_basis_permute[,1])), 4) # 1.1900 1.1923

table(ex2_mat_first_d_0.75_n1000_bandconst_1_const_tuning_0.4_diag_no_basis_permute[,2]) # 269 731
table(ex2_mat_first_d_1.0_n1000_bandconst_1_const_tuning_0.4_diag_no_basis_permute[,2])  # 158 838 4
table(ex2_mat_first_d_1.25_n1000_bandconst_1_const_tuning_0.4_diag_no_basis_permute[,2]) # 256 732 12


# A_structure = "band"

ex2_first_d_0.75_n1000_bandconst_1_const_tuning_0.4_no_basis_permute = foreach(iwk = 1:1000) %dopar% H_estimate(fun_dat = sim_ex2(n = 1000, d = 0.75, seed_number = iwk, A_structure = "band", basis_permute = "FALSE"),
                                                                                             band_const = 1, tuning_const = 0.4, delta_val = 1e-2,
                                                                                             choice_score = "first", taper_choice = "none")

ex2_first_d_1.0_n1000_bandconst_1_const_tuning_0.4_no_basis_permute = foreach(iwk = 1:1000) %dopar% H_estimate(fun_dat = sim_ex2(n = 1000, d = 1.0, seed_number = iwk, A_structure = "band", basis_permute = "FALSE"),
                                                                                             band_const = 1, tuning_const = 0.4, delta_val = 1e-2, 
                                                                                             choice_score = "first", taper_choice = "none")

ex2_first_d_1.25_n1000_bandconst_1_const_tuning_0.4_no_basis_permute = foreach(iwk = 1:1000) %dopar% H_estimate(fun_dat = sim_ex2(n = 1000, d = 1.25, seed_number = iwk, A_structure = "band", basis_permute = "FALSE"),
                                                                                             band_const = 1, tuning_const = 0.4, delta_val = 1e-2, 
                                                                                             choice_score = "first", taper_choice = "none")

ex2_mat_first_d_0.75_n1000_bandconst_1_const_tuning_0.4_no_basis_permute = ex2_mat_first_d_1.0_n1000_bandconst_1_const_tuning_0.4_no_basis_permute = ex2_mat_first_d_1.25_n1000_bandconst_1_const_tuning_0.4_no_basis_permute = matrix(NA, 1000, 3)
for(ik in 1:1000)
{
    dum = ex2_first_d_0.75_n1000_bandconst_1_const_tuning_0.4_no_basis_permute[[ik]]
    ex2_mat_first_d_0.75_n1000_bandconst_1_const_tuning_0.4_no_basis_permute[ik,] = c(dum$d_est, dum$ncomp, dum$diff_order)
    rm(dum) 
    
    dum = ex2_first_d_1.0_n1000_bandconst_1_const_tuning_0.4_no_basis_permute[[ik]]
    ex2_mat_first_d_1.0_n1000_bandconst_1_const_tuning_0.4_no_basis_permute[ik,] = c(dum$d_est, dum$ncomp, dum$diff_order)
    rm(dum) 
    
    dum = ex2_first_d_1.25_n1000_bandconst_1_const_tuning_0.4_no_basis_permute[[ik]]
    ex2_mat_first_d_1.25_n1000_bandconst_1_const_tuning_0.4_no_basis_permute[ik,] = c(dum$d_est, dum$ncomp, dum$diff_order)
    rm(dum); print(ik); rm(ik)
}

round(c(mean(ex2_mat_first_d_0.75_n1000_bandconst_1_const_tuning_0.4_no_basis_permute[,1]), median(ex2_mat_first_d_0.75_n1000_bandconst_1_const_tuning_0.4_no_basis_permute[,1])), 4) # 0.7945 0.7979
round(c(mean(ex2_mat_first_d_1.0_n1000_bandconst_1_const_tuning_0.4_no_basis_permute[,1]), median(ex2_mat_first_d_1.0_n1000_bandconst_1_const_tuning_0.4_no_basis_permute[,1])), 4)   # 0.9701 0.9728
round(c(mean(ex2_mat_first_d_1.25_n1000_bandconst_1_const_tuning_0.4_no_basis_permute[,1]), median(ex2_mat_first_d_1.25_n1000_bandconst_1_const_tuning_0.4_no_basis_permute[,1])), 4) # 1.1898 1.1923

table(ex2_mat_first_d_0.75_n1000_bandconst_1_const_tuning_0.4_no_basis_permute[,2]) # 289 710 1
table(ex2_mat_first_d_1.0_n1000_bandconst_1_const_tuning_0.4_no_basis_permute[,2])  # 166 832 2
table(ex2_mat_first_d_1.25_n1000_bandconst_1_const_tuning_0.4_no_basis_permute[,2]) # 255 733 12


###########
# n = 2000
###########

# A_structure = "diag"

registerDoMC(24)
ex2_first_d_0.75_n2000_bandconst_1_const_tuning_0.4_diag_no_basis_permute = foreach(iwk = 1:1000) %dopar% H_estimate(fun_dat = sim_ex2(n = 2000, d = 0.75, seed_number = iwk, A_structure = "diag", basis_permute = "FALSE"),
                                                                                              band_const = 1, tuning_const = 0.4, delta_val = 1e-2, 
                                                                                              choice_score = "first", taper_choice = "none")

ex2_first_d_1.0_n2000_bandconst_1_const_tuning_0.4_diag_no_basis_permute = foreach(iwk = 1:1000) %dopar% H_estimate(fun_dat = sim_ex2(n = 2000, d = 1.0, seed_number = iwk, A_structure = "diag", basis_permute = "FALSE"),
                                                                                              band_const = 1, tuning_const = 0.4, delta_val = 1e-2, 
                                                                                              choice_score = "first", taper_choice = "none")

ex2_first_d_1.25_n2000_bandconst_1_const_tuning_0.4_diag_no_basis_permute = foreach(iwk = 1:1000) %dopar% H_estimate(fun_dat = sim_ex2(n = 2000, d = 1.25, seed_number = iwk, A_structure = "diag", basis_permute = "FALSE"),
                                                                                              band_const = 1, tuning_const = 0.4, delta_val = 1e-2, 
                                                                                              choice_score = "first", taper_choice = "none")

ex2_mat_first_d_0.75_n2000_bandconst_1_const_tuning_0.4_diag_no_basis_permute = ex2_mat_first_d_1.0_n2000_bandconst_1_const_tuning_0.4_diag_no_basis_permute = ex2_mat_first_d_1.25_n2000_bandconst_1_const_tuning_0.4_diag_no_basis_permute = matrix(NA, 1000, 3)
for(ik in 1:1000)
{
    dum = ex2_first_d_0.75_n2000_bandconst_1_const_tuning_0.4_diag_no_basis_permute[[ik]]
    ex2_mat_first_d_0.75_n2000_bandconst_1_const_tuning_0.4_diag_no_basis_permute[ik,] = c(dum$d_est, dum$ncomp, dum$diff_order)
    rm(dum) 
    
    dum = ex2_first_d_1.0_n2000_bandconst_1_const_tuning_0.4_diag_no_basis_permute[[ik]]
    ex2_mat_first_d_1.0_n2000_bandconst_1_const_tuning_0.4_diag_no_basis_permute[ik,] = c(dum$d_est, dum$ncomp, dum$diff_order)
    rm(dum)
    
    dum = ex2_first_d_1.25_n2000_bandconst_1_const_tuning_0.4_diag_no_basis_permute[[ik]]
    ex2_mat_first_d_1.25_n2000_bandconst_1_const_tuning_0.4_diag_no_basis_permute[ik,] = c(dum$d_est, dum$ncomp, dum$diff_order)
    rm(dum); print(ik); rm(ik)
}

round(c(mean(ex2_mat_first_d_0.75_n2000_bandconst_1_const_tuning_0.4_diag_no_basis_permute[,1]), median(ex2_mat_first_d_0.75_n2000_bandconst_1_const_tuning_0.4_diag_no_basis_permute[,1])), 4) # 0.7918 0.7933
round(c(mean(ex2_mat_first_d_1.0_n2000_bandconst_1_const_tuning_0.4_diag_no_basis_permute[,1]), median(ex2_mat_first_d_1.0_n2000_bandconst_1_const_tuning_0.4_diag_no_basis_permute[,1])), 4)   # 0.9812 0.9842
round(c(mean(ex2_mat_first_d_1.25_n2000_bandconst_1_const_tuning_0.4_diag_no_basis_permute[,1]), median(ex2_mat_first_d_1.25_n2000_bandconst_1_const_tuning_0.4_diag_no_basis_permute[,1])), 4) # 1.2075 1.2099

table(ex2_mat_first_d_0.75_n2000_bandconst_1_const_tuning_0.4_diag_no_basis_permute[,2]) # 121 879
table(ex2_mat_first_d_1.0_n2000_bandconst_1_const_tuning_0.4_diag_no_basis_permute[,2])  # 117 876 7
table(ex2_mat_first_d_1.25_n2000_bandconst_1_const_tuning_0.4_diag_no_basis_permute[,2]) # 218 766 16


# A_structure = "band"

ex2_first_d_0.75_n2000_bandconst_1_const_tuning_0.4_no_basis_permute = foreach(iwk = 1:1000) %dopar% H_estimate(fun_dat = sim_ex2(n = 2000, d = 0.75, seed_number = iwk, A_structure = "band", basis_permute = "FALSE"),
                                                                                              band_const = 1, tuning_const = 0.4, delta_val = 1e-2, 
                                                                                              choice_score = "first", taper_choice = "none")

ex2_first_d_1.0_n2000_bandconst_1_const_tuning_0.4_no_basis_permute = foreach(iwk = 1:1000) %dopar% H_estimate(fun_dat = sim_ex2(n = 2000, d = 1.0, seed_number = iwk, A_structure = "band", basis_permute = "FALSE"),
                                                                                              band_const = 1, tuning_const = 0.4, delta_val = 1e-2, 
                                                                                              choice_score = "first", taper_choice = "none")

ex2_first_d_1.25_n2000_bandconst_1_const_tuning_0.4_no_basis_permute = foreach(iwk = 1:1000) %dopar% H_estimate(fun_dat = sim_ex2(n = 2000, d = 1.25, seed_number = iwk, A_structure = "band", basis_permute = "FALSE"),
                                                                                              band_const = 1, tuning_const = 0.4, delta_val = 1e-2, 
                                                                                              choice_score = "first", taper_choice = "none")

ex2_mat_first_d_0.75_n2000_bandconst_1_const_tuning_0.4_no_basis_permute = ex2_mat_first_d_1.0_n2000_bandconst_1_const_tuning_0.4_no_basis_permute = ex2_mat_first_d_1.25_n2000_bandconst_1_const_tuning_0.4_no_basis_permute = matrix(NA, 1000, 3)
for(ik in 1:1000)
{
    dum = ex2_first_d_0.75_n2000_bandconst_1_const_tuning_0.4_no_basis_permute[[ik]]
    ex2_mat_first_d_0.75_n2000_bandconst_1_const_tuning_0.4_no_basis_permute[ik,] = c(dum$d_est, dum$ncomp, dum$diff_order)
    rm(dum) 
    
    dum = ex2_first_d_1.0_n2000_bandconst_1_const_tuning_0.4_no_basis_permute[[ik]]
    ex2_mat_first_d_1.0_n2000_bandconst_1_const_tuning_0.4_no_basis_permute[ik,] = c(dum$d_est, dum$ncomp, dum$diff_order)
    rm(dum)
    
    dum = ex2_first_d_1.25_n2000_bandconst_1_const_tuning_0.4_no_basis_permute[[ik]]
    ex2_mat_first_d_1.25_n2000_bandconst_1_const_tuning_0.4_no_basis_permute[ik,] = c(dum$d_est, dum$ncomp, dum$diff_order)
    rm(dum); print(ik); rm(ik)
}

round(c(mean(ex2_mat_first_d_0.75_n2000_bandconst_1_const_tuning_0.4_no_basis_permute[,1]), median(ex2_mat_first_d_0.75_n2000_bandconst_1_const_tuning_0.4_no_basis_permute[,1])), 4) # 0.7919 0.7928
round(c(mean(ex2_mat_first_d_1.0_n2000_bandconst_1_const_tuning_0.4_no_basis_permute[,1]), median(ex2_mat_first_d_1.0_n2000_bandconst_1_const_tuning_0.4_no_basis_permute[,1])), 4)   # 0.9811 0.9842
round(c(mean(ex2_mat_first_d_1.25_n2000_bandconst_1_const_tuning_0.4_no_basis_permute[,1]), median(ex2_mat_first_d_1.25_n2000_bandconst_1_const_tuning_0.4_no_basis_permute[,1])), 4) # 1.2075 1.2094

table(ex2_mat_first_d_0.75_n2000_bandconst_1_const_tuning_0.4_no_basis_permute[,2]) # 132 868
table(ex2_mat_first_d_1.0_n2000_bandconst_1_const_tuning_0.4_no_basis_permute[,2])  # 119 875 6
table(ex2_mat_first_d_1.25_n2000_bandconst_1_const_tuning_0.4_no_basis_permute[,2]) # 218 766 16

