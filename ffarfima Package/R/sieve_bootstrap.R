sieve_bootstrap <- function(fun_dat, CPV_percent = 0.95, VAR_type = "none", B, burn_in)
{

    n_row = nrow(fun_dat)
    n_col = ncol(fun_dat)

    center_fun_dat = scale(fun_dat, center = TRUE, scale = FALSE)

    X_bar = colMeans(fun_dat)

    if(is.null(colnames(fun_dat)))
    {
      grid = seq(1, n_col, length.out = n_col)
    }
    else
    {
      grid = seq(min(as.numeric(colnames(fun_dat))),
                 max(as.numeric(colnames(fun_dat))),
                 length.out = n_col)
    }

    colnames(center_fun_dat) = grid
    rownames(center_fun_dat) = 1:n_row
    fun_dat_object = rainbow::fts(grid, t(center_fun_dat))

    ncomp = utils::head(which(cumsum(ftsa::ftsm(y = fun_dat_object, order = min(nrow(center_fun_dat), ncol(center_fun_dat)))$varprop)>= CPV_percent), 1)

    ftsm_object_scores = ftsa::ftsm(y = fun_dat_object, order = ncomp, mean = FALSE)$coeff

    VAR_max_order = min(9, nrow(ftsm_object_scores) - 2)
    VAR_AICC = vector("numeric", VAR_max_order)
    if(ncomp == 1)
    {
        for(VAR_order in 1:VAR_max_order)
        {
            ftsm_object_resi = stats::na.omit(stats::residuals(stats::ar(x = ftsm_object_scores, aic = FALSE, order.max = VAR_order)))
            sigma = det(crossprod(ftsm_object_resi)/n_row)
            VAR_AICC[VAR_order] = n_row * log(sigma) + n_row * (n_row * ncomp + VAR_order * (ncomp^2))/(n_row - ncomp * (VAR_order + 1) - 1)
        }
    }
    if(ncomp > 1)
    {
        for(VAR_order in 1:VAR_max_order)
        {
            ftsm_object_resi = stats::residuals(vars::VAR(y = ftsm_object_scores, p = VAR_order, type = "none"))
            sigma = det(crossprod(ftsm_object_resi)/n_row)
            VAR_AICC[VAR_order] = n_row * log(sigma) + n_row * (n_row * ncomp + VAR_order * (ncomp^2))/(n_row - ncomp * (VAR_order + 1) - 1)
        }
    }

    p_order = which.min(VAR_AICC[is.finite(VAR_AICC)])
    m_component = ncomp

    variance_est = stats::var(center_fun_dat)

    eigen_decomp = eigen(variance_est)
    eigen_decomp_vector = as.matrix(eigen_decomp$vectors[,1:m_component])
    eigen_decomp_score  = center_fun_dat %*% eigen_decomp_vector
    colnames(eigen_decomp_score) = 1:ncol(eigen_decomp_score)

    reconstruction = eigen_decomp_score %*% t(eigen_decomp_vector[,1:m_component])
    reconstruction_mean = (sweep(reconstruction, 2, X_bar, FUN = "+"))
    reconstruction_err = fun_dat - reconstruction_mean

    if(m_component == 1)
    {
      A_est = matrix(stats::ar(eigen_decomp_score, aic = FALSE, order.max = p_order)$ar, p_order, 1)
      residual_e = matrix(stats::ar(eigen_decomp_score, aic = FALSE, order.max = p_order)$resid[-(1:p_order)],,1)
    }
    if(m_component > 1)
    {
      VAR_forward = vars::VAR(eigen_decomp_score, p = p_order, type = VAR_type)
      A_est = matrix(NA, p_order * m_component, m_component)
      for(iw in 1:m_component)
      {
        A_est[,iw] = VAR_forward$varresult[[iw]]$coef
      }
      residual_e = stats::resid(VAR_forward)
      rownames(residual_e) = (p_order + 1):n_row
    }
    colnames(A_est) = colnames(eigen_decomp_score) = colnames(residual_e) = paste("score", 1:m_component, sep="_")
    A_est_list = list()
    for(j in 1:p_order)
    {
      A_est_list[[j]] = A_est[((j-1)*m_component+1):(j*m_component),]
    }

    residual_e_centered = scale(residual_e, center = TRUE, scale = FALSE)

    eigen_decomp_boot = array(NA, dim = c(n_row + burn_in, m_component, B))
    for(iwk in 1:B)
    {
      eigen_decomp_boot[1:p_order,,iwk] = eigen_decomp_score[1:p_order,]
      for(iw in p_order:(n_row + (burn_in - 1)))
      {
        fitted_val = rep(0, m_component)
        for(ij in 1:p_order)
        {
          fitted_val = fitted_val + eigen_decomp_boot[rev(1:iw)[ij],,iwk] %*% A_est_list[[ij]]
        }
        eigen_decomp_boot[(iw+1),,iwk] = fitted_val + residual_e_centered[sample(1:nrow(residual_e_centered), size = 1),]
      }
    }

    U_centered = scale(reconstruction_err, center = TRUE, scale = FALSE)
    X_boot = array(NA, dim = c(n_col, (n_row + burn_in), B))
    for(ik in 1:B)
    {
      X_boot[,,ik] = (eigen_decomp_vector %*% t(eigen_decomp_boot[,,ik])) + t(U_centered[sample(1:n_row, size = (n_row + burn_in), replace = TRUE),]) + matrix(rep(X_bar, (n_row + burn_in)), n_col, (n_row + burn_in))
    }
    X_boot_sample = X_boot[,(burn_in+1):(n_row+burn_in),,drop = FALSE]
    return(list(X_boot = X_boot_sample, p_order = p_order, m_component = m_component))
}
