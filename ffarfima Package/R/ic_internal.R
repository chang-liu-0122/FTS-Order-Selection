residual_autocov <- function(E, h)
{
  T <- ncol(E)

  A <- E[, (h + 1):T, drop = FALSE]
  B <- E[, 1:(T - h), drop = FALSE]

  Gamma <- A %*% t(B) / T

  return(Gamma)
}

temporal_dependence <- function(E, H)
{
  fit <- 0

  for(h in 1:H)
  {
    G <- residual_autocov(E, h)
    fit <- fit + sum(G^2)
  }

  return(fit / H)
}

