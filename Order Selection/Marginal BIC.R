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



compute_FT_far <- function(X, p, q, Npc = 1)
{
  T_full <- ncol(X)

  if(p == 0 && q == 0) {
    return(list(F_T = NA, converged = FALSE))
  }

  fit <- tryCatch(
    farma_estimation(X, p = p, q = q, m = NULL, Npc = Npc),
    error = function(e) NULL
  )

  if(is.null(fit)) {
    return(list(F_T = Inf, converged = FALSE))
  }

  Lmax <- max(p, q)

  E <- fit$residuals[, (Lmax + 1):T_full, drop = FALSE]

  # Center residual curves pointwise.
  E <- t(scale(t(E), center = TRUE, scale = FALSE))

  if(any(!is.finite(E))) {
    return(list(F_T = Inf, converged = FALSE))
  }

  T_eff <- ncol(E)

  H <- floor(log(T_eff))

  F_T <- temporal_dependence(E, H)

  return(list(F_T = F_T, converged = TRUE))
}



select_order_far_full <- function(X, p_max = 5, q_max = 5, Npc = 1)
{
  out <- list()
  counter <- 1

  for(p in 0:p_max)
  {
    for(q in 0:q_max)
    {
      if(p == 0 && q == 0) next

      res <- compute_FT_far(X, p, q, Npc = Npc)

      out[[counter]] <- data.frame(
        p = p,
        q = q,
        F_T = res$F_T,
        converged = res$converged
      )

      counter <- counter + 1
    }
  }

  ic_results <- dplyr::bind_rows(out)

  # Estimate adaptive noise level
  # delta_hat = min F_T over finite candidate set M
  delta_hat <- ic_results %>%
    dplyr::filter(is.finite(F_T)) %>%
    dplyr::filter(converged == TRUE) %>%
    dplyr::summarise(delta_hat = min(F_T, na.rm = TRUE)) %>%
    dplyr::pull(delta_hat)

  # If all models fail
  if(length(delta_hat) == 0 || !is.finite(delta_hat)) {
    ic_results <- ic_results %>%
      dplyr::mutate(
        delta_hat = NA_real_,
        penalty = NA_real_,
        IC = NA_real_
      )

    return(list(
      p_hat = NA,
      q_hat = NA,
      delta_hat = NA,
      IC_results = ic_results
    ))
  }

  T <- ncol(X)

  ic_results <- ic_results %>%
    dplyr::mutate(
      delta_hat = delta_hat,
      penalty = log(log(T)) * delta_hat * (p + q),
      IC = F_T + penalty
    )

  best <- ic_results %>%
    dplyr::filter(is.finite(IC)) %>%
    dplyr::slice_min(IC, n = 1, with_ties = FALSE)

  list(
    p_hat = best$p,
    q_hat = best$q,
    delta_hat = delta_hat,
    IC_results = ic_results
  )
}




select_order_far_auto <- function(X, p_max = 5, q_max = 5, Npc = 1)
{
  T <- ncol(X)

  # Step 1: approximate delta_hat using the largest model
  big_fit <- compute_FT_far(X, p = p_max, q = q_max, Npc = Npc)

  if(!isTRUE(big_fit$converged) || !is.finite(big_fit$F_T)) {
    return(list(
      p_hat = NA,
      q_hat = NA,
      delta_hat = NA_real_,
      IC_results = data.frame(),
      search_path = data.frame()
    ))
  }

  delta_hat <- big_fit$F_T

  # This stores all models tried by the automatic search
  tried_models <- list()
  counter <- 1

  # Small helper only to avoid repeating code
  get_IC <- function(p, q)
  {
    res <- compute_FT_far(X, p = p, q = q, Npc = Npc)

    penalty <- log(log(T)) * delta_hat * (p + q)

    IC <- res$F_T + penalty

    if(!isTRUE(res$converged) || !is.finite(res$F_T)) {
      IC <- Inf
    }

    data.frame(
      p = p,
      q = q,
      F_T = res$F_T,
      converged = res$converged,
      delta_hat = delta_hat,
      penalty = penalty,
      IC = IC
    )
  }

  # Step 2: starting models
  starts <- data.frame(
    p = c(2, 1, 0),
    q = c(2, 0, 1)
  ) %>%
    dplyr::filter(
      p <= p_max,
      q <= q_max,
      !(p == 0 & q == 0)
    )

  start_results <- list()

  for(i in 1:nrow(starts))
  {
    tmp <- get_IC(starts$p[i], starts$q[i])
    start_results[[i]] <- tmp

    tried_models[[counter]] <- tmp
    counter <- counter + 1
  }

  start_results <- dplyr::bind_rows(start_results) %>%
    dplyr::filter(is.finite(IC), converged == TRUE)

  if(nrow(start_results) == 0) {
    return(list(
      p_hat = NA,
      q_hat = NA,
      delta_hat = delta_hat,
      IC_results = dplyr::bind_rows(tried_models),
      search_path = data.frame()
    ))
  }

  current <- start_results %>%
    dplyr::slice_min(IC, n = 1, with_ties = FALSE)

  current_p <- current$p
  current_q <- current$q
  current_IC <- current$IC

  search_path <- data.frame(
    step = 0,
    p = current_p,
    q = current_q,
    IC = current_IC
  )

  # Step 3: local search
  improved <- TRUE
  step <- 0

  while(improved)
  {
    improved <- FALSE

    neighbours <- expand.grid(
      dp = c(-1, 0, 1),
      dq = c(-1, 0, 1)
    ) %>%
      dplyr::filter(!(dp == 0 & dq == 0)) %>%
      dplyr::mutate(
        p = current_p + dp,
        q = current_q + dq
      ) %>%
      dplyr::filter(
        p >= 0,
        q >= 0,
        p <= p_max,
        q <= q_max,
        !(p == 0 & q == 0)
      ) %>%
      dplyr::select(p, q) %>%
      dplyr::distinct()

    neighbour_results <- list()

    for(i in 1:nrow(neighbours))
    {
      tmp <- get_IC(neighbours$p[i], neighbours$q[i])
      neighbour_results[[i]] <- tmp

      tried_models[[counter]] <- tmp
      counter <- counter + 1
    }

    neighbour_results <- dplyr::bind_rows(neighbour_results) %>%
      dplyr::filter(is.finite(IC), converged == TRUE)

    if(nrow(neighbour_results) > 0)
    {
      best_neighbour <- neighbour_results %>%
        dplyr::slice_min(IC, n = 1, with_ties = FALSE)

      if(best_neighbour$IC < current_IC)
      {
        current_p <- best_neighbour$p
        current_q <- best_neighbour$q
        current_IC <- best_neighbour$IC

        step <- step + 1

        search_path <- dplyr::bind_rows(
          search_path,
          data.frame(
            step = step,
            p = current_p,
            q = current_q,
            IC = current_IC
          )
        )

        improved <- TRUE
      }
    }
  }

  IC_results <- dplyr::bind_rows(tried_models) %>%
    dplyr::distinct(p, q, .keep_all = TRUE) %>%
    dplyr::arrange(p, q)

  return(list(
    p_hat = current_p,
    q_hat = current_q,
    delta_hat = delta_hat,
    IC_results = IC_results,
    search_path = search_path
  ))
}







