parallel_sim_order_selection_FARMA_compare <- function(
    n_rep = 100,
    grid,
    p_max = 5,
    q_max = 5,
    Npc = 1,
    no_grid = 101,
    lag_test = NULL,
    alpha = 0.05
) {
  with_progress({

    seeds <- 1:n_rep

    p_prog <- progressor(
      steps = length(seeds) * nrow(grid)
    )

    results <- future_lapply(seeds, function(s) {

      simple_list <- list()
      IC_detail_list <- list()
      WN_detail_list <- list()

      counter <- 1

      for (i in 1:nrow(grid)) {

        out <- tryCatch({

          n_i  <- grid$T[i]
          ar_i <- grid$ar_strength[i]
          ma_i <- grid$ma_strength[i]
          p_i  <- grid$p[i]
          q_i  <- grid$q[i]

          sim <- sim_FARMA(
            n = n_i,
            d = 0,
            no_grid = no_grid,
            seed_number = s * 1000 + i,
            p = p_i,
            q = q_i,
            ar_strength = ar_i,
            ma_strength = ma_i,
            lag_strength = "equal"
          )

          X <- sim$X

          res_ic <- select_order_far_full(
            X = X,
            p_max = p_max,
            q_max = q_max,
            Npc = Npc
          )

          simple_ic <- data.frame(
            n = n_i,
            ar_strength = ar_i,
            ma_strength = ma_i,
            p_true = p_i,
            q_true = q_i,
            seed = s,
            method = "Adaptive_IC",
            p_hat = res_ic$p_hat,
            q_hat = res_ic$q_hat,
            delta_hat = res_ic$delta_hat,
            p_value = NA_real_,
            lag_test = NA_integer_,
            correct = (res_ic$p_hat == p_i & res_ic$q_hat == q_i),
            error_message = NA_character_
          )

          IC_detail <- res_ic$IC_results %>%
            dplyr::mutate(
              n = n_i,
              ar_strength = ar_i,
              ma_strength = ma_i,
              p_true = p_i,
              q_true = q_i,
              seed = s,
              p_hat = res_ic$p_hat,
              q_hat = res_ic$q_hat,
              correct_selected = (res_ic$p_hat == p_i & res_ic$q_hat == q_i),
              .before = 1
            )

          res_wn <- select_order_residual_test(
            X = X,
            Npc = Npc,
            p_max = p_max,
            q_max = q_max,
            lag_test = lag_test,
            alpha = alpha,
            verbose = FALSE
          )

          simple_wn <- data.frame(
            n = n_i,
            ar_strength = ar_i,
            ma_strength = ma_i,
            p_true = p_i,
            q_true = q_i,
            seed = s,
            method = "WN_Test",
            p_hat = res_wn$p_hat,
            q_hat = res_wn$q_hat,
            delta_hat = NA_real_,
            p_value = res_wn$p_value,
            lag_test = res_wn$lag_test,
            correct = (res_wn$p_hat == p_i & res_wn$q_hat == q_i),
            error_message = NA_character_
          )

          WN_detail <- res_wn$detail %>%
            dplyr::mutate(
              n = n_i,
              ar_strength = ar_i,
              ma_strength = ma_i,
              p_true = p_i,
              q_true = q_i,
              seed = s,
              p_hat = res_wn$p_hat,
              q_hat = res_wn$q_hat,
              correct_selected = (res_wn$p_hat == p_i & res_wn$q_hat == q_i),
              .before = 1
            )

          list(
            simple = dplyr::bind_rows(simple_ic, simple_wn),
            IC_detail = IC_detail,
            WN_detail = WN_detail
          )

        }, error = function(e) {

          simple_error <- data.frame(
            n = grid$T[i],
            ar_strength = grid$ar_strength[i],
            ma_strength = grid$ma_strength[i],
            p_true = grid$p[i],
            q_true = grid$q[i],
            seed = s,
            method = c("Adaptive_IC", "WN_Test"),
            p_hat = NA_integer_,
            q_hat = NA_integer_,
            delta_hat = NA_real_,
            p_value = NA_real_,
            lag_test = NA_integer_,
            correct = NA,
            error_message = e$message
          )

          list(
            simple = simple_error,
            IC_detail = data.frame(),
            WN_detail = data.frame()
          )
        })

        simple_list[[counter]] <- out$simple
        IC_detail_list[[counter]] <- out$IC_detail
        WN_detail_list[[counter]] <- out$WN_detail

        counter <- counter + 1

        p_prog()
      }

      list(
        simple = dplyr::bind_rows(simple_list),
        IC_detail = dplyr::bind_rows(IC_detail_list),
        WN_detail = dplyr::bind_rows(WN_detail_list)
      )

    }, future.seed = TRUE)

    simple_result <- dplyr::bind_rows(
      lapply(results, function(x) x$simple)
    )

    IC_detail <- dplyr::bind_rows(
      lapply(results, function(x) x$IC_detail)
    )

    WN_detail <- dplyr::bind_rows(
      lapply(results, function(x) x$WN_detail)
    )

    list(
      simple_result = simple_result,
      IC_detail = IC_detail,
      WN_detail = WN_detail
    )
  })
}
