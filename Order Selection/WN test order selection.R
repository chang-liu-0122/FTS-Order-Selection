select_order_residual_test <- function(X,
                                       Npc,
                                       p_max = 5,
                                       q_max = 5,
                                       lag_test = NULL,
                                       alpha = 0.05,
                                       verbose = TRUE) {
  T <- ncol(X)

  if (is.null(lag_test)) {
    lag_test = 20
    # lag_test <- floor(log(T))
  }

  lag_test <- max(1, min(lag_test, T - 1))

  detail_list <- list()
  counter <- 1

  for (total_order in 1:(p_max + q_max)) {

    for (q in 0:min(q_max, total_order)) {

      p <- total_order - q

      if (p < 0 || p > p_max) next
      if (p == 0 && q == 0) next

      if (verbose) {
        cat("\n-----------------------------------\n")
        cat("Trying FARMA(", p, ",", q, ")\n", sep = "")
        cat("lag_test = ", lag_test, "\n", sep = "")
      }

      out_row <- tryCatch({

        fit <- farma_estimation(
          X,
          p = p,
          q = q,
          m = NULL,
          Npc = Npc
        )

        Lmax <- max(p, q)

        E <- fit$residuals[, (Lmax + 1):T, drop = FALSE]

        E <- t(scale(t(E), center = TRUE, scale = FALSE))

        if (any(!is.finite(E))) {

          if (verbose) {
            cat("Residuals contain NA/Inf -> skip\n")
          }

          data.frame(
            p = p,
            q = q,
            total_order = total_order,
            lag_test = lag_test,
            p_value = NA_real_,
            passed = FALSE,
            converged = FALSE,
            selected = FALSE,
            error_message = "Residuals contain NA/Inf"
          )

        } else {

          wn_test <- fport_test(
            f_data = E,
            test = "multi-lag",
            lag = lag_test,
            iid = FALSE,
            suppress_print_output = TRUE
          )

          p_value <- wn_test$p_value

          passed <- is.finite(p_value) && p_value > alpha

          if (verbose) {
            cat("White-noise test p-value: ",
                round(p_value, 4),
                "\n",
                sep = "")
          }

          data.frame(
            p = p,
            q = q,
            total_order = total_order,
            lag_test = lag_test,
            p_value = p_value,
            passed = passed,
            converged = TRUE,
            selected = FALSE,
            error_message = NA_character_
          )
        }

      }, error = function(e) {

        data.frame(
          p = p,
          q = q,
          total_order = total_order,
          lag_test = lag_test,
          p_value = NA_real_,
          passed = FALSE,
          converged = FALSE,
          selected = FALSE,
          error_message = e$message
        )
      })

      detail_list[[counter]] <- out_row
      counter <- counter + 1

      if (isTRUE(out_row$passed)) {

        detail <- dplyr::bind_rows(detail_list)

        detail <- detail %>%
          dplyr::mutate(
            selected = (p == out_row$p & q == out_row$q)
          )

        if (verbose) {
          cat("\nSelected FARMA(", p, ",", q, ")\n", sep = "")
        }

        return(list(
          p_hat = p,
          q_hat = q,
          p_value = out_row$p_value,
          lag_test = lag_test,
          detail = detail
        ))
      }
    }
  }

  if (verbose) {
    cat("\nNo model passed the white-noise test.\n")
  }

  detail <- dplyr::bind_rows(detail_list)

  return(list(
    p_hat = NA_integer_,
    q_hat = NA_integer_,
    p_value = NA_real_,
    lag_test = lag_test,
    detail = detail
  ))
}
