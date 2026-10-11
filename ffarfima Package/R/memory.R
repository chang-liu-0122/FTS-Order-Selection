#' Estimate the memory parameter of a functional time series
#'
#' The curves are centred, projected onto the leading eigenfunction of the
#' chosen covariance-type operator, and the memory parameter is estimated from
#' the resulting score series. When the estimate is close to one, the score is
#' differenced and the estimate is shifted back by the number of differences.
#'
#' @param x Numeric matrix of curves, with grid points in rows and time in columns.
#' @param method Estimator:
#'   * `"bc_lpwn2"`, `"bc_lpwn1"`: sieve-bootstrap bias-corrected LPWN with two
#'     or one short-memory polynomial terms;
#'   * `"lpwn2"`, `"lpwn1"`: local polynomial Whittle with noise;
#'   * `"lw"`: local Whittle;
#'   * `"dfa"`: detrended fluctuation analysis.
#' @param cov_type Operator for the leading eigenfunction: `"sample"`
#'   (empirical covariance), `"lrc_RS"` or `"lrc_LRS"` (long-run covariance).
#' @param m_const Bandwidth exponent, giving `floor(1 + T^m_const)` Fourier
#'   frequencies. Defaults to 0.75 for the LPWN estimators and 0.65 for LW and DFA.
#' @param B,burn_in Number of sieve-bootstrap samples and their burn-in length
#'   (bias-corrected estimators only).
#' @param R_noise Number of polynomial terms for the noise in LPWN.
#' @param tuning_const Constant \eqn{c} in the threshold \eqn{1 - c / \log m}
#'   above which the score is differenced.
#' @param cov_band_const Bandwidth constant for `cov_type = "lrc_LRS"`.
#' @param dfa_levels,dfa_minnpts,dfa_cutoff Block-size grid for DFA: number of
#'   block sizes, minimum number of points per block, and the range of block
#'   sizes used in the regression.
#' @param seed Optional seed for the bootstrap.
#' @param verbose Whether to report progress.
#' @return An object of class `ffarfima_memory` with the estimate `d`, the
#'   uncorrected estimate `d_initial`, the bootstrap `bias`, `sd` and
#'   `bootstrap` draws (bias-corrected estimators only), and the leading
#'   `scores` and `eigenfunction`. Printing it adds notes when an LPWN
#'   estimate sits at the edge of its search range.
#' @examples
#' x <- sim_ffarfima(n = 300, d = 0.3, p = 1, q = 0, seed = 1)$X
#' estimate_memory(x, method = "lw")
#' estimate_memory(x, method = "lpwn2")
#' @export
estimate_memory <- function(x,
                            method = c("bc_lpwn2", "bc_lpwn1", "lpwn2", "lpwn1", "lw", "dfa"),
                            cov_type = c("sample", "lrc_RS", "lrc_LRS"),
                            m_const = NULL,
                            B = 100,
                            burn_in = 200,
                            R_noise = 0,
                            tuning_const = 1,
                            cov_band_const = 1,
                            dfa_levels = 50,
                            dfa_minnpts = 3,
                            dfa_cutoff = 10^c(0.7, 2.5),
                            seed = NULL,
                            verbose = FALSE) {
  method <- match.arg(method)
  cov_type <- match.arg(cov_type)
  x <- as_curves(x, "x")

  if (is.null(m_const)) {
    m_const <- if (method %in% c("lw", "dfa")) 0.65 else 0.75
  }

  if (!is.null(seed)) {
    withr::local_seed(seed)
  }

  xc <- sweep(x, 1, rowMeans(x), "-")
  lead <- leading_pc(xc, get_covariance(xc, type = cov_type, band_const = cov_band_const))
  R_short <- if (method %in% c("lpwn1", "bc_lpwn1")) 1 else 2

  d_initial <- as.numeric(switch(
    method,
    lw = local_whittle_general(lead$scores, m_const = m_const, tuning_const = tuning_const),
    dfa = peng_general(lead$scores, levels = dfa_levels, minnpts = dfa_minnpts,
                       cut.off = dfa_cutoff, tuning_const = tuning_const, m_const = m_const),
    LPWN_general(lead$scores, m_const = m_const, R_short = R_short,
                 R_noise = R_noise, tuning_const = tuning_const)
  ))

  edge <- function(v) {
    r <- v %% 1
    is.finite(v) & (abs(r - 0.01) < 1e-6 | abs(r - 0.99) < 1e-6)
  }
  at_bound <- method %in% c("lpwn1", "lpwn2", "bc_lpwn1", "bc_lpwn2") && edge(d_initial)
  boot_at_bound <- NA_real_
  d <- d_initial
  bias <- NA_real_
  boot_sd <- NA_real_
  boot <- NULL

  if (method %in% c("bc_lpwn1", "bc_lpwn2")) {
    x_diff <- t(FCVAR::FracDiff(t(xc), d = d_initial))
    x_boot <- sieve_bootstrap(fun_dat = t(x_diff), CPV_percent = 0.95,
                              VAR_type = "none", B = B, burn_in = burn_in)$X_boot
    boot <- rep(NA_real_, B)
    bar <- progress_start(verbose, B, sprintf("Estimating d by BC-LPWN%d with %d bootstrap samples", R_short, B))

    for (b in seq_len(B)) {
      xb <- t(FCVAR::FracDiff(t(x_boot[, , b]), d = -d_initial))
      lead_b <- leading_pc(xb, get_covariance(xb, type = cov_type, band_const = cov_band_const))
      boot[b] <- tryCatch(
        as.numeric(LPWN_general(lead_b$scores, m_const = m_const, R_short = R_short,
                                R_noise = R_noise, tuning_const = tuning_const)),
        error = function(e) NA_real_
      )
      progress_step(bar, b)
    }

    progress_end(bar)
    boot_at_bound <- mean(edge(boot[is.finite(boot)]))
    bias <- mean(boot, na.rm = TRUE) - d_initial
    boot_sd <- sd(boot - d_initial, na.rm = TRUE)
    d <- d_initial - bias
  }

  admissible <- is.finite(d) && (!method %in% c("bc_lpwn1", "bc_lpwn2") || d >= 0)

  structure(
    list(
      d = d,
      method = method,
      d_initial = d_initial,
      bias = bias,
      sd = boot_sd,
      bootstrap = boot,
      m_const = m_const,
      cov_type = cov_type,
      scores = lead$scores,
      eigenfunction = lead$psi1
    ),
    class = "ffarfima_memory",
    diagnostics = list(admissible = admissible, at_bound = at_bound, boot_at_bound = boot_at_bound)
  )
}

#' @export
print.ffarfima_memory <- function(x, ...) {
  cat("Memory estimate d = ", format(x$d, digits = 4), "\n", sep = "")
  cat("  ", memory_label(x), "\n", sep = "")
  if (is.finite(x$bias)) {
    cat("  uncorrected d = ", format(x$d_initial, digits = 4),
        ", bootstrap bias = ", format(x$bias, digits = 3),
        ", bootstrap sd = ", format(x$sd, digits = 3), "\n", sep = "")
  }
  if (!isTRUE(attr(x, "diagnostics")$admissible)) {
    cat("  inadmissible: the bias-corrected estimate is below zero\n")
  }
  cat(memory_notes(x, "  "), sep = "")
  invisible(x)
}

memory_notes <- function(m, indent) {
  notes <- character(0)
  dg <- attr(m, "diagnostics")
  if (isTRUE(dg$at_bound)) {
    notes <- c(notes, paste0(indent, "note: LPWN stopped at the edge of its search range (0.01 or 0.99 per difference),\n",
                             indent, "      so d is imprecise; compare with method = \"lw\" or supply d\n"))
  }
  if (isTRUE(dg$boot_at_bound > 0.25)) {
    notes <- c(notes, paste0(indent, "note: ", round(100 * dg$boot_at_bound), "% of bootstrap estimates are at that edge,",
                             " so the bias correction is unreliable\n"))
  }
  notes
}
