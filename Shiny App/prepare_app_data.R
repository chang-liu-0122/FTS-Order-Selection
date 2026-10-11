suppressMessages(library(tidyverse))

root <- "/Users/liuchang/Desktop/FTS Order Selection Thesis"
app_dir <- file.path(root, "Shiny App")

strength_label <- function(p, q, ar, ma) {
  case_when(
    p > 0 & q > 0 ~ sprintf("\u03c8 = %.1f, \u03bd = %.1f", ar, ma),
    p > 0 ~ sprintf("\u03c8 = %.1f", ar),
    q > 0 ~ sprintf("\u03bd = %.1f", ma),
    TRUE ~ "none"
  )
}

memory <- bind_rows(
  readRDS(file.path(root, "Fractional Differencing/covariance_matrix_comparison_FAR.rds")) |>
    mutate(model = "FARFIMA(1,d,0)", ma_strength = 0),
  readRDS(file.path(root, "Fractional Differencing/covariance_matrix_comparison_FARMA.rds")) |>
    rename(ar_strength = strength_ar, ma_strength = strength_ma) |>
    mutate(model = "FARFIMA(1,d,1)")
) |>
  pivot_longer(c(sample_whittle, sample_peng, lrc_RS_whittle, lrc_RS_peng, lrc_LRS_whittle, lrc_LRS_peng), names_to = "key", values_to = "d_hat") |>
  mutate(
    operator = factor(case_when(startsWith(key, "sample") ~ "Covar", grepl("_RS_", key) ~ "RS", TRUE ~ "LRS"), levels = c("Covar", "RS", "LRS")),
    estimator = factor(ifelse(endsWith(key, "whittle"), "LW", "DFA"), levels = c("LW", "DFA"))
  ) |>
  transmute(model, T, d_true, ar_strength, ma_strength, seed, operator, estimator, d_hat)

order_sel <- readRDS(file.path(root, "Order Selection/sim_order_compare.rds"))$simple_result |>
  transmute(
    T = n, p_true, q_true, ar_strength, ma_strength, seed,
    Method = factor(recode(method, Adaptive_IC = "Functional IC", WN_Test = "Functional Box\u2013Jenkins"), levels = c("Functional IC", "Functional Box\u2013Jenkins")),
    p_hat, q_hat, correct = as.logical(correct)
  ) |>
  mutate(strength = strength_label(p_true, q_true, ar_strength, ma_strength))

sim_long <- function(res, study) {
  has_os <- "order_selection_method" %in% names(res$ratio_summary)
  fix <- function(df) {
    if (!has_os) df$order_selection_method <- "ideal"
    if (!"d_true" %in% names(df)) df$d_true <- 1
    df
  }
  keys <- c("n", "window_type", "d_true", "p_true", "q_true", "ar_strength", "ma_strength", "h")
  farma_label <- function(os) ifelse(os == "ideal", "Plug-in", ifelse(os == "adaptive_ic", "Functional IC", "CV"))
  rs <- fix(res$ratio_summary)
  point <- bind_rows(
    rs |> group_by(across(all_of(keys))) |>
      summarise(RMSPE = mean(RMSE_far), MAPE = mean(MAE_far), `RMSPE ratio` = 1, `MAPE ratio` = 1, .groups = "drop") |>
      mutate(Method = "FARFIMA(1,1,0)"),
    rs |> transmute(across(all_of(keys)), Method = farma_label(order_selection_method), RMSPE = RMSE_farma, MAPE = MAE_farma, `RMSPE ratio` = RMSE_ratio, `MAPE ratio` = MAE_ratio)
  ) |>
    pivot_longer(c(RMSPE, MAPE, `RMSPE ratio`, `MAPE ratio`), names_to = "metric", values_to = "value") |>
    mutate(Rule = "none", alpha = NA_real_)
  iv <- fix(res$interval_summary)
  interval <- bind_rows(
    iv |> filter(method == "far") |> group_by(across(all_of(c(keys, "alpha")))) |>
      summarise(across(c(ECP, CPD, Score, ECP_scores, CPD_scores, Score_scores, ECP_aci, CPD_aci, Score_aci), mean), .groups = "drop") |>
      mutate(Method = "FARFIMA(1,1,0)"),
    iv |> filter(method == "farma") |> mutate(Method = farma_label(order_selection_method)) |>
      dplyr::select(all_of(c(keys, "alpha")), Method, ECP, CPD, Score, ECP_scores, CPD_scores, Score_scores, ECP_aci, CPD_aci, Score_aci)
  ) |>
    pivot_longer(c(ECP, CPD, Score, ECP_scores, CPD_scores, Score_scores, ECP_aci, CPD_aci, Score_aci), names_to = "key", values_to = "value") |>
    mutate(
      Rule = case_when(endsWith(key, "_scores") ~ "NS", endsWith(key, "_aci") ~ "AC", TRUE ~ "SR"),
      metric = recode(sub("_scores$|_aci$", "", key), Score = "IS")
    ) |>
    dplyr::select(-key)
  ws <- fix(res$wis_summary)
  wis <- bind_rows(
    ws |> filter(method == "far") |> group_by(across(all_of(keys))) |>
      summarise(across(c(WIS, WIS_scores, WIS_aci), mean), .groups = "drop") |> mutate(Method = "FARFIMA(1,1,0)"),
    ws |> filter(method == "farma") |> mutate(Method = farma_label(order_selection_method)) |>
      dplyr::select(all_of(keys), Method, WIS, WIS_scores, WIS_aci)
  ) |>
    pivot_longer(c(WIS, WIS_scores, WIS_aci), names_to = "key", values_to = "value") |>
    mutate(Rule = case_when(endsWith(key, "_scores") ~ "NS", endsWith(key, "_aci") ~ "AC", TRUE ~ "SR"), metric = "WIS", alpha = NA_real_) |>
    dplyr::select(-key)
  bind_rows(point, interval, wis) |>
    rename(T = n, window = window_type) |>
    mutate(study = study, strength = strength_label(p_true, q_true, ar_strength, ma_strength))
}

ideal <- readRDS(file.path(root, "Forecasting/ideal_summary.rds"))
unified <- readRDS(file.path(root, "Forecasting/unified_summary.rds"))

sim_fc <- bind_rows(sim_long(ideal, "Ideal"), sim_long(unified, "Unified")) |>
  mutate(across(c(study, window, Method, Rule, metric, strength), factor))

unified_d <- unified$estimate_detail |>
  filter(order_selection_method == "adaptive_ic") |>
  transmute(T, d_true, p_true, q_true, ar_strength, ma_strength, d_hat, strength = strength_label(p_true, q_true, ar_strength, ma_strength))

unified_order <- unified$estimate_detail |>
  transmute(T, d_true, p_true, q_true, ar_strength, ma_strength,
            Method = factor(ifelse(order_selection_method == "adaptive_ic", "Functional IC", "CV"), levels = c("Functional IC", "CV")),
            p_hat = p_est, q_hat = q_est, correct = order_correct,
            strength = strength_label(p_true, q_true, ar_strength, ma_strength))

rm(ideal, unified)
gc()

setwd(root)
purled <- knitr::purl("Empirical Study/Empirical analysis.Rmd", output = file.path(tempdir(), "emp_purled.R"), quiet = TRUE)
lines <- readLines(purled)
starts <- grep("^## ----", lines)
labels <- sub("(,.*|-{2,}.*)$", "", sub("^## ----", "", lines[starts]), perl = TRUE)
ends <- c(starts[-1] - 1, length(lines))
params <- list(rerun_models = FALSE)
for (k in c("setup", "data-functions", "run-models", "metric-functions", "mortality-load-data", "mortality-run-models", "mortality-metrics")) {
  i <- which(labels == k)
  eval(parse(text = lines[starts[i]:ends[i]]), envir = globalenv())
}

emp_series <- list(
  `England and Wales fertility` = list(result = forecast_results[["England and Wales"]], fitted = fertility_data[["England and Wales"]]$smooth, observed = fertility_data[["England and Wales"]]$raw, scale = "Box-Cox transformed fertility rate"),
  `French female mortality` = list(result = mortality_forecast_results$Female, fitted = mortality_series$Female, observed = mortality_raw_series$Female, scale = "Log mortality rate"),
  `French male mortality` = list(result = mortality_forecast_results$Male, fitted = mortality_series$Male, observed = mortality_raw_series$Male, scale = "Log mortality rate")
)

emp_fc <- list()
emp_forecasts <- list()
emp_mult <- list()
emp_gamma <- list()
emp_models <- list()
for (s in names(emp_series)) {
  es <- emp_series[[s]]
  point <- summarise_point_metrics(es$result, es$fitted, es$observed)
  interval <- summarise_interval_metrics_by_rule(es$result, es$fitted, es$observed) |>
    mutate(Rule = recode(as.character(Rule), grid_sd = "SR", scores = "NS", aci = "AC"))
  wis <- bind_rows(lapply(c("SR", "NS", "AC"), function(r) summarise_wis(point, interval |> filter(Rule == r)) |> mutate(Rule = r)))
  emp_fc[[s]] <- bind_rows(
    point |> transmute(Method = as.character(Method), h, Rule = "none", alpha = NA_real_, RMSPE = RMSE, MAPE = MAE) |>
      pivot_longer(c(RMSPE, MAPE), names_to = "metric", values_to = "value"),
    interval |> transmute(Method = as.character(Method), h, Rule, alpha, ECP, CPD, IS = Score) |>
      pivot_longer(c(ECP, CPD, IS), names_to = "metric", values_to = "value"),
    wis |> transmute(Method = as.character(Method), h, Rule, alpha = NA_real_, metric = "WIS", value = WIS)
  ) |> mutate(series = s)

  ed <- es$result$error_detail |> attach_method()
  cal <- es$result$calibration_detail |> attach_method()
  fc_rows <- list()
  for (i in seq_len(nrow(ed))) {
    E <- ed$error[[i]]
    cols <- ed$t0[i] + seq_len(ncol(E))
    f <- es$fitted[, cols, drop = FALSE] + E
    o <- es$observed[, cols, drop = FALSE]
    fc_rows[[i]] <- tibble(
      Method = as.character(ed$Method[i]), t0 = ed$t0[i],
      h = rep(seq_len(ncol(E)), each = nrow(E)),
      year = as.integer(colnames(es$observed)[cols])[rep(seq_len(ncol(E)), each = nrow(E))],
      age = rep(parse_number(rownames(es$observed)), ncol(E)),
      forecast = as.vector(f), observed = as.vector(o)
    )
  }
  emp_forecasts[[s]] <- bind_rows(fc_rows) |> mutate(series = s)
  mult_rows <- list()
  gamma_rows <- list()
  for (j in seq_len(nrow(cal))) {
    H <- length(cal$zeta[[j]])
    gamma_rows[[j]] <- bind_rows(lapply(seq_len(H), function(hh) tibble(
      Method = as.character(cal$Method[j]), alpha = cal$alpha[j], h = hh,
      age = parse_number(rownames(es$observed)), gamma = as.numeric(cal$gamma[[j]][[hh]])
    )))
    aq <- cal$aci_Q[[j]]
    t0s <- unique(ed$t0[ed$Method == cal$Method[j]])
    mult_rows[[j]] <- bind_rows(
      expand_grid(t0 = t0s, h = seq_len(H)) |> mutate(Rule = "SR", mult = cal$zeta[[j]][h]),
      expand_grid(t0 = t0s, h = seq_len(H)) |> mutate(Rule = "NS", mult = 1 + 2 * cal$Q_hat[[j]][h]),
      expand_grid(t0 = t0s, h = seq_len(H)) |> mutate(Rule = "AC", mult = 1 + 2 * aq[cbind(match(as.character(t0), rownames(aq)), h)])
    ) |> mutate(Method = as.character(cal$Method[j]), alpha = cal$alpha[j], eta = cal$aci_eta[j])
  }
  emp_gamma[[s]] <- bind_rows(gamma_rows) |> mutate(series = s)
  emp_mult[[s]] <- bind_rows(mult_rows) |> mutate(series = s)
  ms <- es$result$model_summary |> filter(order_selection_method %in% c("adaptive_ic", "cv_expanding"))
  eta <- cal |> distinct(Method, alpha, aci_eta) |> mutate(Method = as.character(Method)) |>
    pivot_wider(names_from = alpha, values_from = aci_eta, names_prefix = "eta_")
  emp_models[[s]] <- tibble(
    Method = c("FARFIMA(1,1,0)", "Functional IC", "CV"),
    d = c(1, ms$d_hat[ms$order_selection_method == "adaptive_ic"], ms$d_hat[ms$order_selection_method == "cv_expanding"]),
    p = c(1L, ms$p_est[ms$order_selection_method == "adaptive_ic"], ms$p_est[ms$order_selection_method == "cv_expanding"]),
    q = c(0L, ms$q_est[ms$order_selection_method == "adaptive_ic"], ms$q_est[ms$order_selection_method == "cv_expanding"])
  ) |>
    left_join(eta |> mutate(Method = recode(Method, `Benchmark (FARFIMA(1,1,0))` = "FARFIMA(1,1,0)")), by = "Method") |>
    mutate(series = s)
}

bench <- function(x) recode(x, `Benchmark (FARFIMA(1,1,0))` = "FARFIMA(1,1,0)")
emp_fc <- bind_rows(emp_fc) |> mutate(Method = bench(Method))
emp_forecasts <- bind_rows(emp_forecasts) |> mutate(Method = bench(Method))
emp_gamma <- bind_rows(emp_gamma) |> mutate(Method = bench(Method))
emp_mult <- bind_rows(emp_mult) |> mutate(Method = bench(Method))
emp_models <- bind_rows(emp_models)

emp_coverage <- emp_forecasts |>
  inner_join(emp_gamma, by = c("series", "Method", "h", "age"), relationship = "many-to-many") |>
  inner_join(emp_mult, by = c("series", "Method", "alpha", "t0", "h"), relationship = "many-to-many") |>
  mutate(half = mult * gamma, covered = abs(observed - forecast) <= half) |>
  group_by(series, Method, Rule, alpha, h, t0, year) |>
  summarise(coverage = mean(covered), .groups = "drop")

check <- emp_coverage |> group_by(series, Method, Rule, alpha, h) |> summarise(ecp_bands = mean(coverage), .groups = "drop") |>
  inner_join(emp_fc |> filter(metric == "ECP") |> dplyr::select(series, Method, Rule, alpha, h, value), by = c("series", "Method", "Rule", "alpha", "h"))
cat("max |ECP from bands - ECP from Rmd functions|:", max(abs(check$ecp_bands - check$value)), "\n")

curve_rows <- list()
for (s in names(emp_series)) {
  es <- emp_series[[s]]
  for (kind in c("Raw", "Smoothed")) {
    m <- if (kind == "Raw") es$observed else es$fitted
    curve_rows[[paste(s, kind)]] <- tibble(
      series = s, kind = kind,
      age = rep(parse_number(rownames(m)), ncol(m)),
      year = rep(as.integer(colnames(m)), each = nrow(m)),
      value = as.vector(m)
    )
  }
}
emp_curves <- bind_rows(curve_rows)

set.seed(123)
emp_facf <- bind_rows(lapply(names(emp_series), function(s) {
  m <- emp_series[[s]]$fitted
  f <- fdaACF::obtain_FACF(t(m), v = parse_number(rownames(m)), nlags = 40, ci = 0.95, estimation = "MC", figure = FALSE)
  tibble(series = s, lag = seq_along(f$rho), rho = f$rho, bound = f$Blueline)
}))

emp_scale <- tibble(series = names(emp_series), scale = vapply(emp_series, function(x) x$scale, ""))

app_data <- list(
  memory = memory, order_sel = order_sel, sim_fc = sim_fc, unified_d = unified_d, unified_order = unified_order,
  emp_fc = emp_fc, emp_forecasts = emp_forecasts, emp_gamma = emp_gamma, emp_mult = emp_mult, emp_models = emp_models,
  emp_curves = emp_curves, emp_facf = emp_facf, emp_scale = emp_scale
)
saveRDS(app_data, file.path(app_dir, "app_data.rds"), compress = "xz")
for (nm in names(app_data)) cat(sprintf("%-15s %8d rows\n", nm, nrow(app_data[[nm]])))
cat("app_data.rds MB:", round(file.size(file.path(app_dir, "app_data.rds")) / 1e6, 2), "\n")
