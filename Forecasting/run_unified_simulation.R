# Directory setup
file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)

project_root <- if (length(file_arg) == 1) {
  # Rscript encodes spaces in the script path as "~+~".
  script_path <- gsub("~+~", " ", sub("^--file=", "", file_arg), fixed = TRUE)
  normalizePath(file.path(dirname(script_path), ".."))
} else if (file.exists("Forecasting/unified_compare_sim.R")) {
  normalizePath(".")
} else {
  normalizePath("..")
}

setwd(project_root)

if (!file.exists("Forecasting/unified_compare_sim.R")) {
  stop("Could not locate the project root; run this script from inside the project.")
}

library(longmemo)
library(matrixStats)
library(far)
library(FCVAR)
library(fda)
library(ftsa)
library(refund)
library(MASS)
library(fdaACF)
library(fracdiff)
library(wwntests)
library(vars)
library(tidyverse)
library(future)
library(future.apply)
library(progressr)

source("Order Selection/Marginal BIC.R")
source("Order Selection/CV order selection.R")
source("DGP non-stationary FTS/sim_FARMA.R")
source("Fractional Differencing/local.W.R")
source("Fractional Differencing/Covariance Matrix.R")
source("Fractional Differencing/LPWN.R")
source("Fractional Differencing/sieve_bootstrap.R")
source("Fractional Differencing/Sieve-bootstrapped corrected LW estimation.R")
source("Operator Estimation/Kuenzer FARMA Estimation.R")
source("Forecasting/Forecasting Helper.R")
source("Forecasting/forecast_compare_sim.R")
source("Forecasting/unified_compare_sim.R")



n_rep <- 100
future::plan(future::multisession, workers = 100)
progressr::handlers("progress")


# TEST CASE
# grid_pipeline <- data.frame(
#   T = 500,
#   d = 0.6,
#   p = 1L,
#   q = 1L,
#   ar_strength = 0.5,
#   ma_strength = 0.5
# )

grid_pipeline <- expand.grid(
  T = c(500, 1000),
  d = c(0.2, 0.4, 0.6, 0.8, 1, 1.2, 1.4),
  ar_strength = c(0.2, 0.5, 0.8),
  ma_strength = c(0.2, 0.5, 0.8),
  model = c("1_0", "0_1", "1_1", "2_0")
) %>%
  dplyr::mutate(
    p = dplyr::case_when(
      model == "1_0" ~ 1L,
      model == "0_1" ~ 0L,
      model == "1_1" ~ 1L,
      model == "2_0" ~ 2L
    ),
    q = dplyr::case_when(
      model == "1_0" ~ 0L,
      model == "0_1" ~ 1L,
      model == "1_1" ~ 1L,
      model == "2_0" ~ 0L
    ),
    ar_strength = dplyr::if_else(p == 0, 0, ar_strength),
    ma_strength = dplyr::if_else(q == 0, 0, ma_strength)
  ) %>%
  dplyr::distinct() %>%
  dplyr::select(T, d, p, q, ar_strength, ma_strength) %>%
  dplyr::arrange(T, d, p, q, ar_strength, ma_strength)

out_ic <- parallel_repeat_full_forecasting_pipeline_compact_order_method(
  n_rep = n_rep,
  grid = grid_pipeline,
  m_const = 0.75,
  fracdiff_trim = 200,
  far_benchmark_d = 1,
  window_types = c("expanding", "fixed"),
  order_selection_method = c("both"),
  cv_order_window_type = c("expanding")
)

saveRDS(out_ic, "Forecasting/unified_raw.rds")

pipeline_metrics <- extract_required_metrics_from_compact_order_method(
  out_ic,
  summary_stat = "both",
  keep_replication_details = TRUE
)

saveRDS(pipeline_metrics, "Forecasting/unified_summary.rds")
