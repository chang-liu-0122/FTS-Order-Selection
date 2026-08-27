suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(tibble)
  library(demography)
  library(ftsa)
  library(fdaACF)
  library(imputeTS)
})

project_root <- normalizePath(file.path(dirname(commandArgs(trailingOnly = FALSE)[1]), ".."))
if (!file.exists(file.path(project_root, "Empirical Study", "demog_list_asfr.rds"))) {
  project_root <- normalizePath(".")
}

method_levels <- c("Benchmark", "Functional IC", "CV")
method_key <- tibble::tribble(
  ~order_selection_method, ~internal_method, ~Method,
  "adaptive_ic",  "far",   "Benchmark",
  "adaptive_ic",  "farma", "Functional IC",
  "cv_expanding", "farma", "CV"
)

attach_method <- function(data) {
  data |>
    rename(internal_method = method) |>
    inner_join(method_key, by = c("order_selection_method", "internal_method")) |>
    mutate(Method = factor(Method, levels = method_levels))
}

fertility_change <- function(x) {
  previous <- x[, -ncol(x), drop = FALSE]
  current <- x[, -1, drop = FALSE]
  denominator <- previous + current
  change <- 2 * (previous - current) / denominator
  change[denominator == 0] <- 0
  change
}

extract_horizon <- function(result, fitted, observed, grid, series, h = 10, alpha = 0.2) {
  forecasts <- attach_method(result$error_detail)
  forecasts$n_horizons <- vapply(forecasts$error, ncol, integer(1))
  forecasts <- forecasts |> filter(n_horizons >= h)

  latest_t0 <- forecasts |>
    distinct(t0, Method) |>
    count(t0) |>
    filter(n == length(method_levels)) |>
    summarise(t0 = max(t0)) |>
    pull(t0)

  calibration <- attach_method(result$calibration_detail) |>
    filter(abs(.data$alpha - alpha) < 1e-10)

  rows <- vector("list", length(method_levels))
  target_column <- latest_t0 + h

  for (i in seq_along(method_levels)) {
    method_now <- method_levels[i]
    forecast_row <- forecasts |>
      filter(t0 == latest_t0, Method == method_now) |>
      slice(1)
    calibration_row <- calibration |>
      filter(Method == method_now) |>
      slice(1)

    error_matrix <- as.matrix(forecast_row$error[[1]])
    forecast <- fitted[, target_column] + error_matrix[, h]
    half_width <- calibration_row$zeta[[1]][h] *
      as.numeric(calibration_row$gamma[[1]][[h]])

    rows[[i]] <- tibble(
      Series = series,
      Method = factor(method_now, levels = method_levels),
      x = grid,
      observed = observed[, target_column],
      forecast = forecast,
      lower = forecast - half_width,
      upper = forecast + half_width
    )
  }

  bind_rows(rows)
}

fer <- readRDS(file.path(project_root, "Empirical Study", "demog_list_asfr.rds"))
fertility_raw <- as.matrix(fer$FRATNP$raw$rate$FRATNP_ASFR)
fertility_smooth <- as.matrix(fer$FRATNP$smooth$rate$FRATNP_ASFR)
fertility_smooth[!is.finite(fertility_smooth)] <- fertility_raw[!is.finite(fertility_smooth)]

fertility_plot_data <- extract_horizon(
  result = readRDS(file.path(
    project_root,
    "Empirical Study/results/france_fertility_change_forecast.rds"
  )),
  fitted = fertility_change(fertility_smooth),
  observed = fertility_change(fertility_raw),
  grid = fer$FRATNP$raw$age,
  series = "Fertility-rate change"
)

mortality <- demography::read.demogdata(
  file = file.path(project_root, "Empirical Study/Mortality Data/FRATNP_Mx_1x1.txt"),
  popfile = file.path(
    project_root,
    "Empirical Study/Mortality Data/FRATNP_Exposures_1x1.txt"
  ),
  type = "mortality",
  label = "France"
) |>
  demography::extract.ages(0:95)

for (sex in c("female", "male")) {
  mortality$rate[[sex]][
    mortality$rate[[sex]] <= 0 | mortality$rate[[sex]] > 1
  ] <- NA_real_
  for (age in seq_len(nrow(mortality$rate[[sex]]))) {
    mortality$rate[[sex]][age, ] <- imputeTS::na_interpolation(
      mortality$rate[[sex]][age, ], option = "linear"
    )
  }
}

mortality_smooth <- demography::smooth.demogdata(
  mortality,
  method = "mspline",
  b = 65,
  weight = TRUE
)
keep_years <- mortality$year >= 1861

mortality_plot_rows <- list()
for (sex in c("female", "male")) {
  mortality_plot_rows[[sex]] <- extract_horizon(
    result = readRDS(file.path(
      project_root,
      paste0("Empirical Study/results/france_", sex, "_mortality_forecast.rds")
    )),
    fitted = log10(mortality_smooth$rate[[sex]][, keep_years, drop = FALSE]),
    observed = log10(mortality$rate[[sex]][, keep_years, drop = FALSE]),
    grid = mortality$age,
    series = paste(tools::toTitleCase(sex), "log mortality")
  )
}

plot_data <- bind_rows(
  fertility_plot_data,
  mortality_plot_rows$female,
  mortality_plot_rows$male
) |>
  mutate(
    Series = factor(
      Series,
      levels = c("Fertility-rate change", "Female log mortality", "Male log mortality")
    )
  )

observed_data <- plot_data |>
  distinct(Series, x, observed)

method_colours <- c(
  "Benchmark" = "#1B9E77",
  "Functional IC" = "#D95F02",
  "CV" = "#7570B3"
)
method_linetypes <- c(
  "Benchmark" = "longdash",
  "Functional IC" = "solid",
  "CV" = "dotdash"
)

make_forecast_plot <- function(data) {
  observed <- data |>
    distinct(Series, x, observed)

  ggplot() +
  geom_ribbon(
    data = data,
    aes(x = x, ymin = lower, ymax = upper, fill = Method, group = Method),
    alpha = 0.70,
    colour = NA
  ) +
  geom_line(
    data = observed,
    aes(x = x, y = observed),
    colour = "#111111",
    linewidth = 0.9
  ) +
  geom_line(
    data = data,
    aes(x = x, y = forecast, colour = Method, linetype = Method, group = Method),
    linewidth = 0.8
  ) +
  facet_wrap(~Series, ncol = 1, scales = "free") +
  scale_colour_manual(values = method_colours) +
  scale_fill_manual(values = method_colours) +
  scale_linetype_manual(values = method_linetypes) +
  labs(
    x = "Age",
    y = NULL,
    colour = "Method",
    fill = "Method",
    linetype = "Method"
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "top",
    legend.box = "horizontal",
    legend.title = element_text(face = "bold"),
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )
}

fertility_forecast_plot <- make_forecast_plot(
  filter(plot_data, Series == "Fertility-rate change")
)

mortality_forecast_plot <- make_forecast_plot(
  filter(plot_data, Series != "Fertility-rate change")
)

output_dir <- file.path(
  project_root,
  "Empirical Study/Thesis Draft/Figures/Empirical"
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

grDevices::cairo_pdf(
  file.path(output_dir, "france_fertility_rates.pdf"),
  width = 12 / 2.54,
  height = 10 / 2.54,
  family = "Helvetica"
)
plot(
  fts(fer$FRATNP$raw$age, fertility_smooth),
  xlab = "Age",
  ylab = "Fertility rate"
)
grDevices::dev.off()

grDevices::cairo_pdf(
  file.path(output_dir, "france_fertility_changes.pdf"),
  width = 12 / 2.54,
  height = 10 / 2.54,
  family = "Helvetica"
)
plot(
  fts(fer$FRATNP$raw$age, fertility_change(fertility_smooth)),
  xlab = "Age",
  ylab = "Fertility-rate change"
)
grDevices::dev.off()

ggsave(
  file.path(output_dir, "france_h10_forecast_example.pdf"),
  make_forecast_plot(plot_data),
  width = 7.2,
  height = 9.4,
  units = "in",
  device = cairo_pdf
)

ggsave(
  file.path(output_dir, "france_fertility_h10_forecast.pdf"),
  fertility_forecast_plot,
  width = 7.2,
  height = 4.1,
  units = "in",
  device = cairo_pdf
)

ggsave(
  file.path(output_dir, "france_mortality_h10_forecast.pdf"),
  mortality_forecast_plot,
  width = 7.2,
  height = 6.6,
  units = "in",
  device = cairo_pdf
)

ggsave(
  file.path(output_dir, "france_h10_forecast_example.png"),
  make_forecast_plot(plot_data),
  width = 7.2,
  height = 9.4,
  units = "in",
  dpi = 220
)

extract_all_horizons <- function(
    result,
    fitted,
    observed,
    grid,
    series,
    alpha_set = c(0.20, 0.05)) {
  forecasts <- attach_method(result$error_detail)
  forecasts$n_horizons <- vapply(forecasts$error, ncol, integer(1))
  forecasts <- forecasts |>
    filter(n_horizons >= 10)

  latest_t0 <- forecasts |>
    distinct(t0, Method) |>
    count(t0) |>
    filter(n == length(method_levels)) |>
    summarise(t0 = max(t0)) |>
    pull(t0)

  calibration <- attach_method(result$calibration_detail) |>
    filter(alpha %in% alpha_set)

  rows <- list()
  row_index <- 0L

  for (method_now in method_levels) {
    forecast_row <- forecasts |>
      filter(t0 == latest_t0, Method == method_now) |>
      slice(1)

    error_matrix <- as.matrix(forecast_row$error[[1]])
    fitted_matrix <- fitted[, latest_t0 + seq_len(10), drop = FALSE]
    observed_matrix <- observed[, latest_t0 + seq_len(10), drop = FALSE]
    forecast_matrix <- fitted_matrix + error_matrix[, seq_len(10), drop = FALSE]

    method_calibration <- calibration |>
      filter(Method == method_now)

    for (calibration_index in seq_len(nrow(method_calibration))) {
      alpha_now <- method_calibration$alpha[calibration_index]
      gamma_now <- method_calibration$gamma[[calibration_index]]
      zeta_now <- method_calibration$zeta[[calibration_index]]

      for (h_now in seq_len(10)) {
        half_width <- zeta_now[h_now] * as.numeric(gamma_now[[h_now]])
        row_index <- row_index + 1L
        rows[[row_index]] <- tibble(
          Series = series,
          Method = factor(method_now, levels = method_levels),
          `Interval level` = scales::percent(1 - alpha_now, accuracy = 1),
          h = h_now,
          panel = factor(paste0("h = ", h_now), levels = paste0("h = ", 1:10)),
          age = grid,
          observed = observed_matrix[, h_now],
          forecast = forecast_matrix[, h_now],
          lower = forecast_matrix[, h_now] - half_width,
          upper = forecast_matrix[, h_now] + half_width
        )
      }
    }
  }

  bind_rows(rows)
}

make_all_horizon_plot <- function(data, interval_level, y_label) {
  band_data <- data |>
    filter(`Interval level` == interval_level)
  observed_data <- band_data |>
    distinct(panel, h, age, observed)
  point_data <- band_data |>
    distinct(panel, h, age, Method, forecast)

  ggplot() +
    geom_ribbon(
      data = band_data,
      aes(age, ymin = lower, ymax = upper, fill = Method, group = Method),
      alpha = 0.70,
      colour = NA
    ) +
    geom_line(
      data = observed_data,
      aes(age, observed),
      colour = "#111111",
      linewidth = 1.1
    ) +
    geom_line(
      data = point_data,
      aes(age, forecast, colour = Method, linetype = Method, group = Method),
      linewidth = 0.9
    ) +
    facet_wrap(~panel, ncol = 5, scales = "free_y") +
    scale_colour_manual(name = "Point forecast", values = method_colours, drop = FALSE) +
    scale_linetype_manual(name = "Point forecast", values = method_linetypes, drop = FALSE) +
    scale_fill_manual(
      name = paste(interval_level, "prediction band"),
      values = method_colours,
      drop = FALSE
    ) +
    guides(
      colour = guide_legend(order = 1, override.aes = list(linewidth = 1.2)),
      linetype = guide_legend(order = 1),
      fill = guide_legend(order = 2, override.aes = list(alpha = 0.7))
    ) +
    labs(
      x = "Age",
      y = y_label
    ) +
    theme_bw(base_size = 15) +
    theme(
      legend.position = "top",
      legend.box = "vertical",
      legend.title = element_text(face = "bold"),
      strip.text = element_text(face = "bold"),
      panel.grid.minor = element_blank(),
      plot.margin = margin(8, 8, 8, 16)
    )
}

empirical_plot_dir <- file.path(project_root, "Empirical Study/Empirical Plots")
dir.create(empirical_plot_dir, recursive = TRUE, showWarnings = FALSE)

fertility_all_horizons <- extract_all_horizons(
  result = readRDS(file.path(
    project_root,
    "Empirical Study/results/france_fertility_change_forecast.rds"
  )),
  fitted = fertility_change(fertility_smooth),
  observed = fertility_change(fertility_raw),
  grid = fer$FRATNP$raw$age,
  series = "France fertility"
)

mortality_all_horizons <- list()
for (sex in c("female", "male")) {
  mortality_all_horizons[[sex]] <- extract_all_horizons(
    result = readRDS(file.path(
      project_root,
      paste0("Empirical Study/results/france_", sex, "_mortality_forecast.rds")
    )),
    fitted = log10(mortality_smooth$rate[[sex]][, keep_years, drop = FALSE]),
    observed = log10(mortality$rate[[sex]][, keep_years, drop = FALSE]),
    grid = mortality$age,
    series = paste("France", sex, "mortality")
  )
}

for (interval_now in c("80%", "95%")) {
  fertility_plot <- make_all_horizon_plot(
    fertility_all_horizons,
    interval_now,
    "Increment fertility rate"
  )

  ggsave(
    filename = file.path(
      empirical_plot_dir,
      paste0("France_fertility_increment_forecasts_", sub("%", "", interval_now), ".pdf")
    ),
    plot = fertility_plot,
    width = 12,
    height = 10,
    dpi = 1200,
    bg = "white"
  )

  for (sex in c("female", "male")) {
    mortality_plot <- make_all_horizon_plot(
      mortality_all_horizons[[sex]],
      interval_now,
      "Log mortality rate"
    )

    ggsave(
      filename = file.path(
        empirical_plot_dir,
        paste0(
          "France_",
          sex,
          "_mortality_forecasts_",
          sub("%", "", interval_now),
          ".pdf"
        )
      ),
      plot = mortality_plot,
      width = 12,
      height = 10,
      dpi = 1200,
      bg = "white"
    )
  }
}

for (sex in c("female", "male")) {
  identification <- fdaACF::obtain_FACF(
    t(log10(mortality_smooth$rate[[sex]][, keep_years, drop = FALSE])),
    v = mortality$age,
    nlags = 20,
    ci = 0.95,
    estimation = "MC",
    figure = FALSE
  )

  grDevices::cairo_pdf(
    file.path(
      empirical_plot_dir,
      paste0("France_", sex, "_mortality_facf.pdf")
    ),
    width = 12 / 2.54,
    height = 10 / 2.54,
    family = "Helvetica"
  )
  par(
    mgp = c(2.2, 0.45, 0),
    tcl = -0.4,
    mar = c(3.2, 3.5, 1.8, 1) + 0.1,
    pch = 1
  )
  plot(
    seq_along(identification$rho),
    identification$rho,
    type = "h",
    lwd = 8,
    lend = 2,
    xlab = "Lag",
    ylab = paste(tools::toTitleCase(sex), "mortality fACF"),
    xlim = c(0, length(identification$rho) + 1),
    ylim = c(0, 1),
    yaxt = "n",
    las = 1,
    yaxs = "i",
    xaxs = "i"
  )
  grid(lty = 1)
  lines(
    seq_along(identification$rho),
    identification$rho,
    type = "h",
    col = "lightgrey",
    lwd = 6,
    lend = 2
  )
  abline(h = identification$Blueline, col = "#0073C2FF", lwd = 4, lty = 2)
  legend(
    "topright",
    legend = "i.i.d. bound (95 % conf.)",
    col = "#0073C2FF",
    lty = 2,
    lwd = 4
  )
  axis(2, at = seq(0, 1, by = 0.2), labels = seq(0, 1, by = 0.2))
  box()
  dev.off()

  message(
    tools::toTitleCase(sex),
    " mortality fACF at lag 20: ",
    sprintf("%.3f", identification$rho[20])
  )
}
