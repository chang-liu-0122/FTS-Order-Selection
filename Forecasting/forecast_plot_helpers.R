ratio_axis_limits <- function(values, reference = 1, padding = 0.06) {
  value_range <- range(c(values, reference), na.rm = TRUE)
  pad <- padding * diff(value_range)

  if (!is.finite(pad) || pad == 0) {
    pad <- 0.05
  }

  c(value_range[1] - pad, value_range[2] + pad)
}

make_interval_plot_summary <- function(df, alpha_show = 0.05) {
  T_name <- if ("T" %in% names(df)) {
    "T"
  } else if ("n" %in% names(df)) {
    "n"
  } else {
    stop("Need a `T` or `n` column.")
  }

  if (!is.null(alpha_show) && !"alpha" %in% names(df)) {
    stop("Cannot find the `alpha` column.")
  }

  if (!"method" %in% names(df)) {
    stop("Cannot find the `method` column.")
  }

  df_use <- df

  if (!is.null(alpha_show)) {
    df_use <- dplyr::filter(
      df_use,
      abs(as.numeric(.data$alpha) - alpha_show) < 1e-10
    )
  }

  dplyr::mutate(
    df_use,
    T_value = as.numeric(.data[[T_name]]),
    T_factor = factor(T_value, levels = sort(unique(T_value))),
    h = as.numeric(.data$h),
    ar_strength = as.numeric(.data$ar_strength),
    ma_strength = as.numeric(.data$ma_strength),
    p_true = as.integer(.data$p_true),
    q_true = as.integer(.data$q_true),
    method = dplyr::case_when(
      stringr::str_detect(tolower(as.character(.data$method)), "farma") ~ "farma",
      stringr::str_detect(
        tolower(as.character(.data$method)),
        "far\\(1\\)|far1|far_1|far$"
      ) ~ "far",
      TRUE ~ tolower(as.character(.data$method))
    ),
    window_type = dplyr::case_when(
      stringr::str_detect(tolower(as.character(.data$window_type)), "exp") ~ "Expanding",
      stringr::str_detect(tolower(as.character(.data$window_type)), "fix") ~ "Fixed",
      TRUE ~ as.character(.data$window_type)
    ),
    window_type = factor(window_type, levels = c("Expanding", "Fixed")),
    model = paste0("(", p_true, ",", q_true, ")"),
    model = factor(model, levels = c("(1,0)", "(0,1)", "(2,0)", "(1,1)")),
    model_label = dplyr::case_when(
      model == "(1,0)" ~ "FARFIMA(1,d,0)",
      model == "(0,1)" ~ "FARFIMA(0,d,1)",
      model == "(2,0)" ~ "FARFIMA(2,d,0)",
      model == "(1,1)" ~ "FARFIMA(1,d,1)",
      TRUE ~ as.character(model)
    ),
    psi_label = factor(
      paste0("psi==", ar_strength),
      levels = paste0("psi==", sort(unique(ar_strength)))
    ),
    nu_label = factor(
      paste0("nu==", ma_strength),
      levels = paste0("nu==", sort(unique(ma_strength)))
    )
  )
}

# Shared visual specification for paper figures.  The long benchmark label is
# retained for joins, while the legend displays the compact model name.
paper_method_colours <- c(
  "Benchmark (FARFIMA(1,1,0))" = "#1B9E77",
  "FARFIMA(1,1,0)" = "#1B9E77",
  "Functional IC" = "#D95F02",
  "CV" = "#7570B3"
)

paper_method_linetypes <- c(
  "Benchmark (FARFIMA(1,1,0))" = "longdash",
  "FARFIMA(1,1,0)" = "longdash",
  "Functional IC" = "solid",
  "CV" = "dashed"
)

paper_method_shapes <- c(
  "Benchmark (FARFIMA(1,1,0))" = 15,
  "FARFIMA(1,1,0)" = 15,
  "Functional IC" = 16,
  "CV" = 17
)

paper_method_labels <- c(
  "Benchmark (FARFIMA(1,1,0))" = "FARFIMA(1,1,0)",
  "FARFIMA(1,1,0)" = "FARFIMA(1,1,0)",
  "Functional IC" = "Functional IC",
  "CV" = "CV"
)

paper_legend_theme <- function(base_size = 15, position = "top") {
  ggplot2::theme(
    legend.position = position,
    legend.direction = "horizontal",
    legend.box = "horizontal",
    legend.justification = "center",
    legend.box.just = "center",
    legend.title = ggplot2::element_text(
      size = base_size * 0.92,
      face = "bold"
    ),
    legend.text = ggplot2::element_text(size = base_size * 0.88),
    legend.key.width = grid::unit(1.65, "lines"),
    legend.key.height = grid::unit(0.85, "lines"),
    legend.spacing.x = grid::unit(0.55, "lines"),
    legend.spacing.y = grid::unit(0.10, "lines"),
    legend.margin = ggplot2::margin(0, 0, 0, 0),
    legend.box.margin = ggplot2::margin(0, 0, 0, 0)
  )
}

paper_plot_theme <- function(base_size = 17, legend_position = "top") {
  ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      plot.title = ggplot2::element_blank(),
      plot.subtitle = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(
        colour = "grey25",
        fill = NA,
        linewidth = 0.5
      ),
      strip.background = ggplot2::element_rect(
        fill = "grey92",
        colour = "grey25",
        linewidth = 0.5
      ),
      strip.text = ggplot2::element_text(
        face = "bold",
        size = base_size * 0.88
      ),
      axis.title = ggplot2::element_text(size = base_size),
      axis.text = ggplot2::element_text(size = base_size * 0.82),
      panel.grid = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(5, 5, 5, 5)
    ) +
    paper_legend_theme(base_size, legend_position)
}

get_visible_legend <- function(plot) {
  plot_grob <- ggplot2::ggplotGrob(plot)
  guide_index <- which(grepl("^guide-box", plot_grob$layout$name))
  guide_index <- guide_index[
    !vapply(
      plot_grob$grobs[guide_index],
      inherits,
      logical(1),
      what = "zeroGrob"
    )
  ]

  if (length(guide_index) == 0L) {
    return(NULL)
  }

  plot_grob$grobs[[guide_index[1]]]
}

paper_method_scales <- function(methods,
                                title = "Forecasting method",
                                drop = FALSE) {
  methods <- as.character(methods)

  list(
    ggplot2::scale_colour_manual(
      name = title,
      values = paper_method_colours,
      breaks = methods,
      labels = paper_method_labels[methods],
      drop = drop
    ),
    ggplot2::scale_linetype_manual(
      name = title,
      values = paper_method_linetypes,
      breaks = methods,
      labels = paper_method_labels[methods],
      drop = drop
    ),
    ggplot2::scale_shape_manual(
      name = title,
      values = paper_method_shapes,
      breaks = methods,
      labels = paper_method_labels[methods],
      drop = drop
    )
  )
}

make_custom_method_legend <- function(methods,
                                      colours,
                                      linetypes,
                                      shapes,
                                      labels = methods,
                                      title = "Forecasting method",
                                      base_size = 21,
                                      point_size = 4.5,
                                      line_width = 1.5) {
  methods <- as.character(methods)
  labels <- stats::setNames(as.character(labels), methods)
  legend_data <- expand.grid(
    method = factor(methods, levels = methods),
    x = c(1, 2),
    KEEP.OUT.ATTRS = FALSE
  )
  legend_data$y <- as.numeric(legend_data$method)

  legend_plot <- ggplot2::ggplot(
    legend_data,
    ggplot2::aes(
      x = x,
      y = y,
      colour = method,
      linetype = method,
      shape = method,
      group = method
    )
  ) +
    ggplot2::geom_line(linewidth = line_width) +
    ggplot2::geom_point(size = point_size) +
    ggplot2::scale_colour_manual(
      name = title,
      values = colours,
      breaks = methods,
      labels = labels[methods],
      drop = FALSE
    ) +
    ggplot2::scale_linetype_manual(
      name = title,
      values = linetypes,
      breaks = methods,
      labels = labels[methods],
      drop = FALSE
    ) +
    ggplot2::scale_shape_manual(
      name = title,
      values = shapes,
      breaks = methods,
      labels = labels[methods],
      drop = FALSE
    ) +
    ggplot2::guides(
      colour = ggplot2::guide_legend(
        nrow = 1,
        byrow = TRUE,
        override.aes = list(linewidth = line_width, size = point_size)
      ),
      linetype = "none",
      shape = "none"
    ) +
    ggplot2::theme_void(base_size = base_size) +
    paper_legend_theme(base_size = base_size, position = "right")

  cowplot::get_legend(legend_plot)
}

make_method_legend <- function(methods = c("Functional IC", "CV"),
                               title = "Forecasting method",
                               base_size = 21,
                               point_size = 4.5,
                               line_width = 1.5) {
  make_custom_method_legend(
    methods = methods,
    colours = paper_method_colours,
    linetypes = paper_method_linetypes,
    shapes = paper_method_shapes,
    labels = paper_method_labels[methods],
    title = title,
    base_size = base_size,
    point_size = point_size,
    line_width = line_width
  )
}

make_two_method_legend <- function(method_cols,
                                   method_linetypes,
                                   method_shapes,
                                   title = "Forecasting method",
                                   title_size = 7.5,
                                   text_size = 7,
                                   point_size = 4.5,
                                   line_width = 1.5,
                                   base_size = 25) {
  make_method_legend(
    methods = c("Functional IC", "CV"),
    title = title,
    base_size = text_size * ggplot2::.pt,
    point_size = point_size,
    line_width = line_width
  )
}

assemble_ratio_panel_grid <- function(df,
                                      models,
                                      titles,
                                      value_col,
                                      plot_model,
                                      legend,
                                      y_axis_label,
                                      file_name = NULL,
                                      width = 16,
                                      height = 12,
                                      y_limits = NULL,
                                      common_y = TRUE,
                                      shared_label_size = 30,
                                      dpi = 1200) {
  stopifnot(
    length(models) == 4L,
    length(titles) == 4L,
    value_col %in% names(df)
  )

  df_use <- df[as.character(df$model) %in% models, , drop = FALSE]
  stopifnot(nrow(df_use) > 0L)

  if (is.null(y_limits) && common_y) {
    y_limits <- ratio_axis_limits(df_use[[value_col]])
  }

  panels <- Map(
    function(model, title) {
      plot_model(
        df = df_use,
        model_show = model,
        title = title,
        value_col = value_col,
        y_limits = y_limits
      )
    },
    models,
    titles
  )

  main_grid <- cowplot::plot_grid(
    plotlist = panels,
    ncol = 2,
    align = "hv",
    axis = "tblr"
  )

  labelled_grid <- cowplot::ggdraw() +
    cowplot::draw_plot(
      main_grid,
      x = 0.078,
      y = 0.085,
      width = 0.90,
      height = 0.89
    ) +
    cowplot::draw_label(
      y_axis_label,
      x = 0.020,
      y = 0.53,
      angle = 90,
      size = shared_label_size
    ) +
    cowplot::draw_label(
      "Forecast horizon h",
      x = 0.53,
      y = 0.022,
      size = shared_label_size
    )

  final_plot <- if (is.null(legend)) {
    labelled_grid
  } else {
    cowplot::plot_grid(
      legend,
      labelled_grid,
      ncol = 1,
      rel_heights = c(0.13, 1)
    )
  }

  if (!is.null(file_name)) {
    dir.create(dirname(file_name), recursive = TRUE, showWarnings = FALSE)
    ggplot2::ggsave(
      filename = file_name,
      plot = final_plot,
      width = width,
      height = height,
      dpi = dpi,
      bg = "white"
    )
  }

  final_plot
}
