ratio_axis_limits <- function(values, reference = 1, padding = 0.06) {
  value_range <- range(c(values, reference), na.rm = TRUE)
  pad <- padding * diff(value_range)

  if (!is.finite(pad) || pad == 0) {
    pad <- 0.05
  }

  c(value_range[1] - pad, value_range[2] + pad)
}

make_two_method_legend <- function(method_cols,
                                   method_linetypes,
                                   method_shapes,
                                   title = "Order selection",
                                   title_size = 7.5,
                                   text_size = 7,
                                   point_size = 4.5,
                                   line_width = 1.5,
                                   base_size = 25) {
  methods <- c("Functional IC", "CV")
  legend_data <- tibble::tibble(
    method = factor(methods, levels = methods),
    x_start = c(3.5, 7.3),
    x_end = c(4.0, 7.8),
    x_text = c(4.3, 8.1),
    y = 1
  )

  ggplot2::ggplot(legend_data) +
    ggplot2::annotate(
      "text", x = 0, y = 1, label = title, hjust = 0,
      size = title_size, fontface = "bold"
    ) +
    ggplot2::geom_segment(
      ggplot2::aes(
        x = x_start, xend = x_end, y = y, yend = y,
        colour = method, linetype = method
      ),
      linewidth = line_width,
      show.legend = FALSE
    ) +
    ggplot2::geom_point(
      ggplot2::aes(x = x_end, y = y, colour = method, shape = method),
      size = point_size,
      show.legend = FALSE
    ) +
    ggplot2::geom_text(
      ggplot2::aes(x = x_text, y = y, label = method),
      hjust = 0,
      size = text_size
    ) +
    ggplot2::scale_colour_manual(values = method_cols) +
    ggplot2::scale_linetype_manual(values = method_linetypes) +
    ggplot2::scale_shape_manual(values = method_shapes) +
    ggplot2::coord_cartesian(
      xlim = c(0, 11),
      ylim = c(0.85, 1.15),
      clip = "off"
    ) +
    ggplot2::theme_void(base_size = base_size) +
    ggplot2::theme(
      legend.position = "none",
      plot.margin = ggplot2::margin(0, 0, 0, 0)
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

  final_plot <- cowplot::plot_grid(
    legend,
    labelled_grid,
    ncol = 1,
    rel_heights = c(0.13, 1)
  )

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
