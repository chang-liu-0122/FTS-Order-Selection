invisible(Sys.setlocale("LC_CTYPE", "en_US.UTF-8"))
library(shiny)
library(bslib)
library(dplyr)
library(tidyr)
library(ggplot2)
library(DT)
library(shinyWidgets)

source("forecast_plot_helpers.R")
source("intro_pages.R")
app <- readRDS("app_data.rds")

mq_red <- "#A6192E"
mq_purple <- "#80225F"
mq_charcoal <- "#373A36"
mq_slate <- "#415364"
mq_sand <- "#EDEBE5"
psi <- "\u03c8"
nu <- "\u03bd"
dot <- "  \u00b7  "
bench <- "FARFIMA(1,1,0)"
bj <- "Functional Box\u2013Jenkins"
strength_col <- paste0("(", psi, ", ", nu, ")")

method_colours <- c("FARFIMA(1,1,0)" = mq_slate, "Functional IC" = mq_red, "CV" = mq_purple, "Plug-in" = mq_red)
method_colours[bj] <- mq_purple
method_shapes <- c(paper_method_shapes, "Plug-in" = 16)
method_shapes[bj] <- 17
method_linetypes <- c(paper_method_linetypes, "Plug-in" = "solid")
method_linetypes[bj] <- "dashed"
three_colours <- c(mq_slate, mq_red, mq_purple)
strength_colours <- c("0.2" = "#D38C97", "0.5" = mq_red, "0.8" = mq_purple)

theme_set(
  paper_plot_theme(base_size = 14) +
    theme(plot.margin = margin(10, 18, 10, 10), legend.title = element_text(face = "bold"), panel.spacing = unit(1, "lines"),
          text = element_text(colour = mq_charcoal), axis.text = element_text(colour = mq_charcoal),
          panel.border = element_rect(colour = mq_charcoal, fill = NA, linewidth = 0.5),
          strip.background = element_rect(fill = mq_sand, colour = mq_charcoal, linewidth = 0.5))
)

measure_groups <- list(
  "RMSPE and MAPE" = c("RMSPE", "MAPE"),
  "ECP, CPD, and IS" = c("ECP", "CPD", "IS"),
  "WIS" = "WIS"
)
ratio_measures <- c("RMSPE", "MAPE", "IS", "WIS")
order_label <- function(p, q) sprintf("(%d,%d)", p, q)
parse_order <- function(x) as.integer(strsplit(gsub("[()]", "", x), ",")[[1]])
avg_label <- "Average over strengths"

fmt <- function(v, digits = 4) ifelse(is.na(v), "\u2013", formatC(v, format = "f", digits = digits))

with_summary <- function(wide, cols, key = "h") {
  bind_rows(
    wide |> mutate(across(all_of(key), as.character)),
    wide |> summarise(across(all_of(cols), \(v) mean(v, na.rm = TRUE))) |> mutate(!!key := "Mean"),
    wide |> summarise(across(all_of(cols), \(v) median(v, na.rm = TRUE))) |> mutate(!!key := "Median")
  ) |> relocate(all_of(key))
}

mark_min <- function(df, groups, highlight = TRUE, digits = 4) {
  out <- df
  highlight <- rep_len(highlight, nrow(df))
  for (g in groups) {
    m <- round(as.matrix(df[g]), digits)
    mins <- suppressWarnings(apply(m, 1, min, na.rm = TRUE))
    for (cn in g) {
      s <- fmt(df[[cn]], digits)
      is_min <- highlight & !is.na(df[[cn]]) & round(df[[cn]], digits) == mins
      out[[cn]] <- ifelse(is_min, paste0("<span class='min-value'>", s, "</span>"), s)
    }
  }
  out
}

dt_table <- function(df, group = NULL, page = FALSE) {
  opts <- list(dom = if (page) "ftip" else "t", paging = page, pageLength = 20, ordering = FALSE, scrollX = TRUE,
               columnDefs = list(list(className = "dt-right", targets = "_all")))
  ext <- character(0)
  if (!is.null(group)) {
    gi <- match(group, names(df)) - 1
    opts$rowGroup <- list(dataSrc = gi)
    opts$columnDefs <- c(opts$columnDefs, list(list(visible = FALSE, targets = gi)))
    ext <- "RowGroup"
  }
  key <- setdiff(names(df), group)[1]
  datatable(df, escape = FALSE, rownames = FALSE, class = "compact hover", extensions = ext, options = opts) |>
    formatStyle(key, target = "row",
                backgroundColor = styleEqual(c("Mean", "Median"), c("#F6F5F1", "#F6F5F1")),
                fontWeight = styleEqual(c("Mean", "Median"), c("600", "600")))
}

section_label <- function(x) div(class = "control-section", x)

results_card <- function(id) {
  ns <- NS(id)
  card(
    full_screen = TRUE,
    card_header(
      class = "results-header",
      div(class = "results-caption", textOutput(ns("caption"), inline = TRUE)),
      radioGroupButtons(ns("view"), NULL, c("Plot", "Table"), size = "sm", status = "outline-primary")
    ),
    card_body(
      conditionalPanel("input.view == 'Plot'", ns = ns, plotOutput(ns("plot"), height = "auto")),
      conditionalPanel("input.view == 'Table'", ns = ns, DTOutput(ns("table")))
    )
  )
}

card_server <- function(id, plot_r, table_r, caption, height = reactive(460)) {
  moduleServer(id, function(input, output, session) {
    output$caption <- renderText(caption())
    output$plot <- renderPlot({
      p <- plot_r()
      if (inherits(p, "ggplot")) p else grid::grid.draw(p)
    }, height = function() height(), res = 110)
    output$table <- renderDT(table_r())
  })
}

method_scales <- function(title = "Forecasting method") {
  list(
    scale_colour_manual(name = title, values = method_colours),
    scale_shape_manual(name = title, values = method_shapes),
    scale_linetype_manual(name = title, values = method_linetypes)
  )
}

metrics_controls <- function(id) {
  ns <- NS(id)
  tagList(
    section_label("Display"),
    selectInput(ns("group"), "Measures", names(measure_groups)),
    conditionalPanel("input.group != 'RMSPE and MAPE'", ns = ns,
                     selectInput(ns("rule"), "Conformal method", c("SR", "NS", "AC", "All three"))),
    conditionalPanel("input.group == 'ECP, CPD, and IS'", ns = ns,
                     radioGroupButtons(ns("level"), "Nominal coverage", c("80%" = "0.2", "95%" = "0.05"), size = "sm", status = "outline-primary")),
    radioGroupButtons(ns("values"), "Values", c("Raw", "Ratio"), size = "sm", status = "outline-primary"),
    conditionalPanel("input.group == 'ECP, CPD, and IS'", ns = ns,
                     div(class = "control-help", "Ratios apply to IS; ECP and CPD are always shown as raw values."))
  )
}

metrics_server <- function(id, data, methods, setting_text) {
  moduleServer(id, function(input, output, session) {
    all_rules <- reactive(input$group != "RMSPE and MAPE" && input$rule == "All three")
    current <- reactive({
      d <- data()
      lvl <- as.numeric(input$level)
      rows <- lapply(measure_groups[[input$group]], function(m) {
        use_ratio <- input$values == "Ratio" && m %in% ratio_measures
        base <- if (use_ratio && any(d$metric == paste(m, "ratio"))) paste(m, "ratio") else m
        x <- d |> filter(metric == base)
        if (m %in% c("ECP", "CPD", "IS")) x <- x |> filter(abs(alpha - lvl) < 1e-9)
        if (m %in% c("ECP", "CPD", "IS", "WIS") && input$rule != "All three") x <- x |> filter(Rule == input$rule)
        if (use_ratio && base == m) x <- x |> group_by(Rule, h) |> mutate(value = value / value[Method == bench]) |> ungroup()
        if (use_ratio) x <- x |> filter(Method != bench)
        x |> mutate(Measure = if (use_ratio) paste(m, "ratio") else m)
      })
      x <- bind_rows(rows)
      x |> mutate(Measure = factor(Measure, levels = unique(x$Measure)), Method = factor(Method, levels = methods),
                  Rule = factor(Rule, levels = c("none", "SR", "NS", "AC"))) |>
        arrange(Measure, Rule, Method, h)
    })
    caption <- reactive({
      g <- input$group
      extra <- c(
        if (g != "RMSPE and MAPE") if (input$rule == "All three") "SR, NS, and AC" else input$rule,
        if (g == "ECP, CPD, and IS") if (input$level == "0.2") "80% bands" else "95% bands",
        if (input$values == "Ratio") paste(if (g == "ECP, CPD, and IS") "IS as ratio to" else "ratios to", bench)
      )
      paste(c(g, extra, setting_text()), collapse = dot)
    })
    height <- reactive(if (all_rules() && input$group == "ECP, CPD, and IS") 760 else 440)
    plot_r <- reactive({
      x <- current()
      lvl <- as.numeric(input$level)
      ms <- levels(x$Measure)
      refs <- bind_rows(
        tibble(Measure = ms[grepl("ratio", ms)], y = 1),
        tibble(Measure = ms[ms == "ECP"], y = 1 - lvl)
      ) |> mutate(Measure = factor(Measure, levels = ms))
      p <- ggplot(x, aes(h, value, colour = Method, shape = Method, linetype = Method, group = Method))
      if (nrow(refs) > 0) p <- p + geom_hline(data = refs, aes(yintercept = y), inherit.aes = FALSE, linetype = "dashed", colour = "#555555")
      p <- p + geom_line(linewidth = 0.9) + geom_point(size = 2.4) + method_scales() +
        scale_x_continuous(breaks = 1:10) +
        labs(x = "Forecast horizon h", y = NULL)
      if (all_rules()) p + facet_grid(Measure ~ Rule, scales = "free_y") else p + facet_wrap(~ Measure, nrow = 1, scales = "free_y")
    })
    table_r <- reactive({
      x <- current() |> mutate(col = if (all_rules()) ifelse(Rule == "none", as.character(Method), paste0(Rule, ": ", Method)) else as.character(Method))
      wide <- x |> select(Measure, h, col, value) |> pivot_wider(names_from = col, values_from = value)
      cols <- setdiff(names(wide), c("Measure", "h"))
      raw <- wide |> group_by(Measure) |> group_modify(\(df, key) with_summary(df, cols)) |> ungroup() |>
        mutate(Measure = as.character(Measure))
      groups <- if (all_rules()) split(cols, sub(":.*", "", cols)) else list(cols)
      shown <- mark_min(raw, groups, highlight = raw$Measure != "ECP")
      dt_table(shown, group = "Measure")
    })
    card_server("card", plot_r, table_r, caption, height)
  })
}

order_choices <- function(df) {
  df |> distinct(p_true, q_true) |> arrange(p_true + q_true, desc(p_true)) |> mutate(lab = order_label(p_true, q_true)) |> pull(lab)
}

strength_choices <- function(df, order) {
  pq <- parse_order(order)
  df |> filter(p_true == pq[1], q_true == pq[2]) |> distinct(strength) |> pull(strength) |> as.character() |> sort()
}

sim_forecast_ui <- function(id, study_name) {
  ns <- NS(id)
  d <- app$sim_fc |> filter(study == study_name)
  layout_sidebar(
    sidebar = sidebar(
      width = 300,
      section_label("Setting"),
      selectInput(ns("T"), "Sample size T", sort(unique(d$T)), selected = max(d$T)),
      radioGroupButtons(ns("window"), "Window", c("expanding", "fixed"), size = "sm", status = "outline-primary"),
      if (study_name == "Unified") selectInput(ns("d"), "Memory parameter d", sort(unique(d$d_true)), selected = 0.2),
      selectInput(ns("order"), "True order (p, q)", order_choices(d)),
      selectInput(ns("strength"), "Dependence strength", avg_label),
      metrics_controls(ns("m"))
    ),
    results_card(NS(ns("m"), "card"))
  )
}

sim_forecast_server <- function(id, study_name, methods) {
  moduleServer(id, function(input, output, session) {
    sim <- app$sim_fc |> filter(study == study_name)
    observeEvent(input$order, updateSelectInput(session, "strength", choices = c(avg_label, strength_choices(sim, input$order))))
    data_r <- reactive({
      req(input$order, input$strength %in% c(avg_label, strength_choices(sim, input$order)))
      pq <- parse_order(input$order)
      x <- sim |> filter(T == as.numeric(input$T), window == input$window, p_true == pq[1], q_true == pq[2])
      if (study_name == "Unified") x <- x |> filter(abs(d_true - as.numeric(input$d)) < 1e-9)
      if (input$strength != avg_label) x <- x |> filter(strength == input$strength)
      x <- x |> mutate(across(c(Method, Rule, metric, strength), as.character))
      r <- x |> filter(metric %in% c("IS", "WIS")) |>
        group_by(strength, Rule, alpha, h, metric) |>
        mutate(value = value / value[Method == bench]) |> ungroup() |>
        mutate(metric = paste(metric, "ratio"))
      bind_rows(x, r) |> group_by(Method, Rule, alpha, h, metric) |> summarise(value = mean(value), .groups = "drop")
    })
    setting_text <- reactive(paste(c(
      paste0("T = ", input$T), paste(input$window, "window"),
      if (study_name == "Unified") paste0("d = ", input$d),
      paste0("(p, q) = ", input$order), if (input$strength == avg_label) paste0("averaged over ", psi, " and ", nu) else input$strength
    ), collapse = dot))
    metrics_server("m", data_r, methods, setting_text)
  })
}

heatmap_plot <- function(s, pq, method) {
  ggplot(s |> filter(Method == method), aes(q_hat, p_hat, fill = Proportion)) +
    geom_tile(colour = "white", linewidth = 0.8) +
    geom_tile(data = s |> filter(Method == method, p_hat == pq[1], q_hat == pq[2]), colour = mq_charcoal, linewidth = 1.2, fill = NA) +
    geom_text(aes(label = sprintf("%.2f", Proportion), colour = Proportion > 0.5), size = 3.8) +
    scale_fill_gradient(low = "#F6F5F1", high = mq_red, limits = c(0, 1), guide = "none") +
    scale_colour_manual(values = c("FALSE" = mq_charcoal, "TRUE" = "white"), guide = "none") +
    scale_x_discrete(drop = FALSE) + scale_y_discrete(drop = FALSE) +
    labs(title = method, x = expression("Selected " * hat(q)), y = expression("Selected " * hat(p))) +
    theme(plot.title = element_text(face = "bold", size = 14, hjust = 0, colour = mq_charcoal))
}

stacked_heatmaps <- function(sel, pq, methods) {
  s <- sel |> filter(!is.na(p_hat), !is.na(q_hat)) |> count(Method, p_hat, q_hat) |>
    group_by(Method) |> mutate(Proportion = n / sum(n)) |> ungroup() |>
    mutate(p_hat = factor(p_hat, levels = 0:5), q_hat = factor(q_hat, levels = 0:5))
  gridExtra::arrangeGrob(grobs = lapply(methods, function(m) ggplotGrob(heatmap_plot(s, pq, m))), ncol = 1)
}

selected_table <- function(sel, first) {
  raw <- sel |> mutate(`Selected (p, q)` = ifelse(is.na(p_hat) | is.na(q_hat), "failed", order_label(p_hat, q_hat))) |>
    count(Method, `Selected (p, q)`) |> group_by(Method) |> mutate(n = n / sum(n)) |> ungroup() |>
    pivot_wider(names_from = Method, values_from = n, values_fill = 0) |>
    arrange(desc(.data[[first]]))
  dt_table(mark_min(raw, list(setdiff(names(raw), "Selected (p, q)")), highlight = FALSE))
}

accuracy_plot <- function(s, x, xlab, breaks) {
  ggplot(s, aes(.data[[x]], Accuracy, colour = Method, shape = Method, linetype = Method)) +
    geom_line(linewidth = 0.9) + geom_point(size = 2.6) + method_scales("Method") +
    scale_x_continuous(breaks = breaks) + coord_cartesian(ylim = c(0, 1)) +
    labs(x = xlab, y = "Exact-order selection accuracy")
}

memory_ui <- layout_sidebar(
  sidebar = sidebar(
    width = 300,
    section_label("Setting"),
    selectInput("mem_model", "Model", c("FARFIMA(1,d,0)", "FARFIMA(1,d,1)")),
    radioGroupButtons("mem_estimator", "Estimator", c("LW", "DFA"), size = "sm", status = "outline-primary"),
    selectInput("mem_T", "Sample size T", c(500, 1000, 2000), selected = 1000),
    conditionalPanel("input.mem_model == 'FARFIMA(1,d,1)'",
                     selectInput("mem_nu", paste("MA strength", nu), setNames(c(0.2, 0.5, 0.8), paste(nu, "=", c(0.2, 0.5, 0.8))), selected = 0.5))
  ),
  results_card("mem_card")
)

order_ui <- layout_sidebar(
  sidebar = sidebar(
    width = 300,
    section_label("Setting"),
    selectInput("os_order", "True order (p, q)", order_choices(app$order_sel)),
    selectInput("os_strength", "Dependence strength", NULL),
    conditionalPanel("input.os_view == 'Heat map'", selectInput("os_T", "Sample size T", c(1000, 2000, 5000))),
    section_label("Display"),
    radioGroupButtons("os_view", NULL, c("Heat map", "Accuracy across T"), size = "sm", status = "outline-primary")
  ),
  results_card("os_card")
)

unified_d_ui <- layout_sidebar(
  sidebar = sidebar(
    width = 300,
    section_label("Setting"),
    selectInput("ud_T", "Sample size T", c(500, 1000), selected = 1000),
    div(class = "control-help", paste0("For FARFIMA(1,d,1), the dependence strengths are balanced, ", psi, " = ", nu, "."))
  ),
  results_card("ud_card")
)

unified_order_ui <- layout_sidebar(
  sidebar = sidebar(
    width = 300,
    section_label("Setting"),
    selectInput("uo_T", "Sample size T", c(500, 1000), selected = 1000),
    selectInput("uo_order", "True order (p, q)", order_choices(app$unified_order)),
    selectInput("uo_strength", "Dependence strength", NULL),
    conditionalPanel("input.uo_view == 'Heat map'", selectInput("uo_d", "Memory parameter d", sort(unique(app$unified_order$d_true)), selected = 0.2)),
    section_label("Display"),
    radioGroupButtons("uo_view", NULL, c("Accuracy across d", "Heat map"), size = "sm", status = "outline-primary")
  ),
  results_card("uo_card")
)

series_names <- unique(app$emp_fc$series)
emp_methods <- c(bench, "Functional IC", "CV")

emp_data_ui <- layout_sidebar(
  sidebar = sidebar(
    width = 300,
    section_label("Setting"),
    selectInput("ed_series", "Series", series_names),
    section_label("Display"),
    radioGroupButtons("ed_view", NULL, c("Curves", "Functional ACF"), size = "sm", status = "outline-primary")
  ),
  results_card("ed_card")
)

emp_models_ui <- div(class = "intro", card(
  card_header(class = "results-header", div(class = "results-caption", "Memory parameters, short-memory orders, and AC learning rates fixed before test-sample evaluation")),
  card_body(DTOutput("emp_models"))
))

emp_accuracy_ui <- layout_sidebar(
  sidebar = sidebar(width = 300, section_label("Setting"), selectInput("ea_series", "Series", series_names), metrics_controls("ea_m")),
  results_card(NS("ea_m", "card"))
)

emp_visual_ui <- layout_sidebar(
  sidebar = sidebar(
    width = 300,
    section_label("Setting"),
    selectInput("eb_series", "Series", series_names),
    selectInput("eb_origin", "Forecast origin (last observed year)", NULL),
    checkboxGroupInput("eb_methods", "Forecasting methods", emp_methods, selected = emp_methods),
    section_label("Band"),
    radioGroupButtons("eb_rule", "Conformal method", c("SR", "NS", "AC"), size = "sm", status = "outline-primary"),
    radioGroupButtons("eb_level", "Nominal coverage", c("80%" = "0.2", "95%" = "0.05"), size = "sm", status = "outline-primary")
  ),
  results_card("eb_card")
)

app_css <- r"(
:root {
  --mq-red: #A6192E;
  --mq-deep-red: #76232F;
  --mq-purple: #80225F;
  --mq-charcoal: #373A36;
  --mq-charcoal-dark: #262826;
  --mq-sand-light: #EDEBE5;
  --mq-sand: #D6D2C4;
  --soft: #F6F5F1;
  --muted: #6B6E68;
  --serif: 'Source Serif 4', 'Source Serif Pro', Georgia, serif;
}

body {
  background: #fff;
  color: var(--mq-charcoal);
}

h1, h2, h3, h4, h5, h6 {
  color: var(--mq-charcoal);
  font-weight: 600;
}

.navbar {
  background: #fff !important;
  border-top: 8px solid var(--mq-charcoal-dark);
  border-bottom: 1px solid var(--mq-sand);
  padding-top: 0;
  padding-bottom: 0;
}

.navbar-brand {
  display: flex;
  align-items: center;
  gap: 0.65rem;
  font-weight: 700;
  font-size: 1.15rem;
  color: var(--mq-charcoal) !important;
}

.navbar-brand::before {
  content: "";
  flex-shrink: 0;
  width: 1.05rem;
  height: 1.05rem;
  background: linear-gradient(-45deg, transparent 50%, var(--mq-red) 50%);
}

.navbar .navbar-nav .nav-link {
  color: var(--mq-charcoal) !important;
  font-size: 0.75rem;
  font-weight: 500;
  letter-spacing: 0.08em;
  text-transform: uppercase;
  white-space: nowrap;
  padding: 1.25rem 0 1rem !important;
  margin: 0 0.65rem;
  border-bottom: 3px solid transparent;
}

.navbar .navbar-nav .nav-link:hover {
  color: var(--mq-red) !important;
}

.navbar .navbar-nav .nav-link.active {
  color: var(--mq-red) !important;
  font-weight: 600;
  border-bottom-color: var(--mq-red);
}

.btn-primary {
  --bs-btn-hover-bg: var(--mq-deep-red);
  --bs-btn-hover-border-color: var(--mq-deep-red);
  --bs-btn-active-bg: var(--mq-deep-red);
  --bs-btn-active-border-color: var(--mq-deep-red);
}

.btn-outline-primary {
  --bs-btn-hover-bg: var(--mq-red);
  --bs-btn-hover-border-color: var(--mq-red);
  --bs-btn-active-bg: var(--mq-red);
  --bs-btn-active-border-color: var(--mq-red);
}

.nav-underline:not(.navbar-nav) {
  gap: 1.8rem;
  margin-bottom: 1rem;
  border-bottom: 1px solid var(--mq-sand);
}

.nav-underline:not(.navbar-nav) .nav-link {
  color: var(--mq-charcoal);
  font-weight: 500;
}

.nav-underline:not(.navbar-nav) .nav-link:hover {
  color: var(--mq-red);
  border-bottom-color: transparent;
}

.nav-underline:not(.navbar-nav) .nav-link.active {
  color: var(--mq-red);
  font-weight: 600;
  border-bottom-color: var(--mq-red);
}

.bslib-sidebar-layout {
  border: 1px solid var(--mq-sand);
  border-radius: 0;
  background: #fff;
}

.bslib-sidebar-layout > .sidebar {
  background: var(--mq-sand-light);
  border-right: 1px solid var(--mq-sand);
}

.control-section {
  font-size: 0.72rem;
  font-weight: 600;
  letter-spacing: 0.12em;
  text-transform: uppercase;
  color: var(--mq-red);
  margin: 0.4rem 0 0.6rem;
  padding-bottom: 0.35rem;
  border-bottom: 1px solid var(--mq-sand);
}

.sidebar .form-label,
.sidebar .control-label {
  font-weight: 600;
  font-size: 0.88rem;
}

.control-help {
  font-size: 0.8rem;
  color: var(--muted);
  line-height: 1.4;
  margin: -0.2rem 0 0.8rem;
}

.card {
  --bs-card-border-radius: 0;
  --bs-card-inner-border-radius: 0;
  border: 1px solid var(--mq-sand);
  border-top: 3px solid var(--mq-red);
  border-radius: 0;
  box-shadow: none;
  background: #fff;
}

.results-header {
  display: flex;
  justify-content: space-between;
  align-items: center;
  gap: 1rem;
  background: #fff;
  border-bottom: 1px solid var(--mq-sand);
}

.results-caption {
  font-family: var(--serif);
  font-size: 1.04rem;
  line-height: 1.4;
  color: var(--mq-charcoal);
}

.results-header .btn-group .btn {
  font-size: 0.8rem;
  padding: 0.2rem 0.9rem;
}

table.dataTable {
  font-size: 0.9rem;
}

table.dataTable td,
table.dataTable th {
  font-variant-numeric: tabular-nums;
  white-space: nowrap;
}

table.dataTable thead th {
  border-bottom: 2px solid var(--mq-charcoal) !important;
  font-weight: 600;
}

table.dataTable tr.dtrg-group td,
table.dataTable tr.dtrg-group th {
  background: var(--mq-sand-light) !important;
  color: var(--mq-red);
  font-weight: 600;
  text-align: left !important;
  letter-spacing: 0.02em;
}

.min-value {
  color: var(--mq-red);
  font-weight: 700;
}

.intro {
  max-width: 980px;
  margin: 0 auto;
}

.intro-wrap {
  max-width: 1240px;
  margin: 0 auto;
}

.hero {
  padding: 2.8rem 2.6rem 2.6rem;
  margin-bottom: 1.8rem;
  color: #fff;
  background: linear-gradient(135deg, var(--mq-red) 62%, var(--mq-deep-red) 62%);
}

.hero h1 {
  color: #fff;
  font-weight: 700;
  font-size: 2.25rem;
  line-height: 1.15;
  max-width: 24ch;
  margin: 0;
}

.hero .lead {
  font-family: var(--serif);
  font-size: 1.12rem;
  font-weight: 500;
  line-height: 1.55;
  max-width: 60ch;
  margin: 1.1rem 0 0;
}

.intro-wrap .nav-pills {
  position: sticky;
  top: 1rem;
  background: #fff;
  border: 1px solid var(--mq-sand);
  border-top: 3px solid var(--mq-red);
  padding: 0;
}

.intro-wrap .nav-pills .nav-link {
  color: var(--mq-charcoal);
  border-radius: 0;
  border-left: 4px solid transparent;
  border-bottom: 1px solid var(--mq-sand-light);
  padding: 0.7rem 1rem;
  font-size: 0.92rem;
}

.intro-wrap .nav-pills .nav-link:hover {
  background: var(--soft);
  color: var(--mq-red);
}

.intro-wrap .nav-pills .nav-link.active {
  background: var(--mq-sand-light);
  color: var(--mq-red);
  font-weight: 600;
  border-left-color: var(--mq-red);
}

.intro-card .card-body {
  padding: 1.8rem 2.2rem;
}

.intro-section .kicker {
  font-size: 0.74rem;
  font-weight: 600;
  letter-spacing: 0.12em;
  text-transform: uppercase;
  color: var(--mq-red);
}

.intro-section h2 {
  font-weight: 600;
  font-size: 1.9rem;
  margin: 0.25rem 0 1rem;
}

.intro-section h2 + p {
  font-family: var(--serif);
  font-size: 1.1rem;
  font-weight: 500;
  line-height: 1.6;
}

.intro-section h4 {
  font-weight: 600;
  margin-top: 1.4rem;
}

.intro-section p,
.intro-section li {
  line-height: 1.7;
}

.intro-section .katex-display {
  overflow-x: auto;
  overflow-y: hidden;
  padding: 0.3rem 0;
}

.intro-actions {
  margin-top: 1.6rem;
}

.intro-actions .btn {
  font-weight: 600;
  padding: 0.55rem 1.4rem;
  border-radius: 8px;
}

.settings-table {
  width: 100%;
  border-collapse: collapse;
  margin: 0.4rem 0 0.8rem;
}

.settings-table th,
.settings-table td {
  padding: 0.5rem 0.75rem;
  border-bottom: 1px solid var(--mq-sand);
  vertical-align: top;
}

.settings-table th {
  width: 32%;
  font-weight: 600;
  background: var(--mq-sand-light);
}

@media (max-width: 1300px) {
  .navbar .navbar-nav .nav-link {
    margin: 0 0.4rem;
    letter-spacing: 0.04em;
  }
}

@media (max-width: 1120px) {
  .navbar .navbar-nav .nav-link {
    font-size: 0.7rem;
    margin: 0 0.3rem;
    letter-spacing: 0.02em;
  }
}

@media (max-width: 900px) {
  .hero {
    padding: 2rem 1.6rem;
  }

  .hero h1 {
    font-size: 1.65rem;
  }
}

@media (max-width: 480px) {
  .navbar-brand {
    font-size: 0.95rem;
    white-space: normal;
    max-width: calc(100% - 4rem);
  }

  .hero h1 {
    font-size: 1.35rem;
  }
}
)"

katex_head <- tagList(
  tags$link(rel = "stylesheet", href = "https://cdn.jsdelivr.net/npm/katex@0.16.11/dist/katex.min.css"),
  tags$script(defer = NA, src = "https://cdn.jsdelivr.net/npm/katex@0.16.11/dist/katex.min.js"),
  tags$script(defer = NA, src = "https://cdn.jsdelivr.net/npm/katex@0.16.11/dist/contrib/auto-render.min.js"),
  tags$script(HTML("document.addEventListener('DOMContentLoaded', function () { renderMathInElement(document.body, {delimiters: [{left: '\\\\[', right: '\\\\]', display: true}, {left: '\\\\(', right: '\\\\)', display: false}], throwOnError: false}); });"))
)

intro_panel <- function(title, content) nav_panel(title, card(class = "intro-card", card_body(content)))

intro_ui <- div(
  class = "intro-wrap",
  div(
    class = "hero",
    h1("Order Selection and Forecasting in Long Memory Functional Time Series"),
    p(class = "lead", "Interactive results for the simulation study and the empirical applications. Choose a section below to read its design, then open its results.")
  ),
  navset_pill_list(
    id = "intro_nav", widths = c(3, 9), well = FALSE,
    intro_panel("Overview", intro_overview),
    intro_panel("Simulation design", intro_dgp),
    intro_panel("Memory estimation", intro_memory),
    intro_panel("Order selection", intro_order),
    intro_panel("Ideal forecasting", intro_ideal),
    intro_panel("Unified framework", intro_unified),
    intro_panel("Empirical study", intro_empirical),
    intro_panel("Forecast evaluation", intro_evaluation)
  )
)

ui <- page_navbar(
  id = "main",
  title = "Long Memory FTS \u00b7 Results",
  theme = bs_theme(
    version = 5, primary = mq_red, bg = "#FFFFFF", fg = mq_charcoal,
    base_font = font_google("Work Sans", local = FALSE, wght = c(400, 500, 600, 700)),
    font_scale = 0.95
  ),
  navbar_options = navbar_options(bg = "#FFFFFF", theme = "light"),
  header = tagList(
    katex_head,
    tags$link(rel = "stylesheet", href = "https://fonts.googleapis.com/css2?family=Source+Serif+4:opsz,wght@8..60,400;8..60,500;8..60,600&display=swap"),
    tags$style(HTML(app_css))
  ),
  nav_panel("Introduction", value = "intro", intro_ui),
  nav_panel("Memory estimation", value = "memory", memory_ui),
  nav_panel("Order selection", value = "order", order_ui),
  nav_panel("Ideal forecasting", value = "ideal", sim_forecast_ui("ideal", "Ideal")),
  nav_panel(
    "Unified framework", value = "unified",
    navset_underline(
      nav_panel("Forecasting", sim_forecast_ui("unified", "Unified")),
      nav_panel("Memory estimate", unified_d_ui),
      nav_panel("Order selection", unified_order_ui)
    )
  ),
  nav_panel(
    "Empirical study", value = "empirical",
    navset_underline(
      nav_panel("Data", emp_data_ui),
      nav_panel("Selected models", emp_models_ui),
      nav_panel("Forecast accuracy", emp_accuracy_ui),
      nav_panel("Forecast visualisation", emp_visual_ui)
    )
  )
)

server <- function(input, output, session) {
  for (target in c("memory", "order", "ideal", "unified", "empirical")) {
    local({
      tg <- target
      observeEvent(input[[paste0("go_", tg)]], nav_select("main", tg))
    })
  }

  sim_forecast_server("ideal", "Ideal", c(bench, "Plug-in"))
  sim_forecast_server("unified", "Unified", c(bench, "Functional IC", "CV"))

  mem_summary <- reactive({
    x <- app$memory |> filter(model == input$mem_model, estimator == input$mem_estimator, T == as.numeric(input$mem_T))
    if (input$mem_model == "FARFIMA(1,d,1)") x <- x |> filter(abs(ma_strength - as.numeric(input$mem_nu)) < 1e-9)
    x |> group_by(ar_strength, d_true, operator) |>
      summarise(Bias = mean(d_hat - d_true), MSE = mean((d_hat - d_true)^2), .groups = "drop")
  })
  card_server(
    "mem_card",
    reactive({
      s <- mem_summary() |> pivot_longer(c(Bias, MSE), names_to = "Measure", values_to = "value") |>
        mutate(Measure = factor(recode(Measure, Bias = "Mean bias"), levels = c("Mean bias", "MSE")),
               panel = factor(paste(psi, "=", ar_strength)))
      ggplot(s, aes(d_true, value, colour = operator, shape = operator, linetype = operator)) +
        geom_hline(data = tibble(Measure = factor("Mean bias", levels = c("Mean bias", "MSE")), y = 0), aes(yintercept = y),
                   inherit.aes = FALSE, linetype = "dashed", colour = "#555555") +
        geom_line(linewidth = 0.8) + geom_point(size = 2.4) +
        scale_colour_manual(name = "Covariance operator", values = setNames(three_colours, c("Covar", "RS", "LRS"))) +
        scale_shape_manual(name = "Covariance operator", values = c(Covar = 16, RS = 17, LRS = 15)) +
        scale_linetype_manual(name = "Covariance operator", values = c(Covar = "solid", RS = "dashed", LRS = "dotted")) +
        scale_x_continuous(breaks = sort(unique(s$d_true))) +
        facet_grid(Measure ~ panel, scales = "free_y") +
        labs(x = "True memory parameter d", y = NULL)
    }),
    reactive({
      raw <- mem_summary() |> pivot_longer(c(Bias, MSE), names_to = "Measure", values_to = "value") |>
        pivot_wider(names_from = operator, values_from = value) |>
        arrange(Measure, ar_strength, d_true) |>
        transmute(Measure, Setting = paste0(psi, " = ", ar_strength, ", d = ", d_true), Covar, RS, LRS)
      raw <- raw |> group_by(Measure) |> group_modify(\(df, key) with_summary(df, c("Covar", "RS", "LRS"), key = "Setting")) |> ungroup()
      dt_table(mark_min(raw, list(c("Covar", "RS", "LRS")), highlight = FALSE), group = "Measure")
    }),
    reactive(paste(c(paste("Mean bias and MSE of the", input$mem_estimator, "estimates"), input$mem_model, paste0("T = ", input$mem_T),
                     if (input$mem_model == "FARFIMA(1,d,1)") paste(nu, "=", input$mem_nu)), collapse = dot)),
    reactive(560)
  )

  observeEvent(input$os_order, updateSelectInput(session, "os_strength", choices = strength_choices(app$order_sel, input$os_order)))
  os_data <- reactive({
    req(input$os_strength %in% strength_choices(app$order_sel, input$os_order))
    pq <- parse_order(input$os_order)
    app$order_sel |> filter(p_true == pq[1], q_true == pq[2], strength == input$os_strength)
  })
  card_server(
    "os_card",
    reactive({
      if (input$os_view == "Accuracy across T") {
        s <- os_data() |> group_by(T, Method) |> summarise(Accuracy = mean(correct %in% TRUE), .groups = "drop")
        return(accuracy_plot(s, "T", "Sample size T", c(1000, 2000, 5000)))
      }
      stacked_heatmaps(os_data() |> filter(T == as.numeric(input$os_T)), parse_order(input$os_order), c("Functional IC", bj))
    }),
    reactive({
      if (input$os_view == "Accuracy across T") {
        raw <- os_data() |> group_by(T, Method) |> summarise(Accuracy = mean(correct %in% TRUE), .groups = "drop") |>
          pivot_wider(names_from = Method, values_from = Accuracy) |> mutate(T = as.character(T))
        return(dt_table(mark_min(raw, list(setdiff(names(raw), "T")), highlight = FALSE)))
      }
      selected_table(os_data() |> filter(T == as.numeric(input$os_T)), "Functional IC")
    }),
    reactive(paste(c(
      if (input$os_view == "Heat map") "Proportion of replications selecting each order, true order outlined" else "Exact-order selection accuracy across sample sizes",
      paste0("(p, q) = ", input$os_order), input$os_strength, if (input$os_view == "Heat map") paste0("T = ", input$os_T)
    ), collapse = dot)),
    reactive(if (input$os_view == "Heat map") 780 else 460)
  )

  ud_summary <- reactive({
    app$unified_d |> filter(T == as.numeric(input$ud_T), !(p_true == 1 & q_true == 1 & abs(ar_strength - ma_strength) > 1e-9)) |>
      mutate(Model = factor(sprintf("FARFIMA(%d,d,%d)", p_true, q_true), levels = c("FARFIMA(1,d,0)", "FARFIMA(2,d,0)", "FARFIMA(0,d,1)", "FARFIMA(1,d,1)")),
             Strength = ifelse(p_true > 0, ar_strength, ma_strength)) |>
      group_by(Model, Strength, ar_strength, ma_strength, d_true) |>
      summarise(Bias = mean(d_hat - d_true), MSE = mean((d_hat - d_true)^2), .groups = "drop")
  })
  card_server(
    "ud_card",
    reactive({
      s <- ud_summary() |> pivot_longer(c(Bias, MSE), names_to = "Measure", values_to = "value") |>
        mutate(Measure = factor(recode(Measure, Bias = "Mean bias"), levels = c("Mean bias", "MSE")), Strength = factor(Strength))
      ggplot(s, aes(d_true, value, colour = Strength, group = Strength)) +
        geom_hline(data = tibble(Measure = factor("Mean bias", levels = c("Mean bias", "MSE")), y = 0), aes(yintercept = y),
                   inherit.aes = FALSE, linetype = "dashed", colour = "#555555") +
        geom_line(linewidth = 0.8) + geom_point(size = 2.2) +
        scale_colour_manual(name = "Dependence strength", values = strength_colours) +
        scale_x_continuous(breaks = sort(unique(s$d_true))) +
        facet_grid(Measure ~ Model, scales = "free_y") +
        labs(x = "True memory parameter d", y = NULL) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
    }),
    reactive({
      raw <- ud_summary() |> pivot_longer(c(Bias, MSE), names_to = "Measure", values_to = "value") |>
        arrange(Measure, Model, ar_strength, ma_strength) |>
        mutate(Model = sub("FARFIMA", "", as.character(Model)), Strengths = sprintf("(%g, %g)", ar_strength, ma_strength), d = as.character(d_true)) |>
        select(Measure, Model, Strengths, d, value) |>
        pivot_wider(names_from = d, values_from = value)
      dcols <- setdiff(names(raw), c("Measure", "Model", "Strengths"))
      raw <- raw |> group_by(Measure) |>
        group_modify(\(df, key) bind_rows(df, df |> summarise(across(all_of(dcols), mean)) |> mutate(Model = "Mean", Strengths = ""),
                                          df |> summarise(across(all_of(dcols), median)) |> mutate(Model = "Median", Strengths = ""))) |>
        ungroup()
      shown <- mark_min(raw, list(dcols), highlight = FALSE, digits = 3)
      names(shown)[names(shown) == "Strengths"] <- strength_col
      dt_table(shown, group = "Measure")
    }),
    reactive(paste(c("Mean bias and MSE of the BC-LPWN estimates", paste0("T = ", input$ud_T),
                     paste0("dependence strengths 0.2, 0.5, 0.8, with ", psi, " = ", nu, " for FARFIMA(1,d,1)")), collapse = dot)),
    reactive(560)
  )

  observeEvent(input$uo_order, updateSelectInput(session, "uo_strength", choices = c(avg_label, strength_choices(app$unified_order, input$uo_order))))
  uo_data <- reactive({
    req(input$uo_strength %in% c(avg_label, strength_choices(app$unified_order, input$uo_order)))
    pq <- parse_order(input$uo_order)
    x <- app$unified_order |> filter(T == as.numeric(input$uo_T), p_true == pq[1], q_true == pq[2])
    if (input$uo_strength != avg_label) x <- x |> filter(strength == input$uo_strength)
    x
  })
  card_server(
    "uo_card",
    reactive({
      if (input$uo_view == "Accuracy across d") {
        s <- uo_data() |> group_by(d_true, Method) |> summarise(Accuracy = mean(correct %in% TRUE), .groups = "drop")
        return(accuracy_plot(s, "d_true", "True memory parameter d", sort(unique(s$d_true))))
      }
      stacked_heatmaps(uo_data() |> filter(abs(d_true - as.numeric(input$uo_d)) < 1e-9), parse_order(input$uo_order), c("Functional IC", "CV"))
    }),
    reactive({
      if (input$uo_view == "Accuracy across d") {
        raw <- uo_data() |> group_by(d = d_true, Method) |> summarise(Accuracy = mean(correct %in% TRUE), .groups = "drop") |>
          pivot_wider(names_from = Method, values_from = Accuracy) |> mutate(d = as.character(d))
        return(dt_table(mark_min(raw, list(setdiff(names(raw), "d")), highlight = FALSE)))
      }
      selected_table(uo_data() |> filter(abs(d_true - as.numeric(input$uo_d)) < 1e-9), "Functional IC")
    }),
    reactive(paste(c(
      if (input$uo_view == "Heat map") "Proportion of replications selecting each order, true order outlined" else "Exact-order selection accuracy across d",
      paste0("T = ", input$uo_T), paste0("(p, q) = ", input$uo_order),
      if (input$uo_strength == avg_label) paste0("averaged over ", psi, " and ", nu) else input$uo_strength,
      if (input$uo_view == "Heat map") paste0("d = ", input$uo_d)
    ), collapse = dot)),
    reactive(if (input$uo_view == "Heat map") 780 else 460)
  )

  scale_of <- function(s) app$emp_scale$scale[app$emp_scale$series == s]

  card_server(
    "ed_card",
    reactive({
      if (input$ed_view == "Functional ACF") {
        f <- app$emp_facf |> filter(series == input$ed_series)
        return(ggplot(f, aes(lag, rho)) +
                 geom_col(fill = "#D6D2C4", colour = "#8C887C", width = 0.6) +
                 geom_hline(aes(yintercept = bound), colour = mq_red, linetype = "dashed", linewidth = 1) +
                 labs(x = "Lag", y = "fACF"))
      }
      ggplot(app$emp_curves |> filter(series == input$ed_series), aes(age, value, colour = year, group = year)) +
        geom_line(linewidth = 0.4) +
        scale_colour_gradientn(name = "Year", colours = rainbow(7, end = 0.83)) +
        guides(colour = guide_colourbar(theme = theme(legend.key.width = unit(16, "lines")))) +
        facet_wrap(~ kind) +
        labs(x = "Age", y = scale_of(input$ed_series))
    }),
    reactive({
      if (input$ed_view == "Functional ACF") {
        raw <- app$emp_facf |> filter(series == input$ed_series) |> transmute(Lag = as.character(lag), fACF = rho, `95% white-noise bound` = bound)
        return(dt_table(mark_min(raw, list(c("fACF", "95% white-noise bound")), highlight = FALSE)))
      }
      raw <- app$emp_curves |> filter(series == input$ed_series) |>
        pivot_wider(names_from = kind, values_from = value) |>
        arrange(year, age) |> transmute(Year = as.character(year), Age = as.character(age), Raw, Smoothed)
      dt_table(mark_min(raw, list(c("Raw", "Smoothed")), highlight = FALSE), page = TRUE)
    }),
    reactive(paste(c(
      if (input$ed_view == "Curves") "Raw and smoothed curves, coloured from the earliest year (red) to the latest (violet)" else "Functional autocorrelation function with the 95% white-noise bound",
      input$ed_series
    ), collapse = dot))
  )

  output$emp_models <- renderDT({
    tab <- app$emp_models |> transmute(Series = series, Method, d = fmt(d), `(p, q)` = order_label(p, q), eta80 = eta_0.2, eta95 = eta_0.05)
    names(tab)[5:6] <- paste0("AC \u03b7 (", c("80%", "95%"), ")")
    dt_table(tab)
  })

  metrics_server(
    "ea_m",
    reactive(app$emp_fc |> filter(series == input$ea_series) |> select(Method, Rule, alpha, h, metric, value)),
    emp_methods, reactive(input$ea_series)
  )

  observeEvent(input$eb_series, {
    o <- app$emp_forecasts |> filter(series == input$eb_series) |> group_by(t0) |>
      summarise(origin = min(year) - 1L, H = max(h), .groups = "drop") |> arrange(desc(origin))
    updateSelectInput(session, "eb_origin", choices = setNames(o$t0, o$origin), selected = o$t0[o$H == 10][1])
  })
  eb_data <- reactive({
    req(input$eb_origin, length(input$eb_methods) > 0)
    req(as.integer(input$eb_origin) %in% app$emp_forecasts$t0[app$emp_forecasts$series == input$eb_series])
    lvl <- as.numeric(input$eb_level)
    app$emp_forecasts |> filter(series == input$eb_series, t0 == as.integer(input$eb_origin), Method %in% input$eb_methods) |>
      inner_join(app$emp_gamma |> filter(series == input$eb_series, abs(alpha - lvl) < 1e-9), by = c("series", "Method", "h", "age")) |>
      inner_join(app$emp_mult |> filter(series == input$eb_series, Rule == input$eb_rule, abs(alpha - lvl) < 1e-9), by = c("series", "Method", "alpha", "t0", "h")) |>
      mutate(half = mult * gamma, lower = forecast - half, upper = forecast + half,
             Method = factor(Method, levels = emp_methods),
             panel = factor(paste0("h = ", h, "  (", year, ")"), levels = unique(paste0("h = ", sort(unique(h)), "  (", year[match(sort(unique(h)), h)], ")"))))
  })
  card_server(
    "eb_card",
    reactive({
      x <- eb_data()
      ggplot(x) +
        geom_ribbon(aes(age, ymin = lower, ymax = upper, fill = Method), alpha = 0.28) +
        geom_line(aes(age, forecast, colour = Method, linetype = Method), linewidth = 0.8) +
        geom_line(data = x |> distinct(panel, age, observed), aes(age, observed), colour = "black", linewidth = 0.7) +
        facet_wrap(~ panel, ncol = 5, scales = "free_y") +
        scale_colour_manual(name = "Forecasting method", values = method_colours) +
        scale_fill_manual(name = "Forecasting method", values = method_colours) +
        scale_linetype_manual(name = "Forecasting method", values = method_linetypes) +
        labs(x = "Age", y = scale_of(input$eb_series)) +
        theme(strip.text = element_text(size = 11), axis.text = element_text(size = 10))
    }),
    reactive({
      a <- as.numeric(input$eb_level)
      raw <- eb_data() |> group_by(h, year, Method) |>
        summarise(Coverage = mean(abs(observed - forecast) <= half),
                  IS = mean(2 * half + (2 / a) * pmax(abs(observed - forecast) - half, 0)), .groups = "drop") |>
        pivot_wider(names_from = Method, values_from = c(Coverage, IS), names_sep = ": ") |>
        mutate(h = as.character(h), year = as.character(year))
      cov_cols <- grep("^Coverage", names(raw), value = TRUE)
      is_cols <- grep("^IS", names(raw), value = TRUE)
      shown <- mark_min(raw, list(is_cols), highlight = TRUE)
      shown[cov_cols] <- lapply(raw[cov_cols], fmt)
      dt_table(shown)
    }),
    reactive({
      req(input$eb_origin, as.integer(input$eb_origin) %in% app$emp_forecasts$t0[app$emp_forecasts$series == input$eb_series])
      o <- app$emp_forecasts |> filter(series == input$eb_series, t0 == as.integer(input$eb_origin)) |> summarise(o = min(year) - 1L) |> pull(o)
      paste(c(paste0("Forecasts and ", input$eb_rule, " ", if (input$eb_level == "0.2") "80%" else "95%", " prediction bands from origin ", o),
              "black lines are the observed curves", input$eb_series), collapse = dot)
    }),
    reactive(640)
  )
}

shinyApp(ui, server)
