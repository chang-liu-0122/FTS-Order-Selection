settings_table <- function(rows) {
  tags$table(class = "settings-table", tags$tbody(lapply(names(rows), function(k) tags$tr(tags$th(HTML(k)), tags$td(HTML(rows[[k]]))))))
}

intro_section <- function(kicker, title, ..., target = NULL) {
  div(
    class = "intro-section",
    div(class = "kicker", kicker),
    h2(title),
    ...,
    if (!is.null(target)) div(class = "intro-actions", actionButton(paste0("go_", target), "View results \u2192", class = "btn-primary"))
  )
}

intro_overview <- intro_section(
  "Overview", "About this site",
  p("This site presents the simulation and empirical results for order selection and forecasting in long-memory functional time series. The framework estimates the memory parameter, fractionally differences the series to obtain a stationary representation, selects the FARMA order in this stationary domain, and produces forecasts with conformal prediction bands on the original scale."),
  p("The simulations first examine each component separately, namely memory estimation, order selection, and forecasting with known memory and order, and then evaluate the unified framework, in which estimation errors propagate through the whole procedure. The empirical study applies the framework to fertility and mortality data."),
  p("Select a section on the left to read its design, then open its results. Every results page shows the selected setting as a plot or as a table.")
)

intro_dgp <- intro_section(
  "Simulation", "Data-generating process",
  p("The simulated curves follow a functional FARFIMA\\((p,d,q)\\) model, \\((1-B)^d X_t = Y_t\\), where \\(Y_t\\) is a functional FARMA\\((p,q)\\) process driven by Brownian motion innovations. Each curve is observed at 101 grid points on \\([0,1]\\)."),
  p("The autoregressive kernels have a Gaussian shape and the moving-average kernels a \\(\\min(u,v)\\) shape. Their overall strengths are controlled by \\(\\psi\\) for the autoregressive part and \\(\\nu\\) for the moving-average part, shared equally across lags. Each setting is replicated 100 times.")
)

intro_memory <- intro_section(
  "Simulation", "Memory parameter estimation",
  p("The curves are projected onto the leading eigenfunction of a covariance-type operator, and the memory parameter \\(d\\) is estimated from the resulting score series. Three operators are compared: the empirical covariance operator (Covar), the long-run covariance operator with an adaptive bandwidth (RS), and the long-run covariance operator in sandwich form (LRS). For each, \\(d\\) is estimated by the local Whittle (LW) and detrended fluctuation analysis (DFA) estimators, and accuracy is summarised by the bias and mean squared error of \\(\\widehat d\\)."),
  settings_table(list(
    "Models" = "FARFIMA\\((1,d,0)\\) and FARFIMA\\((1,d,1)\\)",
    "Sample size \\(T\\)" = "500, 1000, 2000",
    "Memory parameter \\(d\\)" = "0, 0.4, 0.8, 1.2",
    "Strengths" = "\\(\\psi, \\nu \\in \\{0.2, 0.5, 0.8\\}\\)"
  )),
  target = "memory"
)

intro_order <- intro_section(
  "Simulation", "FARMA order selection",
  p("Order selection is examined separately from memory estimation by generating short-memory FARMA data directly, so that the selected order can be compared with the true order. Candidate orders \\((p,q)\\) range over \\(0 \\le p, q \\le 5\\), excluding \\((0,0)\\). Functional IC applies the adaptive functional Hannan\u2013Quinn criterion to every candidate, while the functional Box\u2013Jenkins procedure fits candidates sequentially until the residuals pass a white-noise test. Performance is the proportion of replications in which the true order is recovered."),
  settings_table(list(
    "True orders" = "FAR: (1,0), (2,0); FMA: (0,1), (0,2); FARMA: (1,1), (2,1), (1,2), (2,2)",
    "Sample size \\(T\\)" = "1000, 2000, 5000",
    "Strengths" = "\\(\\psi, \\nu \\in \\{0.2, 0.5, 0.8\\}\\)"
  )),
  target = "order"
)

intro_ideal <- intro_section(
  "Simulation", "Forecasting with known memory and order",
  p("Forecasting is first compared in an ideal setting where \\(d\\) and \\((p,q)\\) are known. The plug-in method fits the correctly specified model, whereas the benchmark always fits FARFIMA\\((1,1,0)\\), which isolates the effect of the short-memory specification."),
  settings_table(list(
    "True orders" = "(1,0), (0,1), (2,0), (1,1)",
    "Memory parameter" = "\\(d = 1\\)",
    "Sample size \\(T\\)" = "500, 1000",
    "Strengths" = "\\(\\psi, \\nu \\in \\{0.2, 0.5, 0.8\\}\\)",
    "Design" = "60/20/20 split, \\(h = 1,\\ldots,10\\), expanding and fixed windows"
  )),
  target = "ideal"
)

intro_unified <- intro_section(
  "Simulation", "Unified FARFIMA framework",
  p("The complete procedure is applied: \\(d\\) is estimated by BC-LPWN, the series is fractionally differenced, the order is selected, and forecasts are produced. Functional IC selects the order by the adaptive information criterion, whereas CV selects the order with the best validation forecast performance. Both are compared with the FARFIMA\\((1,1,0)\\) benchmark."),
  settings_table(list(
    "True orders" = "(1,0), (0,1), (2,0), (1,1)",
    "Memory parameter \\(d\\)" = "0.2, 0.4, 0.6, 0.8, 1, 1.2, 1.4",
    "Sample size \\(T\\)" = "500, 1000",
    "Strengths" = "\\(\\psi, \\nu \\in \\{0.2, 0.5, 0.8\\}\\)",
    "Design" = "As in the ideal setting"
  )),
  target = "unified"
)

intro_empirical <- intro_section(
  "Empirical study", "Fertility and mortality",
  p("The framework is applied to age-specific fertility rates for England and Wales and to age-specific mortality rates for French females and males, with each year treated as a curve over age. The fertility memory estimate lies well above one, whereas the mortality estimates are close to one."),
  p("Fertility is modelled on the Box\u2013Cox scale with \\(\\lambda = 0.4\\) and mortality on the \\(\\log_{10}\\) scale. The models are fitted to smoothed curves, and calibration and evaluation use the transformed raw curves."),
  settings_table(list(
    "Fertility" = "England and Wales, 1938\u20132022, ages 15\u201349",
    "Mortality" = "France, 1901\u20132023, ages 0\u201395+, by sex",
    "Design" = "60/20/20 split, expanding window, \\(h = 1,\\ldots,10\\)",
    "Methods" = "Functional IC, CV, and the FARFIMA\\((1,1,0)\\) benchmark"
  )),
  target = "empirical"
)

intro_evaluation <- intro_section(
  "Forecast evaluation", "Criteria and prediction bands",
  settings_table(list(
    "RMSPE, MAPE" = "Root mean squared and mean absolute prediction errors over test origins and grid points",
    "ECP" = "Empirical coverage probability of the prediction band",
    "CPD" = "Absolute deviation of the ECP from the nominal coverage",
    "IS" = "Interval score: band width plus a penalty for observations outside the band",
    "WIS" = "Weighted interval score combining MAPE with the IS of the 80% and 95% bands"
  )),
  p("Ratios are taken relative to FARFIMA\\((1,1,0)\\) under the same conformal method, and values below one favour the corresponding method."),
  h4("Conformal prediction bands"),
  settings_table(list(
    "SR" = "Scaled residual band: a common multiplier of the validation error scale",
    "NS" = "Non-conformity score band: an empirical quantile of width-standardised scores",
    "AC" = "Adaptive conformal band: the NS quantile level is updated over the test period according to recent coverage, with a learning rate chosen on the validation sample"
  )),
  p("All three bands are calibrated on the validation sample. AC reduces to NS when its learning rate is zero.")
)
