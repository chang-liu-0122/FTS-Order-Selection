# Thesis figure audit

## Recommended main-text figures

### Monte Carlo simulation study

1. **Memory estimation:** use the vertically combined LW--DFA figure for FARFIMA\((1,d,0)\). If a FARMA comparison is needed, use the corresponding combined FARFIMA\((1,d,1)\), \(\nu=0.5\), figure. Each combined figure has one legend and one pair of shared axes. The separate estimator figures and the remaining moving-average strengths belong in the appendix.
2. **Order selection:** retain both order-accuracy figures:
   - `order_accuracy_far_fma_2x2.pdf`
   - `order_accuracy_farma_2x2.pdf`
3. **Ideal forecasting:** retain:
   - `forecast_rmse_ratio_2x2.pdf`
   - `forecast_ecp_80_expanding_2x2.pdf`
   - `forecast_score_ratio_80_2x2.pdf`
   - `forecast_wis_ratio_2x2.pdf`
4. **Unified forecasting:** retain:
   - `unified_d_estimation_T1000.pdf`
   - `unified_order_accuracy_T1000_2x2.pdf`
   - `unified_rmse_ratio_balanced_expanding_T1000_2x2.pdf`
   - `unified_ecp_80_expanding_2x2.pdf`
   - `unified_score_ratio_80_expanding_2x2.pdf`
   - `unified_wis_ratio_expanding_2x2.pdf`

Use only \(T=1000\) and the expanding-window scheme in the main ideal and unified forecasting figures. This makes the two studies directly comparable and keeps each retained figure focused on estimation of \(d\), recovery of \((p,q)\), point accuracy, coverage, or overall probabilistic accuracy. The other sample sizes and window schemes remain useful as robustness results.

## Interactive supplementary figures

The following supporting settings remain reproducible from the simulation results and R Markdown files, but are not exported to the thesis figure folders. They may later be exposed through a Shiny application:

- separate LW and DFA memory-estimation figures after the combined main-text figure is created;
- ideal-study figures for \(T=500\) and fixed-window forecasting;
- FARFIMA\((1,d,1)\) memory-estimation figures for \(\nu=0.2\) and \(\nu=0.8\);
- all \(T=500\) unified figures;
- fixed-window ideal coverage figures;
- all 95% empirical-coverage-probability figures;
- 95% interval-score-ratio figures;
- full raw-value tables if chapter length is a concern. Their formatting can remain unchanged.

## Empirical chapter

For each application, retain the data display, dependence diagnostic, point-ratio figure, ECP figure, WIS-ratio figure, and one ten-horizon forecast illustration in the main text. Move the separate interval-score-ratio figure and detailed raw-value tables to the appendix when space is limited.

Recommended main-text set:

- fertility rainbow plot;
- fertility fACF;
- fertility RMSE/MAE ratios;
- fertility ECP;
- fertility WIS ratio;
- fertility ten-horizon forecast illustration;
- mortality rainbow plot;
- female/male mortality fACFs;
- female/male RMSE/MAE ratios;
- female/male ECP;
- female/male WIS ratios;
- one female and one male ten-horizon forecast illustration, or a single combined figure if the panels remain readable.

## Plot specification applied in the code

- no R-generated plot title or subtitle;
- legend at the top and centred;
- equal legend key widths and horizontal spacing;
- benchmark: green square with a long-dashed line;
- Functional IC: orange circle with a solid line;
- CV: purple triangle with a dashed line;
- compact benchmark label `FARFIMA(1,1,0)`;
- dark reference lines for ratio one and nominal coverage;
- shared axes for combined panels;
- separate export sizes for paper figures and slides;
- white background, grey facet strips, and a consistent classic paper theme.

## Additional figure recommendation

No additional performance plot is necessary. The only useful new figure is the combined 2-by-2 memory-estimation comparison described above. A compact simulation-design table is preferable to another results figure for summarising the parameter grid, sample sizes, horizons, and repetitions.
