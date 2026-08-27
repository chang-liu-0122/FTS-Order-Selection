# Empirical chapter draft

This folder is a drop-in draft for the thesis project in
`Order_Selection_and_Forecasting_in_Long_Memory_Functional_Time_Series.zip`.

- Replace the thesis stub `Chapters/05-Real Data Application.tex` with the file
  in `Chapters/` after reviewing the prose.
- Copy `Figures/Empirical/` into the thesis `Figures/` directory.
- Merge the entries in `Bibliography/empirical-additions.bib` into the thesis
  bibliography file.
- Re-run `scripts/make_forecast_example.R` from the repository root whenever
  the saved empirical results change.

The draft intentionally reports forecasts against raw observations while the
models are fitted to smoothed curves. It also distinguishes stationary long
memory in the fertility changes from nonstationary fractional integration in
the mortality series.

Before final submission, record the exact data-download dates in the HFD and
HMD bibliography entries. The corrected fertility figures in this folder use
ages 15--49; the older copies in `Empirical Study/Empirical Plots/` use curve
indices 1--35 on the horizontal axis.
