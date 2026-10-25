<br />
<div align="center">
  <a href="https://www.mq.edu.au">
    <img src="Macquarie_University_logo.png" alt="Macquarie University" width="360">
  </a>

<h3 align="center">
Order Selection and Forecasting in Long-Memory Functional Time Series
</h3>

<p align="center">
Chang Liu
<br />
Supervisors: Han Lin Shang · Samuel Müller · Garth Tarr
</p>

</div>

---

# Abstract

Functional time series (FTS) are increasingly used to model data observed over continuous domains such as time, wavelength, and age. Functional autoregressive moving-average (FARMA) models provide a flexible framework for capturing temporal dependence; however, reliable model specification, particularly consistent order selection, remains a fundamental challenge. This issue is further complicated by the presence of long-range dependence or non-stationarity, both of which violate standard modelling assumptions.

This thesis studies order selection for FARMA models within a framework that accommodates long-memory and non-stationary functional time series, where the functional autoregressive fractionally integrated moving-average (FARFIMA) model is used to capture stationary and non-stationary long memory. The proposed framework consists of three main steps. First, the memory parameter is estimated and fractional differencing is used to obtain a stationary representation of the process. Second, order selection and FARMA model fitting are conducted in this stationary domain. Finally, forecasts with conformal prediction intervals are obtained on the original non-stationary scale.

This thesis makes two main contributions. First, two order selection approaches for FARFIMA models in long-memory functional time series are developed: a functional Box–Jenkins procedure based on residual diagnostics and an adaptive functional Hannan–Quinn information criterion. Second, a unified framework is proposed that integrates memory estimation, order selection, model estimation, and forecasting. Simulation studies show that the functional information criterion provides the most stable order selection, and that both memory estimation and short-memory order selection improve forecasts relative to a fixed benchmark, with the former giving larger gains at longer horizons and the latter mainly at short horizons. Applications to England and Wales fertility and French female and male mortality rates illustrate the practical usefulness of the proposed methodology.

---

# Repository Structure

The folders follow the steps of the framework:

| Section | Folder |
|---|---|
| 1. Simulation setting | `DGP non-stationary FTS` |
| 2. Memory parameter estimation | `Fractional Differencing` |
| 3. FARMA order selection | `Order Selection` |
| 4. Model estimation | `Operator Estimation` |
| 5. Forecasting | `Forecasting` |
| 6. Empirical data analysis | `Empirical Study` |
| 7. Shiny app | `Shiny App` |

---

# 1. Simulation Setting

| File | Description |
|---|---|
| `sim_FARMA.R` | Simulation of FARFIMA(p,d,q) processes with Brownian motion innovations, used in all simulation studies |

---

# 2. Memory Parameter Estimation

## Main Files

| File | Description |
|---|---|
| `local.W.R` | Local Whittle estimator (LW) and detrended fluctuation analysis (DFA) |
| `LPWN.R` | Local polynomial Whittle with noise estimators: LPWN<sub>1</sub> and LPWN<sub>2</sub> |
| `sieve_bootstrap.R` | Sieve bootstrap for functional time series |
| `Sieve-bootstrapped corrected LW estimation.R` | Bootstrap bias-corrected estimators: BC-LPWN<sub>1</sub> and BC-LPWN<sub>2</sub> |
| `Covariance Matrix.R` | Covariance-type operators for the leading eigenfunction: the empirical covariance operator (Covar) and the long-run covariance operators (RS and LRS) |
| `cov_compare_sim.Rmd` | Simulation comparing LW and DFA under the three covariance-type operators |
| `demo differencing.Rmd` | Illustration of fractional differencing for scalar and functional time series |
| `hurst.R` | Hurst exponent estimators from the fArma package, kept for reference |

## Simulation Result Files

| File | Description |
|---|---|
| `covariance_matrix_comparison_FAR.rds`, `covariance_matrix_comparison_FAR_complmentary.rds` | Results for FARFIMA(1,d,0) |
| `covariance_matrix_comparison_FARMA.rds`, `covariance_matrix_comparison_FARMA_supplementary.rds` | Results for FARFIMA(1,d,1) |
| `LPWN_LW_DFA_sim_comparison_FAR.rds`, `LPWN_LW_DFA_sim_comparison_FARMA.rds` | Comparison of LW, DFA and LPWN |

Figures are in `Sim Plots`.

---

# 3. FARMA Order Selection

## Main Files

| File | Description |
|---|---|
| `Marginal BIC.R` | Adaptive functional Hannan–Quinn information criterion, with exhaustive search (functional IC) and automatic search (auto functional IC) |
| `WN test order selection.R` | Functional Box–Jenkins procedure based on the functional portmanteau test |
| `CV order selection.R` | Forecast-based cross-validation (CV) |
| `Sim Order Selection.R` | Simulation comparing the order selection methods |
| `Order selection AIC simulation.Rmd` | Order selection simulation study and its figures |
| `Order Selection WN simulation.Rmd` | Checks of the functional white-noise test |

## Simulation Result Files

| File | Description |
|---|---|
| `sim_order_compare.rds` | Order selection results for all simulation settings |

Figures are in `Sim Plots`.

---

# 4. Model Estimation

| File | Description |
|---|---|
| `Kuenzer FARMA Estimation.R` | Functional Hannan–Rissanen estimation of the FARMA operators (Kuenzer, 2024), used in all other steps |
| `EstimationFARMA.R` | Comparison of the FARMA estimators of Kuenzer and Kühnert |

---

# 5. Forecasting

## Main Files

| File | Description |
|---|---|
| `forecast_helpers.R` | Forecasts on the original scale and the scaled residual (SR), non-conformity score (NS) and adaptive conformal (AC) prediction intervals |
| `ideal_pipeline.R` | Rolling forecasts with known memory parameter and orders |
| `unified_pipeline.R` | Unified procedure: memory estimation, order selection, model estimation and forecasting |
| `forecast_plot_helpers.R` | Functions for figures and tables |
| `ideal_analysis.Rmd` | Simulation study: forecasting with known memory and order |
| `unified_analysis.Rmd` | Simulation study: unified FARFIMA framework |

## Simulation Result Files

| File | Description |
|---|---|
| `ideal_raw.rds`, `ideal_summary.rds` | Results of forecasting with known memory and order |

The results of the unified framework exceed the GitHub file size limit and are not included. They are produced by `unified_analysis.Rmd`.

## Plot Folders

| Folder | Description |
|---|---|
| `Sim Plots Ideal` | Figures for forecasting with known memory and order |
| `Sim Plots Unified` | Figures for the unified framework |
| `Plots Slides` | Figures for presentations |

---

# 6. Empirical Data Analysis

## Main Files

| File | Description |
|---|---|
| `Empirical analysis.Rmd` | Data processing, FACF, forecasts and their evaluation for England and Wales fertility and French female and male mortality |
| `save_function.R` | Functions for saving figures |

## Datasets

| File | Description |
|---|---|
| `Fertility Data/asfrRR.txt`, `Fertility Data/exposRR.txt` | Age-specific fertility rates and exposures from the Human Fertility Database (England and Wales, ages 15 to 49) |
| `Mortality Data/FRATNP_Mx_1x1.txt`, `Mortality Data/FRATNP_Exposures_1x1.txt` | Age-specific death rates and exposures from the Human Mortality Database (France) |

## Results and Figures

| Folder | Description |
|---|---|
| `results` | Forecast results for each dataset |
| `Empirical Plots` | Figures used in the thesis |
| `Shiny Screenshots` | Screenshots of the Shiny app |

---

# 7. Shiny App

The simulation and empirical results can be explored in the Shiny app: <https://changliu0122.shinyapps.io/Order_selection_shiny/>

| File | Description |
|---|---|
| `app.R` | The Shiny app |
| `intro_pages.R` | Introduction pages of the app |
| `prepare_app_data.R` | Prepares the app data from the simulation results |
| `app_data.rds` | Data used by the app |

---

# Suggested Order of Running Files

Open `FTS Order Selection.Rproj` first, as all files use paths relative to the project folder. Each notebook loads the scripts it needs.

1. `Fractional Differencing/cov_compare_sim.Rmd`
2. `Order Selection/Order selection AIC simulation.Rmd`
3. `Forecasting/ideal_analysis.Rmd`
4. `Forecasting/unified_analysis.Rmd`
5. `Empirical Study/Empirical analysis.Rmd`
6. `Shiny App/prepare_app_data.R`, then `Shiny App/app.R`

---

# Contact

### Chang Liu

Department of Actuarial Studies and Business Analytics

Macquarie University

Email: chang.liu45@students.mq.edu.au

---

### Han Lin Shang

Department of Actuarial Studies and Business Analytics

Macquarie University

Email: hanlin.shang@mq.edu.au

---

### Samuel Müller

Macquarie University

Email: samuel.muller@mq.edu.au

---

### Garth Tarr

University of Sydney

Email: garth.tarr@sydney.edu.au
