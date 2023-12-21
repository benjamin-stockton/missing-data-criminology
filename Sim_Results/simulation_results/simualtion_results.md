---
title: "Criminology Simulation Results"
author: "Ben Stockton"
number-sections: true 
df-print: paged
pdf-engine: pdflatex
toc: true
format: 
  html: 
    code-fold: false
  pdf:
    geometry: 
      - top=30mm
      - left=30mm
    shift-heading-level-by: -1
  docx:
    shift-heading-level-by: -1
editor: visual
fig-width: 7
fig-height: 5
fig-dpi: 300
keep-tex: true
keep-md: true
execute: 
  cache: true
bibliography: ../../Literature/Criminology.bib
---





## Overview

I ran the simulations under a variety of settings, each for 1000 iterations.

| Sample Size ($n$) | Number of Imputations ($M$) |
|-------------------|-----------------------------|
| 500               | 3                           |
| 500               | 5                           |
| 500               | 8                           |
| 1000              | 3                           |
| 1000              | 5                           |
| 1000              | 8                           |
| 2500              | 3                           |
| 2500              | 5                           |
| 2500              | 8                           |

From the estimated rate of missing information, I determined that there should be between five and eight imputations to achieve a relative efficiency of \>0.99. For a relative efficiency of \~0.95, three imputations is sufficient.

In previous simulations with $n = 1000$ and $M = 3$ , I was able to find that on average, depending on my parameter settings for the missingness model, I was able to create data sets that over- and under-estimated the race effect on average for each data set. In other words, for a given data set, if you take a sample from the population with race effect $\gamma = \exp(\beta)$ and estimate the race effect with the complete data to get $\hat{\gamma} = \exp(\hat{\beta})$, then estimates made with complete case analysis after imposing missing values for that same data set will be over or under the complete data estimate $\hat{\gamma}$ on average, depending on what parameters are used in the missingness model.

In this document, I will demonstrate that for appropriate numbers of imputations to achieve 95% relative efficiency, 99% relative efficiency, and 99+% relative efficiency, we see the same behavior in the sampling distributions and medians of the differences $\hat{\gamma}^{miss} - \hat{\gamma}$.

## Results

### n = 500


::: {.cell hash='simualtion_results_cache/pdf/n = 500 Results Tables_2d5d76a861e2e022ae8717bf95aa6175'}
::: {.cell-output .cell-output-stderr}
```
Warning: `aes_string()` was deprecated in ggplot2 3.0.0.
i Please use tidy evaluation idioms with `aes()`.
i See also `vignette("ggplot2-in-packages")` for more information.
```
:::

::: {.cell-output-display}
![](simualtion_results_files/figure-pdf/n = 500 Results Tables-1.pdf)
:::

::: {.cell-output-display}
\begin{tabular}[t]{llllrrr}
\toprule
Num of Imp & Miss Type & Direction & Analysis & OR Difference & 0.025 Quantile & 0.975 Quantile\\
\midrule
3 imps & MAR & OVER & CCA & 0.015 & -0.063 & 0.107\\
3 imps & MAR & OVER & MI & 0.000 & -0.080 & 0.074\\
3 imps & MAR & UNDR & CCA & -0.045 & -0.141 & 0.040\\
3 imps & MAR & UNDR & MI & 0.002 & -0.065 & 0.073\\
3 imps & MAR & No Missing Data & COMP & 0.000 & 0.000 & 0.000\\
\addlinespace
5 imps & MAR & OVER & CCA & 0.014 & -0.074 & 0.108\\
5 imps & MAR & OVER & MI & 0.001 & -0.080 & 0.073\\
5 imps & MAR & UNDR & CCA & -0.043 & -0.147 & 0.036\\
5 imps & MAR & UNDR & MI & 0.001 & -0.059 & 0.075\\
5 imps & MAR & No Missing Data & COMP & 0.000 & 0.000 & 0.000\\
\addlinespace
8 imps & MAR & OVER & CCA & 0.015 & -0.067 & 0.111\\
8 imps & MAR & OVER & MI & 0.000 & -0.076 & 0.076\\
8 imps & MAR & UNDR & CCA & -0.043 & -0.150 & 0.036\\
8 imps & MAR & UNDR & MI & 0.002 & -0.060 & 0.073\\
8 imps & MAR & No Missing Data & COMP & 0.000 & 0.000 & 0.000\\
\bottomrule
\end{tabular}
:::

::: {.cell-output-display}
\begin{tabular}[t]{llrrr}
\toprule
Miss Type & Direction & Prop. of Inc. Cases & 0.025 Quantile & 0.975 Quantile\\
\midrule
MAR & OVER & 0.007 & 0.006 & 0.008\\
MAR & UNDR & 0.008 & 0.008 & 0.008\\
\bottomrule
\end{tabular}
:::

::: {.cell-output-display}
![](simualtion_results_files/figure-pdf/n = 500 Results Tables-2.pdf)
:::

::: {.cell-output-display}
\begin{tabular}[t]{llllrrr}
\toprule
Num of Imp & Miss Type & Direction & Analysis & Bias & 0.025 Quantile & 0.975 Quantile\\
\midrule
3 imps & MAR & OVER & CCA & 0.155 & -0.562 & 1.170\\
3 imps & MAR & OVER & MI & 0.140 & -0.563 & 1.128\\
3 imps & MAR & UNDR & CCA & 0.095 & -0.586 & 1.085\\
3 imps & MAR & UNDR & MI & 0.142 & -0.564 & 1.153\\
3 imps & MAR & No Missing Data & COMP & 0.140 & -0.568 & 1.125\\
\addlinespace
5 imps & MAR & OVER & CCA & 0.147 & -0.522 & 1.225\\
5 imps & MAR & OVER & MI & 0.133 & -0.524 & 1.177\\
5 imps & MAR & UNDR & CCA & 0.089 & -0.544 & 1.082\\
5 imps & MAR & UNDR & MI & 0.134 & -0.528 & 1.182\\
5 imps & MAR & No Missing Data & COMP & 0.133 & -0.524 & 1.177\\
\addlinespace
8 imps & MAR & OVER & CCA & 0.133 & -0.562 & 1.225\\
8 imps & MAR & OVER & MI & 0.118 & -0.565 & 1.189\\
8 imps & MAR & UNDR & CCA & 0.074 & -0.591 & 1.078\\
8 imps & MAR & UNDR & MI & 0.120 & -0.558 & 1.183\\
8 imps & MAR & No Missing Data & COMP & 0.118 & -0.559 & 1.163\\
\bottomrule
\end{tabular}
:::
:::


### n = 1000


::: {.cell hash='simualtion_results_cache/pdf/n = 1000 Results Tables_2f603f30bda828dca25850daf70ab8e8'}
::: {.cell-output-display}
![](simualtion_results_files/figure-pdf/n = 1000 Results Tables-1.pdf)
:::

::: {.cell-output-display}
\begin{tabular}[t]{llllrrr}
\toprule
Num of Imp & Miss Type & Direction & Analysis & OR Difference & 0.025 Quantile & 0.975 Quantile\\
\midrule
3 imps & MAR & OVER & CCA & 0.014 & -0.031 & 0.060\\
3 imps & MAR & OVER & MI & -0.001 & -0.047 & 0.040\\
3 imps & MAR & UNDR & CCA & -0.049 & -0.100 & -0.002\\
3 imps & MAR & UNDR & MI & 0.001 & -0.035 & 0.040\\
3 imps & MAR & No Missing Data & COMP & 0.000 & 0.000 & 0.000\\
\addlinespace
3 imps & MNAR & OVER & CCA & 0.031 & -0.002 & 0.072\\
3 imps & MNAR & OVER & MI & 0.030 & -0.008 & 0.077\\
3 imps & MNAR & UNDR & CCA & -0.047 & -0.105 & 0.007\\
3 imps & MNAR & UNDR & MI & -0.047 & -0.118 & 0.010\\
3 imps & MNAR & No Missing Data & COMP & 0.000 & 0.000 & 0.000\\
\addlinespace
5 imps & MAR & OVER & CCA & 0.014 & -0.036 & 0.059\\
5 imps & MAR & OVER & MI & -0.001 & -0.045 & 0.038\\
5 imps & MAR & UNDR & CCA & -0.049 & -0.102 & -0.001\\
5 imps & MAR & UNDR & MI & 0.001 & -0.034 & 0.037\\
5 imps & MAR & No Missing Data & COMP & 0.000 & 0.000 & 0.000\\
\addlinespace
5 imps & MNAR & OVER & CCA & 0.031 & 0.003 & 0.070\\
5 imps & MNAR & OVER & MI & 0.031 & -0.001 & 0.080\\
5 imps & MNAR & UNDR & CCA & -0.046 & -0.109 & 0.007\\
5 imps & MNAR & UNDR & MI & -0.045 & -0.112 & 0.007\\
5 imps & MNAR & No Missing Data & COMP & 0.000 & 0.000 & 0.000\\
\addlinespace
8 imps & MAR & OVER & CCA & 0.013 & -0.036 & 0.060\\
8 imps & MAR & OVER & MI & -0.001 & -0.046 & 0.036\\
8 imps & MAR & UNDR & CCA & -0.049 & -0.106 & -0.003\\
8 imps & MAR & UNDR & MI & 0.002 & -0.034 & 0.039\\
8 imps & MAR & No Missing Data & COMP & 0.000 & 0.000 & 0.000\\
\addlinespace
8 imps & MNAR & OVER & CCA & 0.031 & 0.002 & 0.071\\
8 imps & MNAR & OVER & MI & 0.031 & 0.000 & 0.071\\
8 imps & MNAR & UNDR & CCA & -0.046 & -0.108 & 0.011\\
8 imps & MNAR & UNDR & MI & -0.046 & -0.108 & 0.010\\
8 imps & MNAR & No Missing Data & COMP & 0.000 & 0.000 & 0.000\\
\bottomrule
\end{tabular}
:::

::: {.cell-output-display}
\begin{tabular}[t]{llrrr}
\toprule
Miss Type & Direction & Prop. of Inc. Cases & 0.025 Quantile & 0.975 Quantile\\
\midrule
MAR & OVER & 0.008 & 0.008 & 0.009\\
MAR & UNDR & 0.009 & 0.009 & 0.009\\
MNAR & OVER & 0.008 & 0.008 & 0.008\\
MNAR & UNDR & 0.008 & 0.008 & 0.008\\
\bottomrule
\end{tabular}
:::

::: {.cell-output-display}
![](simualtion_results_files/figure-pdf/n = 1000 Results Tables-2.pdf)
:::

::: {.cell-output-display}
\begin{tabular}[t]{llllrrr}
\toprule
Num of Imp & Miss Type & Direction & Analysis & Bias & 0.025 Quantile & 0.975 Quantile\\
\midrule
3 imps & MAR & OVER & CCA & 0.057 & -0.375 & 0.579\\
3 imps & MAR & OVER & MI & 0.043 & -0.386 & 0.547\\
3 imps & MAR & UNDR & CCA & -0.006 & -0.415 & 0.492\\
3 imps & MAR & UNDR & MI & 0.045 & -0.380 & 0.534\\
3 imps & MAR & No Missing Data & COMP & 0.044 & -0.385 & 0.546\\
\addlinespace
3 imps & MNAR & OVER & CCA & 0.092 & -0.369 & 0.680\\
3 imps & MNAR & OVER & MI & 0.091 & -0.370 & 0.669\\
3 imps & MNAR & UNDR & CCA & 0.014 & -0.429 & 0.576\\
3 imps & MNAR & UNDR & MI & 0.014 & -0.434 & 0.576\\
3 imps & MNAR & No Missing Data & COMP & 0.062 & -0.389 & 0.631\\
\addlinespace
5 imps & MAR & OVER & CCA & 0.066 & -0.354 & 0.629\\
5 imps & MAR & OVER & MI & 0.052 & -0.357 & 0.597\\
5 imps & MAR & UNDR & CCA & 0.003 & -0.396 & 0.543\\
5 imps & MAR & UNDR & MI & 0.054 & -0.348 & 0.620\\
5 imps & MAR & No Missing Data & COMP & 0.053 & -0.354 & 0.611\\
\addlinespace
5 imps & MNAR & OVER & CCA & 0.079 & -0.362 & 0.658\\
5 imps & MNAR & OVER & MI & 0.079 & -0.363 & 0.657\\
5 imps & MNAR & UNDR & CCA & 0.002 & -0.420 & 0.541\\
5 imps & MNAR & UNDR & MI & 0.003 & -0.416 & 0.543\\
5 imps & MNAR & No Missing Data & COMP & 0.048 & -0.380 & 0.614\\
\addlinespace
8 imps & MAR & OVER & CCA & 0.057 & -0.375 & 0.656\\
8 imps & MAR & OVER & MI & 0.042 & -0.375 & 0.623\\
8 imps & MAR & UNDR & CCA & -0.005 & -0.419 & 0.553\\
8 imps & MAR & UNDR & MI & 0.046 & -0.376 & 0.625\\
8 imps & MAR & No Missing Data & COMP & 0.044 & -0.383 & 0.632\\
\addlinespace
8 imps & MNAR & OVER & CCA & 0.078 & -0.356 & 0.617\\
8 imps & MNAR & OVER & MI & 0.078 & -0.356 & 0.622\\
8 imps & MNAR & UNDR & CCA & 0.001 & -0.409 & 0.515\\
8 imps & MNAR & UNDR & MI & 0.001 & -0.410 & 0.515\\
8 imps & MNAR & No Missing Data & COMP & 0.047 & -0.374 & 0.576\\
\bottomrule
\end{tabular}
:::
:::


### n = 2500


::: {.cell hash='simualtion_results_cache/pdf/n = 2500 Results Tables_125b14c6c79c70cbb3d60fb7a2aa4498'}
::: {.cell-output-display}
![](simualtion_results_files/figure-pdf/n = 2500 Results Tables-1.pdf)
:::

::: {.cell-output-display}
\begin{tabular}[t]{llllrrr}
\toprule
Num of Imp & Miss Type & Direction & Analysis & OR Difference & 0.025 Quantile & 0.975 Quantile\\
\midrule
3 imps & MAR & OVER & CCA & 0.014 & -0.012 & 0.040\\
3 imps & MAR & OVER & MI & -0.002 & -0.030 & 0.023\\
3 imps & MAR & UNDR & CCA & -0.044 & -0.074 & -0.018\\
3 imps & MAR & UNDR & MI & 0.002 & -0.019 & 0.024\\
3 imps & MAR & No Missing Data & COMP & 0.000 & 0.000 & 0.000\\
\addlinespace
3 imps & MNAR & OVER & CCA & 0.033 & 0.015 & 0.054\\
3 imps & MNAR & OVER & MI & 0.030 & 0.010 & 0.056\\
3 imps & MNAR & UNDR & CCA & -0.052 & -0.088 & -0.019\\
3 imps & MNAR & UNDR & MI & -0.048 & -0.087 & -0.014\\
3 imps & MNAR & No Missing Data & COMP & 0.000 & 0.000 & 0.000\\
\addlinespace
5 imps & MAR & OVER & CCA & 0.014 & -0.012 & 0.039\\
5 imps & MAR & OVER & MI & -0.001 & -0.026 & 0.021\\
5 imps & MAR & UNDR & CCA & -0.044 & -0.074 & -0.019\\
5 imps & MAR & UNDR & MI & 0.001 & -0.019 & 0.021\\
5 imps & MAR & No Missing Data & COMP & 0.000 & 0.000 & 0.000\\
\addlinespace
5 imps & MNAR & OVER & CCA & 0.033 & 0.015 & 0.053\\
5 imps & MNAR & OVER & MI & 0.030 & 0.011 & 0.053\\
5 imps & MNAR & UNDR & CCA & -0.052 & -0.088 & -0.019\\
5 imps & MNAR & UNDR & MI & -0.049 & -0.084 & -0.015\\
5 imps & MNAR & No Missing Data & COMP & 0.000 & 0.000 & 0.000\\
\addlinespace
8 imps & MAR & OVER & CCA & 0.014 & -0.011 & 0.041\\
8 imps & MAR & OVER & MI & -0.002 & -0.024 & 0.019\\
8 imps & MAR & UNDR & CCA & -0.045 & -0.071 & -0.018\\
8 imps & MAR & UNDR & MI & 0.001 & -0.018 & 0.021\\
8 imps & MAR & No Missing Data & COMP & 0.000 & 0.000 & 0.000\\
\addlinespace
8 imps & MNAR & OVER & CCA & 0.033 & 0.016 & 0.054\\
8 imps & MNAR & OVER & MI & 0.030 & 0.012 & 0.051\\
8 imps & MNAR & UNDR & CCA & -0.053 & -0.086 & -0.019\\
8 imps & MNAR & UNDR & MI & -0.049 & -0.084 & -0.017\\
8 imps & MNAR & No Missing Data & COMP & 0.000 & 0.000 & 0.000\\
\bottomrule
\end{tabular}
:::

::: {.cell-output-display}
\begin{tabular}[t]{llrrr}
\toprule
Miss Type & Direction & Prop. of Inc. Cases & 0.025 Quantile & 0.975 Quantile\\
\midrule
MAR & OVER & 0.009 & 0.009 & 0.01\\
MAR & UNDR & 0.009 & 0.009 & 0.01\\
MNAR & OVER & 0.010 & 0.010 & 0.01\\
MNAR & UNDR & 0.010 & 0.010 & 0.01\\
\bottomrule
\end{tabular}
:::

::: {.cell-output-display}
![](simualtion_results_files/figure-pdf/n = 2500 Results Tables-2.pdf)
:::

::: {.cell-output-display}
\begin{tabular}[t]{llllrrr}
\toprule
Num of Imp & Miss Type & Direction & Analysis & Bias & 0.025 Quantile & 0.975 Quantile\\
\midrule
3 imps & MAR & OVER & CCA & 0.027 & -0.216 & 0.327\\
3 imps & MAR & OVER & MI & 0.011 & -0.231 & 0.310\\
3 imps & MAR & UNDR & CCA & -0.031 & -0.268 & 0.249\\
3 imps & MAR & UNDR & MI & 0.015 & -0.229 & 0.311\\
3 imps & MAR & No Missing Data & COMP & 0.013 & -0.226 & 0.306\\
\addlinespace
3 imps & MNAR & OVER & CCA & 0.056 & -0.200 & 0.362\\
3 imps & MNAR & OVER & MI & 0.052 & -0.204 & 0.350\\
3 imps & MNAR & UNDR & CCA & -0.030 & -0.278 & 0.262\\
3 imps & MNAR & UNDR & MI & -0.026 & -0.270 & 0.271\\
3 imps & MNAR & No Missing Data & COMP & 0.022 & -0.227 & 0.320\\
\addlinespace
5 imps & MAR & OVER & CCA & 0.024 & -0.233 & 0.323\\
5 imps & MAR & OVER & MI & 0.009 & -0.243 & 0.303\\
5 imps & MAR & UNDR & CCA & -0.034 & -0.282 & 0.254\\
5 imps & MAR & UNDR & MI & 0.012 & -0.246 & 0.303\\
5 imps & MAR & No Missing Data & COMP & 0.010 & -0.237 & 0.298\\
\addlinespace
5 imps & MNAR & OVER & CCA & 0.050 & -0.210 & 0.345\\
5 imps & MNAR & OVER & MI & 0.047 & -0.214 & 0.337\\
5 imps & MNAR & UNDR & CCA & -0.035 & -0.286 & 0.247\\
5 imps & MNAR & UNDR & MI & -0.031 & -0.283 & 0.258\\
5 imps & MNAR & No Missing Data & COMP & 0.017 & -0.238 & 0.317\\
\addlinespace
8 imps & MAR & OVER & CCA & 0.022 & -0.254 & 0.327\\
8 imps & MAR & OVER & MI & 0.006 & -0.257 & 0.307\\
8 imps & MAR & UNDR & CCA & -0.037 & -0.296 & 0.247\\
8 imps & MAR & UNDR & MI & 0.008 & -0.256 & 0.308\\
8 imps & MAR & No Missing Data & COMP & 0.007 & -0.265 & 0.303\\
\addlinespace
8 imps & MNAR & OVER & CCA & 0.059 & -0.207 & 0.397\\
8 imps & MNAR & OVER & MI & 0.056 & -0.212 & 0.391\\
8 imps & MNAR & UNDR & CCA & -0.027 & -0.281 & 0.282\\
8 imps & MNAR & UNDR & MI & -0.023 & -0.277 & 0.282\\
8 imps & MNAR & No Missing Data & COMP & 0.026 & -0.233 & 0.356\\
\bottomrule
\end{tabular}
:::
:::


## MAR Overall Results Table


::: {.cell hash='simualtion_results_cache/pdf/unnamed-chunk-5_919c9008e6d66bb8eab801d8682d88ee'}

```{.r .cell-code}
sim_settings <- readr::read_csv("../sim_settings.csv")
```

::: {.cell-output .cell-output-stderr}
```
Rows: 36 Columns: 3
-- Column specification --------------------------------------------------------
Delimiter: ","
dbl (3): n, M, simQ

i Use `spec()` to retrieve the full column specification for this data.
i Specify the column types or set `show_col_types = FALSE` to quiet this message.
```
:::

```{.r .cell-code}
sim_settings[sim_settings$n < 5000, "simQ"] <- 1000

sim_settings <- sim_settings[-c(22, 23, 24, 35, 36),]

sim_res_list <- list()

for (i in 1:nrow(sim_settings)) {
    Q <- sim_settings[i, "simQ"]
    n <- sim_settings[i, "n"]
    m <- sim_settings[i, "M"]
    fname <- paste0("../simulation_results_Q", Q, "_n_", n, "_m_", m, "_p_miss_0.01.csv")
    tmp <- read.csv(fname, header = T)
    tmp$IMPUTATION <- paste0(m, " imps")
    tmp$N <- as.numeric(n)
    sim_res_list[[i]] <- tmp
    
}

sim_res_list %>% purrr::reduce(bind_rows) -> sim_res

sim_res[which(sim_res$DIRECTION == "COMP"), "ANALYSIS"] <- "COMP"
sim_res[which(sim_res$DIRECTION == "COMP"), "DIRECTION"] <- "No Missing Data"
sim_res$DIRECTION %<>% fct_relevel(c("OVER", "UNDR", "No Missing Data"))
sim_res$ANALYSIS %<>% fct_relevel(c("CCA", "COMP", "MI"))
sim_res$IMPUTATION %<>% fct_relevel(paste0(c(3, 5, 8, 25, 50, 100), " imps"))
```
:::

::: {.cell hash='simualtion_results_cache/pdf/unnamed-chunk-6_395b117ca45438e281a920d0dd36eb5f'}

```{.r .cell-code}
options(dplyr.summarise.inform = FALSE)
sim_res |> 
  filter(IMPUTATION == "3 imps", 
         MISS_TYPE == "MAR") |>
  group_by(N, DIRECTION, ANALYSIS) |>
  summarize(
    OR_DIFF = mean(OR_DIFF),
    OR_BIAS = mean(OR_BIAS),
    PROP_MISS = mean(PROP_MISS)
  ) |> 
  kableExtra::kbl(format = "latex", 
                  digits = 3, 
                  booktabs = TRUE,
                  caption = "MAR simulation results with 3 imputations for MI.")
```

::: {.cell-output-display}
\begin{table}

\caption{\label{tab:unnamed-chunk-6}MAR simulation results with 3 imputations for MI.}
\centering
\begin{tabular}[t]{rllrrr}
\toprule
N & DIRECTION & ANALYSIS & OR\_DIFF & OR\_BIAS & PROP\_MISS\\
\midrule
500 & OVER & CCA & 0.015 & 0.155 & 0.007\\
500 & OVER & MI & 0.000 & 0.140 & 0.007\\
500 & UNDR & CCA & -0.045 & 0.095 & 0.008\\
500 & UNDR & MI & 0.002 & 0.142 & 0.008\\
500 & No Missing Data & COMP & 0.000 & 0.140 & 0.000\\
\addlinespace
1000 & OVER & CCA & 0.014 & 0.057 & 0.008\\
1000 & OVER & MI & -0.001 & 0.043 & 0.008\\
1000 & UNDR & CCA & -0.049 & -0.006 & 0.009\\
1000 & UNDR & MI & 0.001 & 0.045 & 0.009\\
1000 & No Missing Data & COMP & 0.000 & 0.044 & 0.000\\
\addlinespace
2500 & OVER & CCA & 0.014 & 0.027 & 0.009\\
2500 & OVER & MI & -0.002 & 0.011 & 0.009\\
2500 & UNDR & CCA & -0.044 & -0.031 & 0.009\\
2500 & UNDR & MI & 0.002 & 0.015 & 0.009\\
2500 & No Missing Data & COMP & 0.000 & 0.013 & 0.000\\
\addlinespace
5000 & OVER & CCA & 0.014 & 0.035 & 0.010\\
5000 & OVER & MI & -0.002 & 0.019 & 0.010\\
5000 & UNDR & CCA & -0.046 & -0.025 & 0.010\\
5000 & UNDR & MI & 0.001 & 0.022 & 0.010\\
5000 & No Missing Data & COMP & 0.000 & 0.021 & 0.000\\
\addlinespace
10000 & OVER & CCA & 0.015 & 0.013 & 0.010\\
10000 & OVER & MI & -0.001 & -0.003 & 0.010\\
10000 & UNDR & CCA & -0.047 & -0.048 & 0.010\\
10000 & UNDR & MI & 0.002 & 0.000 & 0.010\\
10000 & No Missing Data & COMP & 0.000 & -0.002 & 0.000\\
\addlinespace
25000 & OVER & CCA & 0.015 & 0.019 & 0.010\\
25000 & OVER & MI & -0.001 & 0.003 & 0.010\\
25000 & UNDR & CCA & -0.047 & -0.043 & 0.010\\
25000 & UNDR & MI & 0.002 & 0.006 & 0.010\\
25000 & No Missing Data & COMP & 0.000 & 0.004 & 0.000\\
\addlinespace
50000 & OVER & CCA & 0.014 & 0.014 & 0.010\\
50000 & OVER & MI & -0.002 & -0.002 & 0.010\\
50000 & UNDR & CCA & -0.046 & -0.046 & 0.010\\
50000 & UNDR & MI & 0.002 & 0.002 & 0.010\\
50000 & No Missing Data & COMP & 0.000 & 0.000 & 0.000\\
\bottomrule
\end{tabular}
\end{table}
:::

```{.r .cell-code}
sim_res |> 
  filter(IMPUTATION == "5 imps", 
         MISS_TYPE == "MAR") |>
  group_by(N, DIRECTION, ANALYSIS) |>
  summarize(
    OR_DIFF = mean(OR_DIFF),
    OR_BIAS = mean(OR_BIAS),
    PROP_MISS = mean(PROP_MISS)
  ) |> 
  kableExtra::kbl(format = "latex", 
                  digits = 3, 
                  booktabs = TRUE,
                  caption = "MAR simulation results with 5 imputations for MI.")
```

::: {.cell-output-display}
\begin{table}

\caption{\label{tab:unnamed-chunk-6}MAR simulation results with 5 imputations for MI.}
\centering
\begin{tabular}[t]{rllrrr}
\toprule
N & DIRECTION & ANALYSIS & OR\_DIFF & OR\_BIAS & PROP\_MISS\\
\midrule
500 & OVER & CCA & 0.014 & 0.147 & 0.007\\
500 & OVER & MI & 0.001 & 0.133 & 0.007\\
500 & UNDR & CCA & -0.043 & 0.089 & 0.008\\
500 & UNDR & MI & 0.001 & 0.134 & 0.008\\
500 & No Missing Data & COMP & 0.000 & 0.133 & 0.000\\
\addlinespace
1000 & OVER & CCA & 0.014 & 0.066 & 0.008\\
1000 & OVER & MI & -0.001 & 0.052 & 0.008\\
1000 & UNDR & CCA & -0.049 & 0.003 & 0.009\\
1000 & UNDR & MI & 0.001 & 0.054 & 0.009\\
1000 & No Missing Data & COMP & 0.000 & 0.053 & 0.000\\
\addlinespace
2500 & OVER & CCA & 0.014 & 0.024 & 0.009\\
2500 & OVER & MI & -0.001 & 0.009 & 0.009\\
2500 & UNDR & CCA & -0.044 & -0.034 & 0.009\\
2500 & UNDR & MI & 0.001 & 0.012 & 0.009\\
2500 & No Missing Data & COMP & 0.000 & 0.010 & 0.000\\
\addlinespace
5000 & OVER & CCA & 0.014 & 0.021 & 0.010\\
5000 & OVER & MI & -0.002 & 0.005 & 0.010\\
5000 & UNDR & CCA & -0.047 & -0.040 & 0.010\\
5000 & UNDR & MI & 0.001 & 0.008 & 0.010\\
5000 & No Missing Data & COMP & 0.000 & 0.007 & 0.000\\
\addlinespace
10000 & OVER & CCA & 0.014 & 0.023 & 0.010\\
10000 & OVER & MI & -0.002 & 0.007 & 0.010\\
10000 & UNDR & CCA & -0.046 & -0.037 & 0.010\\
10000 & UNDR & MI & 0.002 & 0.011 & 0.010\\
10000 & No Missing Data & COMP & 0.000 & 0.009 & 0.000\\
\addlinespace
25000 & OVER & CCA & 0.014 & 0.013 & 0.010\\
25000 & OVER & MI & -0.002 & -0.003 & 0.010\\
25000 & UNDR & CCA & -0.047 & -0.048 & 0.010\\
25000 & UNDR & MI & 0.001 & 0.000 & 0.010\\
25000 & No Missing Data & COMP & 0.000 & -0.001 & 0.000\\
\addlinespace
50000 & OVER & CCA & 0.014 & 0.017 & 0.010\\
50000 & OVER & MI & -0.002 & 0.002 & 0.010\\
50000 & UNDR & CCA & -0.047 & -0.044 & 0.010\\
50000 & UNDR & MI & 0.001 & 0.005 & 0.010\\
50000 & No Missing Data & COMP & 0.000 & 0.003 & 0.000\\
\bottomrule
\end{tabular}
\end{table}
:::

```{.r .cell-code}
sim_res |> 
  filter(IMPUTATION == "8 imps", 
         MISS_TYPE == "MAR") |>
  group_by(N, DIRECTION, ANALYSIS) |>
  summarize(
    OR_DIFF = mean(OR_DIFF),
    OR_BIAS = mean(OR_BIAS),
    PROP_MISS = mean(PROP_MISS)
  ) |> 
  kableExtra::kbl(format = "latex", 
                  digits = 3, 
                  booktabs = TRUE,
                  caption = "MAR simulation results with 8 imputations for MI.")
```

::: {.cell-output-display}
\begin{table}

\caption{\label{tab:unnamed-chunk-6}MAR simulation results with 8 imputations for MI.}
\centering
\begin{tabular}[t]{rllrrr}
\toprule
N & DIRECTION & ANALYSIS & OR\_DIFF & OR\_BIAS & PROP\_MISS\\
\midrule
500 & OVER & CCA & 0.015 & 0.133 & 0.007\\
500 & OVER & MI & 0.000 & 0.118 & 0.007\\
500 & UNDR & CCA & -0.043 & 0.074 & 0.008\\
500 & UNDR & MI & 0.002 & 0.120 & 0.008\\
500 & No Missing Data & COMP & 0.000 & 0.118 & 0.000\\
\addlinespace
1000 & OVER & CCA & 0.013 & 0.057 & 0.008\\
1000 & OVER & MI & -0.001 & 0.042 & 0.008\\
1000 & UNDR & CCA & -0.049 & -0.005 & 0.009\\
1000 & UNDR & MI & 0.002 & 0.046 & 0.009\\
1000 & No Missing Data & COMP & 0.000 & 0.044 & 0.000\\
\addlinespace
2500 & OVER & CCA & 0.014 & 0.022 & 0.009\\
2500 & OVER & MI & -0.002 & 0.006 & 0.009\\
2500 & UNDR & CCA & -0.045 & -0.037 & 0.009\\
2500 & UNDR & MI & 0.001 & 0.008 & 0.009\\
2500 & No Missing Data & COMP & 0.000 & 0.007 & 0.000\\
\addlinespace
5000 & OVER & CCA & 0.014 & 0.023 & 0.010\\
5000 & OVER & MI & -0.002 & 0.008 & 0.010\\
5000 & UNDR & CCA & -0.046 & -0.037 & 0.010\\
5000 & UNDR & MI & 0.001 & 0.011 & 0.010\\
5000 & No Missing Data & COMP & 0.000 & 0.010 & 0.000\\
\addlinespace
10000 & OVER & CCA & 0.013 & 0.010 & 0.010\\
10000 & OVER & MI & -0.003 & -0.006 & 0.010\\
10000 & UNDR & CCA & -0.046 & -0.049 & 0.010\\
10000 & UNDR & MI & 0.002 & -0.001 & 0.010\\
10000 & No Missing Data & COMP & 0.000 & -0.003 & 0.000\\
\addlinespace
25000 & OVER & CCA & 0.014 & 0.013 & 0.010\\
25000 & OVER & MI & -0.002 & -0.003 & 0.010\\
25000 & UNDR & CCA & -0.047 & -0.048 & 0.010\\
25000 & UNDR & MI & 0.001 & 0.000 & 0.010\\
25000 & No Missing Data & COMP & 0.000 & -0.001 & 0.000\\
\bottomrule
\end{tabular}
\end{table}
:::
:::

::: {.cell hash='simualtion_results_cache/pdf/fig-mar-diff-pmiss-001_9e0bd433c6dae1444cec8350aa54e96c'}

```{.r .cell-code}
sim_res |> 
  filter(MISS_TYPE == "MAR",
         ANALYSIS != "COMP") |>
  mutate(
    ANALYSIS = case_when(
      ANALYSIS == "MI" ~ paste0("MI ", IMPUTATION),
      ANALYSIS == "CCA" ~ "CCA"
    )
  ) -> sim_res_sum
sim_res_sum$ANALYSIS %<>% fct_relevel(c("CCA", paste0("MI ", c(3, 5, 8, 25, 50, 100), " imps")))
sim_res_sum |>
  ggplot(aes(DIRECTION, OR_DIFF, color = as.factor(N))) +
        geom_hline(yintercept = 0, color = "gray") +
        geom_boxplot(fill = "white", alpha = .75, outlier.shape = 1) +
        facet_wrap(ANALYSIS ~ .) +
        theme_bw() +
        ggthemes::scale_fill_colorblind() +
        ggthemes::scale_color_colorblind(name = 'Sample Size') +
        theme(legend.position = "right") +
        labs(title = "Difference in Race Effect Estimates \n from CCA/MI to the Complete Data Estimate",
             x = "Intended Direction of Bias due to Missing Data",
             y = "Difference in the OR")
```

::: {.cell-output-display}
![Distributions of the difference between the complete data estimate and the incomplete data estimates under complete case analysis (CCA), or multiple imputation (MI) with 3, 5, 8, 25, 50, and 100 imputations (imps) from simulations at sample sizes from N = 500 to N = 50000.](simualtion_results_files/figure-pdf/fig-mar-diff-pmiss-001-1.pdf){#fig-mar-diff-pmiss-001 fig-pos='H'}
:::
:::

::: {.cell hash='simualtion_results_cache/pdf/fig-mar-bias-pmiss-001_a91f6162687731f0880240600a82b552'}

```{.r .cell-code}
sim_res |> 
  filter(MISS_TYPE == "MAR") |>
  mutate(
    ANALYSIS = case_when(
      ANALYSIS == "MI" ~ paste0("MI ", IMPUTATION),
      ANALYSIS == "CCA" ~ "CCA",
      ANALYSIS == "COMP" ~ "COMP"
    )
  ) -> sim_res_sum
sim_res_sum$ANALYSIS %<>% fct_relevel(c("COMP", "CCA", paste0("MI ", c(3, 5, 8, 25, 50, 100), " imps")))
sim_res_sum |>
  ggplot(aes(DIRECTION, OR_BIAS, color = as.factor(N))) +
        geom_hline(yintercept = 0, color = "gray") +
        geom_boxplot(fill = "white", alpha = .75, outlier.shape = 1) +
        facet_wrap(. ~ ANALYSIS, scales = "free_x") +
        theme_bw() +
        ggthemes::scale_fill_colorblind() +
        ggthemes::scale_color_colorblind(name = 'Sample Size') +
        theme(legend.position = "right") +
        labs(title = "Statistical Bias of the Race Effect Estimates", 
             x = "Intended Direction of Bias due to Missing Data",
             y = "Statistical Bias")
```

::: {.cell-output-display}
![Distributions of the statistical bias between the true race effect parameter and the estimates from the complete data (COMP), complete case analysis (CCA) estimate, and multiple imputation (MI) with 3, 5, 8, 25, 50, and 100 imputations (imps) from simulations at sample sizes from N = 500 to N = 50000.](simualtion_results_files/figure-pdf/fig-mar-bias-pmiss-001-1.pdf){#fig-mar-bias-pmiss-001 fig-pos='H'}
:::
:::
