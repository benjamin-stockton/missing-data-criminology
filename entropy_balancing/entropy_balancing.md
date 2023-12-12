---
title: "Entropy Balancing"
author: "Benjamin Stockton"
number-sections: true 
df-print: paged
pdf-engine: pdflatex
keep-tex: true
keep-md: true
toc: true
tbl-cap-location: bottom
format: 
  pdf:
    geometry: 
      - top=30mm
      - left=30mm
    shift-heading-level-by: -1
  docx:
    shift-heading-level-by: -1
editor: visual
execute: 
  cache: true
  keep-md: true
  keep-tex: true
bibliography: ../Literature/Criminology.bib
csl: '../Drafts/JQC_Manuscript/journal-of-quantitative-criminology.csl'
---



## Entropy Balancing

To make causal inferences in observational data, the observed Y(1)\|D=1 and counterfactual Y(0)\|D=1 are compared to get the Population Average Treatment Effect on Treated (PATT) defined as $\tau = E(Y(1)|D=1) - E(Y(0) | D=1)$. In experimental studies where treatment assignment is independent of the potential outcomes, approximate the second expectation by $E(Y(0) | D=0)$ i.e. the mean of the control group. In observational studies, if we can "assume ignorable treatment assignment and overlap", then $Y(0) \perp D|X$ and $P(D = 1| X=x) \leq 1$ for all $x$ in the support of $f_{X|D=1}$. That means that if the confounding covariates are similar on both the treatment and control groups, we can estimate the counterfactual average outcome by $\tau = E(Y|D=1) = \int E(Y|X=x, D=0) f_{X|D=1}(x) dx.$

To make causal inferences in observational data, the observed $Y(1)|D=1$ and counterfactual $Y(0)|D=1$ are compared to get the Population Average Treatment Effect on Treated (PATT) defined as $\tau = E(Y(1)|D=1) - E(Y(0) | D=1)$. In experimental studies where treatment assignment is independent of the potential outcomes, approximate the second expectation by $E(Y(0) | D=0)$ i.e. the mean of the control group. In observational studies, if we can "assume ignorable treatment assignment and overlap", then $Y(0) \perp D|X$ and $P(D = 1| X=x) \leq 1$ for all $x$ in the support of $f_{X|D=1}$. That means that if the confounding covariates are similar on both the treatment and control groups, we can estimate the counterfactual average outcome by $\tau = E(Y|D=1) = \int E(Y|X=x, D=0) f_{X|D=1}(x) dx.$

> Notice that the last term in this expression is equal to the covariate adjusted mean, that is, the estimated mean of $Y$ in the source population if its covariates were distributed as in the target population.

-   p\. 28 [@hainmueller2012]

> \[Rosenbaum and Rubin (1983)\] showed that the multivariate matching pre-processing problem\] can be reduced to a single dimension if the counterfactual mean can be identified as $E(Y(0)|D=1) = \int E(Y|p(X) = \rho, D=0) f_{p|D=1}(\rho)d\rho$ where $f_{p|D=1}$ is the dist. of the propensity score $p(x) = P(D=1 | X=x)$ in the target population.

-   p\. 28 [@hainmueller2012]

Propensity score weighting is performed by using a binary response (logit/probit) regression to estimate a probability $p_i$ to be in the treatment given the covariates. These are converted to weights $d_i$ by the inverse of the link function $d_i = p_i / (1-p_i)$ so that the counterfactual mean is estimated as $\widehat{E(Y(0)|D=1)} = \frac{\sum_{i|D=0}Y_i d_i}{\sum_{i|D=0}d_i}.$

Drawbacks of PSW:

-   The true propensity score is valuable because it balances the covariate distributions of the two groups, but is unknown and difficult to accurately estimate.

-   Mis-specified scores can lead to biased estimation of TE

-   Practitioners often iterate between weighting and matching, modeling the PS, and then evaluating the balance until a suitable balance is achieved. Imai, King and Stuart (2008) call this the "propensity score tautology".

    -   Despite this balance often isn't achieved and can make balance worse among confounders.

### Entropy Balancing Procedure:

**Goal:** Estimate $\tau = E(Y(1) |D=1) - E(Y(0) | D=1).$

The counterfactual mean could be estimated by $\widehat{E(Y(0)|D=1)}=\frac{\sum_{i|D=0} Y_i w_i}{\sum_{i|D=0} w_i}.$ The weights are found by minimizing $H(w) = \sum_{i|D=0} h(w_i)$ subject to $\sum_{i|D=0} w_i c_{ri}(X_i) = m_r$ for $r =1,…,R$ and $\sum_{i|D=0} w_i = 1$ and $w_i \geq 0$ for all $i$ given $D_i = 0$ where $h(.)$ is a distance metric and $c_{ri}(X_i)=m_r$ describes the set of $R$ balance constraints imposed by the covariate moments of the re-weighted control group.

Authors choose to use $h(w_i)=w_i\log(w_i/q_i)$ the Kullback entropy divergence.

In conventional PSW, the researcher (1) estimates $d_i$ then (2) checks if the weights balance the covariate distribution. Entropy balancing reverses the approach by obtaining weights by minimizing a linear equation with constraints that guarantee balance while remaining close to the uniform weights to guarantee efficiency.

Optimization is performed using Lagrangian Multipliers.

#### 3 Main Issues with Entropy Balancing:

1.  "No weighting solution exists if the balance constraints are inconsistent." Easy to avoid as constraints are researcher imposed.
2.  "Balance constraints are consistent but there exists no set of positive weights to actually satisfy the constraints." For example, there's heavy imbalance in a binary categorical variable between the treatment and control groups. "If there aren't enough controls that look anything like the treated units, then the existing data do not contain sufficient information to reliably infer the counterfactual."
3.  "A solution exists, but due to limited overlap, the solution involves an extreme adjustment to the weights of some control units." Only a few units receive relatively large weights and all others are set near 0, then the variance will increase and effective sample size is small.

Issues 2 and 3 are also relevant to most pre-processing that involves balancing including propensity score weighting.

### Entropy Balancing in Sentencing Research

These notes are largely based on and quotions from [@macdonaldEvaluatingRoleRace2019]. The authors note that a cetnral goal of contemporary criminal justice reform is to "reduce racial disparities in prisons" in part by building an "understanding of the sources of these inequalities." Regression has been the primary tool used to estimate the race effect (note that the authors make causal claims without specifying the use of causal methods). MacDonald and Donnelly motivate entropy balancing by noting that

> "\[M\]ultivariate regression approaches to estimating the effect of race on sentencing may produce biased estimates if there are important subgroup differences across covariates that are not adequately removed through mean adjustment."

Given that is the extent of the motivation and the lack of theoretical basis for entropy balancing beyond making the covariate distributions identical on the moments, this doesn't seem to be the best path forward for causal methods in criminology or other social sciences. MacDonald and Donnelly engage with propensity score matching and weighting, but dismiss it for not producing well-balanced covariate data sets in every case. The balancing seems to me to be a side-effect of constructing suitable counterfactuals, and not the direct goal of propensity score methods. In other words, propensity scores (and resulting weights) have the theoretical meaning of being the likelihood of belonging to the treatment group given other characterisitcs. Entropy weights have no such meaning because they work backward from the assumption that the moments are matched.

The authors also claim entropy balancing as doubly robust, while other papers have called the doubly robust claims into question [@freedman2008] (notably the co-author, Berk, is a criminology/statistics professor who has been critical of the way statistical methods have been used in criminology in other papers as well [@berkWhatYouCan2010]). MacDonald and Donnelly seem to define doubly robust estimation as "assum\[ing\] consistent measurement of treatment effects if either the weights from the **balancing approaches** or **regression model** are correctly specified [@bang2005]." [@macdonaldEvaluatingRoleRace2019, p. 662].[^1] This does not seem to be the correct definition.

[^1]: When actually checking Bang and Robins' @bang2005 definition they say "In a missing data model, an estimator is DR if it remains consistent when either (but not necessarily both) a model for the missingness mechanism or a model for the distribution of the complete data is correctly specified." which is consistent with Funk et al's definition @funk2011 . These are markedly different definitions in that Bang and Robins are discussing the model specification as being the source of biasedness while MacDonald and Donnelly seem to be paraphrasing the definition as caring only about whether the estimates are biased regardless of the cause.

The central claim by the authors is that

> "Although the coefficients produced from the linear regression and the other approaches are not statistically different from each other (e.g. all coefficients are within 1 standard deviation of each other), the change in the size of the estimated mean sentence length difference between Blacks and Whites suggests that the regression model is downwardly biased. Specifically, race coefficients have similar precision, but the parameter estimated from the linear regression is one-half the size of the coefficients estimated from entropy or propensity score weighting." [@macdonaldEvaluatingRoleRace2019, p. 675]

But there's no justification for any of these claims. No engagement with the uncertainty measures; in fact saying they're all within one standard error of each other should be the end of the conversation. That means there are no real distinctions between results so trying to draw any is pointless.

A central issue with entropy weights is the lack of a clear causal estimand resulting from the weighting. IPW results in ATE (total population), "treated" results in ATT (target population is only the treated individuals), and "overlap" weights result in ATO (overlap of treatment and control in population) [@zhou2022]. Generally, entropy weighting will also make comparisons on the portion of the total population where there's significant overlap between the treatment and control groups, but without theoretical optimality properties.

## Implemented by the ebal Package

### Simulated Data

First, we'll demonstrate entropy balancing in a logistic regression setting to estimate the race effect on the simulated data. This data contains the same number of observations as the PCS data on the same variables which have been simulated independently, so a priori we should expect the marginal means of the variables to be the same for "Black" and "Non-Black" simulated defendants. The simulated data set is complete.


::: {.cell hash='entropy_balancing_cache/pdf/set-up_1f1507fabc9f179afc7e26b46f1ad572'}

```{.r .cell-code}
library(ebal)
```

::: {.cell-output .cell-output-stderr}
```
##
## ebal Package: Implements Entropy Balancing.
```
:::

::: {.cell-output .cell-output-stderr}
```
## See http://www.stanford.edu/~jhain/ for additional information.
```
:::

```{.r .cell-code}
library(dplyr)
```

::: {.cell-output .cell-output-stderr}
```

Attaching package: 'dplyr'
```
:::

::: {.cell-output .cell-output-stderr}
```
The following objects are masked from 'package:stats':

    filter, lag
```
:::

::: {.cell-output .cell-output-stderr}
```
The following objects are masked from 'package:base':

    intersect, setdiff, setequal, union
```
:::

```{.r .cell-code}
dat <- read.csv("../Data/simulated_data.csv")

dat$YEAR <- factor(dat$YEAR)
dat$COUNTY <- factor(dat$COUNTY)

dat$OFF_RACER <- factor(dat$OFF_RACER, levels = c("WHITE", "BLACK", "LATINO", "OTHER"))
```
:::


Entropy balancing seems to be performed very quickly. Next we check the marginal means for each of the covariates.


::: {.cell hash='entropy_balancing_cache/pdf/ebal-sim_045d105dbdd87cfe893b42ad2a5ba4c6'}

```{.r .cell-code}
dat %>%
    select(-c(INCAR, YEAR, COUNTY, OFF_RACER)) -> X
X <- model.matrix(~., data = X)[,-1]

treatment <- ifelse(dat$OFF_RACER == "BLACK", TRUE, FALSE)

eb.out <- ebalance(Treatment = treatment, 
                   X = X, print.level = 2)
```

::: {.cell-output .cell-output-stdout}
```
Iteration 1 maximum deviation is = 211656 
Iteration 2 maximum deviation is = 9074 
Converged within tolerance 
Converged within tolerance 
```
:::

```{.r .cell-code}
c1 <- apply(X[treatment,], 2, mean)


c2 <- apply(X[!treatment,], 2, weighted.mean, w = eb.out$w)


c3 <- apply(X[!treatment,], 2, mean)

X.means <- bind_rows(c1, c2, c3) %>% as.data.frame()
rownames(X.means) <- c("Black", "Non-Black - EB", "Non-Black - Unbalanced")
round(t(X.means), 4) |> 
    kableExtra::kbl(format = "latex",
                    booktabs = TRUE,
                    digits = 3) |> 
    print()
```

::: {.cell-output .cell-output-stdout}
```

\begin{tabular}[t]{lrrr}
\toprule
  & Black & Non-Black - EB & Non-Black - Unbalanced\\
\midrule
CRIMETYPEDUI & 0.237 & 0.237 & 0.238\\
CRIMETYPEOTHER & 0.124 & 0.124 & 0.123\\
CRIMETYPEPERSONS & 0.156 & 0.156 & 0.157\\
CRIMETYPEPROPERTY & 0.256 & 0.256 & 0.255\\
OGS & 3.428 & 3.428 & 3.418\\
\addlinespace
TRIALYes & 0.978 & 0.978 & 0.979\\
MALEMale & 0.773 & 0.773 & 0.774\\
PRS4/5 & 0.181 & 0.181 & 0.182\\
PRSNone & 0.486 & 0.486 & 0.486\\
PRSREVOC/RFEL & 0.030 & 0.030 & 0.030\\
\addlinespace
DOSAGE & 34.267 & 34.267 & 34.276\\
RECMINYes & 0.670 & 0.670 & 0.671\\
OGSQ & 17.944 & 17.944 & 17.850\\
DOSAGEQ & 1306.290 & 1306.290 & 1307.200\\
\bottomrule
\end{tabular}
```
:::
:::


The means are essentially the same between the covariates before balancing, so it won't have any effect on the actual analysis. This is a key point regarding use of the modal approach in the simulations.


::: {.cell hash='entropy_balancing_cache/pdf/ebal-sim-weights_860ab475ef8cf0b1be29cf82c3d5f59b'}

```{.r .cell-code}
c(summary(eb.out$w), "Std. Dev" = sd(eb.out$w))
```

::: {.cell-output .cell-output-stdout}
```
       Min.     1st Qu.      Median        Mean     3rd Qu.        Max. 
0.353880254 0.366262566 0.368087992 0.368239110 0.369965186 0.385803692 
   Std. Dev 
0.002779261 
```
:::
:::


The weights are essentially uniform still indicating little needs to be done to balance the covariates.


::: {.cell hash='entropy_balancing_cache/pdf/ebal-sim-analysis_7b7047ba5bf0d753141834485763f29b'}

```{.r .cell-code}
weights <- numeric(nrow(dat))
for (i in 1:nrow(dat)) {
    if (treatment[i]) {
        weights[i] <- 1
    }
    else {
        weights[i] <- eb.out$w[i]
    }
}

ptime <- system.time({
    fit_sim_unweighted <- summary(glm(INCAR ~ ., 
                              data = dat, 
                              family = binomial(link = "logit")))
    
    fit_sim_weighted <- summary(glm(INCAR ~ ., 
                              data = dat, 
                              weights = weights,
                              family = binomial(link = "logit")))
})[3]
```

::: {.cell-output .cell-output-stderr}
```
Warning in eval(family$initialize): non-integer #successes in a binomial glm!
```
:::

```{.r .cell-code}
print(ptime)
```

::: {.cell-output .cell-output-stdout}
```
elapsed 
 41.589 
```
:::

```{.r .cell-code}
print("Unweighted Race Effect Estimate:")
```

::: {.cell-output .cell-output-stdout}
```
[1] "Unweighted Race Effect Estimate:"
```
:::

```{.r .cell-code}
fit_sim_unweighted$coefficients["OFF_RACERBLACK",]
```

::: {.cell-output .cell-output-stdout}
```
    Estimate   Std. Error      z value     Pr(>|z|) 
 0.225326688  0.005659322 39.815138611  0.000000000 
```
:::

```{.r .cell-code}
print("Weighted Race Effect Estimate:")
```

::: {.cell-output .cell-output-stdout}
```
[1] "Weighted Race Effect Estimate:"
```
:::

```{.r .cell-code}
fit_sim_weighted$coefficients["OFF_RACERBLACK",]
```

::: {.cell-output .cell-output-stdout}
```
     Estimate    Std. Error       z value      Pr(>|z|) 
 2.246423e-01  7.484070e-03  3.001606e+01 6.057134e-198 
```
:::
:::


### Real PCS Data


::: {.cell hash='entropy_balancing_cache/pdf/PCS-set-up_533865441e535412de83cb1f80e8dcd5'}

```{.r .cell-code}
pcs <- read.csv("../Data/most_serious_sentence_2010-2019_slim.csv")

set.seed(7822)
pcs %>%
    group_by(COUNTY) %>%
    # sample_n(size = 300) %>%
    ungroup() %>%
    mutate(
        YEAR = as.factor(YEAR),
        COUNTY = as.factor(COUNTY),
        OFF_RACER = factor(OFF_RACER, levels = c("WHITE", "BLACK", "LATINO", "OTHER"))
        ) -> pcs.sub

mice::md.pattern(pcs.sub, rotate.names = T)
```

::: {.cell-output-display}
![](entropy_balancing_files/figure-pdf/PCS-set-up-1.pdf){fig-pos='H'}
:::

::: {.cell-output .cell-output-stdout}
```
       INCAR CRIMETYPE OGS OGSQ TRIAL MALE COUNTY YEAR PRS DOSAGE DOSAGEQ
834546     1         1   1    1     1    1      1    1   1      1       1
26206      1         1   1    1     1    1      1    1   1      1       1
2076       1         1   1    1     1    1      1    1   1      1       1
34         1         1   1    1     1    1      1    1   1      1       1
1465       1         1   1    1     1    1      1    1   1      0       0
92         1         1   1    1     1    1      1    1   1      0       0
3          1         1   1    1     1    1      1    1   0      1       1
           0         0   0    0     0    0      0    0   3   1557    1557
       RECMIN OFF_RACER      
834546      1         1     0
26206       1         0     1
2076        0         1     1
34          0         0     2
1465        1         1     2
92          1         0     3
3           0         1     2
         2113     26332 31562
```
:::
:::

::: {.cell hash='entropy_balancing_cache/pdf/pcs-ebal_16a99a0ec641e44e703561d5402df6b3'}

```{.r .cell-code}
pcs.sub %>% 
    filter(!is.na(OFF_RACER),
           !is.na(DOSAGE),
           !is.na(RECMIN),
           !is.na(PRS)) -> pcs.cc
pcs.cc %>%
    select(-c(INCAR, YEAR, COUNTY, OFF_RACER)) -> X
X <- model.matrix(~., data = X)[,-1]

treatment <- ifelse(pcs.cc$OFF_RACER == "BLACK", TRUE, FALSE)

eb.pcs <- ebalance(Treatment = treatment, 
                   X = X, print.level = 2)
```

::: {.cell-output .cell-output-stdout}
```
Iteration 1 maximum deviation is = 28887568 
Iteration 2 maximum deviation is = 26506872 
Iteration 3 maximum deviation is = 23179050 
Iteration 4 maximum deviation is = 18192459 
Iteration 5 maximum deviation is = 8447031 
Iteration 6 maximum deviation is = 222757 
Iteration 7 maximum deviation is = 315.2 
Converged within tolerance 
Converged within tolerance 
```
:::

```{.r .cell-code}
c1 <- apply(X[treatment,], 2, mean)


c2 <- apply(X[!treatment,], 2, weighted.mean, w = eb.pcs$w)


c3 <- apply(X[!treatment,], 2, mean)

X.means <- bind_rows(c1, c2, c3) %>% as.data.frame()
rownames(X.means) <- c("Black", "Non-Black - EB", "Non-Black - Unbalanced")
round(t(X.means), 4) |> 
    kableExtra::kbl(format = "latex",
                    booktabs = TRUE,
                    digits = 3) |> 
    print()
```

::: {.cell-output .cell-output-stdout}
```

\begin{tabular}[t]{lrrr}
\toprule
  & Black & Non-Black - EB & Non-Black - Unbalanced\\
\midrule
CRIMETYPEDUI & 0.127 & 0.127 & 0.279\\
CRIMETYPEOTHER & 0.163 & 0.163 & 0.108\\
CRIMETYPEPERSONS & 0.184 & 0.184 & 0.146\\
CRIMETYPEPROPERTY & 0.244 & 0.244 & 0.259\\
OGS & 4.060 & 4.060 & 3.183\\
\addlinespace
OGSQ & 24.727 & 24.727 & 15.321\\
RECMINYes & 0.530 & 0.530 & 0.720\\
TRIALYes & 0.959 & 0.959 & 0.986\\
PRS4/5 & 0.272 & 0.272 & 0.152\\
PRSNone & 0.386 & 0.386 & 0.515\\
\addlinespace
PRSREVOC/RFEL & 0.051 & 0.051 & 0.023\\
MALEMale & 0.827 & 0.827 & 0.754\\
DOSAGE & 32.912 & 32.912 & 34.789\\
DOSAGEQ & 1213.208 & 1213.208 & 1342.119\\
\bottomrule
\end{tabular}
```
:::

```{.r .cell-code}
c(summary(eb.pcs$w), "Std. Dev" = sd(eb.pcs$w))
```

::: {.cell-output .cell-output-stdout}
```
       Min.     1st Qu.      Median        Mean     3rd Qu.        Max. 
 0.07609185  0.18832187  0.29698135  0.36708400  0.45848952 16.65560610 
   Std. Dev 
 0.28137031 
```
:::
:::

::: {.cell hash='entropy_balancing_cache/pdf/pcs-ebal-analysis_741366f401823cb44fed2bed6ef1fdf1'}

```{.r .cell-code}
weights <- numeric(nrow(pcs.cc))
for (i in 1:nrow(pcs.cc)) {
    if (treatment[i]) {
        weights[i] <- 1
    }
    else {
        weights[i] <- eb.pcs$w[i]
    }
}

fit_pcs_unweighted <- summary(glm(INCAR ~ ., 
                          data = pcs.cc, 
                          family = binomial(link = "logit"),
                          model = FALSE, 
                          y = FALSE))

fit_pcs_weighted <- summary(glm(INCAR ~ ., 
                          data = pcs.cc, 
                          weights = weights,
                          family = binomial(link = "logit"),
                          model = FALSE, 
                          y = FALSE))
```

::: {.cell-output .cell-output-stderr}
```
Warning in eval(family$initialize): non-integer #successes in a binomial glm!
```
:::

```{.r .cell-code}
print("Unweighted Race Effect Estimate:")
```

::: {.cell-output .cell-output-stdout}
```
[1] "Unweighted Race Effect Estimate:"
```
:::

```{.r .cell-code}
fit_pcs_unweighted$coefficients["OFF_RACERBLACK",]
```

::: {.cell-output .cell-output-stdout}
```
    Estimate   Std. Error      z value     Pr(>|z|) 
3.332389e-01 5.358802e-02 6.218533e+00 5.018254e-10 
```
:::

```{.r .cell-code}
print("Weighted Race Effect Estimate:")
```

::: {.cell-output .cell-output-stdout}
```
[1] "Weighted Race Effect Estimate:"
```
:::

```{.r .cell-code}
fit_pcs_weighted$coefficients["OFF_RACERBLACK",]
```

::: {.cell-output .cell-output-stdout}
```
    Estimate   Std. Error      z value     Pr(>|z|) 
0.3464378239 0.0799957463 4.3307030664 0.0000148634 
```
:::
:::


### Incomplete Data

The impact of incomplete data is the next topic to tackle. We know that CCA can be biased under several (untestable) assumptions for the missingness mechanisms for logistic regression. Need to show that missingness also has an impact on weights doubly impacting the inferential results.


::: {.cell hash='entropy_balancing_cache/pdf/missing-model_23f3b0d7245fb82c778a787cb4a52eac'}

```{.r .cell-code}
p <- 1 - (1 + exp(-(-9 + .5 * pcs.cc$DOSAGE + .5 * pcs.cc$OGS)))^-1
boxplot(p)
```

::: {.cell-output-display}
![](entropy_balancing_files/figure-pdf/missing-model-1.pdf){fig-pos='H'}
:::

```{.r .cell-code}
set.seed(12948)
mis <- sample(nrow(pcs.cc), size = ceiling(0.01 * nrow(pcs.cc)), prob = p)
str(mis)
```

::: {.cell-output .cell-output-stdout}
```
 int [1:194] 11771 6050 2603 15054 2490 11909 13617 2250 14735 12330 ...
```
:::
:::

::: {.cell hash='entropy_balancing_cache/pdf/pcs-inc-ebal_9b92200edb3bbd04250308ddf5593a7f'}

```{.r .cell-code}
pcs.inc <- pcs.cc
pcs.inc[mis, "OFF_RACER"] <- NA


pcs.inc %>% 
    filter(!is.na(OFF_RACER),
           !is.na(DOSAGE),
           !is.na(RECMIN),
           !is.na(PRS)) -> pcs.cc2
pcs.cc2 %>%
    select(-c(INCAR, YEAR, COUNTY, OFF_RACER)) -> X
X <- model.matrix(~., data = X)[,-1]

treatment <- ifelse(pcs.cc2$OFF_RACER == "BLACK", TRUE, FALSE)

eb.pcs2 <- ebalance(Treatment = treatment, 
                   X = X, print.level = 2)
```

::: {.cell-output .cell-output-stdout}
```
Iteration 1 maximum deviation is = 28904124 
Iteration 2 maximum deviation is = 26520140 
Iteration 3 maximum deviation is = 23186669 
Iteration 4 maximum deviation is = 18189093 
Iteration 5 maximum deviation is = 8411645 
Iteration 6 maximum deviation is = 221142 
Iteration 7 maximum deviation is = 310.9 
Converged within tolerance 
Converged within tolerance 
```
:::

```{.r .cell-code}
c1 <- apply(X[treatment,], 2, mean)


c2 <- apply(X[!treatment,], 2, weighted.mean, w = eb.pcs2$w)


c3 <- apply(X[!treatment,], 2, mean)

X.means <- bind_rows(c1, c2, c3) %>% as.data.frame()
rownames(X.means) <- c("Black", "Non-Black - EB", "Non-Black - Unbalanced")
round(t(X.means), 4) |> 
    kableExtra::kbl(format = "latex",
                    booktabs = TRUE,
                    digits = 3) |> 
    print()
```

::: {.cell-output .cell-output-stdout}
```

\begin{tabular}[t]{lrrr}
\toprule
  & Black & Non-Black - EB & Non-Black - Unbalanced\\
\midrule
CRIMETYPEDUI & 0.127 & 0.127 & 0.279\\
CRIMETYPEOTHER & 0.163 & 0.163 & 0.108\\
CRIMETYPEPERSONS & 0.184 & 0.184 & 0.146\\
CRIMETYPEPROPERTY & 0.244 & 0.244 & 0.259\\
OGS & 4.060 & 4.060 & 3.183\\
\addlinespace
OGSQ & 24.726 & 24.726 & 15.321\\
RECMINYes & 0.530 & 0.530 & 0.720\\
TRIALYes & 0.959 & 0.959 & 0.986\\
PRS4/5 & 0.272 & 0.272 & 0.152\\
PRSNone & 0.386 & 0.386 & 0.515\\
\addlinespace
PRSREVOC/RFEL & 0.051 & 0.051 & 0.023\\
MALEMale & 0.827 & 0.827 & 0.754\\
DOSAGE & 32.909 & 32.909 & 34.787\\
DOSAGEQ & 1212.969 & 1212.969 & 1341.995\\
\bottomrule
\end{tabular}
```
:::

```{.r .cell-code}
c(summary(eb.pcs2$w), "Std. Dev" = sd(eb.pcs2$w))
```

::: {.cell-output .cell-output-stdout}
```
       Min.     1st Qu.      Median        Mean     3rd Qu.        Max. 
 0.07608166  0.18827935  0.29694217  0.36704165  0.45845299 16.66036196 
   Std. Dev 
 0.28137440 
```
:::
:::

::: {.cell hash='entropy_balancing_cache/pdf/pcs-inc-analysis_3fa790fb1d281fd0ff3318f3579afae6'}

```{.r .cell-code}
weights2 <- numeric(nrow(pcs.cc2))
for (i in 1:nrow(pcs.cc2)) {
    if (treatment[i]) {
        weights2[i] <- 1
    }
    else {
        weights2[i] <- eb.pcs2$w[i]
    }
}

fit_pcs_unweighted2 <- summary(glm(INCAR ~ ., 
                          data = pcs.cc2, 
                          family = binomial(link = "logit"),
                          model = FALSE, 
                          y = FALSE))

fit_pcs_weighted2 <- summary(glm(INCAR ~ ., 
                          data = pcs.cc2, 
                          weights = weights2,
                          family = binomial(link = "logit"),
                          model = FALSE, 
                          y = FALSE))
```

::: {.cell-output .cell-output-stderr}
```
Warning in eval(family$initialize): non-integer #successes in a binomial glm!
```
:::

```{.r .cell-code}
print("Unweighted Race Effect Estimate with Incomplete Data:")
```

::: {.cell-output .cell-output-stdout}
```
[1] "Unweighted Race Effect Estimate with Incomplete Data:"
```
:::

```{.r .cell-code}
fit_pcs_unweighted2$coefficients["OFF_RACERBLACK",]
```

::: {.cell-output .cell-output-stdout}
```
    Estimate   Std. Error      z value     Pr(>|z|) 
3.259340e-01 5.390277e-02 6.046702e+00 1.478410e-09 
```
:::

```{.r .cell-code}
print("Weighted Race Effect Estimate with Incomplete Data:")
```

::: {.cell-output .cell-output-stdout}
```
[1] "Weighted Race Effect Estimate with Incomplete Data:"
```
:::

```{.r .cell-code}
fit_pcs_weighted2$coefficients["OFF_RACERBLACK",]
```

::: {.cell-output .cell-output-stdout}
```
    Estimate   Std. Error      z value     Pr(>|z|) 
3.640781e-01 8.093915e-02 4.498171e+00 6.854065e-06 
```
:::
:::


## Implemented by the PSweight Package

Following the [vignette](https://cran.r-project.org/web/packages/PSweight/vignettes/vignette.pdf) from the PSweight package, I will implement the estimation of the ATE, ATO, and ATT by using overlap weights (OW) and entropy balancing (EB). The vignette includes formulae for how to calculate causal estimands like the ATE, ATO, ATT, ATEN and how those estimands relate to the choice of weights. They also note that the entropy balancing has no theoretical foundation is not as frequently used as inverse propensity weighting (IPW) or OW.

Following the steps in this framework, I will obtain causal estimands and be able to demonstrate that these inferences are also subject to the influence of missingness.

### Binary Treatment

As an initial step, I'll reduce the comparison to just White and Black defendants being Black set to be the "treatment". Later I will consider the same analysis with multiple treatment levels where White is the control and Black, Latino, and Other are treatments.


::: {.cell hash='entropy_balancing_cache/pdf/psweight-pcs-design_881c3fcd8498ea556869645819b6d0f2'}

```{.r .cell-code}
# library(PSweight)
# 
# pcs.cc %>% 
#     filter(OFF_RACER %in% c("WHITE", "BLACK")) %>% 
#     mutate(BLACK = ifelse(OFF_RACER == "BLACK", 1, 0)) %>% 
#     select(-c("OFF_RACER")) -> pcs.bw
# 
# ps.any <- BLACK ~ OGS + OGSQ + as.factor(CRIMETYPE) + as.factor(RECMIN) + as.factor(TRIAL) + as.factor(PRS) + as.factor(MALE) + DOSAGE + DOSAGEQ
# 
# bal.any <- SumStat(ps.formula = ps.any, data = pcs.bw,
#                    weight = c("IPW", "overlap", "treated", "entropy"))
# bal.any
```
:::

::: {.cell hash='entropy_balancing_cache/pdf/psweight-pcs-design-plots_58b7f38ed86856b4ff8338d6397eb4ea'}

```{.r .cell-code}
# plot(bal.any, type = "density")
# plot(bal.any, type = "balance", metric = "PSD")
```
:::


Now we continue to the analysis step, first without including confounders.


::: {.cell hash='entropy_balancing_cache/pdf/psweight-pcs-analysis_47214e627675abe309384b883fa10c14'}

```{.r .cell-code}
# ate.any <- PSweight(ps.formula = ps.any, 
#                     yname = "INCAR", 
#                     data = pcs.bw,
#                     weight = "IPW",
#                     family = "binomial")
# 
# att.any <- PSweight(ps.formula = ps.any, 
#                     yname = "INCAR", 
#                     data = pcs.bw,
#                     weight = "treated",
#                     family = "binomial")
# 
# ato.any <- PSweight(ps.formula = ps.any, 
#                     yname = "INCAR", 
#                     data = pcs.bw,
#                     weight = "overlap",
#                     family = "binomial")
# 
# aten.any <- PSweight(ps.formula = ps.any, 
#                     yname = "INCAR", 
#                     data = pcs.bw,
#                     weight = "entropy",
#                     family = "binomial")
# 
# summary(ate.any, type = "OR")
# summary(att.any, type = "OR")
# summary(ato.any, type = "OR")
# summary(aten.any, type = "OR")
```
:::


And including confounders.


::: {.cell hash='entropy_balancing_cache/pdf/psweight-pcs-aug-analysis_d34f1b2ee4d4e69b357149d53000f0e5'}

```{.r .cell-code}
# out.incar <- INCAR ~ as.factor(CRIMETYPE) + OGS + OGSQ + as.factor(RECMIN) + as.factor(TRIAL) + as.factor(PRS) + as.factor(MALE) + DOSAGE + DOSAGEQ + as.factor(YEAR) + as.factor(COUNTY)
# 
# ate.any.aug <- PSweight(ps.formula = ps.any, 
#                         out.formula = out.incar,
#                         augmentation = TRUE,
#                         yname = "INCAR", 
#                         data = pcs.sub,
#                         weight = "IPW",
#                         family = "binomial")
# 
# summary(ate.any.aug, type = "OR")
```
:::

::: {.cell hash='entropy_balancing_cache/pdf/unnamed-chunk-15_77e0d23421b47393698f2b982723e821'}

```{.r .cell-code}
# ato.any.aug <- PSweight(ps.formula = ps.any, 
#                         out.formula = out.incar,
#                         augmentation = TRUE,
#                         yname = "INCAR", 
#                         data = pcs.bw[s,],
#                         weight = "overlap",
#                         family = "binomial")
# 
# summary(ato.any.aug, type = "OR")
```
:::

::: {.cell hash='entropy_balancing_cache/pdf/unnamed-chunk-16_481c7d5062476e967dfc701c5c441cd6'}

```{.r .cell-code}
# aten.any.aug <- PSweight(ps.formula = ps.any, 
#                         out.formula = out.incar,
#                         augmentation = TRUE,
#                         yname = "INCAR", 
#                         data = pcs.bw[s,],
#                         weight = "entropy",
#                         family = "binomial")
# 
# summary(aten.any.aug, type = "OR")
```
:::


## References
