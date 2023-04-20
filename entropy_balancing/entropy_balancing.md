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

## Implementation by ebal

### Simulated Data

First, we'll demonstrate entropy balancing in a logistic regression setting to estimate the race effect on the simulated data. This data contains the same number of observations as the PCS data on the same variables which have been simulated independently, so a priori we should expect the marginal means of the variables to be the same for "Black" and "Non-Black" simulated defendants. The simulated data set is complete.


::: {.cell hash='entropy_balancing_cache/pdf/unnamed-chunk-1_169f333dca97614f92442d3be76b7d1a'}

```{.r .cell-code}
library(ebal)
library(dplyr)

dat <- read.csv("../Data/simulated_data.csv")

dat$YEAR <- factor(dat$YEAR)
dat$COUNTY <- factor(dat$COUNTY)

dat$OFF_RACER <- factor(dat$OFF_RACER, levels = c("WHITE", "BLACK", "LATINO", "OTHER"))

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
:::


Entropy balancing seems to be performed very quickly. Next we check the marginal means for each of the covariates


::: {.cell hash='entropy_balancing_cache/pdf/unnamed-chunk-2_19e250440965f62c87baf1c7ef919db7'}

```{.r .cell-code}
c1 <- apply(X[treatment,], 2, mean)


c2 <- apply(X[!treatment,], 2, weighted.mean, w = eb.out$w)


c3 <- apply(X[!treatment,], 2, mean)

X.means <- bind_rows(c1, c2, c3) %>% as.data.frame()
rownames(X.means) <- c("Black", "Non-Black - EB", "Non-Black - Unbalanced")
round(t(X.means), 4)
```

::: {.cell-output .cell-output-stdout}
```
                      Black Non-Black - EB Non-Black - Unbalanced
CRIMETYPEDUI         0.2367         0.2367                 0.2379
CRIMETYPEOTHER       0.1241         0.1241                 0.1230
CRIMETYPEPERSONS     0.1560         0.1560                 0.1568
CRIMETYPEPROPERTY    0.2556         0.2556                 0.2548
OGS                  3.4283         3.4283                 3.4184
TRIALYes             0.9782         0.9782                 0.9786
MALEMale             0.7731         0.7731                 0.7739
PRS4/5               0.1813         0.1813                 0.1820
PRSNone              0.4856         0.4856                 0.4859
PRSREVOC/RFEL        0.0300         0.0300                 0.0300
DOSAGE              34.2666        34.2666                34.2765
RECMINYes            0.6701         0.6701                 0.6707
OGSQ                17.9444        17.9444                17.8505
DOSAGEQ           1306.2898      1306.2898              1307.1996
```
:::
:::


The means are essentially the same between the covariates before balancing, so it won't have any effect on the actual analysis. This is a key point regarding use of the modal approach in the simulations.


::: {.cell hash='entropy_balancing_cache/pdf/unnamed-chunk-3_08c7aaa514f673837809b09d7cc42741'}

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


::: {.cell hash='entropy_balancing_cache/pdf/unnamed-chunk-4_b4adc9563c757f1e33e09d287981d799'}

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

fit_sim_unweighted <- summary(glm(INCAR ~ ., 
                          data = dat, 
                          family = binomial(link = "logit")))

fit_sim_weighted <- summary(glm(INCAR ~ ., 
                          data = dat, 
                          weights = weights,
                          family = binomial(link = "logit")))
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


::: {.cell hash='entropy_balancing_cache/pdf/unnamed-chunk-5_d130609669e5cb1951f8288534548ad8'}

```{.r .cell-code}
pcs <- read.csv("../Data/most_serious_sentence_2010-2019_slim.csv")

pcs$YEAR <- factor(pcs$YEAR)
pcs$COUNTY <- factor(pcs$COUNTY)
pcs$OFF_RACER <- factor(pcs$OFF_RACER, levels = c("WHITE", "BLACK", "LATINO", "OTHER"))

mice::md.pattern(pcs, rotate.names = T)
```

::: {.cell-output-display}
![](entropy_balancing_files/figure-pdf/unnamed-chunk-5-1.pdf){fig-pos='H'}
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

::: {.cell hash='entropy_balancing_cache/pdf/unnamed-chunk-6_c5fe6aaa78eb161eddf7f28eec89cc53'}

```{.r .cell-code}
pcs %>% 
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
round(t(X.means), 4)
```

::: {.cell-output .cell-output-stdout}
```
                      Black Non-Black - EB Non-Black - Unbalanced
CRIMETYPEDUI         0.1274         0.1274                 0.2790
CRIMETYPEOTHER       0.1632         0.1632                 0.1084
CRIMETYPEPERSONS     0.1838         0.1838                 0.1458
CRIMETYPEPROPERTY    0.2445         0.2445                 0.2588
OGS                  4.0601         4.0601                 3.1826
OGSQ                24.7268        24.7268                15.3210
RECMINYes            0.5303         0.5303                 0.7202
TRIALYes             0.9591         0.9591                 0.9858
PRS4/5               0.2720         0.2720                 0.1517
PRSNone              0.3863         0.3863                 0.5152
PRSREVOC/RFEL        0.0506         0.0506                 0.0228
MALEMale             0.8272         0.8272                 0.7545
DOSAGE              32.9122        32.9122                34.7886
DOSAGEQ           1213.2078      1213.2078              1342.1189
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

::: {.cell hash='entropy_balancing_cache/pdf/unnamed-chunk-7_d960e7f744dc6bba311b1a3eed34296e'}

```{.r .cell-code}
weights <- numeric(nrow(pcs.cc))
for (i in 1:nrow(pcs.cc)) {
    if (treatment[i]) {
        weights[i] <- 1
    }
    else {
        weights[i] <- eb.out$w[i]
    }
}

fit_pcs_unweighted <- summary(glm(INCAR ~ ., 
                          data = pcs, 
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
     Estimate    Std. Error       z value      Pr(>|z|) 
 2.268135e-01  6.511253e-03  3.483408e+01 7.418640e-266 
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
     Estimate    Std. Error       z value      Pr(>|z|) 
 2.578608e-01  9.128484e-03  2.824794e+01 1.508926e-175 
```
:::
:::


## Incomplete Data

The impact of incomplete data is the next topic to tackle. We know that CCA can be biased under several (untestable) assumptions for the missingness mechanisms for logistic regression. Need to show that missingness also has an impact on weights doubly impacting the inferential results.


::: {.cell hash='entropy_balancing_cache/pdf/unnamed-chunk-8_5637070ba21c43bf9fa7db5f0cb9b05a'}

```{.r .cell-code}
p <- 1 - (1 + exp(-(-9 + .5 * pcs.cc$DOSAGE + .5 * pcs.cc$OGS)))^-1
boxplot(p)
```

::: {.cell-output-display}
![](entropy_balancing_files/figure-pdf/unnamed-chunk-8-1.pdf){fig-pos='H'}
:::

```{.r .cell-code}
mis <- sample(nrow(pcs.cc), size = ceiling(0.25 * nrow(pcs.cc)), prob = p)
str(mis)
```

::: {.cell-output .cell-output-stdout}
```
 int [1:208637] 238658 215529 148383 356638 752797 68920 88571 434253 722658 330656 ...
```
:::
:::

::: {.cell hash='entropy_balancing_cache/pdf/unnamed-chunk-9_285d3d5426e3ab350e1930bf0d8886e7'}

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
Iteration 1 maximum deviation is = 19573597 
Iteration 2 maximum deviation is = 18358262 
Iteration 3 maximum deviation is = 16777100 
Iteration 4 maximum deviation is = 14647799 
Iteration 5 maximum deviation is = 11521672 
Iteration 6 maximum deviation is = 5521037 
Iteration 7 maximum deviation is = 108551 
Iteration 8 maximum deviation is = 84.07 
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
round(t(X.means), 4)
```

::: {.cell-output .cell-output-stdout}
```
                      Black Non-Black - EB Non-Black - Unbalanced
CRIMETYPEDUI         0.1341         0.1341                 0.2904
CRIMETYPEOTHER       0.1618         0.1618                 0.1096
CRIMETYPEPERSONS     0.1933         0.1933                 0.1531
CRIMETYPEPROPERTY    0.2357         0.2357                 0.2541
OGS                  4.4561         4.4561                 3.4369
OGSQ                28.9857        28.9857                17.5288
RECMINYes            0.4369         0.4369                 0.6647
TRIALYes             0.9513         0.9513                 0.9834
PRS4/5               0.3346         0.3346                 0.1822
PRSNone              0.3054         0.3054                 0.4560
PRSREVOC/RFEL        0.0654         0.0654                 0.0279
MALEMale             0.8267         0.8267                 0.7489
DOSAGE              36.8090        36.8090                38.4811
DOSAGEQ           1475.2901      1475.2901              1595.4856
```
:::

```{.r .cell-code}
c(summary(eb.pcs2$w), "Std. Dev" = sd(eb.pcs2$w))
```

::: {.cell-output .cell-output-stdout}
```
       Min.     1st Qu.      Median        Mean     3rd Qu.        Max. 
 0.08111835  0.17454630  0.26378869  0.35167721  0.43921475 11.91469178 
   Std. Dev 
 0.28622420 
```
:::
:::

::: {.cell hash='entropy_balancing_cache/pdf/unnamed-chunk-10_259ec7e9c7e2211d7f8a280559e2b91d'}

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
     Estimate    Std. Error       z value      Pr(>|z|) 
 2.088674e-01  7.604425e-03  2.746657e+01 4.405814e-166 
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
     Estimate    Std. Error       z value      Pr(>|z|) 
 2.432750e-01  1.084716e-02  2.242754e+01 2.120633e-111 
```
:::
:::


## References
