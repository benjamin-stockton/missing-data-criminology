# Data and Missingness SImulation

----
## Simulating Predictors
There are several options for simulating the covariates:
1. Bootstrap sampling from the observed data which maintains relationships between covariates (least control).
2. Sample from observed data within each covariate independent of other columns which breaks relationships while maintaining observed distributions.
3. Sample from distributions estimated by each covariate independently of other covariates (most control).
	1. Could keep some relationships by choice; ie distribution of race by county and year
4. Sampling from regressions with each covariate as the outcome and all other covariates and the response as the predictors.

I'll continue with option 3 for now with race distributed by county. 

## Simulating the Outcome
Regardless of the choice for simulation of the predictors, the response is simulated by drawing $Y_i \sim Bern(\hat{p}_i)$ where $\hat{p}_i = (1 + \exp(-\mathbf{x}_i '\hat{\beta}))^{-1}$ and $\hat{\beta}$ comes from the logistic regression fitted with CCA on the observed data.

## Simulations to Determine the Missingness Strategy
In the observed data, the estimate for the race effect is $\hat{\beta} = 0.227$ which gives an odds ratio of 1.255 which I'll use for the basis of my simulations to determine how to create missingness. In those simulations I'm creating a "County" variable with four different prevalences of category 2 in X. X is then determined by sampling by county from $Bern(p_j)$ for $j = 1,2,3,4$. The outcome Y is determined by sampling $Y_i \sim Bern(\hat{p}_i)$ where $\hat{p}_i = (1 + \exp{-(.05 + .227 x_i)})^{-1}$. A new data set and new missingness is generated for each iteration of the simulation.

Those simulations demonstrate that I can get a distinct under or over estimation of the effect in MNAR and MAR (with flipping the effect under MNAR) if I'm interested in the difference between $OR_{miss} - OR_{comp}$ with similar amounts of missingness to what's in the reported criminology data sets (< 5% missingness in total; primarily in Race).
![[Pasted image 20220622111435.png]]

![[Pasted image 20220622111402.png]]

Since I can't get the distinct under/over estimation with MAR missingness in X based solely on Y, I'll need to include the "County" grouping in this simulation and in the simulated criminology data.
![[Pasted image 20220622111616.png]]

## First Pass at Criminology Simulations

Running the same simulations but now with the simulated criminology data set. The regression is now on all twelve predictors instead of just race/X. If a defendant doesn't meet the qualifications for missingness, we assume they are completely observed. Probability of missingness is set to get roughly 3-5% of defendants missing Race observations in the sample.

**MNAR:**

Overstate the Race effect: Missingness in Race is conditional on Race being non-White & not being incarcerated.
$$P(\mathrm{RACE ~is ~Miss} ~|~ \mathrm{RACE ~is~not~WHITE}, ~\mathrm{INCAR} = 0) = .25$$

Understate the Race effect: Missingness in Race is conditional on Race being non-White & being incarcerated.
$$P(\mathrm{OFF\_RACER ~is ~Miss} ~|~ \mathrm{OFF\_RACER ~is~not~WHITE}, ~\mathrm{INCAR} = 1) = .3$$

**MAR:**

Overstate the Race effect: Missingness in Race conditional on being in a county with  <50% White defendants & not being incarcerated OR being in county with >90% White defendants & being incarcerated.
$$P[\mathrm{RACE ~is ~Miss} ~|~ (\mathrm{County ~is~<50\% ~ WHITE}, ~\mathrm{INCAR} = 0) \cup (\mathrm{County ~is~>90\% ~ WHITE}, ~\mathrm{INCAR} = 1)] = .2$$

Understate the Race effect: Missingness in Race is conditional on being in a county with <50% White defendants & being incarcerated OR being in a county with >90% White defendants & not being incarcerated.

$$P[\mathrm{RACE ~is ~Miss} ~|~ (\mathrm{County ~is~<50\% ~ WHITE}, ~\mathrm{INCAR} = 1) \cup (\mathrm{County ~is~>90\% ~ WHITE}, ~\mathrm{INCAR} = 0)] = .2$$

![[Pasted image 20220623142725.png]]

County is not strongly enough associated with Race to affect the Race effect even with 90% of the defendents from the selected counties missing Race. 

## Second Pass at Criminology Simulations

Simulations were good for the MNAR condition, so I'll focus on the MAR condition from here on.

Whiteness is predicted by INCAR and COUNTY, nothing else. That's as it should be.

>Call:
glm(formula = RACE_WHITE ~ INCAR + COUNTY + YEAR + RECMIN + TRIAL + 
    PRS + OGS + OGSQ + DOSAGE + DOSAGEQ + CRIMETYPE + MALE, family = binomial(link = "logit"), 
    data = dat.norace)
>
Deviance Residuals: 
    Min       1Q   Median       3Q      Max  
-2.6634  -0.9384   0.5979   0.7642   1.7350  
>
Coefficients:
                   Estimate Std. Error z value Pr(>|z|)    
(Intercept)       3.825e-01  2.185e-01   1.750 0.080055 .  
INCAR            -2.911e-01  6.019e-02  -4.836 1.32e-06 ***
COUNTY52          2.319e+00  6.097e-01   3.804 0.000143 ***
COUNTY46          1.556e-01  1.013e-01   1.537 0.124308    
COUNTY65          1.396e+00  1.619e-01   8.619  < 2e-16 ***
COUNTY9           1.358e+00  1.319e-01  10.300  < 2e-16 ***
COUNTY51         -1.013e+00  1.089e-01  -9.296  < 2e-16 ***
COUNTY23         -1.130e-01  1.027e-01  -1.100 0.271267    
COUNTY67          1.098e+00  1.208e-01   9.086  < 2e-16 ***
COUNTY28          1.423e+00  2.075e-01   6.859 6.95e-12 ***
COUNTY6           1.111e+00  1.520e-01   7.312 2.64e-13 ***
COUNTY22         -1.770e-01  1.301e-01  -1.361 0.173669    
COUNTY21          1.469e+00  2.001e-01   7.341 2.12e-13 ***
COUNTY39          9.333e-01  1.460e-01   6.395 1.61e-10 ***
COUNTY15          1.155e+00  1.602e-01   7.212 5.52e-13 ***
COUNTY54          2.498e+00  3.727e-01   6.702 2.06e-11 ***
COUNTY13          2.920e+00  5.173e-01   5.644 1.66e-08 ***
COUNTY14          1.725e+00  3.312e-01   5.209 1.90e-07 ***
COUNTY40          1.308e+00  1.704e-01   7.678 1.62e-14 ***
COUNTY10          2.416e+00  3.190e-01   7.576 3.57e-14 ***
COUNTY18          2.789e+00  7.319e-01   3.810 0.000139 ***
COUNTY26          1.333e+00  2.083e-01   6.398 1.57e-10 ***
COUNTY25          7.354e-01  1.621e-01   4.537 5.71e-06 ***
COUNTY36          1.458e+00  1.742e-01   8.367  < 2e-16 ***
COUNTY48          1.154e+00  1.636e-01   7.051 1.78e-12 ***
COUNTY33          1.646e+01  3.568e+02   0.046 0.963201    
COUNTY35          9.700e-01  1.846e-01   5.254 1.49e-07 ***
COUNTY38          1.645e+00  2.878e-01   5.715 1.10e-08 ***
COUNTY31          1.785e+00  4.840e-01   3.689 0.000225 ***
COUNTY45          9.074e-01  2.134e-01   4.252 2.12e-05 ***
COUNTY8           1.651e+01  3.936e+02   0.042 0.966531    
COUNTY20          1.669e+00  3.213e-01   5.193 2.07e-07 ***
COUNTY56          1.732e+00  4.447e-01   3.895 9.81e-05 ***
COUNTY55          2.265e+00  6.078e-01   3.726 0.000195 ***
COUNTY44          3.300e+00  7.230e-01   4.564 5.03e-06 ***
COUNTY7           1.624e+00  2.496e-01   6.508 7.61e-11 ***
COUNTY30          2.619e+00  1.039e+00   2.519 0.011755 *  
COUNTY41          9.027e-01  1.902e-01   4.746 2.07e-06 ***
COUNTY43          9.294e-01  2.341e-01   3.970 7.18e-05 ***
COUNTY32          1.552e+00  3.219e-01   4.820 1.44e-06 ***
COUNTY1           1.813e+00  3.029e-01   5.987 2.14e-09 ***
COUNTY49          1.933e+00  3.592e-01   5.381 7.42e-08 ***
COUNTY16          2.269e+00  6.071e-01   3.737 0.000186 ***
COUNTY11          1.332e+00  2.511e-01   5.306 1.12e-07 ***
COUNTY64          1.655e+01  4.690e+02   0.035 0.971848    
COUNTY61          1.651e+01  3.315e+02   0.050 0.960283    
COUNTY63          9.095e-01  2.396e-01   3.796 0.000147 ***
COUNTY4           5.314e-01  1.821e-01   2.919 0.003513 ** 
COUNTY42          3.097e+00  7.246e-01   4.274 1.92e-05 ***
COUNTY34          2.947e+00  1.027e+00   2.870 0.004098 ** 
COUNTY62          3.091e+00  1.023e+00   3.021 0.002516 ** 
COUNTY58          2.347e+00  7.395e-01   3.174 0.001502 ** 
COUNTY3           3.129e+00  7.248e-01   4.316 1.59e-05 ***
COUNTY37          8.410e-01  3.001e-01   2.803 0.005067 ** 
COUNTY59          2.150e+00  7.450e-01   2.886 0.003896 ** 
COUNTY66          3.131e+00  1.023e+00   3.061 0.002205 ** 
COUNTY47          2.738e+00  1.037e+00   2.639 0.008320 ** 
COUNTY29          1.642e+01  5.812e+02   0.028 0.977459    
COUNTY19          2.912e+00  7.274e-01   4.003 6.26e-05 ***
COUNTY5           2.900e+00  7.294e-01   3.975 7.03e-05 ***
COUNTY50          2.776e+00  7.314e-01   3.796 0.000147 ***
COUNTY57          1.660e+01  1.196e+03   0.014 0.988923    
COUNTY60          2.351e+00  7.412e-01   3.171 0.001519 ** 
COUNTY53          1.648e+01  4.992e+02   0.033 0.973658    
COUNTY17          1.654e+01  4.104e+02   0.040 0.967853    
COUNTY24          1.646e+01  6.169e+02   0.027 0.978714    
COUNTY27          1.663e+01  1.693e+03   0.010 0.992164    
COUNTY12          1.328e+00  1.122e+00   1.183 0.236701    
YEAR2010          2.995e-02  1.082e-01   0.277 0.781904    
YEAR2015         -2.434e-03  1.088e-01  -0.022 0.982156    
YEAR2018         -1.186e-01  1.088e-01  -1.091 0.275419    
YEAR2013         -3.624e-02  1.061e-01  -0.342 0.732634    
YEAR2016         -1.576e-01  1.075e-01  -1.466 0.142765    
YEAR2011         -1.674e-02  1.089e-01  -0.154 0.877884    
YEAR2014         -1.971e-02  1.061e-01  -0.186 0.852598    
YEAR2012          2.291e-02  1.092e-01   0.210 0.833778    
YEAR2019         -2.046e-02  1.068e-01  -0.192 0.848056    
RECMINNo          1.177e-01  5.271e-02   2.233 0.025579 *  
TRIALPlea        -2.225e-01  1.683e-01  -1.322 0.186267    
PRSNone          -1.062e-01  6.726e-02  -1.580 0.114209    
PRS1/2/3         -9.006e-02  7.153e-02  -1.259 0.207973    
PRSREVOC/RFEL    -1.395e-01  1.484e-01  -0.940 0.347225    
OGS              -4.871e-03  2.873e-02  -0.170 0.865349    
OGSQ              2.010e-03  2.684e-03   0.749 0.453926    
DOSAGE           -2.558e-03  9.226e-03  -0.277 0.781611    
DOSAGEQ           8.248e-06  1.123e-04   0.073 0.941458    
CRIMETYPEDRUGNEW -6.160e-02  6.946e-02  -0.887 0.375169    
CRIMETYPEDUI      1.602e-01  7.103e-02   2.256 0.024102 *  
CRIMETYPEOTHER   -9.367e-02  8.237e-02  -1.137 0.255457    
CRIMETYPEPERSONS  6.661e-02  7.772e-02   0.857 0.391405    
MALEFemale       -4.348e-02  5.754e-02  -0.756 0.449859    
>
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
>
(Dispersion parameter for binomial family taken to be 1)
>
    Null deviance: 12047  on 9999  degrees of freedom
Residual deviance: 10440  on 9909  degrees of freedom
AIC: 10622
>
Number of Fisher Scoring iterations: 15

Incarceration is predicted by all predictors except OGS, YEAR, and DOSAGE (DOSAGEQ).

>Call:
glm(formula = INCAR ~ ., family = binomial(link = "logit"), data = dat)
>
Deviance Residuals: 
    Min       1Q   Median       3Q      Max  
-4.0968  -0.7419  -0.2815   0.7552   3.1552  
>
Coefficients:
                   Estimate Std. Error z value Pr(>|z|)    
(Intercept)      -1.878e+00  2.481e-01  -7.571 3.71e-14 ***
CRIMETYPEDRUGNEW -6.516e-01  7.610e-02  -8.563  < 2e-16 ***
CRIMETYPEDUI      1.675e+00  7.456e-02  22.466  < 2e-16 ***
CRIMETYPEOTHER   -2.960e-01  8.834e-02  -3.351 0.000805 ***
CRIMETYPEPERSONS  4.923e-01  8.012e-02   6.145 8.02e-10 ***
OGS              -6.469e-02  3.883e-02  -1.666 0.095716 .  
TRIALPlea         9.926e-01  1.863e-01   5.329 9.89e-08 ***
MALEFemale       -5.473e-01  6.213e-02  -8.809  < 2e-16 ***
COUNTY52          4.279e+00  5.254e-01   8.143 3.86e-16 ***
COUNTY46          2.543e+00  1.340e-01  18.984  < 2e-16 ***
COUNTY65          1.502e+00  1.720e-01   8.734  < 2e-16 ***
COUNTY9           2.320e+00  1.449e-01  16.011  < 2e-16 ***
COUNTY51          1.203e+00  1.441e-01   8.350  < 2e-16 ***
COUNTY23          2.975e+00  1.373e-01  21.672  < 2e-16 ***
COUNTY67          1.428e+00  1.448e-01   9.861  < 2e-16 ***
COUNTY28          1.771e+00  2.111e-01   8.391  < 2e-16 ***
COUNTY6           2.241e+00  1.700e-01  13.182  < 2e-16 ***
COUNTY22          8.863e-01  1.798e-01   4.928 8.31e-07 ***
COUNTY21          2.775e+00  1.991e-01  13.941  < 2e-16 ***
COUNTY39          2.242e+00  1.717e-01  13.055  < 2e-16 ***
COUNTY15          3.966e+00  2.037e-01  19.470  < 2e-16 ***
COUNTY54          2.840e+00  2.428e-01  11.699  < 2e-16 ***
COUNTY13          3.222e+00  2.958e-01  10.891  < 2e-16 ***
COUNTY14          2.069e+00  2.802e-01   7.384 1.54e-13 ***
COUNTY40          1.727e+00  1.811e-01   9.537  < 2e-16 ***
COUNTY10          8.200e-01  2.379e-01   3.447 0.000567 ***
COUNTY18          2.592e+00  4.314e-01   6.007 1.89e-09 ***
COUNTY26          1.210e+00  2.176e-01   5.560 2.70e-08 ***
COUNTY25          1.716e+00  1.921e-01   8.932  < 2e-16 ***
COUNTY36          2.282e+00  1.754e-01  13.012  < 2e-16 ***
COUNTY48          2.564e+00  1.829e-01  14.014  < 2e-16 ***
COUNTY33          2.415e+00  3.574e-01   6.758 1.40e-11 ***
COUNTY35          2.464e+00  2.070e-01  11.900  < 2e-16 ***
COUNTY38          2.665e+00  2.625e-01  10.153  < 2e-16 ***
COUNTY31          2.216e+00  3.993e-01   5.550 2.86e-08 ***
COUNTY45          2.426e+00  2.406e-01  10.083  < 2e-16 ***
COUNTY8           3.292e+00  4.083e-01   8.064 7.41e-16 ***
COUNTY20          3.697e+00  3.260e-01  11.341  < 2e-16 ***
COUNTY56          2.307e+00  4.010e-01   5.754 8.73e-09 ***
COUNTY55          2.443e+00  4.291e-01   5.694 1.24e-08 ***
COUNTY44          3.365e+00  3.460e-01   9.725  < 2e-16 ***
COUNTY7           1.042e+00  2.444e-01   4.263 2.02e-05 ***
COUNTY30          2.565e+00  6.154e-01   4.168 3.07e-05 ***
COUNTY41          1.634e+00  2.194e-01   7.448 9.46e-14 ***
COUNTY43          3.924e+00  2.812e-01  13.953  < 2e-16 ***
COUNTY32          2.278e+00  2.887e-01   7.892 2.97e-15 ***
COUNTY1           1.379e+00  2.616e-01   5.271 1.36e-07 ***
COUNTY49          1.409e+00  3.087e-01   4.562 5.06e-06 ***
COUNTY16          2.157e+00  3.895e-01   5.538 3.06e-08 ***
COUNTY11          2.161e+00  2.536e-01   8.521  < 2e-16 ***
COUNTY64          3.950e+00  5.158e-01   7.658 1.88e-14 ***
COUNTY61          2.752e+00  3.446e-01   7.986 1.40e-15 ***
COUNTY63          1.074e+00  2.826e-01   3.801 0.000144 ***
COUNTY4           7.047e-01  2.501e-01   2.818 0.004838 ** 
COUNTY42          2.276e+00  3.532e-01   6.445 1.15e-10 ***
COUNTY34          2.089e+00  4.901e-01   4.261 2.03e-05 ***
COUNTY62          2.646e+00  4.924e-01   5.375 7.66e-08 ***
COUNTY58          1.860e+00  5.313e-01   3.501 0.000463 ***
COUNTY3           3.333e+00  3.688e-01   9.039  < 2e-16 ***
COUNTY37          1.659e+00  3.414e-01   4.860 1.17e-06 ***
COUNTY59          6.498e-01  5.895e-01   1.102 0.270269    
COUNTY66          3.257e+00  5.123e-01   6.358 2.04e-10 ***
COUNTY47          3.265e+00  6.206e-01   5.261 1.43e-07 ***
COUNTY29          1.779e+00  6.164e-01   2.885 0.003911 ** 
COUNTY19          2.596e+00  3.905e-01   6.648 2.97e-11 ***
COUNTY5           3.314e+00  3.959e-01   8.371  < 2e-16 ***
COUNTY50          2.805e+00  3.978e-01   7.052 1.77e-12 ***
COUNTY57          4.221e+00  1.195e+00   3.532 0.000412 ***
COUNTY60          2.551e+00  5.078e-01   5.025 5.04e-07 ***
COUNTY53          2.679e+00  5.136e-01   5.217 1.82e-07 ***
COUNTY17          4.365e+00  4.704e-01   9.279  < 2e-16 ***
COUNTY24          2.851e+00  6.802e-01   4.192 2.77e-05 ***
COUNTY27          1.397e+01  1.257e+02   0.111 0.911506    
COUNTY12          2.676e+00  1.123e+00   2.384 0.017134 *  
OFF_RACERBLACK    3.219e-01  6.305e-02   5.105 3.31e-07 ***
OFF_RACERLATINO   6.304e-01  2.757e-01   2.286 0.022228 *  
OFF_RACEROTHER   -3.064e-01  3.559e-01  -0.861 0.389316    
YEAR2010          4.192e-01  1.160e-01   3.613 0.000302 ***
YEAR2015          5.809e-02  1.166e-01   0.498 0.618351    
YEAR2018          6.065e-02  1.185e-01   0.512 0.608725    
YEAR2013          2.035e-01  1.148e-01   1.773 0.076232 .  
YEAR2016          2.802e-02  1.160e-01   0.241 0.809194    
YEAR2011          2.400e-01  1.173e-01   2.046 0.040755 *  
YEAR2014          1.060e-01  1.134e-01   0.934 0.350092    
YEAR2012          2.302e-01  1.165e-01   1.977 0.048053 *  
YEAR2019         -1.429e-01  1.148e-01  -1.245 0.213202    
PRSNone          -1.669e+00  7.158e-02 -23.317  < 2e-16 ***
PRS1/2/3         -8.559e-01  7.478e-02 -11.445  < 2e-16 ***
PRSREVOC/RFEL     6.850e-01  1.694e-01   4.044 5.25e-05 ***
DOSAGE           -1.255e-02  1.014e-02  -1.238 0.215739    
RECMINNo          1.268e+00  5.610e-02  22.598  < 2e-16 ***
OGSQ              4.473e-02  4.372e-03  10.231  < 2e-16 ***
DOSAGEQ          -1.168e-04  1.251e-04  -0.933 0.350627    
>
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
>
(Dispersion parameter for binomial family taken to be 1)
>
    Null deviance: 13780.0  on 9999  degrees of freedom
Residual deviance:  9368.5  on 9907  degrees of freedom
AIC: 9554.5
>
Number of Fisher Scoring iterations: 10

In my small simulation, I set it so that the predicted probabilities for the logistic regression of $(Y == 1 \& X == 1)$ or $(Y == 0 \& X == 1)$ on $CTY$ so that the missingness is a function of $X$ and $Y$. But this doesn't change the results at all even though this is now MNAR. 