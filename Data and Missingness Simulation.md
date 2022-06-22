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

