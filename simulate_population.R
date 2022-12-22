# simulate_population.R

# A minimal script that reimplements the data generation in "data_cleaning_prep.Rmd". 
# This file removes the plots, tables, and marginal summaries that aren't used 
# in the sampling procedure.

library(dplyr, warn.conflicts = FALSE)
# Suppress summarise info
options(dplyr.summarise.inform = FALSE)
library(magrittr)
library(haven)
library(readr)

## THE COMMENTED OUT BLOCK REQUIRES THE PCS DATA.
# dat <- read_dta("Data/Most Serious Sentence (2010-2019, limited).dta")
# colnames(dat)
# 
# 
# # Data cleaning and pivoting
# 
# dat$TRIAL <- factor(dat$TRIAL, labels = c("Yes", "Plea"))
# dat$MALE <- factor(dat$MALE, labels = c("Female", "Male"))
# dat$RECMIN <- factor(dat$RECMIN, labels = c("Yes", "No"))
# 
# dat[dat$UNKNOWN == 1, "OFF_RACER"] <- NA
# dat$OFF_RACER <- factor(dat$OFF_RACER,
#                         labels = c("WHITE", "BLACK", "LATINO", "OTHER"))
# 
# # Previous Criminal History
# 
# dat$PRS <- factor(dat$PRS_0 + 2 * dat$PRS_123 + 3 * dat$PRS_45 + 4 * dat$PRS_RR,
#                   labels = c("None", "1/2/3", "4/5", "REVOC/RFEL"))
# 
# # Crime Type
# 
# dat$CRIMETYPE <- factor(dat$PERSONS + 2 * dat$PROPERTY + 3 * dat$DRUGNEW 
#                         + 4 * dat$DUI + 5 * dat$OTHER,
#                         labels = c("PERSONS", "PROPERTY", "DRUGNEW",
#                                    "DUI", "OTHER"))
# 
# # Year
# 
# dat$YEAR <- factor(dat$YEAR1 + 2 * dat$YEAR2 + 3 * dat$YEAR3 + 4 * dat$YEAR4 
#                    + 5 * dat$YEAR5 + 6 * dat$YEAR6 + 7 * dat$YEAR7 
#                    + 8 * dat$YEAR8 + 9 * dat$YEAR9 + 10 * dat$YEAR10 + 2009)
# 
# # County
# 
# dat$COUNTY <- rep(0, nrow(dat))
# for (k in 1:67) {
#     dat[dat[,paste0("COUNTY", k)] == 1, "COUNTY"] <- k
# }
# dat$COUNTY <- factor(dat$COUNTY)
# 
# 
# 
# 
# dat.sm <- dat %>% select(INCAR, CRIMETYPE, OGS, OGSQ, RECMIN, TRIAL, PRS, 
#                          MALE, DOSAGE, DOSAGEQ, OFF_RACER, COUNTY, YEAR)
# 
# write_csv(dat.sm, "Data/most_serious_sentence_2010-2019_slim.csv")
# rm(dat)
# 
# head(dat.sm)
# 
# # CCA Logistic Regression of PCS data set
# # R automatically converts factors to dummy variables
# # INCAR ~ CRIMETYPE + OGS + OGSQ + RECMIN + TRIAL + PRS + MALE + DOSAGE + DOSAGEQ + OFF_RACER + COUNTY + YEAR
# fit <- glm(INCAR ~ INCAR ~ CRIMETYPE + OGS + OGSQ + RECMIN + TRIAL +
#                PRS + MALE + DOSAGE + DOSAGEQ + OFF_RACER + COUNTY + YEAR,
#            data = dat.sm,
#            family = binomial(link = "logit"), model = F)
# # summary(fit)
# coef(fit) %>% round(3)
# beta <- coef(fit)
# 
# fit$data <- NULL
# fit$data
# fit$y <- NULL
# 
# 
# # save(fit, file = "full_data_cca_log_reg_no_data.rda")
# save(beta, file = "full_data_cca_log_reg_coef.rda")
# 
# rm(fit)

## THE REST OF THE SCRIPT CAN BE RUN WITHOUT THE PCS DATA

load("full_data_cca_log_reg_coef.rda")
load("list_of_variable_summaries.rda")

# The odds ratio of the Race Effect (comparing Black defendants to White defendants) is:
round(exp(beta)["OFF_RACERBLACK"], 3)

## Simulating Data 

N.sim <- list_of_var_smrys$N
sim.dat <- data.frame(INCAR = rep(0, N.sim))

# list_of_var_smrys[["N"]] <- N.sim

### Simulating predictors

#### Crimetype (Categorical, 5 levels) - Person, Property, Drug-new, DUI, Other

# list_of_var_smrys[["CRIMETYPE"]] <- prop.table(table(dat.sm$CRIMETYPE))
pt <- list_of_var_smrys$CRIMETYPE

sim.dat$CRIMETYPE <- factor(sample(1:5, size = N.sim, prob = pt, replace = T),
                            labels = names(list_of_var_smrys$CRIMETYPE))

#### Offense Gravity Score (ordinal, 1-14)
# list_of_var_smrys[["OGS"]] <- prop.table(table(dat.sm$OGS))
pt <- list_of_var_smrys$OGS

sim.dat$OGS <- sample(1:14, size = N.sim, prob = pt, replace = T)

#### Trial (Binary) - Yes = 1 or Plea Negotiation = 0
# list_of_var_smrys[["TRIAL"]] <- prop.table(table(dat.sm$TRIAL))
pt <- list_of_var_smrys$TRIAL

sim.dat$TRIAL <- factor(sample(0:1, size = N.sim, prob = pt, replace = T), 
                        labels = c("Yes", "Plea"))

#### Gender (Binary) - Male = 1, Female = 0
# list_of_var_smrys[["MALE"]] <- prop.table(table(dat.sm$MALE))
pt <- list_of_var_smrys$MALE

sim.dat$MALE <- factor(sample(0:1, size = N.sim, prob = pt, replace = T), 
                       labels = c("Female", "Male"))

#### COUNTY (Categorical, 67 levels) - Counties in Pennsylvania
# list_of_var_smrys[["COUNTY"]] <- prop.table(table(dat.sm$COUNTY))
pt <- list_of_var_smrys$COUNTY

sim.dat$COUNTY <- factor(sample(1:67, size = N.sim, prob = pt, replace = T), 
                        labels = names(list_of_var_smrys$COUNTY))

#### Defendent Race (Categorical, 5 levels) - White (ref), Black, Latino, AAPI/Other, Unknown/NA 
# list_of_var_smrys[["OFF_RACER"]] <- prop.table(table(dat.sm$OFF_RACER))
pt <- list_of_var_smrys$OFF_RACER

sim.dat$OFF_RACER <- factor(sample(1:4, size = N.sim, prob = pt, replace = T),
                            labels = c("WHITE", "BLACK", "LATINO", "OTHER"))

#### YEAR (Categorical, 10 levels) Indicators for the year from 2010 to 2019
# list_of_var_smrys[["YEAR"]] <- prop.table(table(dat.sm$YEAR))
pt <- list_of_var_smrys$YEAR

sim.dat$YEAR <- factor(sample(2010:2019, size = N.sim, prob = pt, replace = T))

#### Previous Criminal History PRS (Categorical, 4 levels) - 0, 1/2/3, 4/5, REVOC/RFEL. REVOC/RFEL is the most serious
# list_of_var_smrys[["PRS"]] <- prop.table(table(dat.sm$PRS))
pt <- list_of_var_smrys$PRS

sim.dat$PRS <- factor(sample(1:4, size = N.sim, prob = pt, replace = T), 
                        labels = names(list_of_var_smrys$PRS))

#### Defendant's Age at Sentencing DOSAGE (Numeric, from 11.71 to 113.88)
# X <- dat.sm$DOSAGE - min(dat.sm$DOSAGE, na.rm = T); N <- nrow(dat.sm)
# alpha.hat <- (N - 1) / N *  mean(X, na.rm = T)^2 / var(X, na.rm = T)
# beta.hat <- N / (N - 1) * var(X, na.rm = T) / mean(X, na.rm = T)
alpha.hat <- list_of_var_smrys$DOSAGE["alpha.hat"]
beta.hat <- list_of_var_smrys$DOSAGE["beta.hat"]
min.age <- list_of_var_smrys$DOSAGE["min"]

# list_of_var_smrys[["DOSAGE"]] <- c("alpha.hat" = alpha.hat, "beta.hat" = beta.hat, "min" = min(dat.sm$DOSAGE, na.rm = T))
sim.dat$DOSAGE <- rgamma(N.sim, alpha.hat, scale = beta.hat) + min.age

#### Recommended Minimum Sentence (Binary) - Yes = 1, No = 0
# list_of_var_smrys[["RECMIN"]] <- prop.table(table(dat.sm$RECMIN))
pt <- list_of_var_smrys$RECMIN

sim.dat$RECMIN <- factor(sample(0:1, size = N.sim, prob = pt, replace = T), 
                         labels = names(list_of_var_smrys$RECMIN))

#### Simulating INCAR

# Simulating Incarceration decision based on the simulated predictors using the estimated regression coefficients from the CCA logistic regression. Using the regression fit, I assign with a draw from a Bernoulli distribution with the predicted probability for each observation.

# creating quadratic terms
sim.dat$OGSQ <- sim.dat$OGS ^ 2
sim.dat$DOSAGEQ <- sim.dat$DOSAGE ^ 2

# Fitted probabilities
X <- model.matrix(~ CRIMETYPE + OGS + OGSQ + RECMIN + TRIAL + PRS + MALE + DOSAGE + DOSAGEQ + OFF_RACER + COUNTY + YEAR, data = sim.dat)
pred_prob <- plogis(beta %*% t(X))

sim.dat$INCAR <- rbinom(N.sim, 1, pred_prob)

write_csv(sim.dat, "Data/simulated_data.csv")


