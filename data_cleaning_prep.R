# data_cleaning_prep.R

library(dplyr)
library(magrittr)
library(ggplot2)
library(RColorBrewer)
library(mice)
library(lattice)
library(VIM)
library(haven)

dat <- read_dta("Data/Most Serious Sentence (2010-2019, limited).dta")
colnames(dat)

table(dat$INCAR, useNA = "ifany")
table(dat$MALE, useNA = "ifany")
table(dat$BLACK, dat$LATINO, useNA = "ifany")
table(dat$OFF_RACER, useNA = "ifany")
table(dat$PERSONS, useNA = "ifany")
table(dat$RECMIN, useNA = "ifany")
table(dat$TRIAL, useNA = "ifany")
table(dat$PRS_0, useNA = "ifany")
table(dat$PRS_123, useNA = "ifany")
table(dat$PRS_45, useNA = "ifany")
table(dat$PRS_RR, useNA = "ifany")

summary(dat$DOSAGE)
table(dat$OGS)
summary(dat$OGS)

dat$PRS <- factor(dat$PRS_0 + 2 * dat$PRS_123 + 3 * dat$PRS_45 + 4 * dat$PRS_RR,
                  labels = c("None", "1/2/3", "4/5", "REVOC/RFEL"))
table(dat$PRS, useNA = "ifany")

dat$CRIMETYPE <- factor(dat$PERSONS + 2 * dat$PROPERTY + 3 * dat$DRUGNEW 
                        + 4 * dat$DUI + 5 * dat$OTHER,
                        labels = c("PERSONS", "PROPERTY", "DRUGNEW",
                                   "DUI", "OTHER"))
table(dat$CRIMETYPE, useNA = "ifany")

dat$YEAR <- factor(dat$YEAR1 + 2 * dat$YEAR2 + 3 * dat$YEAR3 + 4 * dat$YEAR4 
                   + 5 * dat$YEAR5 + 6 * dat$YEAR6 + 7 * dat$YEAR7 
                   + 8 * dat$YEAR8 + 9 * dat$YEAR9 + 10 * dat$YEAR10 + 2009)
table(dat$YEAR, useNA = "ifany")

dat$COUNTY <- rep(0, nrow(dat))
for (k in 1:67) {
    dat[dat[,paste0("COUNTY", k)] == 1, "COUNTY"] <- k
}
dat$COUNTY <- factor(dat$COUNTY)
table(dat$COUNTY, useNA = "ifany")
table(dat$COUNTY32, useNA = "ifany")

dat.sm <- dat %>% select(INCAR, CRIMETYPE, OGS, OGSQ, RECMIN, TRIAL, PRS, 
                         MALE, DOSAGE, DOSAGEQ, OFF_RACER, COUNTY, YEAR)

# write_csv(dat.sm, "Data/most_serious_sentence_2010-2019_reduced.csv")
rm(dat)

md.pattern(dat.sm, rotate.names = T)
prop_mis <- 3673 / nrow(dat.sm)
prop_mis

marginplot(dat.sm[,c("RECMIN", "PRS")])
pbox(dat.sm, pos = 3, main = "Offense Gravity Score by Missingness")
pbox(dat.sm, pos = 9, main = "Offender Age by Missingness")

# CCA Logistic Regression INCAR ~ .
fit <- glm(INCAR ~ ., data = dat.sm,
           family = binomial(link = "logit"))
summary(fit)
coef(fit) %>% round(3)
beta <- coef(fit)
# rm(fit)

# MICE Logistic Regression INCAR ~ .
# imp1 <- mice(dat.sm, m = 3)
# imp1
# 
# imp1$imp$PRS

# Data Simulation

# Crimetype

# OGS

# TRIAL

# MALE

# OFF_RACER

# COUNTY

# YEAR

# PRS

# AGE

# RECMIN

