library(palmerpenguins)
library(mice)
library(VIM)
library(lattice)

head(penguins)

md.pattern(penguins)

table(penguins$species)

imp <- mice(penguins, m = 10)
imp

imp$imp$sex
imp_tot2 <- complete(imp, "long", inc = TRUE)

## labels observed data in blue and imputed data in red for y1
col <- rep(c("blue", "red")[1 + as.numeric(is.na(imp$data$sex))], 6)
## plots data for y1 by imputation
stripplot(sex ~ .imp, data = imp_tot2, jit = TRUE, col = col, xlab = "imputation Number")

fitm <- with(imp, lm(bill_length_mm ~ bill_depth_mm + species + island + flipper_length_mm + sex))
summary(fitm)

p.est <- pool(fitm)
cof <- which(p.est$pooled[,1] == "bill_depth_mm")
p.est$pooled[cof, "estimate"]

summary(p.est)
