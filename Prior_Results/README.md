## Previous Results

These results come from the summer '22 simulation runs which used the default settings for MICE:
```{r}
imp <- mice(data, m = m)
```
which used the "logreg" imputation for race and "polyreg" imputation for other variables. I thought it was using PMM by default. That is now remedied and the results in the new Plots and Sim_Results Directory in the parent to this Prior_Results directory are "correct".

New imputation call is:
```{r}
ini <- mice(data, m = 1, maxit = 0, print = F, method = "pmm")
        
methd <- ini$method
methd["DOSAGEQ"] <- "~I(DOSAGE^2)"
print(methd)
pred_mat <- ini$predictorMatrix
pred_mat["DOSAGE", "DOSAGEQ"] <- 0

imp <- mice(data, m = m, predictorMatrix = pred_mat, method = methd)
```

