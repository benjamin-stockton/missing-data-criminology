#' Complete Case Analysis Simulation
#'
#' @param data A data frame.
#' @param sim_size A string indicating whether the simulation is for testing ("small") or results ("full").
#'
#' @return A named numeric vector with the odds ratio and proportion of misingness.
sim_log_reg <- function(data, sim_size = "small") {
    df <- data %>% filter(complete.cases(.))
    # if (sim_size == "full") {
    #     Z1 <- ifelse(df$OFF_RACER == "BLACK" & df$INCAR == 1, 1, 0)
    #     Z2 <- ifelse(df$OFF_RACER == "WHITE" & df$INCAR == 1, 1, 0)
    #     Z3 <- ifelse(df$OFF_RACER == "BLACK" & df$INCAR == 0, 1, 0)
    #     Z4 <- ifelse(df$OFF_RACER == "WHITE" & df$INCAR == 0, 1, 0)
    #     
    #     # print(paste0("P(RACE == BLACK & INCAR == 1) = ", mean(Z1)))
    #     # p_miss <- (mean(Z1) / mean(Z3)) / (mean(Z2) / mean(Z4))
    # }
    # else {
    #     p_miss <- 1 - nrow(df) / nrow(data)
    # }
    p_miss <- 1 - nrow(df) / nrow(data)
    
    
    if (sim_size == "small") {
        fit <- glm(Y ~ ., data = df, family = binomial(link = "logit"))
        OR <- exp(coef(fit))["X"]
    }
    else if (sim_size == "full") {
        # print(prop.table(table(df$INCAR, df$OFF_RACER), margin = 2))
        fit <- glm(INCAR ~ ., data = df, family = binomial(link = "logit"))
        OR <- exp(coef(fit))["OFF_RACERBLACK"]
    }
    return(c("OR" = OR, "Prop_Miss" = p_miss))
}



#' Title
#'
#' @param data A data frame.
#' @param sim_size A string indicating whether the simulation is for testing ("small") or results ("full").
#' @param m An integer for how many imputations to run
#'
#' @return A numeric. Odds Ratio.
mi_log_reg <- function(data, sim_size = "small", m = 3) {
    # Create the multiply imputed data sets
    # Have mice passively impute age^2 by squaring the imputed age
    # Stack the data sets together
    # imp_tot <- complete(imp, "long", inc = TRUE)
    # Run the analyses
    if (sim_size == "small") {
        imp <- mice(data, m = m)
        fitm <- with(imp, glm(Y ~ X + Z1 + Z1Q + Z2 + CTY, family = binomial(link = "logit")))
        pool_fit <- pool(fitm)
        x_index <- which(pool_fit$pooled[,1] == "X")
        OR <- exp(pool_fit$pooled[x_index, "estimate"])
    }
    
    else if (sim_size == "full") {
        ini <- mice(data, m = 1, maxit = 0, print = F, method = "pmm")
        
        methd <- ini$method
        methd["DOSAGEQ"] <- "~I(DOSAGE^2)"
        print(methd)
        pred_mat <- ini$predictorMatrix
        pred_mat["DOSAGE", "DOSAGEQ"] <- 0
        
        imp <- mice(data, m = m, predictorMatrix = pred_mat, method = methd)
        fitm <- with(imp, glm(INCAR ~ CRIMETYPE + OGS + OGSQ + RECMIN + TRIAL + PRS + 
                                  MALE + DOSAGE + DOSAGEQ + OFF_RACER + COUNTY + YEAR, 
                              family = binomial(link = "logit")))
        pool_fit <- pool(fitm)
        race_index <- which(pool_fit$pooled[,1] == "OFF_RACERBLACK")
        OR <- exp(pool_fit$pooled[race_index, "estimate"])
    }
    return(OR)
}