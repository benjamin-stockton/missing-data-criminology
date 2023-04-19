assign_missing <- function(data, miss_pattern, sim_size = "small") {
    df <- data
    
    if (sim_size == "small") {
        df[which(miss_pattern == 2), "X"] <- NA
        df[which(miss_pattern == 3), "Z1"] <- NA
        df[which(miss_pattern == 4), c("X", "Z1")] <- NA
        df[which(miss_pattern == 5), "Z2"] <- NA
        df[which(miss_pattern == 6), c("Z2", "X")] <- NA
        df[which(miss_pattern == 7), c("Z2", "Z1")] <- NA
        df[which(miss_pattern == 8), c("Z1", "Z2", "X")] <- NA
        
        df$Z1Q <- df$Z1^2
    } 
    else if (sim_size == "full") {
        df[which(miss_pattern == 2), "OFF_RACER"] <- NA
        df[which(miss_pattern == 3), "RECMIN"] <- NA
        df[which(miss_pattern == 4), c("OFF_RACER", "RECMIN")] <- NA
        df[which(miss_pattern == 5), "DOSAGE"] <- NA
        df[which(miss_pattern == 6), c("DOSAGE", "OFF_RACER")] <- NA
        df[which(miss_pattern == 7), c("DOSAGE", "RECMIN")] <- NA
        df[which(miss_pattern == 8), c("DOSAGE", "RECMIN", "OFF_RACER")] <- NA
        
        # Induce missingness in dosage^2
        df$DOSAGEQ <- df$DOSAGE^2
    }
    
    return(df)
}

generate_missing_x <- function(data, miss_pars_over, miss_pars_undr, sim_size = "small") {
    p1 <- missingness_model(data, miss_pars = miss_pars_over, sim_size = sim_size)
    p2 <- missingness_model(data, miss_pars = miss_pars_undr, sim_size = sim_size)
    data.over <- assign_missing(data, p1, sim_size = sim_size)
    data.under <- assign_missing(data, p2, sim_size = sim_size)
    
    return(list("COMP" = data, "OVER" = data.over, "UNDR" = data.under))
}

