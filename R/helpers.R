# helpers.R

sample_from_pop <- function(pop_dat, N = 1000, replace = FALSE) {
    s <- sample(1:nrow(pop_dat), size = N, replace = replace)
    samp.dat <- pop_dat[s,]
    return(samp.dat)
}

dummify_data_matrix <- function(data, sim_size = "small") {
    if (sim_size == "small") {
        data %<>% mutate(as.data.frame(model.matrix(~ 0 + CTY, data = cur_data()))) %>% 
            select(-c(CTY, CTY1))
        return(data)
        
    }
    else if (sim_size == "full") {
        data %<>%
            mutate(as.data.frame(model.matrix(~ 0 + OFF_RACER + PRS + CRIMETYPE + YEAR + COUNTY, data = cur_data()))) %>%
            select(-c(COUNTY, OFF_RACER, PRS, YEAR, CRIMETYPE, OFF_RACERWHITE)) 
        data$TRIAL <- ifelse(data$TRIAL == "Yes", 1, 0)
        data$MALE <- ifelse(data$MALE == "Male", 1, 0)
        data$RECMIN <- ifelse(data$RECMIN == "Yes", 1, 0)
        
        return(data)
    }
}



find_label_y_pos <- function(X) {
    d <- density(X)
    return(.5 * max(d$y))
}
