missingness_model <- function(data, miss_pars, sim_size = "small", p_miss_target = 0.01) {
    # data: N x p+1 matrix (Y, X)
    # miss_pars: a matrix of parameters for the missingness model
    
    data <- dummify_data_matrix(data, sim_size = sim_size)
    
    # print(str(data))
    if (sim_size == "small") {
        V <- data %>% select(-c(X, Z1, Z1Q, Z2)) %>% as.matrix()
        L <- data %>% select(X, Z1, Z1Q, Z2) %>% as.matrix()
    }
    else if (sim_size == "full") {
        V <- data %>% 
            select(-c(OFF_RACERBLACK, 
                    OFF_RACERLATINO, 
                    OFF_RACEROTHER, 
                    DOSAGE, DOSAGEQ,
                    RECMIN)) %>% 
            as.matrix()
        
        L <- data %>% 
            select(c(INCAR, OFF_RACERBLACK, 
                   OFF_RACERLATINO, 
                   OFF_RACEROTHER, 
                   DOSAGE, DOSAGEQ,
                   RECMIN)) %>% 
            mutate(Z1 = OFF_RACERBLACK * (1 - INCAR),
                   Z2 = OFF_RACERBLACK * INCAR) %>% 
            select(-INCAR) %>% 
            as.matrix()
    }
    # print(colnames(V[,1:6]))
    # print(colnames(L))
    nbeta <- ncol(V); ngamma <- ncol(L)
    
    alpha <- miss_pars[,1]
    beta <- miss_pars[,2:(nbeta+1)]
    gamma <- miss_pars[,(nbeta+2):(nbeta+ngamma+1)]
    
    EV <- apply(V, 2, mean); EL <- apply(L, 2, mean)
    prob.mat <- calc_pattern_probs(V, L, beta, gamma, EV, EL)
    
    p_true <- apply(prob.mat, 2, mean)
    # print(p_true)
    # prob.mat <- prob.mat * p_miss_target / sum(p_true[-1])
    
    # Create probabilities of missingness in Race/X
    # assigned_pattern <- apply(matrix(1:nrow(data), ncol = 1), 1,
    #                           function(r) {
    #                               return(sample(1:8, 1, prob = prob.mat[r,]))
    #                           })
    id_set <- 1:nrow(data)
    assigned_pattern <- rep(1, nrow(data))
    for (r in 2:8) {
        nr <- nrow(data) * p_true[r] * p_miss_target / sum(p_true[-1]) 
        # print(paste0("N_", r, " = ", nr))
        # print(paste0("len(id_set): ", length(id_set)))
        sr <- sample(id_set, size = floor(nr), replace = F, prob = prob.mat[id_set,r])
        assigned_pattern[sr] <- r
        id_set <- setdiff(id_set, sr)
    }
    
    # print(prop.table(table(assigned_pattern)))
    return(assigned_pattern)
}
