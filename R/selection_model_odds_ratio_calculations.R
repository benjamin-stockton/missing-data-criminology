p_miss_func <- function(beta, gamma, EV, EL, V, L, pat_i) {
    pat.probs <- c(834548, 26206, 2078, 34, 1465, 92, 1, 1) / 864422
    alpha <- numeric(8)
    prob_mat <- matrix(numeric(8 * nrow(V)), ncol = 8)
    
    for (i in 1:8) {
        alpha[i] <- -log(1 / pat.probs[i] - 1) - (EV %*% beta + EL %*% gamma)[,1]
        lp_i <- alpha[i] + (V %*% beta + L %*% gamma)[,1]
        p_i <- exp(lp_i)
        prob_mat[,i] <- p_i
    }
    prob_mat[,1] <- rep(1, nrow(prob_mat))
    rsum <- apply(prob_mat, 1, sum)
    
    nmat <- t(apply(matrix(1:nrow(prob_mat), ncol = 1), 1, function(r) {
        return(prob_mat[r,] / rsum[r])
    }))
    # nmat[,1] <- sapply(p1, function(p) {return(ifelse(p > 0, p, 0))})
    
    return(nmat)
    
    return(prob_mat)
}

V <- rbind(c(1, rep(0, length(EV) - 1)),
           c(1, rep(0, length(EV) - 1)),
           c(0, rep(0, length(EV) - 1)),
           c(0, rep(0, length(EV) - 1)))
L <- rbind(c(1, rep(0, length(EL) - 1)),
           c(0, rep(0, length(EL) - 1)),
           c(1, rep(0, length(EL) - 1)),
           c(0, rep(0, length(EL) - 1)))
L[,7:8] <- c(
    L[,1] * (1 - V[,1]),
    L[,1] * V[,1]
)
V[,1:10]
L

beta
# > beta["OFF_RACERBLACK"]
# OFF_RACERBLACK 
# 0.2268135

prop.table(table(dat$INCAR, dat$OFF_RACER))
prop.table(table(dat$INCAR))
prop.table(table(dat$OFF_RACER))

p_y <- 0.46585
p_x1 <- 0.269

a_y <- -log(1/p_y - 1) - 0.2268 * p_x1

# exp(0.2268)
# [1] 1.254579

lp_y_1 <- a_y + 0.2268
lp_y_0 <- a_y
p_y_1 <- (1 + exp(-lp_y_1))^{-1}
p_y_0 <- (1 + exp(-lp_y_0))^{-1}

calc_odds_ratio <- function(p_r) {
    p_11 <- p_r[1,1] * p_y_1
    p_10 <- p_r[3,1] * (1 - p_y_1)
    p_01 <- p_r[2,1] * p_y_0
    p_00 <- p_r[4,1] * (1 - p_y_0)
    
    odds1 <- (p_11 / p_10)
    odds0 <- (p_01 / p_00)
    # * exp(-0.2268)
    return(odds1 / odds0)
}

# No missing data
p_r <- matrix(rep(1,4), nrow = 4)
calc_odds_ratio(p_r)

# Over estimate
beta <- log(c(1, rep(1, length(EV) - 1)))
gamma <- log(c(1, rep(1, length(EL) - 3), 10, 1))

p_r <- p_miss_func(beta, gamma, EV, EL, V, L, pat_i)
p_r
calc_odds_ratio(p_r)

# Under Estimate
beta <- log(c(1, rep(1, length(EV) - 1)))
gamma <- log(c(1, rep(1, length(EL) - 3), 1, 10))

p_r <- p_miss_func(beta, gamma, EV, EL, V, L, pat_i)
p_r
calc_odds_ratio(p_r)

# 
beta <- log(c(10, rep(1, length(EV) - 1)))
gamma <- log(c(.1, rep(1, length(EL) - 3), 1, 1))

p_r <- p_miss_func(beta, gamma, EV, EL, V, L, pat_i)
p_r
calc_odds_ratio(p_r)

# 
beta <- log(c(10, rep(1, length(EV) - 1)))
gamma <- log(c(10, rep(1, length(EL) - 3), 1, 1))

p_r <- p_miss_func(beta, gamma, EV, EL, V, L, pat_i)
p_r
calc_odds_ratio(p_r)

