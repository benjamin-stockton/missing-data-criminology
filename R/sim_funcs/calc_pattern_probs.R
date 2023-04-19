calc_pattern_probs <- function(V, L, beta, gamma, EV, EL) {
    pat.probs <- c(834548, 26206, 2078, 34, 1465, 92, 1, 1) / 864422
    alpha <- numeric(8)
    prob.mat <- matrix(numeric(8 * nrow(V)), ncol = 8)
    
    for (i in 1:8) {
        alpha[i] <- -log(1 / pat.probs[i] - 1) - (EV %*% beta[i,] + EL %*% gamma[i,])[,1]
        lp_i <- alpha[i] + (V %*% beta[i,] + L %*% gamma[i,])[,1]
        p_i <- exp(lp_i)
        prob.mat[,i] <- p_i
    }
    prob.mat[,1] <- rep(1, nrow(prob.mat))
    r.sum <- apply(prob.mat, 1, sum)
    
    nmat <- t(apply(matrix(1:nrow(prob.mat), ncol = 1), 1, function(r) {
        return(prob.mat[r,] / r.sum[r])
    }))
    # nmat[,1] <- sapply(p1, function(p) {return(ifelse(p > 0, p, 0))})
    
    return(nmat)
}