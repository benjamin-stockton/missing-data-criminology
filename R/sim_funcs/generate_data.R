#' Generate Data for Small Simulated Data Set
#'
#' @param beta A numeric vector of size 7.
#' @param N An integer for sample size.
#'
#' @return A data frame of size N x 11.
#'
#' @examples
#' generate_data(beta = c(.05, .227, 0, -.01, .9, -.2, 0, -.05),
#'               N = 1000)
generate_data <- function(
        beta = c(.05, .227, 0, -.01, .9, -.2, 0, -.05),
        N = 1000) {
    CTY <- sample(1:4, size = N, replace = TRUE, prob = c(.05, .25, .25, .45))
    X <- c(rbinom(N / 4, 1, .05),
           rbinom(N / 4, 1, .25), 
           rbinom(N / 4, 1, .35),
           rbinom(N / 4, 1, .65))
    C1 <- ifelse(CTY == 1, 1, 0)
    C2 <- ifelse(CTY == 2, 1, 0)
    C4 <- ifelse(CTY == 4, 1, 0)
    Z1 <- rnorm(N)
    Z1Q <- Z1^2
    Z2 <- rbinom(N, 1, .45)
    XC <- cbind(rep(1, N), X, Z1, Z1Q, Z2, C1, C2, C4)
    
    EY <- (1 + exp(-XC%*%beta))^(-1)
    Y <- rbinom(N, 1, EY)
    test.df <- data.frame(X = X, Y = Y, Z1 = Z1, Z1Q = Z1Q, Z2 = Z2, CTY = CTY)
    test.df$CTY <- factor(test.df$CTY)
    return(test.df)
}