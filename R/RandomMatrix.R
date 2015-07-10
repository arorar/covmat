library(Matrix)
library(xts)
library(ggplot2)

#'Plots the eigenvalues of the correlation matrix and overlays the Marchenko Pastur density
#' 
#' @details
#' There is a shap cutoff for the density. We are concerned with eigenvalues beyond
#' this cutoff
#' 
#' @author Rohit Arora
#' 
#' 
eigen.plot <- function(lambda, Q, sigma.sq){
    
    lambda.max <- sigma.sq*(1 + 1/Q + 2*sqrt(1/Q)) 
    
    p <- ggplot(data=data.frame(lambda)) + 
        geom_histogram( aes(x = lambda, y=..density..),
                        breaks=seq(min(lambda)-1,1+max(lambda),0.5), 
                        colour="black", fill="white") +
        stat_function(fun = dmp, args=list(svr = Q, var=sigma.sq), 
                      aes(colour = 'MP density')) + xlab("Eigenvalues") +
        labs(title="Actual vs Fitted Marchenko-Pastur") + ylim(0,1.5) + 
        theme(plot.title = element_text(size = 20, face = "bold", vjust = 1),
              axis.title=element_text(size=14,face="bold"))
    
    print(p)
    p
}

#' Implement Denoising of Covariance matrix using Random matrix theory
#' 
#' @details
#' This method takes in data as a matrix or an xts object. It then
#' fits a marchenko pastur density to eigenvalues of the correlation matrix. All
#' eigenvalues above the cutoff are retained and ones below the cutoff are
#' retained such that the trace of the correlation matrix is 1. Finally, 
#' correlation matrix is converted to covariance matrix.
#' 
#' @param  R xts or matrix of asset returns
#' @author Rohit Arora
#' 
#' @export
#' 
rmt.est <- function(R) {
    .data <- as.matrix(R)
    T <- nrow(.data); M <- ncol(.data) 
    if (T < M) stop("Does not work when T < M")
    
    #eigenvalues can be negative. To avoid this e need a positive-definite matrix 
    S <- cov(.data); S <- nearPD(S)$mat 
    D <- diag(S); C <- cov2cor(S); 
    
    # Marchenko Pastur density is defined for eigenvalues of correlation matrix
    eigen.C <- eigen(C,symmetric=T)
    lambda <- eigen.C$values; sigma.sq <- mean(lambda)
    
    #minimize log-likelihood. 
    loglik.marpas <- function(theta) {
        
        Q <- theta[1]; sigma.sq <- theta[2]
        
        lambda.max <- sigma.sq*(1 + 1/Q + 2*sqrt(1/Q)) 
        lambda.min <- sigma.sq*(1 + 1/Q - 2*sqrt(1/Q))
        
        lambda.tmp <- lambda[lambda < lambda.max & lambda > lambda.min]
        val <- sapply(lambda.tmp,     
                      function(x) dmp(x,svr = Q, var=sigma.sq))
        
        -sum(log(val))        
    }
    
    # these paramters for optimization are slightly arbitrary
    start <- c(T/M,1); lb <- c(1,1); ub <- c(5, var(lambda))
    fit.marpas <- optim(par = start, fn = loglik.marpas, method = "L-BFGS-B", 
                        lower = lb, upper = ub)
    
    Q <- fit.marpas$par[1]; sigma.sq <- fit.marpas$par[2]
    
    lambda.max <- sigma.sq*(1 + 1/Q + 2*sqrt(1/Q))  
    lambda.min <- sigma.sq*(1 + 1/Q - 2*sqrt(1/Q))
    eigen.plot(lambda, Q, sigma.sq)
    
    # now that we have a fit. lets denoise eigenvalues below the cutoff
    idx <- which(lambda > lambda.max)
    if (length(idx) == 0) return(S)
    
    val <- eigen.C$values[idx]; vec <- eigen.C$vectors[,idx,drop=FALSE]
    sum <- 0; for (i in idx) sum <- sum + val[i]*vec[,i] %*% t(vec[,i])
    
    # trace of correlation matrix is 1. Use this to determine all the remaining
    # eigenvalues
    clean.C <- sum + sum(eigen.C$values[-idx])/M * diag(rep(1,M))
    
    # convert correlation to covariance matrix and return
    clean.S <- diag(D)^0.5 %*% clean.C %*% diag(D)^0.5
    clean.S
}

#test
data <- read.zoo("./data/all_sp500_price_data.csv", header=TRUE)
data <- diff(log(na.omit(data)))[,1:80]
denoised.cov.mat <- rmt.est(data)
