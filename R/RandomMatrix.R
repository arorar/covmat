library(Matrix)
library(xts)
library(ggplot2)
library(RMTstat)

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
              axis.title=element_text(size=14,face="bold")) + 
        annotate('text', x = 10, y = 0.9, 
                 label = paste("sigma^{2} == ", round(sigma.sq,3)), 
                 parse=TRUE, size = 8) +
        annotate('text', x = 10, y = 1, 
                 label = paste("Q == ", round(Q,3)), parse=TRUE, , size = 8) + 
        annotate('text', x = 10, y = 0.78, 
                 label = paste("lambda[max] ==", round(lambda.max,3)), 
                 parse=TRUE, , size = 8) + 
        scale_colour_manual("", values = c("red"))
    
    options(warn = -1)
    print(p)
    options(warn = 0)
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
rmt.est <- function(R, numEig=1) {
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
    loglik.marpas <- function(theta, sigma.sq) {
        
        Q <- theta
        val <- sapply(lambda,     
                      function(x) dmp(x,svr = Q, var=sigma.sq))
        
        val <- val[val > 0]
        ifelse(is.infinite(-sum(log(val))), .Machine$double.xmax, -sum(log(val)))        
    }
    
    sigma.sq <- 1 - sum(head(lambda,numEig))/M
    
    lb <- 1; ub <- T/M
    cl <- makeCluster(detectCores())
    registerDoSNOW(cl)
    clusterEvalQ(cl, library(RMTstat))
    
    starts <- seq(lb, ub, length.out = 50)
    fit.marpas <- foreach(start = starts, .combine = rbind) %dopar% 
        optim(par = start, fn = loglik.marpas, method = "L-BFGS-B", 
                            lower = lb, upper = ub, sigma.sq = sigma.sq)    
    stopCluster(cl)

    idx <- grep("CONVERGENCE",unlist(fit.marpas[,"message"]))
    vals <- fit.marpas[idx,c("par","value")]
    Q <- unlist(vals[which.min(vals[,"value"]),"par"])
    
    lambda.max <- qmp(1, svr=Q, var = sigma.sq)  
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

cov.mat <- rmt.est(data)
