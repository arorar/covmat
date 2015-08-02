#' Likelihood of Marchenko–Pastur distribution distribution
#' 
#' @details
#' This method calculates the negative likelihood of the eigenvalues given 
#' that they follow a Marchenko–Pastur distribution. The distribution is assumed
#' to have unit variance.
#' 
.neg.mpLogLik <- function(theta, lambdas) {
  
  gamma <- theta
  
  lambda.max <- (1 + sqrt(gamma))^2 
  lambda.min <- (1 - sqrt(gamma))^2
  
  lambdas <- lambdas[(lambdas < lambda.max) & (lambdas > lambda.min)]
  if(length(lambdas) == 0) return(.Machine$double.xmax)
  
  val <- sapply(lambdas,     
                function(x) dmp(x,svr = 1/gamma))
  
  ifelse(is.infinite(-sum(log(val))), .Machine$double.xmax, -sum(log(val)))        
}

#' Fitting an MP distribution to the data
#' 
#' @details
#' This method takes in the eigenvalues of the sample covariance matrix and fits 
#' a Marchenko–Pastur distribution. The eigenvalues are assumed to have unit
#' varianve
#' 
#' @param lambdas eigenvalues of the sample covariance matrix
#' 
.getMPfit <- function(lambdas) {
  
  lambda.max <- lambdas[which.max(lambdas <= 1) - 1]
  lower <- 0; upper <- (sqrt(lambda.max) - 1)^2
  
  fit <- DEoptim(fn = .neg.mpLogLik, 
               lower=lower, upper = upper,
               control = list(itermax = 500),
               lambdas = lambdas)  
    
  gamma <- fit$optim$bestmem
  lambda.max <- (1 + sqrt(gamma))^2
  spikes <- length(lambdas[lambdas > lambda.max])
  list(spikes = spikes, gamma = gamma, lambda.max = lambda.max)
} 

#' Eigenvalue shrinkage using Spiked Covariance Model
#' 
#' @details
#' This method takes in data as an xts object. It calculates a sample covariance
#' matrix and shrinks the eigenvalues based on the procedure listed in
#' (Donoho, Gavish, and Johnstone, 2013)
#' 
#' @param R xts object of asset returns
#' 
#' @author Rohit Arora
#' 
#' @export
#' 
#' 
estSpikedCovariance <- function(R, ...) {
  
  .data <- if(is.xts(R)) coredata(R) else as.matrix(R)
  T <- nrow(.data); M <- ncol(.data) 
  if (T < M) stop("Does not work when T < M")
  
  S <- cov(.data)
  eigen <- eigen(S, symmetric=T)
  lambdas <- eigen$values; scale.factor <- sd(lambdas)
  lambdas <- lambdas/scale.factor
  
  fit <- .getMPfit(lambdas)
  gamma <- fit$gamma; lambda.max <- fit$lambda.max; spikes <- fit$spikes
  
}