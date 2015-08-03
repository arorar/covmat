#' Table 2
#' 
#' @details
#' Table 2
#' 
Frobenius.1 <- function(lambda, gamma) {
  c <- c.lambda(lambda, gamma)
  s <- s.lambda(lambda, gamma)
  lambda*c^2 + s^2
}

#' Equation 4.4
#' 
#' @details
#' Equation 4.4
#' 
ell.lambda <- function(lambda, gamma) {
  temp <- lambda + 1 - gamma
  0.5*(temp + sqrt(temp^2 - 4*lambda))
}

#' Equation 1.2
#' 
#' @details
#' Equation 1.2
#' 
c.lambda <- function(lambda, gamma) {
  temp <- ell.lambda(lambda, gamma)
  sqrt((1 - gamma/(temp - 1)^2)/(1 + gamma/(temp - 1)))
}

#' Equation 6.3
#' 
#' @details
#' Equation 1.2
#' 
s.lambda <- function(lambda, gamma) {
  sqrt(1 - (c.lambda(lambda, gamma))^2)
}

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
#' @param norm Type of matrix norm that must be calculated. Defaults to Frobenius
#' @param pivot takes values from 1...7. Details can be found in the paper
#' @param statistical Stein/Entropy/Bhattarcharya/Frechet. Default is set to NA.
#'        when a valid value is set norm and pivot values are ignored
#' 
#' @author Rohit Arora
#' 
#' @export
#' 
#' 
estSpikedCovariance <- function(R, norm = c("Frobenius", "Operator", "Nuclear"),
                                pivot = 1, statistical = NA) {
  
  .data <- if(is.xts(R)) coredata(R) else as.matrix(R)
  T <- nrow(.data); M <- ncol(.data) 

  if (T < M) stop("Does not work when T < M")
  
  norm <- norm[1]
  if (!norm %in% c("Frobenius", "Operator", "Nuclear"))
    stop("Invalid norm value")
  
  if (pivot < 0 || pivot > 7) stop("Invalid pivot selected")
  if(!is.na(statistical) && 
     !statistical %in% c("Stein","Entropy","Bhattarcharya","Frechet"))
    stop("Invalid statistical parameter selected")
  
  S <- cov(.data)
  eigen <- eigen(S, symmetric=T)
  lambdas <- eigen$values; scale.factor <- sd(lambdas)
  lambdas <- lambdas/scale.factor
  
  fit <- .getMPfit(lambdas)
  gamma <- fit$gamma; lambda.max <- fit$lambda.max; spikes <- fit$spikes
  
  spiked.lambdas <- lambdas[lambdas > lambda.max]
  
  if(is.na(statistical)) {
    fun.name <- paste(norm,".",pivot,sep="")
    loss.func <- match.fun(fun.name)
    formals(loss.func)$gamma <- gamma
    spiked.lambdas <- sapply(spiked.lambdas, loss.func)
    
    lambdas[lambdas > lambda.max] <- spiked.lambdas
    lambdas[lambdas <= lambda.max] <- 1
  }
  
  C <- eigen$vectors %*% diag(lambdas) %*% t(eigen$vectors)
  scale.factor*C
}