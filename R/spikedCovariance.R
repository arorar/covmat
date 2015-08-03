#' Optimization problem 6.2
#' 
#' @details
#' Some of the shrinkers do not have an analytical solution and the non-lineartiy
#' must be computed numerically. This function implements
#' Optimization problem 6.2 as described in (Donoho, Gavish, and Johnstone, 2013)
#' 
#' @param lambda sample eigenvalue
#' @param gamma  fitted value of variables/observations
#' @param type   combination of (pivot+norm)/(statistical loss)
#' 
.shrink.eigen2 <- function(lambda, c, s, type) {
  
  temp <- unlist(strsplit(type,"\\."))
  
  isStatistical <- FALSE; mnorm <- ""; pivot <- -1
  if (length(temp) == 2) {
    mnorm <- temp[1]
    pivot <- as.numeric(temp[2])
  } else { 
    isStatistical <- TRUE
  }
  
  .obj <- function(eta) {
    
    A <- matrix(c(lambda, 0, 0, 1), nrow=2, byrow = TRUE)
    B <- matrix(c(1 + (eta -1)*c^2, (eta -1)*c*s, 
                  (eta -1)*c*s, 1 + (eta -1)*s^2), nrow=2, byrow = TRUE)
  
    if(isStatistical) {
      0.5*log(0.5*det(A + B)/sqrt(det(A)*det(B)))
    } else {
    
      M <- 
        if (pivot == 3) solve(A) %*% B - diag(rep(1,2))
        else if (pivot == 4) solve(B) %*% A - diag(rep(1,2))
        else if (pivot == 5) solve(A) %*% B + solve(B) %*% A - 2*diag(rep(1,2))
        else if (pivot == 7) log(sqrt(solve(A)) %*% B %*% sqrt(solve(A)))
      
      if (mnorm == "Frobenius") norm(M, type = "F")
      else if (mnorm == "Operator") sqrt(eigen(M%*%t(M))$values[1])
      else if (mnorm == "Nuclear") sum(diag(sqrt(t(M) %*% M)))  
    }
  }
  
  fit <- DEoptim(fn = .obj, 
                 lower=1, upper = 500,
                 control = list(itermax = 500, trace = 0))  
  
  fit$optim$bestmem
}

#' Table 2
#' 
#' @details
#' Table 2 as described in (Donoho, Gavish, and Johnstone, 2013)
#' 
#' @param lambda sample eigenvalue
#' @param gamma  fitted value of variables/observations
#' @param type   combination of (pivot+norm)/(statistical loss)
#' 
.shrink.eigen <- function(lambda, gamma, type) {
  
  c <- .c.lambda(lambda, gamma)
  s <- .s.lambda(lambda, gamma)
  
  if (type == "Frobenius.1") lambda*c^2 + s^2
  else if (type == "Frobenius.2") lambda/(c^2 + lambda*s^2)
  else if (type == "Frobenius.3") (lambda*c^2 + lambda^2*s^2)/(c^2 + lambda^2*s^2)
  else if (type == "Frobenius.4") (lambda^2*c^2 + s^2)/(lambda*c^2 + s^2)
  else if (type == "Frobenius.5") .shrink.eigen2(lambda, c, s, type)
  else if (type == "Frobenius.6") c^2*(lambda-1)/(c^2 + lambda*s^2)^2
  else if (type == "Frobenius.7") .shrink.eigen2(lambda, c, s, type)
  
  else if (type == "Operator.1") lambda
  else if (type == "Operator.2") lambda
  else if (type == "Operator.3") .shrink.eigen2(lambda, c, s, type)
  else if (type == "Operator.4") .shrink.eigen2(lambda, c, s, type)
  else if (type == "Operator.5") .shrink.eigen2(lambda, c, s, type)
  else if (type == "Operator.6") (lambda-1)/(c^2 + lambda*s^2)
  else if (type == "Operator.7") .shrink.eigen2(lambda, c, s, type)
  
  else if (type == "Nuclear.1") max(1, 1 + (lambda -1)*(1 - 2*s^2))
  else if (type == "Nuclear.2") max(1, lambda/(c^2 + (2*lambda -1)*s^2))
  else if (type == "Nuclear.3") max(1, lambda/(c^2 + lambda^2*s^2))
  else if (type == "Nuclear.4") max(1, (lambda^2*c^2 + s^2)/lambda)
  else if (type == "Nuclear.5") .shrink.eigen2(lambda, c, s, type)
  else if (type == "Nuclear.6") max(1, (lambda-(lambda-1)^2*c^2*s^2)/(c^2 + lambda*s^2)^2)
  else if (type == "Nuclear.7") .shrink.eigen2(lambda, c, s, type)

  else if (type == "Stein") lambda/(c^2 + lambda*s^2)
  else if (type == "Entropy") lambda*c^2 + s^2
  else if (type == "Divergence") sqrt((lambda^2*c^2 + lambda*s^2)/(c^2 + lambda*s^2))
  else if (type == "Affinity") .shrink.eigen2(lambda, c, s, type)
  else if (type == "Frechet") (sqrt(lambda)*c^2 + s^2)^2
  
  else NA
}

#' Equation 4.4
#' 
#' @details
#' Equation 4.4 as described in (Donoho, Gavish, and Johnstone, 2013)
#' 
#' @param lambda sample eigenvalue
#' @param gamma  fitted value of variables/observations
#' 
#' 
.ell.lambda <- function(lambda, gamma) {
  temp <- lambda + 1 - gamma
  0.5*(temp + sqrt(temp^2 - 4*lambda))
}

#' Equation 1.2
#' 
#' @details
#' Equation 1.2 as described in (Donoho, Gavish, and Johnstone, 2013)
#' 
#' @param lambda sample eigenvalue
#' @param gamma  fitted value of variables/observations
#' 
#' 
.c.lambda <- function(lambda, gamma) {
  temp <- .ell.lambda(lambda, gamma)
  sqrt((1 - gamma/(temp - 1)^2)/(1 + gamma/(temp - 1)))
}

#' Equation 6.3
#' 
#' @details
#' Equation 6.3 as described in (Donoho, Gavish, and Johnstone, 2013)
#' 
#' @param lambda sample eigenvalue
#' @param gamma  fitted value of variables/observations
#' 
.s.lambda <- function(lambda, gamma) {
  sqrt(1 - (c.lambda(lambda, gamma))^2)
}

#' Likelihood of Marchenko–Pastur distribution distribution
#' 
#' @details
#' This method calculates the negative likelihood of the eigenvalues given 
#' that they follow a Marchenko–Pastur distribution. The distribution is assumed
#' to have unit variance.
#' 
#' @param theta parameter to be optimized. In this case it is 
#'          gamma, (variables/observations)
#' @param lambdas eigenvalues to which the distribution is fitted
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
#' variance
#' 
#' @param lambdas eigenvalues of the sample covariance matrix
#' 
.getMPfit <- function(lambdas) {
  
  lambda.max <- lambdas[which.max(lambdas <= 1) - 1]
  lower <- 0; upper <- (sqrt(lambda.max) - 1)^2
  
  fit <- DEoptim(fn = .neg.mpLogLik, 
               lower=lower, upper = upper,
               control = list(itermax = 500, trace = 0),
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
#' @param statistical Stein/Entropy/Divergence/Affinity/Frechet. Default is set to NA.
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
     !statistical %in% c("Stein","Entropy","Divergence","Affinity","Frechet"))
    stop("Invalid statistical parameter selected")
  
  S <- cov(.data)
  eigen <- eigen(S, symmetric=T)
  lambdas <- eigen$values; scale.factor <- sd(lambdas)
  lambdas <- lambdas/scale.factor
  
  fit <- .getMPfit(lambdas)
  gamma <- fit$gamma; lambda.max <- fit$lambda.max; spikes <- fit$spikes
  
  spiked.lambdas <- lambdas[lambdas > lambda.max]
  
  type <- ifelse(is.na(statistical), paste(norm,".",pivot,sep=""), statistical)

  spiked.lambdas <- sapply(spiked.lambdas, 
                           function(lambda) 
                             .shrink.eigen(lambda, gamma, type))
  
  lambdas[lambdas > lambda.max] <- spiked.lambdas
  lambdas[lambdas <= lambda.max] <- 1
  
  C <- eigen$vectors %*% diag(lambdas) %*% t(eigen$vectors)
  list(cov = scale.factor*C)
}