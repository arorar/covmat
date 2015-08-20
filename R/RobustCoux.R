#' Calculation of constant required for consistency
#' 
#' @details
#' The biweight rho function needs a constant whose value is determined based on
#' the consistency equation. Solving the integral and assuming multivariate
#' normality we can get an analytical expression for the constant.
#' 
#' @param p number of variables
#' @param c  square root of the chi-square cutoff
#' 
.gamma_cp <- function(p, c) {
  
  IncGamma <- function(a,z) {
    pgamma(z, a, lower.tail = FALSE)*gamma(a)
  }
  
  p*(c^6*gamma(p/2))/(gamma(p/2 + 1)*(6*c^4 -6*c^2*(p+2) + 2*(p+2)*(p+4)) - 
                        6*c^4*IncGamma(1 + p/2, c^2/2) + 
                        12*c^2*IncGamma(2 + p/2, c^2/2) - 
                        8*IncGamma(3 + p/2, c^2/2) + c^6*IncGamma(p/2, c^2/2))
}

#' Givens rotation matrix
#' 
#' @details
#' Givens rotation matrix is an orthogonal matrix with certain structure.
#' 
#' @param i row number
#' @param j column number
#' @param theta angle in the range [-pi/2, pi/2]
#' @param d number of columns
#' 
.givens.rotation <- function(i, j, theta, d) {
  
  G <- diag(rep(1, d))
  G[i,i] <-  G[j,j] <- cos(theta)
  G[j,i] <- sin(theta); G[i,j] <- -G[j,i]
  
  G
}

#' Product of Givens rotoation matrices
#' 
#' @details
#' An orthogonal matrix of size d can be parameterized by d(d-1)/2 angles
#' 
#' @param d number of columns
#' @param theta angle in the range [-pi/2, pi/2]
#' 
.orthogonal.matrix <- function(d, angles) {
  
  O <- diag(rep(1, d)); k <- 1
  
  for(i in 2:d) 
    for(j in 1:(i-1)) {
      O <- O %*% .givens.rotation(i, j, theta = angles[k], d)
      k <- k + 1
    }
  O
}


#' Smoothing matrix
#' 
#' @details
#' The smooting matrix can be parametrized using the d eigenvalues and d(d-1)/2
#' angles to represent the orthogonal matrix when carrying out a spectral
#' decomposition
#' 
#' @param d number of columns
#' @param params a vector of d(d+1)/2 parameters. First d are the eigenvalues and
#'          remaining d(d-1)/2 are the angles.
#' 
.smoothing.matrix <- function(params, d) {
  end <- d*(d+1)/2
  theta.lambda <- params[1:d]; theta.angles <- params[(d+1):end]
 
  lambdas <- theta.lambda; angles <- theta.angles 
  W <- .orthogonal.matrix(d, angles)
  N <- W %*% diag(lambdas) %*% t(W)
  
  N
}

#' Huber function
#' 
#' @details
#' Calculate the output of a huber function given the value and the number of 
#' variables
#' 
#' @param x the value
#' @param p number of columns
.huber <- function(x, p){
  k <- sqrt(qchisq(0.95, df=p))
  min(k, max(x, -k))
}

#' Biweight Rho function
#' 
#' @details
#' Calculate the output of a biweight rho function given the value and the number of 
#' variables
#' 
#' @param x the value
#' @param p number of columns
#' 
.biweight <- function(x, p){
  
  c <- sqrt(qchisq(0.95, df=p))
  g <- .gamma_cp(p, c)
  
  if(abs(x) <= c) g*(1 - (1 - (x/c)^2)^3)
  else g
}

#' Objective to calculate the optimal smoothing matrix
#' 
#' @details
#' The objective calculates the determinant of the one step ahead forecast errors.
#' The objective is a very noisy function and hence the optimization is time 
#' consuming and replicating paramaters across runs is difficult
#' 
#' @param params d(d+1)/2 paramaters to be optimzed
#' @param R data
#' @param y.hat fitted yalues for the data
#' @param Sigma.hat value of the covarianve matrix
#' @param startup_period length of samples required to calculate initial values
#' @param training_period length of samples required to calculate forecast errors
#'                for evalualating the objective
#' @param lambda known constant as described in the paper
#' 
.obj <- function(params, R, y.hat, Sigma.hat, startup_period, training_period, 
                 lambda = 0.2) {
  
  d <- ncol(R); I <- diag(rep(1, d))
  smoothing.matrix <- .smoothing.matrix(params, d);
  
  forecast.error <- matrix(NA, nrow = training_period - startup_period, 
                           ncol = d)
  
  for(t in (startup_period + 1):training_period) {
    
    forecast.error[t - startup_period,] <- R[t,] - y.hat[t-1,]
    ferr <- t(forecast.error[t - startup_period,,drop=FALSE])
    
    I.Sigma.hat <- solve(Sigma.hat)
    temp <- as.numeric(t(ferr) %*% I.Sigma.hat %*% ferr)
    Sigma.hat <- lambda*.biweight(sqrt(temp), d)/temp * ferr %*% t(ferr) + 
                  (1 - lambda)*Sigma.hat
    
    temp <- as.numeric(sqrt(t(ferr) %*% solve(Sigma.hat) %*% ferr))
    cleaned.val <- .huber(temp, d)/temp*ferr + t(y.hat[t-1,,drop=FALSE])
      
    y.hat[t,] <- smoothing.matrix %*% cleaned.val + 
      (I - smoothing.matrix) %*% t(y.hat[t-1,,drop=FALSE])
  }
  
  h <- floor(0.75 * (training_period - startup_period))
  Sigma.hat <- covMcd(forecast.error, nsamp = h)$cov

  det(Sigma.hat)
}


#' Optimal Smoothing Matrix
#' 
#' @details
#' Calcuation of smoothing matrix is done by assuming that the smoothing matrix
#' is symmetrix and has a spectral decomposition. The orthogonal matrix in the 
#' decomposition is calculated using the product of givens rotation matrices and
#' requires d(d-1)/2 angles for a d dimensional matrix. The eigenvalues are 
#' restricted to lie in [0,1].
#' 
#' @param R data
#' @param startup_period length of samples required to calculate initial values
#' @param training_period length of samples required to calculate forecast errors
#'                for evalualating the objective
#' @param seed random seed to replicate the starting values for optimization
#' @param trials number of strarting values to try for any optimization. 
#'            Large number of trials for high dimensions can be time consuming
#' @param method optimization method to use to evaluate an estimate of 
#'            smoothing matrix. Default is L-BFGS-B
#' 
#' @export
#' @author Rohit Arora
#' 
smoothing.matrix <- function(R, startup_period = 10, training_period = 60 , 
                             seed = 9999, trials = 50, method = "L-BFGS-B") {
  
  M <- nrow(R); d <- ncol(R)
  if(M < 4*d) stop("Not enough data for estimation")
  
  startup_period <- startup_period[1]
  if(is.na(startup_period) || startup_period < 2*d) startup_period <- 2*d
  
  training_period <- training_period[1]
  if(is.na(training_period) || training_period < (startup_period + 2*d)) 
    training_period <- startup_period + 2*d
  
  if ( M < (startup_period + training_period)) 
    stop("Insufficienct data. Reset correct startup & training periods")
  
  startup.fit <- lapply(1:d, function(i) {
    lmRob(coredata(R[1:startup_period,i]) ~ as.matrix(1:startup_period))
  })
  
  y.hat <- matrix(NA, nrow = training_period, ncol = ncol(R))
  y.hat[1:startup_period,] <- do.call(cbind, lapply(startup.fit, fitted))
  
  res <- do.call(cbind, lapply(startup.fit, residuals))  
  Sigma.hat <- covMcd(res)$cov
  
  set.seed(seed)
  
  lower <- c(rep(0,d), rep(-pi/2, d*(d-1)/2))
  upper <- c(rep(1,d), rep(pi/2, d*(d-1)/2))
  nlower <- length(lower); width <- upper - lower
  
  Umin <- matrix(rep.int(lower, trials), nrow = trials, ncol=nlower, byrow=T)
  start <- (Umin + matrix(rep.int(width, trials), nrow = trials, 
                      ncol=nlower, byrow=T)*maximinLHS(n = trials, k = nlower))
  
  cl <- makeCluster(detectCores())
  clusterExport(cl, list = c("lower", "upper", ".obj", "optimx","R"), 
                envir = environment())
  registerDoSNOW(cl)  
  
  objmin <- parRapply(cl, start, function (x)
    try(optimx(x, .obj, lower = lower, upper = upper, method = method, R = R , 
               y.hat = y.hat, Sigma.hat = Sigma.hat, startup_period = startup_period, 
             training_period = training_period), silent=TRUE))
  
  fit <- objmin[[unique(which.min(parSapply(cl, objmin, '[[', "value")))]]
  
  stopCluster(cl)
  registerDoSEQ()
  
  params <- unlist(fit[1:(d*(d+1)/2)])
  N <- .smoothing.matrix(params, d)
  list(smooth.mat = N, fit = fit, all.fit = objmin)
}