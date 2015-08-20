#' @export
gamma_cp <- function(p, c) {
  
  IncGamma <- function(a,z) {
    pgamma(z, a, lower=FALSE)*gamma(a)
  }
  
  p*(c^6*gamma(p/2))/(gamma(p/2 + 1)*(6*c^4 -6*c^2*(p+2) + 2*(p+2)*(p+4)) - 
                        6*c^4*IncGamma(1 + p/2, c^2/2) + 
                        12*c^2*IncGamma(2 + p/2, c^2/2) - 
                        8*IncGamma(3 + p/2, c^2/2) + c^6*IncGamma(p/2, c^2/2))
}

#' @export
givens.rotation <- function(i, j, theta, d) {
  
  G <- diag(rep(1, d))
  G[i,i] <-  G[j,j] <- cos(theta)
  G[j,i] <- sin(theta); G[i,j] <- -G[j,i]
  
  G
}

#' @export
orthogonal.matrix <- function(d, angles) {
  
  O <- diag(rep(1, d)); k <- 1
  
  for(i in 2:d) 
    for(j in 1:(i-1)) {
      O <- O %*% givens.rotation(i, j, theta = angles[k], d)
      k <- k + 1
    }
  O
}

#' @export
.smoothing.matrix <- function(params, d) {
  end <- d*(d+1)/2
  theta.lambda <- params[1:d]; theta.angles <- params[(d+1):end]
 
  lambdas <- theta.lambda; angles <- theta.angles 
#  lambdas <- (1 + tanh(theta.lambda))/2
#  angles <-  atan(theta.angles)
  
  W <- orthogonal.matrix(d, angles)
  N <- W %*% diag(lambdas) %*% t(W)
  
  N
}

#' @export
.huber <- function(x, p){
  k <- sqrt(qchisq(0.95, df=p))
  min(k, max(x, -k))
}

#' @export
.biweight <- function(x, p){
  
  c <- sqrt(qchisq(0.95, df=p))
  g <- gamma_cp(p, c)
  
  if(abs(x) <= c) g*(1 - (1 - (x/c)^2)^3)
  else g
}

#' @export
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


#' @export
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