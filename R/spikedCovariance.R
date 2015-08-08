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
#' @param lambda sample eigenvalue to be shrunk
#' 
.shrink.eigen2 <- function(ell, c, s, type, lambda) {
  
  temp <- unlist(strsplit(type,"\\."))
  
  isStatistical <- FALSE; mnorm <- ""; pivot <- -1
  if (length(temp) == 2) {
    mnorm <- temp[1]
    pivot <- as.numeric(temp[2])
  } else { 
    isStatistical <- TRUE
  }
  
  .obj <- function(eta) {
    
    A <- matrix(c(ell, 0, 0, 1), nrow=2, byrow = TRUE)
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
      else if (mnorm == "Nuclear") sum(sqrt(diag(t(M) %*% M)))  
    }
  }
  
  fit <- DEoptim(fn = .obj, 
                 lower=1, upper = lambda,
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
  
  ell <- .ell.lambda(lambda, gamma)
  c <- .c.ell(lambda, gamma)
  s <- .s.ell(lambda, gamma)
  
  if (type == "Frobenius.1") ell*c^2 + s^2
  else if (type == "Frobenius.2") ell/(c^2 + ell*s^2)
  else if (type == "Frobenius.3") (ell*c^2 + ell^2*s^2)/(c^2 + ell^2*s^2)
  else if (type == "Frobenius.4") (ell^2*c^2 + s^2)/(ell*c^2 + s^2)
  else if (type == "Frobenius.5") .shrink.eigen2(ell, c, s, type, lambda)
  else if (type == "Frobenius.6") c^2*(ell-1)/(c^2 + ell*s^2)^2
  else if (type == "Frobenius.7") .shrink.eigen2(ell, c, s, type, lambda)
  
  else if (type == "Operator.1") ell
  else if (type == "Operator.2") ell
  else if (type == "Operator.3") .shrink.eigen2(ell, c, s, type, lambda)
  else if (type == "Operator.4") .shrink.eigen2(ell, c, s, type, lambda)
  else if (type == "Operator.5") .shrink.eigen2(ell, c, s, type, lambda)
  else if (type == "Operator.6") (ell-1)/(c^2 + ell*s^2)
  else if (type == "Operator.7") .shrink.eigen2(ell, c, s, type, lambda)
  
  else if (type == "Nuclear.1") max(1, 1 + (ell -1)*(1 - 2*s^2))
  else if (type == "Nuclear.2") max(1, ell/(c^2 + (2*ell -1)*s^2))
  else if (type == "Nuclear.3") max(1, ell/(c^2 + ell^2*s^2))
  else if (type == "Nuclear.4") max(1, (ell^2*c^2 + s^2)/ell)
  else if (type == "Nuclear.5") .shrink.eigen2(ell, c, s, type, lambda)
  else if (type == "Nuclear.6") max(1, (ell-(ell-1)^2*c^2*s^2)/(c^2 + ell*s^2)^2)
  else if (type == "Nuclear.7") .shrink.eigen2(ell, c, s, type, lambda)

  else if (type == "Stein") ell/(c^2 + ell*s^2)
  else if (type == "Entropy") ell*c^2 + s^2
  else if (type == "Divergence") sqrt((ell^2*c^2 + ell*s^2)/(c^2 + ell*s^2))
  else if (type == "Affinity") .shrink.eigen2(ell, c, s, type, lambda)
  else if (type == "Frechet") (sqrt(ell)*c^2 + s^2)^2
  
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
.c.ell <- function(lambda, gamma) {
  ell <- .ell.lambda(lambda, gamma)
  sqrt((1 - gamma/(ell - 1)^2)/(1 + gamma/(ell - 1)))
}

#' Equation 6.3
#' 
#' @details
#' Equation 6.3 as described in (Donoho, Gavish, and Johnstone, 2013)
#' 
#' @param lambda sample eigenvalue
#' @param gamma  fitted value of variables/observations
#' 
.s.ell <- function(lambda, gamma) {
  sqrt(1 - (.c.ell(lambda, gamma))^2)
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
.neg.mpLogLik <- function(theta, lambdas, numOfSpikes) {
  
  gamma <- theta
  scale.factor <- lambdas[numOfSpikes + 1]/(1 + sqrt(gamma))^2

  lambda.max <- scale.factor*(1 + sqrt(gamma))^2
  lambda.min <- scale.factor*(1 - sqrt(gamma))^2
  
  spikes <- length(lambdas[lambdas > lambda.max])
  if(spikes != numOfSpikes) return(.Machine$double.xmax)
    
  lambdas <- lambdas[(lambdas <= lambda.max) & (lambdas >= lambda.min)]
  if(length(lambdas) == 0) return(.Machine$double.xmax)
  
  val <- sapply(lambdas,     
                function(x) dmp(x,svr = 1/gamma, var = scale.factor))
  
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
#' @param gamma   ratio of varibales/observations
#' @param numOfSpikes number of spikes in the spike covariance model
#' 
.getMPfit <- function(lambdas, gamma, numOfSpikes) {
  
  if(is.na(gamma)) {
    lambda.min <- min(lambdas)
    lower <- 0; upper <- 1
    
    fit <- DEoptim(fn = .neg.mpLogLik, 
                   lower=lower, upper = upper,
                   control = list(itermax = 1000, trace = 0),
                   lambdas = lambdas, numOfSpikes = numOfSpikes)  
    
    gamma <- fit$optim$bestmem  
  }
  
  scale.factor <- lambdas[numOfSpikes + 1]/(1 + sqrt(gamma))^2
  lambdas <- lambdas/scale.factor
  
  lambda.max <- (1 + sqrt(gamma))^2
  spikes <- length(lambdas[lambdas > lambda.max])
  if(spikes != numOfSpikes)
    browser()
  
  list(lambda.max = lambda.max, lambdas = lambdas, gamma = gamma, 
       scale.factor = scale.factor, numOfSpikes = numOfSpikes)
} 

#' (Donoho, Gavish, and Johnstone, 2013)
#' 
#' @param R xts object of asset returns
#' @param gamma  ratio of varibales/observations. If NA it will be estimated from
#'                 the data. One can set it to variables/observations
#' @param numOfSpikes number of spikes in the spike covariance model. 
#'                    It is defaulted to 1
#' @param norm Type of matrix norm that must be calculated. Defaults to Frobenius
#' @param pivot takes values from 1...7. Details can be found in the paper
#' @param statistical Stein/Entropy/Divergence/Affinity/Frechet. Default is set to NA.
#'                    when a valid value is set norm and pivot values are ignored
#' @param fit list with 5 elements, cutoff for the bulk of MP distribution, 
#'            scaled lambdas, fitted gamma and fitted scaling constant,
#'            numOfSpikes
#' 
#' @author Rohit Arora
#' 
#' @export
#' 
#' 
estSpikedCovariance <- function(R, gamma = NA, numOfSpikes = 1,
                                norm = c("Frobenius", "Operator", "Nuclear"),
                                pivot = 1, statistical = NA,
                                fit = NA) {
  
  .data <- if(is.xts(R)) coredata(R) else as.matrix(R)
  T <- nrow(.data); M <- ncol(.data) 

  if (T < M) stop("Does not work when T < M")
  
  if((!is.na(gamma)) && (gamma > 1 || gamma < 0)) stop("Invalid gamma")
  if(is.na(numOfSpikes) || (numOfSpikes > M)) stop("Invalid numOfSpikes")
  
  if ("gamma" %in% names(fit) && (!is.na(gamma)) && (fit$gamma != gamma))
    stop("Conflicting gamma values")
  
  if ("numOfSpikes" %in% names(fit) && (!is.na(numOfSpikes)) && 
      (fit$numOfSpikes != numOfSpikes))
    stop("Conflicting numOfSpikes values")
  
  norm <- norm[1]
  if (!norm %in% c("Frobenius", "Operator", "Nuclear"))
    stop("Invalid norm value")
  
  if (pivot < 0 || pivot > 7) stop("Invalid pivot selected")
  if(!is.na(statistical) && 
     !statistical %in% c("Stein","Entropy","Divergence","Affinity","Frechet"))
    stop("Invalid statistical parameter selected")
  
  S <- cov(.data)
  eigen <- eigen(S, symmetric=T)
  lambdas <- eigen$values
  if(all(is.na(fit))) fit <- .getMPfit(lambdas, gamma, numOfSpikes)
  
  lambda.max <- fit$lambda.max; lambdas <- fit$lambdas; 
  gamma <- fit$gamma; scale.factor <- fit$scale.factor

  spiked.lambdas <- lambdas[lambdas > lambda.max]
  
  type <- ifelse(is.na(statistical), paste(norm,".",pivot,sep=""), statistical)

  spiked.lambdas <- sapply(spiked.lambdas, 
                           function(lambda) 
                             .shrink.eigen(lambda, gamma, type))
  
  lambdas[lambdas > lambda.max] <- spiked.lambdas
  lambdas[lambdas <= lambda.max] <- 1
  
  rescaled.lambdas <- scale.factor * lambdas
  
  C <- eigen$vectors %*% diag(rescaled.lambdas) %*% t(eigen$vectors)
  model <- list(cov = C, data = R, orig.lambdas = eigen$values, pivot = pivot,
                numOfSpikes = numOfSpikes, norm = norm, new.lambdas = rescaled.lambdas, 
                dist.fit = fit)
  
  class(model) <- "spikedCovariance"
  model
}

#' internal function for Eigenvalue plot
#' 
#' @details
#' Compares shrunk eigenvalues against sample eigenvalues
#' 
#' @param x model of the type spikedCovariance
#' @param norm Type of matrix norm that must be calculated. If missing value from
#'              the model will be used
#' @param statistical Stein/Entropy/Divergence/Affinity/Frechet. Default is set to NA.
#'                    when a valid value is set norm and pivot values are ignored

.plot.spikedCovariance <- function(x, norm = NA, statistical = FALSE) {
  
  data <- x$data; gamma <- x$dist.fit$gamma; numOfSpikes <- x$numOfSpikes 
  if(is.na(norm)) norm <- x$norm

  if (!norm %in% c("Frobenius", "Operator", "Nuclear"))
    stop("Invalid norm value")

  if(statistical) norm <- "Statistical"
  
  loss <- c("Stein","Entropy","Divergence","Affinity","Frechet")
  
  numLosses <- if(statistical) length(loss) else 7
  E <- matrix(NA, nrow = ncol(data), ncol = numLosses)
  
  for(i in 1:numLosses) {
    C <- 
      if (!statistical) {
        estSpikedCovariance(data, gamma = gamma, numOfSpikes = numOfSpikes, 
                               norm = norm, pivot = i, fit = x$dist.fit)
      } else {
        estSpikedCovariance(data, gamma = gamma, numOfSpikes = numOfSpikes, 
                               statistical = loss[i], fit = x$dist.fit)
      }
    
    E[,i] <- C$new.lambdas
  }
  
  colnames(E) <- if (!statistical) paste("pivot",1:7,sep = ".") else loss
  
  df <- data.frame(sample = x$orig.lambdas, E)
  df <- melt(df, id = "sample")
  colnames(df) <- c("samplee" ,"typeofloss" , "shrunke")
  
  colors  <- c("red" , "green" , "blue" , "cyan" , "hotpink" , "gold3" , 
               "black", "purple")
  
  colors <- c(head(colors, numLosses), tail(colors,1))
  names(colors) <- c(colnames(E), "y=x")
  
  p <- ggplot(data = df, aes(x = samplee, y = shrunke )) + 
    geom_point(aes(color = typeofloss)) + geom_line(aes(color = typeofloss)) + 
    geom_abline(aes(color = "y=x"), slope=1, intercept=0) + 
    scale_color_manual("", values = colors) + 
    xlab("Sample Eigenvalues") + ylab("Shrunk Eigenvalues") + 
    scale_x_continuous(breaks=pretty_breaks(n=10)) +
    scale_y_continuous(breaks=pretty_breaks(n=10)) +
    ggtitle(paste("Norm = ",norm, ", NumofSpikes = ", numOfSpikes, 
                  ", gamma = ", round(gamma, digits = 6), sep = "")) + 
    theme_bw() +
    theme(plot.title = element_text(size = 16, face = "bold", vjust = 1),
          axis.title=element_text(size=12), legend.key = element_blank())
  
  p
}

#' Eigenvalue plot
#' 
#' @details
#' Compares shrunk eigenvalues against sample eigenvalues for all norms and losses
#' 
#' @importFrom scales pretty_breaks
#' 
#' @param R xts object of asset returns
#' @param numberOfSpikes model of the type spikedCovariance
#' @param gamma atio of varibales/observations. If NA it will be estimated from
#'                 the data. One can set it to variables/observations
#' @param ... additional arguments unused
#' @author Rohit Arora
#' @examples 
#' \dontrun{
#'  data("rmtdata")
#'  model <- estSpikedCovariance(rmtdata, numOfSpikes=10)
#'  plot(model)
#' }
#' 
#' @export
#' 
plot.spikedCovariance <- function(R, gamma = NA, numOfSpikes = 1,  ...) {
  
  dummyModel <- estSpikedCovariance(R, gamma = gamma, numOfSpikes = numOfSpikes, 
                            norm = "Frobenius", pivot = 1)
  
  p.frob <- .plot.spikedCovariance(dummyModel, norm = "Frobenius")
  p.oper <- .plot.spikedCovariance(dummyModel, norm =  "Operator")
  p.nuc  <- .plot.spikedCovariance(dummyModel, norm =  "Nuclear")
  p.stat <- .plot.spikedCovariance(dummyModel, statistical = TRUE)
  
  grid.arrange(p.frob, p.nuc, p.oper, p.stat, ncol = 2)
}