#' This is an implementation of covariance matrix estimation using Factor Model 
#' Monte Carlo method as described in Jiang and Martin (2013).
#' 
#' @details
#' This method takes in returns data and factors as time series. The returns 
#' and factor time-series must be aligned. The method is applicable 
#' when we have longer factor histories and shorter asset return histories. In 
#' this case past returns can be backfilled using the FMMC methodology. We can 
#' apply FMMC to all the assets in the portfolio. Once the returns are backfilled
#' we can construct a more accurate measure of classical and robust covariance 
#' using backfilled returns. Backfilled returns for different assets must be 
#' aligned by truncating date in the begining or end is returns have unequal
#' available return histories.
#' 
#' @importFrom factorAnalytics fmmc
#' @param  R vector of asset returns in xts format
#' @param  factors matrix of factor returns in xts format
#' @param  robust boolean to indicate if robust methods must be used for fitting
#'          factor models and constucting the covariance matrix. By default this
#'          option is turned off.
#' @param  parallel boolean to indicate if all cores on the system must be used
#' @param  align string to indicate where the longer backfilled return histories
#'         must be truncated when available returns are unequal. The default value
#'         is to align them at the end
#'          truncated to align the returns
#' @param  ... allows passing paramters to factorAnalytics to butild the factor
#'          model or control parameters to covrob for constructing the robust 
#'          covariance matrix
#' @author Rohit Arora
#' @export
#' 
fmmc.cov <-function(R, factors, robust = FALSE, parallel = TRUE, 
                    align=c("end", "begin"), ...) {
    
  #default is end
  align <- align[1]
  if(!align %in% c("end", "begin"))  stop("Invalid align parameter")
  
  fit.method <- ifelse(robust,"Robust","LS")
    
    add.args <- list(...)
    if(!"fit.method" %in% names(add.args)) 
        add.args[["fit.method"]] <- fit.method
    else if(robust && add.args[["fit.method"]] == "LS") 
        stop("fit.method does not agree with robust flag")
        
    if(!"variable.selection" %in% names(add.args)) 
        add.args[["variable.selection"]] <- "subsets"

    # Create an fmmc object that can be used to calc risk and performance estimates
    objs <- fmmc(R, factors, parallel=parallel, add.args)
    
    m <- min(do.call(c,lapply(lapply(objs, function(x) x$bootdist$returns),nrow)))
    
    rets <- if(align == "begin")
        do.call(cbind, lapply(objs, function(x) x$bootdist$returns[1:m]))
    else
        do.call(cbind, lapply(objs, function(x) tail(x$bootdist$returns,m)))
    
    colnames(rets) <- unlist(lapply(objs, function(x) dimnames(x$bootdist$returns)))
    
    cov.estim <- "mcd"; cov.control <- covRob.control("mcd")
    if("estim" %in% names(add.args)) cov.estim <- add.args[["estim"]]
    if("cov.control" %in% names(add.args)) cov.control <- add.args[["cov.control"]]
    
    cov.mat <- if(robust) 
        covRob(rets, estim = cov.estim, control = cov.control)$cov 
    else cov(rets)
    
    rownames(cov.mat) <- colnames(rets)
    colnames(cov.mat) <- colnames(rets)
    
    cov.mat
}