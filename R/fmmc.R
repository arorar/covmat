# Calculate fmmc estimates and their standard errors for risk and performance measures

rm(list = ls())

library(factorAnalytics)
library(quantmod)
library(robust)

symbols <- c('LNKD', 'V', 'LAZ')
getSymbols(symbols, adjust=TRUE, from ="2007-04-01", to = "2014-12-31")

for(symbol in symbols) {
    temp0  <- temp <- adjustOHLC(get(symbol), symbol.name=symbol)
    temp    <- to.monthly(temp0, indexAt='endof', drop.time=FALSE)
    colnames(temp) <- colnames(temp0)    
    assign(x=symbol, value=temp)
}

symbolData <- do.call(merge, lapply(symbols, function(x) Cl(get(x))))
colnames(symbolData) <- symbols
retdata <- diff(log(symbolData))[-1,]

load("./data/factordata.RData")
dates.ret.monthly <- format(index(retdata), "%Y%m")
dates.factors.monthly <- format(index(factor.data), "%Y%m")
index(retdata) <- index(factor.data)[which(dates.factors.monthly %in% dates.ret.monthly)]

# -------------------------------------------
# Part 2. Use FMMC
# -------------------------------------------

fit.method <- "Robust"; robust <- ifelse(fit.method=="Robust",TRUE,FALSE)
# Create an fmmc object that can be used to calc risk and performance estimates
objs <- fmmc(retdata, factor.data, parallel=TRUE, 
    variable.selection="subsets", fit.method=fit.method)

m <- min(do.call(c,lapply(lapply(objs, function(x) x$bootdist$returns),nrow)))
rets <- do.call(cbind, lapply(objs, function(x) x$bootdist$returns[1:m,1]))
colnames(rets) <- unlist(lapply(objs, function(x) dimnames(x$bootdist$returns)))

# We are done. Now we can calculate the covariance in a usual way
cov.mat <- if(robust) covRob(rets)$cov else cov(rets)
rownames(cov.mat) <- colnames(rets)
colnames(cov.mat) <- colnames(rets)

compare.cov <- function(cov1, cov2, labels, corr=FALSE) {
    if (length(labels) != 2) stop("There must be two labels")
    
    complist <- list(list(corr=corr, cov=cov1), list(corr=corr, cov=cov2))
    names(complist) <- labels
    ellipsesPlot.covfm(complist)
}

compare.cov(cov(rets),cov(na.omit(retdata)),c("FMMC","Truncated"))

