library(PortfolioAnalytics)
library(tawny)
library(covmat)
library(ROI)
library(ROI.plugin.quadprog)
library(ROI.plugin.glpk)

custom.portfolio.moments <- function(R, portfolio) {
  momentargs <- list()
  momentargs$mu  <-  matrix(as.vector(apply(R,2, "mean")), ncol = 1)
  momentargs$sigma  <-  estRMT(R)$cov
  momentargs$m3 <- matrix(0, nrow=ncol(R), ncol=ncol(R)^2)
  momentargs$m4 <- matrix(0, nrow=ncol(R), ncol=ncol(R)^3)

  return(momentargs)
}

data("largereturn")
datap <- largesymdata

pspec.lo <- portfolio.spec(assets = colnames(datap))

#long-only
pspec.lo <- add.constraint(pspec.lo, type="full_investment")
pspec.lo <- add.constraint(pspec.lo, type="long_only")

pspec.lo <- add.objective(portfolio=pspec.lo, type="return", name="mean")
pspec.lo <- add.objective(portfolio=pspec.lo, type="risk", name="var")

opt.ordinary <- optimize.portfolio.rebalancing(datap, pspec.lo, 
                              optimize_method="quadprog",
                               rebalance_on="months",
                               training_period=120,
                               trailing_periods=120)
ordinary.wts <- na.omit(extractWeights(opt.ordinary))
ordinary <- Return.rebalancing(R=datap, weights=ordinary.wts)

opt.rmt <- optimize.portfolio.rebalancing(datap, pspec.lo, 
                                               optimize_method="quadprog",
                                               momentFUN = "custom.portfolio.moments",
                                               rebalance_on="months",
                                               training_period=120,
                                               trailing_periods=120)
rmt.wts <- na.omit(extractWeights(opt.rmt))
rmt <- Return.rebalancing(R=datap, weights=rmt.wts)
gmv <- merge.zoo(ordinary,rmt)
colnames(gmv) <- c("ordinary", "rmt")
charts.PerformanceSummary(gmv,wealth.index = T,
                          colorset = c("red","blue"), 
                          cex.legend = 1.3,cex.axis = 1.3,
                          main="Comparison of Portflio Performance using two different covariance matrices",
                          legend.loc = "topleft") 
