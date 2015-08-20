symbols <- c('LNKD', 'V', 'LAZ')
symCount <- length(symbols)

data("missingdata")
data("factordata")

symdata <- missingdata["2007-04-01/2014-12-31",symbols]

dates.ret.monthly <- format(index(symdata), "%Y%m")
dates.factors.monthly <- format(index(factor.data), "%Y%m")
index(symdata) <- index(factor.data)[which(dates.factors.monthly %in% 
                                             dates.ret.monthly)]

rets <- fmmc.cov(symdata, factor.data, parallel = FALSE)

#expect 3 returns
expect_equal(ncol(rets), symCount)

#invalid align paramter
expect_error(fmmc.cov(symdata, factor.data, align = "middle"))
