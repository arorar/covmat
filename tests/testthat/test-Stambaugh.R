library(xts)
library(robust)
library(quantmod)

symbols <- c('BABA', 'TWTR', 'LNKD', 'YHOO', 'GE')
symCount <- length(symbols)

getSymbols(symbols, adjust=TRUE, from ="2010-04-01", to = "2015-05-31")

for(symbol in symbols) {
  temp0  <- temp <- adjustOHLC(get(symbol), symbol.name=symbol)
  temp    <- to.monthly(temp0, indexAt='endof', drop.time=FALSE)
  colnames(temp) <- colnames(temp0)    
  assign(x=symbol, value=temp)
}

data <- do.call(merge, lapply(symbols, function(x) Cl(get(x))))
colnames(data) <- symbols

models <- stambaugh.fit(data, style = c("classic", "truncated"))

#expect 2 models
expect_equal(length(models), 2)

#class check
expect_is(models, "stambaugh")

#dimension check for location
expect_equal(nrow(models$models$Classical$center), symCount)
expect_equal(nrow(models$models$Truncated$center), symCount)

#dimension check for covariance
expect_equal(nrow(models$models$Classical$cov), symCount)
expect_equal(ncol(models$models$Classical$cov), symCount)
expect_equal(nrow(models$models$Classical$cov), symCount)
expect_equal(ncol(models$models$Classical$cov), symCount)

#model type check
expect_equal(models$models$Truncated$type, "Classical")
expect_equal(models$models$Truncated$type, "Classical")

# need to select styles. This defaults to 3
expect_error(stambaugh.fit(data))

# invalid data object
expect_error(stambaugh.fit(c(), style = c("classic", "truncated")))

# invalid style
expect_error(stambaugh.fit(data, style = c("class", "truncated")))

# invalid plot option
expect_error(plot.stambaugh(models,5))

# ellipse plot needs two models
expect_error(plot.stambaugh(models$models$Classical,1))


