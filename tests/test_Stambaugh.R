library(covmat)
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
expect_that(length(models), equals(2))

#class check
expect_that(models, is_a("stambaugh"))

#dimension check for location
expect_that(nrow(models$models$Classical$center), equals(symCount))
expect_that(nrow(models$models$Truncated$center), equals(symCount))

#dimension check for covariance
expect_that(nrow(models$models$Classical$cov), equals(symCount))
expect_that(ncol(models$models$Classical$cov), equals(symCount))
expect_that(nrow(models$models$Classical$cov), equals(symCount))
expect_that(ncol(models$models$Classical$cov), equals(symCount))

#model type check
expect_that(models$models$Truncated$type, equals("Classical"))
expect_that(models$models$Truncated$type, equals("Classical"))

# need to select styles. This defaults to 3
expect_that(stambaugh.fit(data),throws_error())

# invalid data object
expect_that(stambaugh.fit(c(), style = c("classic", "truncated")),throws_error())

# invalid style
expect_that(stambaugh.fit(data, style = c("class", "truncated")),throws_error())

# invalid plot option
expect_that(plot.stambaugh(models,5), throws_error())

# ellipse plot needs two models
expect_that(plot.stambaugh(models$models$Classical,1), throws_error())


