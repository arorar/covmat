symbols <- c('BABA', 'TWTR', 'LNKD', 'YHOO', 'GE')
symCount <- length(symbols)

data("missingdata")

symdata <- missingdata["2007-04-01/", symbols]
models  <- stambaugh.fit(symdata, method = c("classic", "truncated"))

#expect 2 models
expect_equal(length(models), 2)

#class check
expect_is(models, "stambaugh")

#dimension check for location
expect_equal(nrow(models$models$Stambaugh$center), symCount)
expect_equal(nrow(models$models$Truncated$center), symCount)

#dimension check for covariance
expect_equal(nrow(models$models$Stambaugh$cov), symCount)
expect_equal(ncol(models$models$Stambaugh$cov), symCount)
expect_equal(nrow(models$models$Stambaugh$cov), symCount)
expect_equal(ncol(models$models$Stambaugh$cov), symCount)

#model type check
expect_equal(models$models$Truncated$type, "Classical")
expect_equal(models$models$Truncated$type, "Classical")

# need to select styles. This defaults to 3
expect_error(stambaugh.fit(data))

# invalid data object
expect_error(stambaugh.fit(c(), method = c("classic", "truncated")))

# invalid method
expect_error(stambaugh.fit(data, method = c("class", "truncated")))

# invalid plot option
expect_error(plot.stambaugh(models,5))

# ellipse plot needs two models
expect_error(plot.stambaugh(models$models$Classical,1))

# truncated not allowed in distance plot
expect_error(plot(models,2))


