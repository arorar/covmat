symbols <- c('BABA', 'TWTR', 'LNKD', 'YHOO', 'GE')
symCount <- length(symbols)

data("returnsdata")

symdata <- symdata["2007-04-01/", symbols]
models  <- stambaugh.fit(symdata, style = c("classic", "truncated"))

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


