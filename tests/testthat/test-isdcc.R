data("etfdata")
numRegimes = 3; T <- nrow(etfdata)
model.isdcc <- isdccfit(etfdata, numRegimes=numRegimes, 
                        parallelType = 0, itermax = 1)

#class check
expect_is(model.isdcc, "isdcc")

#check size
expect_equal(nrow(model.isdcc$filtProb), T)
expect_equal(ncol(model.isdcc$filtProb), numRegimes)
expect_equal(length(model.isdcc$cov), T)

#check contents
expect_true(all(is.numeric(model.isdcc$filtProb)))
expect_true(all(unlist(lapply(model.isdcc$cov, 
                                      function(cov) {
                                           lapply(cov, 
                                                  function(c) is.numeric(c))
                                      }))))

#missing data
expect_error(isdccfit(numRegimes=NA, parallelType = 0, itermax = 1))

#regimes must be an integer
expect_error(isdccfit(etfdata, numRegimes=3.5, parallelType = 0, itermax = 1))

#regimes must be an integer
expect_error(isdccfit(etfdata, numRegimes=NA, parallelType = 0, itermax = 1))

#transMatbounds multiplicity
expect_error(isdccfit(etfdata, numRegimes=NA, transMatbounds = c(2, 10, 20),
                      parallelType = 0, itermax = 1))

#transMatbounds values
expect_error(isdccfit(etfdata, numRegimes=NA, transMatbounds = c(2, 10, NA),
                      parallelType = 0, itermax = 1))

#transMatbounds order
expect_error(isdccfit(etfdata, numRegimes=NA, transMatbounds = c(10, 2),
                      parallelType = 0, itermax = 1))

#dccBounds multiplicity
expect_error(isdccfit(etfdata, numRegimes=NA, dccBounds = c(2, 10, 20),
                      parallelType = 0, itermax = 1))

#dccBounds values
expect_error(isdccfit(etfdata, numRegimes=NA, dccBounds = c(2, 10, NA),
                      parallelType = 0, itermax = 1))

#dccBounds order
expect_error(isdccfit(etfdata, numRegimes=NA, dccBounds = c(1, 0),
                      parallelType = 0, itermax = 1))

#dccBounds sign
expect_error(isdccfit(etfdata, numRegimes=NA, dccBounds = c(-1, 0),
                      parallelType = 0, itermax = 1))

