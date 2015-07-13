data("largereturn")
symCount <- ncol(largesymdata)

model <- estRMT(largesymdata, cutoff = "each", parallel = FALSE)

#class check
expect_is(model, "RMT")

#basic count check
expect_equal(nrow(model$cov), symCount)
expect_equal(ncol(model$cov), symCount)

# Q check
expect_error(estRMT(largesymdata, Q=0.5))

# numEig check
expect_error(estRMT(largesymdata, numEig=-1))

#Invalid cutoff option
expect_error(estRMT(largesymdata, cutoff="garbage"))

#Invalid eigneTreat option
expect_error(estRMT(largesymdata, eigenTreat="garbage"))

#output checking
expect_more_than(model$Q, 0.999)
expect_more_than(model$var, 0)
expect_more_than(length(model$eigVals), 0)
