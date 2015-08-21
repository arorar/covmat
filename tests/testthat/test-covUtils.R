data("etfdata")

symdata <- coredata(etfdata[,1:5])

# compare covarianve matrices
cov1 <- cov(symdata)
cov2 <- cov(symdata)

# if plot is produced then ellipse plot will return the orgiinal list
y <- compareCov(cov1 = cov1, cov2 = cov2, labels = c("classical", "robust"))
expect_equal(length(y), 2)
