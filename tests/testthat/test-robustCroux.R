data("dow30data")
symbols <- c('AAPL', 'MSFT')

R <- dow30data[,which(colnames(dow30data) %in% symbols)]
symCount <- length(colnames(R))

fit <- robustMultExpSmoothing(R, smoothMat = cov(R)) 

#dimension check for location
expect_equal(nrow(fit$covMat), symCount)
expect_equal(ncol(fit$covMat), symCount)

#smoothing matrix must be symmetric
smat <- matrix(1:4, nrow = 2)
expect_error(robustMultExpSmoothing(R, smoothMat = smat))

#eigenvalues of the smoothing matrix must be in [0,1]
smat <- matrix(c(1,3,3,1), nrow = 2)
expect_error(robustMultExpSmoothing(R, smoothMat = smat))

