data("rmtdata")

# Rows must be greater than columns
expect_error(estSpikedCovariance(rmtdata[1:10,], gamma=2))

# gamma must be between 0 and 1
expect_error(estSpikedCovariance(rmtdata, gamma=2))
expect_error(estSpikedCovariance(rmtdata, gamma=-2))

# numOfSpikes must be positive
expect_error(estSpikedCovariance(rmtdata, numOfSpikes=-2))

# norm must be in "Frobenius", "Operator", "Nuclear"
expect_error(estSpikedCovariance(rmtdata, norm="junk"))

# pivot must be in 1...7
expect_error(estSpikedCovariance(rmtdata, pivot=-1))
expect_error(estSpikedCovariance(rmtdata, pivot=8))

# method must be median-fitting or KNTest
expect_error(estSpikedCovariance(rmtdata, method="junk"))