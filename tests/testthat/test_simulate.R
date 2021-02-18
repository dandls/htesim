
test_that("simulate returns correct output", {
  sx <- simulate(dgp(m = mF_x1, model = "normal", xmodel = "unif"), nsim = 500, dim = 2)
  expect_s3_class(sx, c("simdgp", "data.frame"))
  expect_identical(nrow(sx), 500L)
  expect_identical(ncol(sx), 4L) # x1, x2, y, trt
  expect_true(all(sx$X1 >= 0 & sx$X1 <= 1))

  sx2 <- simulate(dgp(), nsim = 5000, dim = 5)
  expect_identical(nrow(sx2), 5000L)
  expect_identical(ncol(sx2), 7L)
  expect_true(summary(lm(y ~ X1, data = sx2))$coefficients[2,4] > 0.1)
  expect_true(summary(lm(y ~ trt, data = sx2))$coefficients[2,4] > 0.1)

})


test_that("simulate returns correct output for weibull model", {
  sw <- simulate(dgp(model = "weibull"), nsim = 500, dim = 2)

})


