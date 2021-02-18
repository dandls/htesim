
test_that("simulate returns correct output for normal model", {
  sx <- simulate(dgp(m = mF_x1, model = "normal", xmodel = "unif"), nsim = 500, dim = 2)
  expect_s3_class(sx, c("simdgp", "data.frame"))
  expect_identical(nrow(sx), 500L)
  expect_identical(ncol(sx), 4L) # x1, x2, y, trt
  expect_true(all(sx$X1 >= 0 & sx$X1 <= 1))

  sx2 <- simulate(dgp(), nsim = 500, dim = 5)
  expect_identical(nrow(sx2), 500L)
  expect_identical(ncol(sx2), 7L)

  sx3 <- simulate(dgp(xmodel = "normal"), nsim = 5000, dim = 5)
  expect_true(any(sx3$X1 >= 1))
  expect_true(any(sx3$X2 >= 1))

})


test_that("simulate returns correct output for diverse models", {
  library(tram)
  library(survival)
  library(mlt)

  # weibull
  sw <- simulate(dgp(model = "weibull"), nsim = 500, dim = 2)
  expect_s3_class(sw$y, "Surv")

  # polr
  sp <- simulate(dgp(model = "polr"), nsim = 500, dim = 2)
  expect_s3_class(sp$y, c("ordered", "factor"))

  # binomial
  sb <- simulate(dgp(model = "binomial"), nsim = 500, dim = 2)
  expect_s3_class(sb$y, "factor")
})


test_that("setting seed works", {
  set.seed(1234L)
  curseed <- get(".Random.seed", envir = .GlobalEnv)
  simulate(dgp())
  nowseed <- get(".Random.seed", envir = .GlobalEnv)
  expect_false(isTRUE(all.equal(curseed, nowseed)))

  set.seed(1234L)
  curseed <- get(".Random.seed", envir = .GlobalEnv)
  a <- simulate(dgp(), seed = 555L)
  b <- simulate(dgp(), seed = 555L)
  nowseed <- get(".Random.seed", envir = .GlobalEnv)
  expect_equal(curseed, nowseed)
  expect_true(all.equal(a, b))
})



