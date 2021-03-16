test_that("throw comprehensible errors for wrong inputs", {
  expect_error(dgp(model = "hweoh"), "Must be element of set")
  expect_error(dgp(xmodel = "hweoh"), "Must be element of set")
})


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


test_that("throw comprehensible errors for wrong inputs", {
  dgp1 <- dgp(t = tF_max_x1_x5)
  expect_error(simulate(dgp1, nsim = 0.57), "integerish")
  expect_error(simulate(dgp1, dim = 0.57), "integerish")
  expect_error(simulate(dgp1, nsimtest = 0.77), "integerish")
  expect_error(simulate(dgp1, seed = "ab"), "number")
  expect_error(simulate(dgp1, dim = 3), "increase dim")
  expect_error(simulate(dgp(p = pF_sin_x3), dim = 2), "increase dim")
  expect_error(simulate(dgp(m = mF_max2_x1_x5), dim = 3), "increase dim")
})


test_that("predict works properly", {
  dg1 <- dgp()
  sim1 <- simulate(dg1, nsim = 100)
  tau <- predict(dg1, newdata = sim1)
  expect_equal(nrow(tau), 100L)
  expect_equal(ncol(tau), 4L)
  expect_true(all(tau[, "pfct"] == 0.5))

  tau2 <- predict(sim1, newdata = sim1)
  expect_equal(tau, tau2)
})

