test_that("throw comprehensible errors for wrong inputs", {
  expect_error(dgp(model = "hweoh"), "Must be element of set")
  expect_error(dgp(xmodel = "hweoh"), "Must be element of set")
})


test_that("simulate returns correct output for normal model", {
  sx <- simulate(dgp(m = mF_x1, model = "normal", xmodel = "unif"), nsim = 500L,
    nsimtest = 100L, dim = 2)
  expect_s3_class(sx, c("simdgp", "data.frame"))
  expect_identical(nrow(sx), 500L)
  expect_identical(ncol(sx), 4L) # x1, x2, y, trt
  expect_true(all(sx$X1 >= 0 & sx$X1 <= 1))
  test <- attributes(sx)$testxdf
  expect_identical(nrow(test), 100L)
  expect_identical(ncol(test), 3L) # x1, x2, trt
  expect_true(all(test$X1 >= 0 & test$X1 <= 1))

  sx2 <- simulate(dgp(), nsim = 500, dim = 5)
  expect_identical(nrow(sx2), 500L)
  expect_identical(ncol(sx2), 7L)

  sx3 <- simulate(dgp(xmodel = "normal"), nsim = 5000, dim = 5)
  expect_true(any(sx3$X1 >= 1))
  expect_true(any(sx3$X2 >= 1))

})


test_that("simulate returns correct output for diverse models", {
  library("tram")

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

test_that("character function call properly handled", {
  sx <- simulate(dgp(p = "0.5", model = "normal", xmodel = "unif"), nsim = 500L,
    nsimtest = 100L, dim = 2)
  tau <- predict(sx, newdata = sx)
  expect_true(all(tau[, "pfct"] == 0.5))
  expect_error(simulate(dgp(p = "ax")))
})

test_that("removing variables after sampling works", {
  test_if_included <- function(sim, x) {
    expect_true(all(!x %in% names(sim)))
    testxdf <- attributes(sim)$testxdf
    expect_true(all(x %in% names(testxdf)))
  }
  # one variable
  dg1 <- dgp(rmvar = "X2")
  sim1 <- simulate(dg1, nsim = 100, dim = 10, nsimtest = 100)
  test_if_included(sim1, "X2")
  # multiple variables
  dg2 <- dgp(rmvar =  c("X2", "X4", "X6"))
  sim2 <- simulate(dg2, nsim = 100, dim = 10, nsimtest = 100)
  test_if_included(sim2, c("X2", "X4", "X6"))

})


test_that("error occurs of name of variable to remove > dim", {
  dg <- dgp(rmvar = "X9")
  expect_error(simulate(dg, nsim = 100, dim = 5, nsimtest = 5))
})


test_that("correlated features are produces", {
  dgp1 <- dgp(p = pF_x1, m = mF_x1, t = 0, model = "normal", xmodel = "correlated")
  sim1 <- simulate(dgp1, nsim = 1000L, d = 2L, nsimtest = 1000L) # in paper d = {2, 5, 10, 15, 20, 30}
  predict(sim1)
  expect_true(all(cor(sim1[, 1:2]) > 0.85))
})
