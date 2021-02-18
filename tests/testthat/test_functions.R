test_that("simulate returns correct output for normal model", {
  tF <- c("tF_exp_x1_x2", "tF_div_x1_x2", "tF_log_x1_x2", "tF_max_x1_x5",
    "tF_max_x1_x5", 0.5, 1)
  pF <- c("pF_x1", "pF_x3", "pF_x4", "pF_sin_x3", "pF_eta_x1_x2",
    "pF_x2_x3", "pF_exp_x1_x2", 0.5)
  mF <- c("mF_x1", "mF_x3", "mF_sin_x1_x5", "mF_max_x1_x5", "mF_log_x1_x3",
    "mF_max2_x1_x5", 0)

  testgrid <- expand.grid(tF, pF, mF)

  apply(testgrid, MARGIN = 1, FUN = function(row) {
    expect_error(simulate(dgp(p = eval(parse(text = pF)),
      m = eval(parse(text = mF)),
      t = eval(parse(text = tF)))), NA)
    invisible()
  })
})
