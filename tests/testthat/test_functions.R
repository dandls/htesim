test_that("simulate returns correct output for normal model", {
  tF <- list(tF_exp_x1_x2, tF_div_x1_x2, tF_log_x1_x2, tF_max_x1_x5,
    tF_max_x1_x5, 0.5, 1)
  pF <- c(pF_x1, pF_x3, pF_x4, pF_sin_x3, pF_eta_x1_x2,
    pF_x2_x3, pF_exp_x1_x2, 0.5)
  mF <- c(mF_x1, mF_x3, mF_sin_x1_x5, mF_max_x1_x5, mF_log_x1_x3,
    mF_max2_x1_x5, 0)

  testgrid <- expand.grid(pF = seq_along(pF), mF = seq_along(mF), tF = seq_along(tF))

  apply(testgrid, MARGIN = 1, FUN = function(row) {
    expect_error(simulate(dgp(
      p = pF[[row[1]]],
      m = mF[[row[2]]],
      t = tF[[row[3]]]), dim = 6), NA)
    invisible()
  })
})
