# tests/testthat/test-synthprediction-smoke.R
test_that("SynthPrediction basic smoke test", {
  skip_if_not_installed("forecast")
  skip_if_not_installed("xts")
  set.seed(1)
  
  T <- 60
  idx <- seq.Date(as.Date("2023-01-01"), by = "day", length.out = T)
  y1 <- xts::xts(rnorm(T), idx)  # target
  y2 <- xts::xts(rnorm(T), idx)  # donor 1
  y3 <- xts::xts(rnorm(T), idx)  # donor 2
  Y  <- list(y1, y2, y3)
  X  <- Y
  
  shock_time <- rep(idx[40], length(Y))
  shock_len  <- rep(3L, length(Y))
  
  out <- SynthPrediction(
    Y_series_list = Y,
    covariates_series_list = X,
    shock_time_vec = shock_time,
    shock_length_vec = shock_len,
    k = 2,
    plots = FALSE
  )
  
  # 结构/性质断言（快且稳定）
  expect_type(out, "list")
  expect_true(all(c("linear_combinations", "predictions", "meta") %in% names(out)))
  expect_length(out$linear_combinations, 2)           # 2 donors
  expect_equal(length(out$predictions$unadjusted), 2) # k = 2
  expect_equal(length(out$predictions$adjusted),   2)
  expect_true(all(is.finite(out$predictions$adjusted)))
})
