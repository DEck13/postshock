test_that("SynthPrediction runs on tiny data", {
  skip_on_cran()
  set.seed(1)
  
  T  <- 40
  idx <- as.Date("2020-01-01") + 0:(T-1)    # build a valid time index
  
  y1 <- rnorm(T); y2 <- rnorm(T); y3 <- rnorm(T)
  Y  <- list(
    xts::xts(y1, order.by = idx),
    xts::xts(y2, order.by = idx),
    xts::xts(y3, order.by = idx)
  )
  
  X  <- list(
    xts::xts(cbind(x1 = rnorm(T), x2 = rnorm(T)), order.by = idx),
    xts::xts(cbind(x1 = rnorm(T), x2 = rnorm(T)), order.by = idx),
    xts::xts(cbind(x1 = rnorm(T), x2 = rnorm(T)), order.by = idx)
  )
  
  shock_t <- rep(idx[30], length(Y))    # use the same index for shock time
  shock_L <- rep(5L, length(Y))
  
  expect_no_error({
    out <- SynthPrediction(
      Y_series_list          = Y,
      covariates_series_list = X,
      shock_time_vec         = shock_t,
      shock_length_vec       = shock_L,
      k = 2,
      covariate_indices = 1:2,
      plots = FALSE
    )
    expect_type(out, "list")
    expect_true(all(c("linear_combinations","predictions","meta") %in% names(out)))
    expect_length(out$predictions$adjusted, 2)
    expect_true(all(is.finite(out$predictions$adjusted)))
  })
})
