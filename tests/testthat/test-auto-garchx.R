test_that("auto_garchx runs and returns a usable garchx fit", {
  skip_on_cran()
  set.seed(123)
  
  y <- as.numeric(scale(rnorm(400)))
  res <- auto_garchx(
    y,
    max.p = 2, max.q = 2, search_o = c(0L, 1L),
    final_refit = TRUE, vcov.type = "robust", verbose = FALSE
  )
  
  expect_true(is.list(res))
  expect_true(all(c("fit","p","q","o","bic","aic") %in% names(res)))
  expect_s3_class(res$fit, "garchx")
  
  vhat <- fitted(res$fit)      
  eh   <- residuals(res$fit)   
  
  n <- length(y)
  expect_true(length(vhat) %in% c(n, n-1))
  expect_true(length(eh)   %in% c(n, n-1))
  
  expect_true(all(is.finite(vhat)))
  expect_true(all(vhat >= 0))
  expect_true(all(is.finite(eh)))
  
  pr <- predict(res$fit, n.ahead = 3)
  expect_type(pr, "double")
  expect_length(pr, 3)
  expect_true(all(is.finite(pr)))
})
