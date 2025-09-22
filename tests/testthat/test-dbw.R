test_that("dbw returns reasonable outputs on toy data", {
  set.seed(1)
  T <- 50; K <- 3
  target <- data.frame(x1=rnorm(T), x2=rnorm(T, .2), x3=rnorm(T,-.1))
  donor1 <- data.frame(x1=rnorm(T, .1), x2=rnorm(T,-.1), x3=rnorm(T, .3))
  donor2 <- data.frame(x1=rnorm(T,-.2), x2=rnorm(T, .4), x3=rnorm(T, 0))
  X_list <- list(target, donor1, donor2)
  shock  <- c(40, 40, 40)
  
  res <- dbw(
    X = X_list, dbw_indices = 1:K, shock_time_vec = shock,
    center = TRUE, scale = TRUE, sum_to_1 = TRUE,
    bounded_below_by = 0, bounded_above_by = 1,
    normchoice = "l2", penalty_normchoice = "l2", penalty_lambda = 1e-2
  )
  
  expect_type(res$opt_params, "double")
  expect_length(res$opt_params, 2)
  expect_true(is.finite(res$loss))
  expect_true(res$convergence %in% c("convergence","failed_convergence"))
})
