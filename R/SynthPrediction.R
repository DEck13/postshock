#' Synthetic Prediction with Donor Balancing and Shock Adjustment
#'
#' @description
#' Builds k-step-ahead forecasts for a target time series by:
#' (1) estimating post-shock fixed effects on each donor via (S)ARIMA with an
#' indicator regressor, (2) combining donor shock effects with donor weights
#' (uniform or from \code{\link{dbw}}), and (3) adjusting the target ARIMA
#' forecast by the weighted shock effect. Designed for aligned time indices.
#'
#' @param Y_series_list List of time series (e.g., \code{xts}). \code{[[1]]} is
#' the target; \code{[[2]]..[[n+1]]} are donors.
#' @param covariates_series_list List of time series aligned with
#' \code{Y_series_list}. Used for DBW and as optional \code{xreg} for the target.
#' You may pass \code{Y_series_list} here as a simple baseline.
#' @param shock_time_vec Vector of shock start times. Each element can be an
#' index value present in the corresponding series or a positive integer row
#' number. Length must equal \code{length(Y_series_list)}.
#' @param shock_length_vec Integer vector (same length as \code{Y_series_list})
#' giving the length (in rows) of the post-shock window used to estimate
#' donor fixed effects.
#' @param k Integer forecast horizon for the target series.
#' @param dbw_scale,dbw_center Logical flags forwarded to \code{dbw} when donor
#' weights are estimated.
#' @param dbw_indices Integer vector of covariate columns used by \code{dbw}.
#' Defaults to all columns of the first covariate matrix when \code{NULL}.
#' @param princ_comp_input Optional PCA dimension for \code{dbw}; forwarded only.
#' @param covariate_indices Integer vector of columns from
#' \code{covariates_series_list[[i]]} to include as \code{xreg} (lagged for the
#' target; contemporaneous for donors) in ARIMA fits. If \code{NULL}, donors use
#' only the post-shock indicator and the target uses no \code{xreg}.
#' @param geometric_sets,days_before_shocktime_vec Reserved (ignored).
#' @param arima_order Optional fixed ARIMA order \code{c(p,d,q)} for all series.
#' If \code{NULL}, \code{forecast::auto.arima} is used.
#' @param seasonal Logical, whether seasonal ARIMA is allowed/requested.
#' @param seasonal_order Optional seasonal order \code{c(P,D,Q)} when
#' \code{seasonal=TRUE} and \code{arima_order} is supplied.
#' @param seasonal_period Optional seasonal period; auto-detected if \code{NULL}
#' and \code{seasonal=TRUE}.
#' @param user_ic_choice Information criterion for \code{auto.arima}
#' (\code{"aicc"}, \code{"aic"}, or \code{"bic"}).
#' @param stepwise,approximation Passed to \code{forecast::auto.arima}.
#' @param plots Logical; if \code{TRUE} and a function named
#' \code{plot_maker_synthprediction} is found via \code{get0()}, it is called.
#' @param display_ground_truth_choice Logical; forwarded to the plotting helper
#' if it exists.
#' @param penalty_lambda Numeric L1/L2 penalty weight forwarded to \code{dbw}.
#' @param penalty_normchoice Character \code{"l1"} or \code{"l2"}; forwarded to
#' \code{dbw}.
#'
#' @return A list with elements:
#' \itemize{
#'   \item \code{linear_combinations}: numeric donor weights \eqn{W}.
#'   \item \code{predictions}: list with numeric vectors
#'     \code{unadjusted}, \code{adjusted}, \code{arithmetic_mean}.
#'   \item \code{loss}: optional list of SSE vs. ground truth (if provided in env).
#'   \item \code{meta}: list with metadata (weights status, shock info, ARIMA orders, IC used).
#' }
#'
#' @examples
#' \donttest{
#' library(xts)
#' set.seed(1)
#' T <- 80
#' idx <- seq.Date(as.Date("2022-01-01"), by="day", length.out=T)
#' y1 <- xts(rnorm(T), idx)  # target
#' y2 <- xts(rnorm(T), idx)  # donor 1
#' y3 <- xts(rnorm(T), idx)  # donor 2
#' Y <- list(y1, y2, y3)
#' X <- list(y1, y2, y3)     # simple baseline
#' shock_time  <- rep(idx[60], length(Y))
#' shock_len   <- rep(5L, length(Y))
#' out <- SynthPrediction(
#'   Y_series_list=Y, covariates_series_list=X,
#'   shock_time_vec=shock_time, shock_length_vec=shock_len,
#'   k=3, plots=FALSE
#' )
#' str(out$predictions)
#' }
#'
#' @export
SynthPrediction <- function(
    Y_series_list,
    covariates_series_list,
    shock_time_vec,
    shock_length_vec,
    k = 1,
    dbw_scale = TRUE,
    dbw_center = TRUE,
    dbw_indices = NULL,
    princ_comp_input = min(length(shock_time_vec), ncol(covariates_series_list[[1]])),
    covariate_indices = NULL,
    geometric_sets = NULL,
    days_before_shocktime_vec = NULL,
    arima_order = NULL,
    seasonal = FALSE,
    seasonal_order = NULL,
    seasonal_period = NULL,
    user_ic_choice = c("aicc","aic","bic")[1],
    stepwise = TRUE,
    approximation = TRUE,
    plots = TRUE,
    display_ground_truth_choice = FALSE,
    penalty_lambda = 0,
    penalty_normchoice = "l1"
) {
  # 0) Basic checks
  stopifnot(is.list(Y_series_list), length(Y_series_list) >= 2)
  stopifnot(is.list(covariates_series_list),
            length(covariates_series_list) == length(Y_series_list))
  stopifnot(length(shock_time_vec)  == length(Y_series_list))
  stopifnot(length(shock_length_vec) == length(Y_series_list))
  
  n <- length(Y_series_list) - 1
  
  if (is.null(dbw_indices)) {
    dbw_indices <- seq_len(ncol(as.matrix(covariates_series_list[[1]])))
  }
  
  # Coerce to xts so index()/lag.xts work consistently
  to_xts <- function(x) if (xts::is.xts(x)) x else xts::as.xts(x)
  Y_series_list <- lapply(Y_series_list, to_xts)
  covariates_series_list <- lapply(covariates_series_list, to_xts)
  
  # Resolve shock times to integer indices
  int_shock_idx <- integer(length(Y_series_list))
  for (i in seq_along(Y_series_list)) {
    if (is.numeric(shock_time_vec[i]) && length(shock_time_vec[i]) == 1) {
      int_shock_idx[i] <- as.integer(shock_time_vec[i])
    } else {
      idx <- zoo::index(Y_series_list[[i]])
      pos <- which(idx == shock_time_vec[i])
      if (length(pos) != 1)
        stop("Shock time not found (or not unique) for series ", i, ".")
      int_shock_idx[i] <- pos
    }
    if (int_shock_idx[i] < 2) stop("Shock time must be >= 2 (series=", i, ").")
  }
  
  # 1) Donor shock effects via (S)ARIMA with an indicator
  omega_star_hat_vec <- numeric(n)
  order_of_arima     <- vector("list", n + 1L)
  
  for (i in 2:(n + 1L)) {
    shock_start <- int_shock_idx[i]
    last_pt     <- shock_start + as.integer(shock_length_vec[i]) - 1L
    last_pt     <- min(last_pt, NROW(Y_series_list[[i]]))
    zeros <- rep(0, shock_start - 1L)
    ones  <- rep(1, last_pt - shock_start + 1L)
    post_ind <- c(zeros, ones)
    
    if (is.null(covariate_indices)) {
      X_i_final <- matrix(post_ind, ncol = 1L)
    } else {
      X_cov <- as.matrix(covariates_series_list[[i]][1:last_pt, covariate_indices, drop = FALSE])
      X_i_final <- cbind(X_cov, post_ind)
      X_i_final <- stats::na.omit(X_i_final)
    }
    
    y_i <- Y_series_list[[i]][1:last_pt]
    
    seasonal_period_i <- seasonal_period
    if (isTRUE(seasonal) && is.null(seasonal_period_i)) {
      seasonal_period_i <- tryCatch(
        forecast::findfrequency(as.numeric(y_i)),
        error = function(e) 4
      )
    }
    
    if (!is.null(arima_order)) {
      if (isTRUE(seasonal) && !is.null(seasonal_order)) {
        donor_model <- forecast::Arima(
          y = y_i,
          order = arima_order,
          seasonal = list(order = seasonal_order, period = seasonal_period_i),
          xreg = X_i_final
        )
      } else {
        donor_model <- forecast::Arima(
          y = y_i,
          order = arima_order,
          seasonal = FALSE,
          xreg = X_i_final
        )
      }
    } else {
      donor_model <- tryCatch(
        forecast::auto.arima(
          y = y_i,
          xreg = X_i_final,
          ic = user_ic_choice,
          seasonal = seasonal,
          stepwise = stepwise,
          approximation = approximation
        ),
        error = function(e) {
          warning(sprintf("Donor %d: auto.arima() failed -> fallback to Arima(1,1,1).", i))
          forecast::Arima(y = y_i, order = c(1,1,1), seasonal = FALSE, xreg = X_i_final)
        }
      )
    }
    
    order_of_arima[[i]] <- donor_model$arma
    cf  <- lmtest::coeftest(donor_model)
    omega_star_hat_vec[i - 1L] <- cf[nrow(cf), "Estimate"]
  }
  
  # 2) Donor weights
  if (is.null(covariate_indices)) {
    w_hat <- rep(1 / n, n)
    omega_star_hat <- sum(w_hat * omega_star_hat_vec)
    dbw_status <- NA
  } else {
    dbw_fit <- dbw(
      X                  = covariates_series_list,
      dbw_indices        = dbw_indices,
      shock_time_vec     = int_shock_idx,
      scale              = dbw_scale,
      center             = dbw_center,
      sum_to_1           = TRUE,
      bounded_below_by   = 0,
      bounded_above_by   = 1,
      Y                  = Y_series_list,
      Y_lookback_indices = NULL,
      X_lookback_indices = rep(list(1), length(dbw_indices)),
      inputted_transformation = base::identity,
      penalty_lambda     = penalty_lambda,
      penalty_normchoice = penalty_normchoice
    )
    w_hat <- as.numeric(dbw_fit$opt_params)
    omega_star_hat <- sum(w_hat * omega_star_hat_vec)
    dbw_status <- dbw_fit$convergence
  }
  
  # 3) Target ARIMA and k-ahead forecast
  target_pre <- Y_series_list[[1]][1:int_shock_idx[1]]
  
  if (is.null(covariate_indices)) {
    xreg_input  <- NULL
    xreg_future <- NULL
  } else {
    X_pre  <- as.matrix(covariates_series_list[[1]][1:int_shock_idx[1], covariate_indices, drop = FALSE])
    X_lag1 <- xts::lag.xts(X_pre, k = 1)
    xreg_input <- stats::na.omit(X_lag1)
    
    if (NROW(xreg_input) < NROW(target_pre)) {
      drop_n <- NROW(target_pre) - NROW(xreg_input)
      target_pre <- target_pre[-seq_len(drop_n)]
    }
    
    X_last <- X_pre[NROW(X_pre), , drop = FALSE]
    xreg_future <- matrix(rep(as.numeric(X_last), each = k), nrow = k, byrow = TRUE)
  }
  
  target_seasonal_period <- seasonal_period
  if (isTRUE(seasonal) && is.null(target_seasonal_period)) {
    target_seasonal_period <- tryCatch(
      forecast::findfrequency(as.numeric(target_pre)),
      error = function(e) 4
    )
  }
  
  if (!is.null(arima_order)) {
    if (isTRUE(seasonal) && !is.null(seasonal_order)) {
      target_model <- forecast::Arima(
        y = target_pre,
        order = arima_order,
        seasonal = list(order = seasonal_order, period = target_seasonal_period),
        xreg = xreg_input
      )
    } else {
      target_model <- forecast::Arima(
        y = target_pre,
        order = arima_order,
        seasonal = FALSE,
        xreg = xreg_input
      )
    }
  } else {
    target_model <- forecast::auto.arima(
      y = target_pre,
      xreg = xreg_input,
      ic = user_ic_choice,
      seasonal = seasonal,
      stepwise = stepwise,
      approximation = approximation
    )
  }
  
  if (is.null(xreg_future)) {
    fc <- stats::predict(target_model, n.ahead = k)
  } else {
    fc <- stats::predict(target_model, n.ahead = k, newxreg = xreg_future)
  }
  raw_pred        <- as.numeric(fc$pred)
  adj_pred        <- raw_pred + omega_star_hat
  arith_mean_pred <- raw_pred + mean(omega_star_hat_vec)
  
  predictions <- list(
    unadjusted      = raw_pred,
    adjusted        = adj_pred,
    arithmetic_mean = arith_mean_pred
  )
  
  loss <- NULL
  if (exists("ground_truth_vec", inherits = TRUE)) {
    gt <- get("ground_truth_vec", inherits = TRUE)
    if (is.numeric(gt)) {
      loss <- list(
        unadjusted      = sum((raw_pred        - gt)^2),
        adjusted        = sum((adj_pred        - gt)^2),
        arithmetic_mean = sum((arith_mean_pred - gt)^2)
      )
    }
  }
  
  meta <- list(
    n_donors       = n,
    shock_time     = shock_time_vec,
    shock_length   = shock_length_vec,
    dbw_status     = dbw_status,
    weights        = w_hat,
    omega_vec      = omega_star_hat_vec,
    combined_omega = as.numeric(omega_star_hat),
    arima_order    = order_of_arima,
    ic_used        = user_ic_choice
  )
  
  out <- list(
    linear_combinations = w_hat,
    predictions         = predictions,
    meta                = meta
  )
  if (!is.null(loss)) out$loss <- loss
  
  # 5) Optional plot hook (lookup via get0 to avoid NOTE)
  if (isTRUE(plots)) {
    f <- get0("plot_maker_synthprediction", mode = "function", inherits = TRUE)
    if (!is.null(f)) {
      try(
        f(
          Y_series_list,
          shock_time_vec,
          int_shock_idx,
          shock_length_vec,
          raw_pred,
          w_hat,
          omega_star_hat_vec,
          adj_pred,
          arith_mean_pred,
          display_ground_truth = display_ground_truth_choice
        ),
        silent = TRUE
      )
    }
  }
  
  out
}
