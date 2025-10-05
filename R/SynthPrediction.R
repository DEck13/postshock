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
# keep same behavior as your scripts that used `id` symbol:
id <- base::identity

### SynthPrediction Function
### Purpose: Synthetic Prediction for Time Series with Shock Analysis
SynthPrediction <- function(Y_series_list,
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
                            user_ic_choice = c('aicc','aic','bic')[1],
                            stepwise = TRUE,
                            approximation = TRUE,
                            plots = TRUE,
                            display_ground_truth_choice = FALSE,
                            penalty_lambda = 0,
                            penalty_normchoice = "l1"
) {
  ### Function Overview:
  n <- length(Y_series_list) - 1
  if (is.null(dbw_indices)) dbw_indices <- 1:ncol(covariates_series_list[[1]])
  
  integer_shock_time_vec <- c()
  integer_shock_time_vec_for_convex_hull_based_optimization <- c()
  for (i in 1:(n+1)){
    if (is.character(shock_time_vec[i])) {
      integer_shock_time_vec[i] <- which(index(Y_series_list[[i]]) == shock_time_vec[i])
      integer_shock_time_vec_for_convex_hull_based_optimization[i] <- which(index(covariates_series_list[[i]]) == shock_time_vec[i])
    } else {
      integer_shock_time_vec[i] <- shock_time_vec[i]
      integer_shock_time_vec_for_convex_hull_based_optimization[i] <- shock_time_vec[i]
    }
  }
  
  omega_star_hat_vec <- c()
  order_of_arima <- list()
  for (i in 2:(n+1)) {
    vec_of_zeros <- rep(0, integer_shock_time_vec[i])
    vec_of_ones <- rep(1, shock_length_vec[i])
    post_shock_indicator <- c(vec_of_zeros, vec_of_ones)
    last_shock_point <- integer_shock_time_vec[i] + shock_length_vec[i]
    
    if (is.null(covariate_indices)) {
      X_i_final <- matrix(post_shock_indicator, ncol = 1L)
      colnames(X_i_final) <- "post_shock_indicator"
    } else {
      X_cov    <- as.matrix(covariates_series_list[[i]][1:last_shock_point, covariate_indices, drop = FALSE])
      X_i_final <- cbind(X_cov, post_shock_indicator = post_shock_indicator)
      X_i_final <- stats::na.omit(X_i_final)
    }
    
    
    # —— Fit (S)ARIMA model to donor series i —— 
    
    seasonal_period_i <- seasonal_period
    if (isTRUE(seasonal) && is.null(seasonal_period_i)) {
      # Auto‐detect seasonal frequency; if that fails, default to 4
      seasonal_period_i <- tryCatch({
        forecast::findfrequency(
          as.numeric(Y_series_list[[i]][1:last_shock_point])
        )
      }, error = function(e) {
        warning(sprintf(
          "Donor %d: seasonal auto-detect failed → default = 4", i
        ))
        4
      })
    }
    
    # ——— Remove any NA rows in X_i_final created by lagging ———
    if (!is.null(covariate_indices)) {
      X_i_final <- na.omit(X_i_final)
    }
    
    # ——— Fit a (Seasonal) ARIMA model for donor i ———
    # Case A: User provided a fixed arima_order
    if (!is.null(arima_order)) {
      if (isTRUE(seasonal) && !is.null(seasonal_order)) {
        # If seasonal = TRUE and a seasonal_order is given, include the seasonal submodel
        donor_model <- forecast::Arima(
          y        = Y_series_list[[i]][1:last_shock_point],
          order    = arima_order,
          seasonal = list(order  = seasonal_order, 
                          period = seasonal_period_i),
          xreg     = X_i_final
        )
      } else {
        # Otherwise, do not include any seasonal component
        donor_model <- forecast::Arima(
          y        = Y_series_list[[i]][1:last_shock_point],
          order    = arima_order,
          seasonal = FALSE,
          xreg     = X_i_final
        )
      }
      
      # Case B: Let auto.arima() select the best (seasonal) ARIMA specification
    } else {
      donor_model <- tryCatch(
        forecast::auto.arima(
          y        = Y_series_list[[i]][1:last_shock_point],
          xreg     = X_i_final,
          ic       = user_ic_choice,
          seasonal = seasonal,
          stepwise = stepwise,
          approximation = approximation
        ),
        error = function(e) {
          # If auto.arima() fails, fall back to a simple ARIMA(1,1,1) with no seasonal
          warning(sprintf(
            "Donor %d: auto.arima() failed → fallback to Arima(1,1,1)", i
          ))
          forecast::Arima(
            y        = Y_series_list[[i]][1:last_shock_point],
            order    = c(1,1,1),
            seasonal = FALSE,
            xreg     = X_i_final
          )
        }
      )
    }
    
    # Record the chosen ARIMA parameters and extract the shock coefficient
    order_of_arima[[i]] <- donor_model$arma
    coef_test <- lmtest::coeftest(donor_model)
    hit <- which(rownames(coef_test) == "post_shock_indicator")
    if (length(hit) != 1L) stop(sprintf("Series %d: cannot find unique 'post_shock_indicator'.", i))
    omega_star_hat_vec[i - 1L] <- coef_test[hit, "Estimate"]
    
  }
  
  if (is.null(covariate_indices)) {
    n_donors <- length(Y_series_list) - 1
    w_hat <- rep(1 / n_donors, n_donors)
    omega_star_hat <- sum(w_hat * omega_star_hat_vec)
    
  } else {
    dbw_output <- dbw(
      X = covariates_series_list,
      dbw_indices = dbw_indices,
      shock_time_vec = integer_shock_time_vec,
      scale = TRUE,
      center = TRUE,
      sum_to_1 = TRUE,
      bounded_below_by = 0,
      bounded_above_by = 1,
      Y = Y_series_list,
      Y_lookback_indices = NULL,
      X_lookback_indices = rep(list(1), length(dbw_indices)),
      inputted_transformation = id,
      penalty_lambda = penalty_lambda,
      penalty_normchoice = penalty_normchoice
    )
    
    w_hat <- dbw_output$opt_params
    omega_star_hat <- as.numeric(w_hat %*% omega_star_hat_vec)
  }
  
  ### Forecast Target Series
  target_series <- Y_series_list[[1]][1:(integer_shock_time_vec[1])]
  forecast_horizon <- k
  
  # Construct covariate matrices for training and forecasting, if applicable
  if (is.null(covariate_indices)) {
    xreg_input  <- NULL
    xreg_future <- NULL
  } else {
    # Create lagged training covariate values (lag 1 before shock)
    X_lagged <- lag.xts(
      covariates_series_list[[1]][1:integer_shock_time_vec[1], covariate_indices]
    )
    xreg_input <- X_lagged
    
    # For forecasting, repeat the last observed covariate vector 'k' times
    X_last <- covariates_series_list[[1]][integer_shock_time_vec[1], covariate_indices]
    xreg_future <- matrix(
      rep(X_last, forecast_horizon),
      nrow = forecast_horizon,
      byrow = TRUE
    )
  }
  
  # Auto‐detect seasonal period for the target if needed
  target_seasonal_period <- seasonal_period
  if (isTRUE(seasonal) && is.null(target_seasonal_period)) {
    target_seasonal_period <- tryCatch({
      forecast::findfrequency(as.numeric(target_series))
    }, error = function(e) {
      warning("Target: seasonal auto‐detect failed → default = 4")
      4
    })
  }
  
  # Fit the ARIMA model for the target series
  if (!is.null(arima_order)) {
    if (isTRUE(seasonal) && !is.null(seasonal_order)) {
      # Include specified seasonal_order
      target_model <- forecast::Arima(
        y        = target_series,
        order    = arima_order,
        seasonal = list(order  = seasonal_order,
                        period = target_seasonal_period),
        xreg     = xreg_input
      )
    } else {
      # No seasonal component
      target_model <- forecast::Arima(
        y        = target_series,
        order    = arima_order,
        seasonal = FALSE,
        xreg     = xreg_input
      )
    }
  } else {
    # Let auto.arima() choose the best model (possibly seasonal)
    target_model <- forecast::auto.arima(
      y        = target_series,
      xreg     = xreg_input,
      ic       = user_ic_choice,
      seasonal = seasonal,
      stepwise = stepwise,
      approximation = approximation
    )
  }
  
  # Generate k‐step ahead forecasts (with or without xreg)
  if (is.null(xreg_future)) {
    unadjusted_pred <- predict(target_model, n.ahead = forecast_horizon)
  } else {
    unadjusted_pred <- predict(
      target_model,
      n.ahead = forecast_horizon,
      newxreg = xreg_future
    )
  }
  
  # Finally, adjust the raw ARIMA forecast by adding the weighted shock effect
  adjusted_pred <- unadjusted_pred$pred + omega_star_hat
  
  # 1. take vectors from objects
  raw_vec       <- as.numeric(unadjusted_pred$pred)          # ARIMA-only
  adjusted_vec  <- as.numeric(adjusted_pred)                 # synthetic + shock
  arithmetic_vec<- raw_vec + mean(omega_star_hat_vec)        # arithmetic mean
  
  # 2. build predictionsn sublist
  predictions <- list(
    unadjusted      = raw_vec,
    adjusted        = adjusted_vec,
    arithmetic_mean = arithmetic_vec
  )
  
  # 3. optional loss sublist
  loss <- NULL
  if (exists("ground_truth_vec") && !is.null(ground_truth_vec)) {
    loss <- list(
      unadjusted      = sum((raw_vec        - ground_truth_vec)^2),
      adjusted        = sum((adjusted_vec   - ground_truth_vec)^2),
      arithmetic_mean = sum((arithmetic_vec - ground_truth_vec)^2)
    )
  }
  
  # 4. meta sublist
  meta <- list(
    n_donors       = n,
    shock_time     = shock_time_vec,
    shock_length   = shock_length_vec,
    dbw_status     = if (exists("dbw_output")) dbw_output[[2]] else NA,
    weights        = w_hat,
    omega_vec      = omega_star_hat_vec,
    combined_omega = as.numeric(omega_star_hat),
    arima_order    = order_of_arima,
    ic_used        = user_ic_choice
  )
  if (exists("omega_star_std_err_hat_vec")) {
    meta$omega_se <- omega_star_std_err_hat_vec
  }
  
  # 5. final output list
  output_list <- list(
    linear_combinations = w_hat,
    predictions         = predictions
  )
  if (!is.null(loss)) {
    output_list$loss <- loss
  }
  output_list$meta <- meta
  
  # 6. opitional
  if (plots) {
    plot_maker_synthprediction(
      Y_series_list,
      shock_time_vec,
      integer_shock_time_vec,
      shock_length_vec,
      raw_vec,
      w_hat,
      omega_star_hat_vec,
      adjusted_vec,
      arithmetic_vec,
      display_ground_truth = display_ground_truth_choice
    )
  }
  
  return(output_list)
}

