#' Synthetic Volatility Forecast (GARCH + Donor Shock Adjustment)
#'
#' @description
#' Produces k-step-ahead **variance** forecasts for a target series using a
#' pre-shock GARCH model (no xreg on the target), then **adjusts** that forecast
#' by a donor-based post-shock fixed effect. Each donor is fit with a GARCH-X
#' including a post-shock indicator; the indicator coefficient is combined with
#' donor weights (uniform or DBW) to yield the final adjustment.
#'
#' @param Y_series_list List of series; element 1 is the target;
#'   elements 2 through n+1 are donors.
#' @param covariates_series_list List of data.frames/matrices for **donors'**
#' covariates (length = number of donors). Each row aligns with the donor series.
#' @param target_covariates data.frame/matrix of **target** covariates (same
#' columns as donor covariates). Used by DBW only (the target GARCH does not use xreg).
#' @param shock_time_vec Integer vector: row index of the shock for each series
#' (length = \code{length(Y_series_list)}).
#' @param shock_length_vec Integer vector: how many rows the post-shock
#' indicator stays at 1 after the shock (per series).
#' @param k Integer forecast horizon (steps after the pre-shock fit end).
#' @param dbw_scale,dbw_center Logical flags forwarded to \code{\link{dbw}}.
#' @param dbw_indices Integer vector: columns used by \code{dbw}. If \code{NULL}
#' uses all columns of \code{target_covariates}.
#' @param princ_comp_input Optional (currently not used here; kept for symmetry).
#' @param covariate_indices Integer vector: which donor covariate columns are
#' allowed into donor GARCH-X \code{xreg}. The target GARCH **never** uses xreg here.
#' @param penalty_lambda Numeric penalty weight forwarded to \code{\link{dbw}}.
#' @param penalty_normchoice Character, one of \code{"l1"} or \code{"l2"};
#' forwarded to \code{\link{dbw}}.
#' @param plots Logical; if \code{TRUE} and a function named
#' \code{plot_maker_garch} is found via \code{get0()}, it is called to visualize
#' results.
#' @param ground_truth_vec Optional numeric vector of realized variances (same
#' length as \code{k}) for loss computation (QLIKE).
#' @param shock_time_labels Optional labels for plotting (passed to plot hook).
#' @param return_fits Logical; if \code{TRUE}, return donor/target fit objects
#' for debugging.
#' @param garch_order Optional fixed GARCH-X order \code{c(q,p,o)}. If
#' \code{NULL}, the order is chosen by \code{\link{auto_garchx}}.
#' @param max.p,max.q Integers: maximum AR/MA orders for \code{auto_garchx}.
#' Ignored when \code{garch_order} is provided.
#' @param backcast.initial Optional numeric scalar backcast for \code{garchx}.
#' Used only when \code{garch_order} is fixed.
#' @param vcov.type Character, \code{"ordinary"} or \code{"robust"}; covariance
#' type used by \code{garchx} and \code{lmtest::coeftest}.
#'
#' @return A list with elements:
#' \itemize{
#'   \item \code{linear_combinations}: donor weights (numeric length n_donors).
#'   \item \code{predictions}: a list with numeric vectors \code{unadjusted},
#'         \code{adjusted}, \code{arithmetic_mean} (all length k).
#'   \item \code{meta}: diagnostic info (weights, per-donor omega, selected orders, etc.).
#'   \item \code{loss}: optional list with QLIKE losses if \code{ground_truth_vec} is provided.
#'   \item \code{donor_fits}, \code{target_fit}: returned only if \code{return_fits=TRUE}.
#' }
#'
#' @details
#' • Target: pre-shock GARCH **without** xreg → no need for future \code{newxreg}.  
#' • Donors: GARCH-X with a 0/1 post-shock indicator (and optional donor xreg).  
#' • DBW uses \code{c(target_covariates, covariates_series_list)} to compute weights.
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' y_tgt <- rnorm(200, 0, 1)
#' y_d1  <- rnorm(200, 0, 1)
#' y_d2  <- rnorm(200, 0, 1)
#' Y <- list(y_tgt, y_d1, y_d2)
#'
#' Xd1 <- cbind(x1 = rnorm(200), x2 = rnorm(200))
#' Xd2 <- cbind(x1 = rnorm(200), x2 = rnorm(200))
#' Xt  <- cbind(x1 = rnorm(200), x2 = rnorm(200))
#'
#' shock_t <- c(150L, 150L, 150L)
#' shock_L <- c(10L,  10L,  10L)
#'
#' out <- SynthVolForecast(
#'   Y_series_list = Y,
#'   covariates_series_list = list(Xd1, Xd2),
#'   target_covariates = Xt,
#'   shock_time_vec = shock_t,
#'   shock_length_vec = shock_L,
#'   k = 3,
#'   covariate_indices = 1:2,  # donors may use xreg
#'   dbw_indices = 1:2,        # DBW uses covariates
#'   plots = FALSE
#' )
#' str(out$predictions)
#' }
#'
#' @importFrom garchx garchx
#' @export
SynthVolForecast <- function(
    Y_series_list,
    covariates_series_list,
    target_covariates,
    shock_time_vec,
    shock_length_vec,
    k = 1,
    dbw_scale = TRUE,
    dbw_center = TRUE,
    dbw_indices = NULL,
    princ_comp_input = NULL,
    covariate_indices = NULL,
    penalty_lambda = 0,
    penalty_normchoice = "l1",
    plots = TRUE,
    ground_truth_vec = NULL,
    shock_time_labels = NULL,
    return_fits = FALSE,
    garch_order = NULL,
    max.p = 3,
    max.q = 3,
    backcast.initial = NULL,
    vcov.type = "robust"
) {
  # ---- sanitize args ---------------------------------------------------------
  k  <- max(1L, as.integer(k))
  vc <- match.arg(tolower(as.character(vcov.type)), c("ordinary","robust"))
  
  n_total  <- length(Y_series_list)
  n_donors <- n_total - 1L
  
  if (!vc %in% c("ordinary","robust"))
    stop("vcov.type must be 'ordinary' or 'robust'.")
  if (!is.null(backcast.initial)) {
    stopifnot(is.numeric(backcast.initial), length(backcast.initial) == 1,
              is.finite(backcast.initial), backcast.initial > 0)
  }
  
  # ---- FAST PATH: no donors --------------------------------------------------
  if (n_donors == 0L) {
    target_Y   <- Y_series_list[[1L]]
    target_end <- as.integer(shock_time_vec[1L])
    
    # Target NEVER uses xreg in this function (to avoid needing future newxreg)
    target_xreg <- NULL
    
    if (is.null(garch_order)) {
      sel_t <- auto_garchx(
        y = target_Y[1:target_end], xreg = target_xreg,
        max.p = max.p, max.q = max.q,
        vcov.type = vc
      )
      target_order <- c(sel_t$q, sel_t$p, sel_t$o)
      fit_target   <- sel_t$fit
    } else {
      target_order <- garch_order
      fit_target   <- garchx::garchx(
        y = target_Y[1:target_end],
        order = target_order,
        xreg  = target_xreg,
        vcov.type = vc,
        backcast.values = if (is.null(backcast.initial)) NULL else backcast.initial
      )
    }
    
    raw_pred <- stats::predict(fit_target, n.ahead = k)
    raw_vec  <- as.numeric(if (is.list(raw_pred)) raw_pred$pred else raw_pred)
    
    return(list(
      linear_combinations = numeric(0),
      predictions = list(
        unadjusted      = raw_vec,
        adjusted        = raw_vec,
        arithmetic_mean = raw_vec
      ),
      meta = list(
        n_donors       = 0L,
        weights        = numeric(0),
        omega_vec      = numeric(0),
        omega_se       = numeric(0),
        combined_omega = 0,
        shock_time     = shock_time_vec,
        shock_length   = shock_length_vec,
        target_order   = target_order
      ),
      target_fit = fit_target
    ))
  }
  
  # ---- donor path checks -----------------------------------------------------
  stopifnot(length(shock_time_vec)   == n_total)
  stopifnot(length(shock_length_vec) == n_total)
  stopifnot(length(covariates_series_list) == n_donors)
  
  # helper to coerce to data.frame
  as_df <- function(x) {
    if (is.numeric(x)) data.frame(V1 = x)
    else if (is.matrix(x)) as.data.frame(x)
    else if (inherits(x, "data.frame")) x
    else stop("covariates must be numeric/matrix/data.frame")
  }
  
  target_covariates <- as_df(target_covariates)
  for (i in seq_len(n_donors)) {
    covariates_series_list[[i]] <- as_df(covariates_series_list[[i]])
  }
  
  if (is.null(dbw_indices))      dbw_indices      <- seq_len(ncol(target_covariates))
  if (is.null(princ_comp_input)) princ_comp_input <- min(length(shock_time_vec),
                                                         ncol(target_covariates))
  
  # QLIKE loss if needed
  if (!exists("QL_loss_function", mode = "function")) {
    QL_loss_function <- function(vhat, y) log(vhat) + (y^2)/vhat
  }
  has_lmtest <- requireNamespace("lmtest", quietly = TRUE)
  
  # ---- per-donor fit and omega extraction -----------------------------------
  omega_star_hat_vec     <- numeric(n_donors)
  omega_star_std_err_vec <- rep(NA_real_, n_donors)
  donor_fits             <- vector("list", n_donors)
  donor_selected_orders  <- vector("list", n_donors)
  
  for (i in seq_len(n_donors)) {
    donor_Y <- Y_series_list[[i+1L]]
    donor_X <- covariates_series_list[[i]]
    
    len_i   <- as.integer(shock_length_vec[i+1L])
    start_i <- as.integer(shock_time_vec[i+1L])
    end_i   <- min(start_i + len_i, length(donor_Y))
    post_shock_indicator <- rep(0L, length(donor_Y))
    if (len_i > 0L && end_i >= start_i + 1L) {
      post_shock_indicator[(start_i+1L):end_i] <- 1L
    }
    
    rows_to_use <- length(donor_Y)
    if (is.null(covariate_indices)) {
      X_i_final <- matrix(
        post_shock_indicator[1:rows_to_use],
        ncol = 1,
        dimnames = list(NULL, "post_shock_indicator")
      )
    } else {
      X_cov     <- as.matrix(donor_X[1:rows_to_use, covariate_indices, drop = FALSE])
      X_i_final <- cbind(X_cov, post_shock_indicator[1:rows_to_use])
      colnames(X_i_final)[ncol(X_i_final)] <- "post_shock_indicator"
      X_i_final <- stats::na.omit(X_i_final)
    }
    
    # fit donor model
    if (is.null(garch_order)) {
      sel <- auto_garchx(
        y     = donor_Y[1:rows_to_use],
        xreg  = X_i_final,
        max.p = max.p,
        max.q = max.q,
        vcov.type = vc
      )
      order_i   <- c(sel$q, sel$p, sel$o)
      donor_fit <- sel$fit
    } else {
      order_i   <- garch_order
      donor_fit <- garchx::garchx(
        y = donor_Y[1:rows_to_use],
        order = order_i,
        xreg  = X_i_final,
        vcov.type = vc,
        backcast.values = if (is.null(backcast.initial)) NULL else backcast.initial
      )
    }
    
    donor_selected_orders[[i]] <- order_i
    donor_fits[[i]]            <- donor_fit
    
    # omega = coefficient on the post-shock indicator
    if (has_lmtest) {
      vcov_fun <- function(obj) stats::vcov(obj, vcov.type = vc)
      ct <- lmtest::coeftest(donor_fit, vcov. = vcov_fun)
      hit <- grep("post_shock_indicator", rownames(ct), fixed = TRUE)
      omega_star_hat_vec[i]     <- ct[hit[1L], "Estimate"]
      if ("Std. Error" %in% colnames(ct)) {
        omega_star_std_err_vec[i] <- ct[hit[1L], "Std. Error"]
      }
    } else {
      cf  <- coef(donor_fit)
      nm  <- names(cf)
      hit <- grep("post_shock_indicator", nm, fixed = TRUE)
      stopifnot(length(hit) >= 1L)
      omega_star_hat_vec[i] <- cf[hit[1L]]
    }
  }
  
  # ---- DBW weights (target + donors) ----------------------------------------
  X_for_dbw <- c(list(target_covariates), covariates_series_list)
  dbw_output <- dbw(
    X                = X_for_dbw,
    dbw_indices      = dbw_indices,
    shock_time_vec   = shock_time_vec,
    scale            = dbw_scale,
    center           = dbw_center,
    sum_to_1         = TRUE,
    bounded_below_by = 0,
    bounded_above_by = 1,
    normchoice       = "l2",
    penalty_normchoice = penalty_normchoice,
    penalty_lambda     = penalty_lambda
  )
  w_hat <- as.numeric(dbw_output$opt_params)
  stopifnot(length(w_hat) == length(omega_star_hat_vec))
  omega_star_hat <- sum(w_hat * omega_star_hat_vec)
  
  # ---- target pre-shock model and forecast ----------------------------------
  target_Y   <- Y_series_list[[1L]]
  target_end <- as.integer(shock_time_vec[1L])
  
  # Target NEVER uses xreg in this function (avoid needing future newxreg)
  target_xreg <- NULL
  
  if (is.null(garch_order)) {
    sel_t <- auto_garchx(
      y     = target_Y[1:target_end],
      xreg  = target_xreg,
      max.p = max.p,
      max.q = max.q,
      vcov.type = vc
    )
    target_order <- c(sel_t$q, sel_t$p, sel_t$o)
    fit_target   <- sel_t$fit
  } else {
    target_order <- garch_order
    fit_target   <- garchx::garchx(
      y = target_Y[1:target_end],
      order = target_order,
      xreg  = target_xreg,
      vcov.type = vc,
      backcast.values = if (is.null(backcast.initial)) NULL else backcast.initial
    )
  }
  
  raw_pred <- stats::predict(fit_target, n.ahead = k)
  raw_vec  <- as.numeric(if (is.list(raw_pred)) raw_pred$pred else raw_pred)
  
  adjusted_vec   <- raw_vec + omega_star_hat
  arithmetic_vec <- raw_vec + mean(omega_star_hat_vec)
  
  predictions <- list(
    unadjusted      = raw_vec,
    adjusted        = adjusted_vec,
    arithmetic_mean = arithmetic_vec
  )
  
  loss <- if (!is.null(ground_truth_vec)) {
    list(
      unadjusted      = sum(QL_loss_function(raw_vec,        ground_truth_vec)),
      adjusted        = sum(QL_loss_function(adjusted_vec,   ground_truth_vec)),
      arithmetic_mean = sum(QL_loss_function(arithmetic_vec, ground_truth_vec))
    )
  } else NULL
  
  meta <- list(
    n_donors       = n_donors,
    shock_time     = shock_time_vec,
    shock_length   = shock_length_vec,
    dbw_status     = dbw_output$convergence,
    weights        = w_hat,
    omega_vec      = omega_star_hat_vec,
    omega_se       = omega_star_std_err_vec,
    combined_omega = omega_star_hat,
    donor_orders   = donor_selected_orders,
    target_order   = target_order
  )
  
  out <- list(
    linear_combinations = w_hat,
    predictions         = predictions,
    meta                = meta
  )
  if (!is.null(loss)) out$loss <- loss
  if (return_fits) {
    out$donor_fits <- donor_fits
    out$target_fit <- fit_target
  }
  
  # ---- optional plot hook ----------------------------------------------------
  if (isTRUE(plots)) {
    f <- get0("plot_maker_garch", mode = "function", inherits = TRUE)
    if (!is.null(f)) {
      try(
        f(
          stats::fitted(fit_target),
          shock_time_labels,
          shock_time_vec,
          shock_length_vec,
          raw_vec,
          w_hat,
          omega_star_hat_vec,
          omega_star_std_err_vec,
          adjusted_vec,
          arithmetic_vec,
          ground_truth_vec
        ),
        silent = TRUE
      )
    }
  }
  
  out
}
