#' Auto-select a GARCHX specification by BIC (with optional backcast refit)
#'
#' Grid-search over \eqn{(q, p, o)} orders for \pkg{garchx} and pick the model
#' with the lowest BIC; optionally refit the selected model using a backcast
#' initialized at the unconditional variance (clamped to a reasonable range).
#'
#' @param y Numeric vector of returns/residuals (no NA/Inf/NaN).
#' @param xreg Optional matrix of exogenous regressors (same number of rows as `y`).
#' @param max.p,max.q Nonnegative integers; AR (`p`) and MA (`q`) orders to search up to.
#' @param search_o Integer vector with values in `{0,1}`; 0 = symmetric, 1 = GJR asymmetry.
#' @param final_refit Logical; if `TRUE`, refit the winner using a backcast equal
#'   to the unconditional variance (with clamping).
#' @param clamp_factor Length-2 numeric vector `c(low, high)`. If the computed
#'   unconditional variance falls outside `[low, high] * var(y)`, fall back to
#'   the sample variance for the backcast.
#' @param vcov.type Character; `"ordinary"` or `"robust"` standard errors.
#' @param verbose Logical; print progress/info/warnings.
#'
#' @details
#' The order used by \pkg{garchx} is `c(q, p, o)` (MA, AR, asymmetry). We skip
#' the trivial `(0,0,0)` model. If no model converges in the grid, the function
#' stops with an error.
#'
#' @return A list with elements:
#' \itemize{
#'   \item `fit` : the selected `garchx` fit object
#'   \item `p`, `q`, `o` : selected orders
#'   \item `bic`, `aic` : information criteria of the selected fit
#'   \item `refit_used` : `TRUE/FALSE` if the final backcast refit replaced the grid one
#' }
#'
#' @examples
#' \donttest{
#'   set.seed(1)
#'   y <- scale(rnorm(500))[,1]
#'   out <- auto_garchx(y, max.p = 2, max.q = 2, search_o = c(0,1),
#'                      final_refit = TRUE, vcov.type = "robust", verbose = FALSE)
#'   c(out$p, out$q, out$o); out$bic
#' }
#'
#' @importFrom stats logLik coef var
#' @export
auto_garchx <- function(
    y,
    xreg         = NULL,
    max.p        = 3,
    max.q        = 3,
    search_o     = c(0L, 1L),      # 0 = symmetric, 1 = GJR
    final_refit  = TRUE,
    clamp_factor = c(0.1, 10),
    vcov.type    = c("ordinary","robust"),
    verbose      = TRUE
) {
  vcov.type <- match.arg(vcov.type)
  
  # ---- Basic input checks ----------------------------------------------------
  stopifnot(is.numeric(y), is.vector(y))
  if (any(!is.finite(y))) stop("auto_garchx: 'y' contains NA/Inf/NaN.")
  n <- length(y)
  if (n < 20 && verbose) warning("auto_garchx: very short series; estimation may be unstable.")
  
  if (!is.null(xreg)) {
    xreg <- as.matrix(xreg)
    if (nrow(xreg) != n) stop("auto_garchx: nrow(xreg) must equal length(y).")
    if (any(!is.finite(xreg))) stop("auto_garchx: 'xreg' contains NA/Inf/NaN.")
  }
  
  # Normalize search list for o (0 or 1 only)
  search_o <- as.integer(unique(search_o))
  search_o <- search_o[search_o %in% c(0L,1L)]
  if (!length(search_o)) search_o <- 0L
  
  # Holder for the best model by BIC
  best <- list(bic = Inf, aic = Inf, fit = NULL,
               p = NA_integer_, o = NA_integer_, q = NA_integer_)
  
  # ---- 1) Grid search over (q, p, o) ----------------------------------------
  # NOTE: garchx::garchx(order = c(q, p, o))
  for (q in 0:max.q) {
    for (p in 0:max.p) {
      for (o in search_o) {
        # skip the trivial model
        if (p == 0L && q == 0L && o == 0L) next
        
        fit0 <- try(
          garchx::garchx(y = y, order = c(q, p, o), xreg = xreg, vcov.type = vcov.type),
          silent = TRUE
        )
        if (inherits(fit0, "try-error")) next
        
        ll   <- as.numeric(logLik(fit0))
        k    <- length(coef(fit0))
        bic0 <- -2 * ll + k * log(n)
        aic0 <- -2 * ll + 2 * k
        
        if (is.finite(bic0) && bic0 < best$bic) {
          best <- list(bic = bic0, aic = aic0, fit = fit0, p = p, o = o, q = q)
        }
      }
    }
  }
  if (is.null(best$fit)) stop("auto_garchx: no converged model in grid search.")
  
  # If we don’t want the refit step, return the grid winner
  if (!isTRUE(final_refit)) return(best)
  
  # ---- 2) Refit using a backcast at unconditional variance ------------------
  # Extract coefficients from the grid winner to compute unconditional variance
  cf     <- coef(best$fit)
  om.idx <- grep("^(intercept|omega)$", names(cf), ignore.case = TRUE)
  if (!length(om.idx)) return(best)  # if we cannot detect omega, give up refit
  
  omega  <- as.numeric(cf[om.idx[1]])
  alpha  <- sum(cf[grep("^arch",  names(cf), ignore.case = TRUE)], na.rm = TRUE)
  beta   <- sum(cf[grep("^garch", names(cf), ignore.case = TRUE)], na.rm = TRUE)
  gamma  <- if (best$o == 1L) sum(cf[grep("asym|gamma|gjr", names(cf), ignore.case = TRUE)], na.rm = TRUE) else 0
  
  # Unconditional variance for GJR: omega / (1 - alpha - beta - 0.5*gamma)
  denom  <- 1 - alpha - beta - 0.5 * gamma
  if (!is.finite(omega) || denom <= 0) return(best)  # no valid refit
  
  uncond_var <- as.numeric(omega / denom)
  
  # Clamp to a reasonable range relative to sample variance
  sv  <- stats::var(y, na.rm = TRUE)
  rng <- range(clamp_factor)
  if (!is.finite(uncond_var) || uncond_var <= 0 ||
      uncond_var < rng[1] * sv || uncond_var > rng[2] * sv) {
    if (verbose) warning("[auto_garchx] backcast out of range; fallback to sample variance.")
    uncond_var <- sv
  }
  if (verbose) {
    message(sprintf("[auto_garchx] refit with backcast=%.6f (q=%d,p=%d,o=%d)",
                    uncond_var, best$q, best$p, best$o))
  }
  
  # Refit with backcast.values set to a single numeric (not a list)
  fit2 <- try(
    garchx::garchx(
      y = y,
      order = c(best$q, best$p, best$o),
      xreg  = xreg,
      backcast.values = uncond_var,
      vcov.type = vcov.type
    ),
    silent = TRUE
  )
  
  # If refit succeeded, keep it only if BIC improves
  if (!inherits(fit2, "try-error")) {
    ll2  <- as.numeric(logLik(fit2))
    k2   <- length(coef(fit2))
    bic2 <- -2 * ll2 + k2 * log(n)
    aic2 <- -2 * ll2 + 2 * k2
    
    if (is.finite(bic2) && bic2 < best$bic) {
      best$fit        <- fit2
      best$bic        <- bic2
      best$aic        <- aic2
      best$refit_used <- TRUE
    } else {
      best$refit_used <- FALSE
      if (verbose) warning("[auto_garchx] refit did not improve BIC; keeping grid-search model.")
    }
  } else if (verbose) {
    warning("[auto_garchx] backcast-refit failed; keeping grid-search model.")
  }
  
  best
}

#' @rdname auto_garchx
#' @export
auto_garch <- auto_garchx
