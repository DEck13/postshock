#' Helper Function: Build control-shock xreg matrix
#'
#' @description
#' Internal helper to convert shock specifications into a design matrix for ARIMA xreg.
#'
#' @param last Integer; length of the training period.
#' @param control_spec List or data.frame defining shock timing, length, and shape.
#'
#' @return A numeric matrix or NULL.
#' @export
build_control_xreg <- function(last, control_spec = NULL) {
  # 1. Basic checks
  if (is.null(last) || length(last) != 1L || !is.finite(last) || last < 1) {
    stop("`last` must be a single positive integer.")
  }
  last <- as.integer(last)
  
  # 2. Return NULL if no controls
  if (is.null(control_spec)) return(NULL)
  
  # 3. Normalize input to a list of specs
  specs <- NULL
  if (is.data.frame(control_spec)) {
    if (!all(c("time") %in% names(control_spec))) {
      stop("If `control_spec` is a data.frame, it must contain a `time` column.")
    }
    # Fill defaults
    if (!("length" %in% names(control_spec))) control_spec$length <- 1L
    if (!("shape"  %in% names(control_spec))) control_spec$shape  <- "point"
    
    specs <- lapply(seq_len(nrow(control_spec)), function(i) {
      list(
        time   = control_spec$time[[i]],
        length = control_spec$length[[i]],
        shape  = control_spec$shape[[i]]
      )
    })
  } else if (is.list(control_spec)) {
    # Check if it's a single spec list (has keys like time/length) or a list of lists
    is_single_spec <- !is.null(control_spec$time)
    if (is_single_spec) {
      specs <- list(control_spec)
    } else {
      specs <- control_spec
    }
  } else {
    stop("`control_spec` must be NULL, a list, or a data.frame.")
  }
  
  # 4. Filter empty specs
  specs <- Filter(function(s) !is.null(s) && length(s) > 0, specs)
  if (length(specs) == 0) return(NULL)
  
  # 5. Build Matrix Columns
  cols <- lapply(seq_along(specs), function(j) {
    s <- specs[[j]]
    
    t0  <- if (!is.null(s$time)) as.integer(s$time) else NA_integer_
    L   <- if (!is.null(s$length)) as.integer(s$length) else 1L
    shp <- if (!is.null(s$shape)) tolower(as.character(s$shape)) else "point"
    
    if (is.na(t0)) stop("Control shock must have a valid `time`.")
    if (is.na(L) || L < 1) stop("Control shock `length` must be >= 1.")
    if (!shp %in% c("point", "window", "step")) {
      stop("Unknown shock shape: ", shp)
    }
    
    v <- rep(0, last)
    
    if (shp == "point") {
      if (t0 >= 1 && t0 <= last) v[t0] <- 1
    } else if (shp == "window") {
      a <- max(1L, t0)
      b <- min(last, t0 + L - 1L)
      if (a <= b) v[a:b] <- 1
    } else if (shp == "step") {
      a <- max(1L, t0)
      if (a <= last) v[a:last] <- 1
    }
    v
  })
  
  C <- do.call(cbind, cols)
  if (is.null(C)) return(NULL) # Safety if cols was empty
  colnames(C) <- paste0("ctrl_shock_", seq_len(ncol(C)))
  storage.mode(C) <- "double"
  return(C)
}