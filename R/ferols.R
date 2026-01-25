# ------------------------------------------------------------------------------
# A fixed effect robust Huber-M estimator
#
# Inspired the 'robtwfe' packge by David Veenman 
# (https://github.com/dveenman/robtwfe)
#
# Heavily draws from the fixest package and calls feols() as a workhose
# in the IRLS process. UI should mimic the one of fixest as close as possible-
# 
# See LICENSE file for license
# ------------------------------------------------------------------------------


# --- Helper functions ---------------------------------------------------------

# These are not be exported by the package

huber_k_from_eff <- function(eff_target) {
  eff_fun <- function(k) {
    # Asymptotic efficiency of Huber M-estimator under N(0,1) relative to OLS
    # See Page 27 in Maronna et al. Robust Statistics 
    D <- 2*(k^2 * (1 - pnorm(k)) + pnorm(k) - 0.5 - k*dnorm(k))
    N <- (pnorm(k) - pnorm(-k))^2
    N/D
  }
  stats::uniroot(function(k) eff_fun(k) - eff_target, lower = 0.1, upper = 10)$root
}

huber_psi <- function(z, k) ifelse(abs(z) <= k, z, k * sign(z))

huber_w   <- function(z, k) {
  abs_z <- pmax(abs(z), sqrt(.Machine$double.eps))
  pmin(1, k / abs_z)
}

huber_phi <- function(z, k) as.numeric(abs(z) <= k)

madn <- function(x) {
  stats::median(
    abs(x - stats::median(x, na.rm = TRUE)), na.rm = TRUE
  ) / stats::qnorm(0.75)
}

# Returns the group means at observation level
gmean_expand <- function(v, g) {
  g <- as.integer(g)
  rs <- rowsum(v, g, reorder = FALSE, na.rm = TRUE)
  n  <- as.numeric(tabulate(g, nbins = nrow(rs)))
  gm <- rs / n
  gm[g, , drop = FALSE]
}

# build FE dummies without intercept
# drop first level per FE to avoid collinearity
fe_dummies <- function(f) {
  f <- factor(f)
  mm <- stats::model.matrix(~ f - 1)
  if (ncol(mm) > 0) mm <- mm[, -1, drop = FALSE]
  mm
}

# ------------------------------------------------------------------------------
# Method of moments quantile estimator as introduced by
# Machado and Santos Silva (2019) 
# https://doi.org/10.1016/j.jeconom.2019.04.009
# Mirrors the approach of Ben Jann in his Stata moremeta code
# https://github.com/benjann/moremata/blob/master/source/mm_aqreg.mata
# ------------------------------------------------------------------------------

lad_mm_ms <- function(fml, data, tau = 0.5) {
  fml_split <- fixest:::fml_split(fml, raw = TRUE)
  if (length(fml_split) == 2) {
    fe_vars <- all.vars(fml_split[[2]])
    n_levels <- vapply(fe_vars, function(v) nlevels(factor(data[,v])), integer(1))
    idx_abs <- which.max(n_levels)
    main_fml <- fml_split[[1]]
    fe_fml <- fml_split[[2]]
    dummy_fes <- fe_vars[-idx_abs]
    if (length(dummy_fes) > 0) {
      add_terms <- paste0("i(", dummy_fes, ")", collapse = " + ")
      main_fml <- stats::update.formula(main_fml, paste(". ~ . +", add_terms))
    }
    absorb_fe <- fe_vars[idx_abs]
    fe_fml <- stats::as.formula(paste("~", absorb_fe))
    fml_loc <- stats::as.formula(
      paste(deparse(main_fml), "|", deparse(fe_fml[[2]]))
    )
  } else {fml_loc <- fml} 
  
  # --- Estimate location ------------------------------------------------------
  fit_loc <- fixest::feols(fml_loc, data, warn = FALSE, notes = FALSE)
  yhat_loc <- as.numeric(stats::fitted(fit_loc, na.rm = FALSE))
  y <- as.numeric(data[, as.character(fml[[2]]), drop = TRUE])
  e <- y - yhat_loc
  Ipos <- as.numeric(e >= 0)
  Ibar <- mean(Ipos, na.rm = TRUE)
  data$r_raw <-  2 * e * (Ipos - Ibar)

  # --- Estimate scale ---------------------------------------------------------
  fml_scl <- fml_loc
  fml_scl[[2]] <- substitute(r_raw)
  fit_scl <- fixest::feols(fml_scl, data, warn = FALSE, notes = FALSE)
  denom <- as.numeric(stats::fitted(fit_scl, na.rm = FALSE))
  
  # Avoiding numerical issues causing scale estimates to become zero 
  eps <- sqrt(.Machine$double.eps)
  denom <- pmax(denom, eps)
  
  u <- e / denom
  qhat <- as.numeric(stats::quantile(
    u, probs = tau, type = 7, names = FALSE, na.rm = TRUE
  ))
  resid_tau <- y - (yhat_loc + qhat * denom)
  df_initial <- fixest::degrees_freedom(fit_loc, type = "resid") 
  n <- fit_loc$nobs
  pprob <- (2 * n - df_initial) / (2 * n)
  
  s <- as.numeric(
    stats::quantile(
      abs(resid_tau), probs = pprob, type = 7, names = FALSE, na.rm = TRUE)
    ) / stats::qnorm(0.75)
  
  bL <- stats::coef(fit_loc)
  bS <- stats::coef(fit_scl)
  bT <- bL + bS * qhat
  
  list(
    b = as.numeric(bT[names(bL) %in% all.vars(fml_split[[1]])]),
    r = resid_tau,
    s = s
  )
}

# ------------------------------------------------------------------------------
# Following methods of moment approach as in 
# Rios-Avila, Siles, and Canavire-Bacarreza  (2024)
# https://docs.iza.org/dp17262.pdf)
# to absorb multiple fixed effects. Generates scale, beta, and residual
# estimates that are very 
# ------------------------------------------------------------------------------

lad_mm_rsc <- function(fml, data, tau = 0.5) {
  # --- Estimate location ------------------------------------------------------
  fit_loc <- fixest::feols(fml, data, warn = FALSE, notes = FALSE)
  yhat_loc <- as.numeric(stats::fitted(fit_loc, na.rm = FALSE))
  y <- as.numeric(data[, as.character(fml[[2]]), drop = TRUE])
  e <- y - yhat_loc
  Ipos <- as.numeric(e >= 0)
  Ibar <- mean(Ipos, na.rm = TRUE)
  data$r_raw <-  2 * e * (Ipos - Ibar)
  
  # --- Estimate scale ---------------------------------------------------------
  fml_scl <- fml
  fml_scl[[2]] <- substitute(r_raw)
  fit_scl <- fixest::feols(fml_scl, data, warn = FALSE, notes = FALSE)
  denom <- as.numeric(stats::fitted(fit_scl, na.rm = FALSE))
  
  # Avoiding numerical issues causing scale estimates to become zero 
  eps <- sqrt(.Machine$double.eps)
  denom <- pmax(denom, eps)

  u <- e / denom
  qhat <- as.numeric(stats::quantile(
    u, probs = tau, type = 7, names = FALSE, na.rm = TRUE
  ))
  resid_tau <- y - (yhat_loc + qhat * denom)
  df_initial <- fixest::degrees_freedom(fit_loc, type = "resid") 
  n <- fit_loc$nobs
  pprob <- (2 * n - df_initial) / (2 * n)
  
  s <- as.numeric(
    stats::quantile(
      abs(resid_tau), probs = pprob, type = 7, names = FALSE, na.rm = TRUE)
    ) / stats::qnorm(0.75)
  
  bL <- stats::coef(fit_loc)
  bS <- stats::coef(fit_scl)
  bT <- bL + bS * qhat
  
  list(
    b = as.numeric(bT),
    r = resid_tau,
    s = s
  )
}


# --- Exported functions -------------------------------------------------------

#' Fixed-effects robust IRLS M regression (Huber loss)
#'
#' @description
#' Estimates linear models with high-dimensional fixed effects using 
#' iteractively reweighted least squares (IRLS) for Huber M-estimation. 
#' Fixed effects are absorbed using `fixest::demean`.
#' Standard errors can be clustered via one dimension using a Huber (psi/phi) 
#' sandwich estimator.
#'
#' @param fml A fixest formula of the form `y ~ x1 + x2 | fe1 + fe2`.
#' @param data A data.frame.
#' @param family Robust loss family. Currently only `"huber"`.
#' @param scale_est Algorithm to estimate the residuals that are used to 
#'  estimate the scale. Defaults to `"lad_mm_rsc"`.
#'    - `"ols"`: Uses `fixest:feols()` with absorbed fixed effects to derive 
#'      the initial residuals.
#'    - `"lad_mm_ms"`: Uses an method of moments algorithm for the median quantile 
#'      regression introduced by Machado and Santos Silva 
#'      (2019, https://doi.org/10.1016/j.jeconom.2019.04.009). 
#'      This should yield similar scale and beta estimates as the approaches of 
#'      the Stata packages `robreg` (https://github.com/benjann/robreg) and 
#'      `robtwfe` (https://github.com/dveenman/robtwfe). 
#'      This algorithm is currently work-in-process.
#'    - `"lad_mm_rsc"`: Method of moment algorithm as extended by
#'      Rios-Avila, Siles, and Canavire-Bacarreza 
#'      (2024, https://docs.iza.org/dp17262.pdf). This algorithm is
#'      faster as it absorbs multi-dimensional fixed effects. Its scale and beta
#'      estimates should be virtually identical (within numerical precision) to 
#'      the ones generated by `"lad_mm_ms"`.
#' @param scale A positive numerical value to override the estimated scale.
#' @param efficiency Target normal-efficiency in (0.68, 1). Default 0.95.
#' @param tol Convergence tolerance on relative coefficient change.
#' @param max_iter Maximum IRLS iterations.
#' @param vcov Variance-covariance specification. Either a string 
#'  (`"iid"` or `"cluster"`) or a formula (`cluster ~ var`) to indicate the 
#'  how it should be estimated. 
#'   -  `"iid"` for non-clustered Huber sandwich estimator (classical vcov, 
#'      default when no fixed effects are present).
#'   -  `"cluster"` to cluster on the first fixed effect (default when fixed 
#'      effects are present).
#'   -  `cluster ~ id` or simply `~id` for one-way clustering by the specified
#'      variable.
#'  All other features provided by `fixest` are not implemented (yet).
#' @param cluster (deprecated) A one-way clustering variable name (string).
#'  Note that this argument is deprecated, you should use `vcov` instead.
#' @param se (deprecated) A string ("standard", or "cluster") to indicated how
#'  the standard errors should be calculated. Note that this argument is 
#'  deprecated, you should use `vcov` instead.
#' @param ssc Small sample correction for constructing the covariance matrix. 
#'  The code uses a default SSC of the type (G / (G - 1)) * ((n - 1) / (n - p)). 
#'  This correction cannot be modified, so this parameter has to stay `NULL`
#'  for the time being.
#' @param ... Further arguments passed to `fixest::feols`. The `weights` 
#'  argument cannot be used as it is set by the IRLS process.
#'
#' @return An object of class `"fixest"` with an additional list
#'   `robust` containing `k`, `scale`, `efficiency`, `converged`, `iter`,
#'   `weights`, and `phi`.
#'
#' @examples
#' df <- generate_panel_data(n_units = 50, n_time = 10, seed = 42)
#' ferols(y ~ x + z | i + t, vcov = ~i, data = df)
#' fixest::feols(y ~ x + z | i + t, vcov = ~i, data = df)
#' 
#' @export
ferols <- function(
    fml, data, family = "huber",  efficiency = 0.95, 
    scale_est = "lad_mm_rsc", scale = NULL, tol = 1e-10, max_iter = 200, 
    vcov = NULL, cluster = NULL, ssc = NULL, se = NULL, ...
) {
  if (!requireNamespace("fixest", quietly = TRUE)) {
    stop("Package 'fixest' is required to run 'ferols()'.")
  }
  
  # --- Argument checks --------------------------------------------------------
  dots <- list(...)
  if ("weights" %in% dots) {
    stop("Argument 'weights' is not supported by 'ferols()'.")
  }
  if (!is.character(family) || family != "huber") {
    stop("Only family = 'huber' is currently supported.")
  }
  if (!is.numeric(efficiency) || efficiency < 0.68 || efficiency > 1) {
    stop("efficiency must be in (0.68, 1).")
  }
  fml_split <- fixest:::fml_split(fml, raw = TRUE)
  if (length(fml_split) > 2) stop(
    "Three part formulae (including IVs) are not supported by ferols()."
  )
  allowed_scale_ests <- c("ols", "lad_mm_ms", "lad_mm_rsc")
  if (! scale_est %in% allowed_scale_ests) {
    stop(sprintf(
      paste(
        "Scale estimation alogorithm '%s' is not supported.", 
        "Supported algorithms: '%s'"
      ), scale_est, paste(allowed_scale_ests, collapse = "', '")
    ))
  } 
  
  # --- Initial fit to estimate scale and starting beta ------------------------
  if (scale_est == "ols")  {
    fit <- fixest::feols(fml, data = data, warn = FALSE, notes = FALSE, ...)
    r <- stats::residuals(fit, na.rm = FALSE)
    if (is.null(scale)) scale <- madn(r)
    beta_old <- stats::coef(fit)
  } else if (scale_est == "lad_mm_rsc") {
    brs <- lad_mm_rsc(fml, data)
    r <- brs$r
    if (is.null(scale)) scale <- brs$s
    beta_old <- brs$b
  } else if (scale_est == "lad_mm_ms") {
    brs <- lad_mm_ms(fml, data)
    r <- brs$r
    if (is.null(scale)) scale <- brs$s
    beta_old <- brs$b
  } else { # this should not happen if allowed_scale_ests is current
    stop(paste(
      sprintf("Unknown scale estimation: '%s'.", scale_est),
      "This should not happen. Please report this error."
    ))
  }

  if (length(r) == 0) stop(
    "Initial fit has no residuals; check formula/data."
  )
  if (anyNA(beta_old)) stop("Initial fit has NA coefficients (collinearity?).")
  if (!is.finite(scale) || scale <= 0) stop(sprintf(
    "Invalid scale %.3f. Scale needs to be strictly positive and finite.",
    scale
  ))

  k <- huber_k_from_eff(efficiency)
  z <- r / scale
  w <- huber_w(z, k)
  phi <- huber_phi(z, k)
  
  # --- IRLS loop --------------------------------------------------------------
  converged <- FALSE
  iter_done <- 0L
  
  for (it in seq_len(max_iter)) {
    iter_done <- it
    fit_new <- fixest::feols(
      fml, data = data, weights = w, warn = FALSE, notes = FALSE, ...
    )
    
    beta_new <- stats::coef(fit_new)
    if (anyNA(beta_new)) stop(
      "IRLS produced NA coefficients; check collinearity/data issues."
    )
    
    # relative max change
    diff <- max(abs(beta_new - beta_old) / pmax(1, abs(beta_old)))
    if (diff <= tol) {
      fit <- fit_new
      converged <- TRUE
      break
    }
    
    fit <- fit_new
    beta_old <- beta_new
    
    r <- stats::residuals(fit, na.rm = FALSE)
    z <- r / scale
    w <- huber_w(z, k)
    phi <- huber_phi(z, k)
  }
  
  if (!converged) stop(
    "IRLS convergence not achieved (increase max_iter or relax tol)."
  )

  #--- Create fit object and call vcov.ferols() --------------------------------
  fit$robust <- list(
    family = family,
    scale_est = scale_est,
    efficiency = efficiency,
    k = k,
    scale = scale,
    converged = converged,
    iter = iter_done,
    weights = w,
    phi = phi
  )
  
  fit$call <- match.call()             
  fit$method <- "ferols"               
  class(fit) <- unique(c("ferols", class(fit)))
  
  V0 <- vcov.ferols(fit, vcov = "iid", cluster = NULL, se = NULL, ssc = NULL)
  V <- vcov(fit, vcov = vcov, cluster = cluster, se = se, ssc = ssc)
  
  # store vcov (names already set in vcov.ferols)
  fit$cov.iid <- V0
  fit$cov.scaled <- V
  
  # compute se + coeftable like fixest does
  se0 <- sqrt(diag(V))
  names(se0) <- names(coef(fit))
  attr(se0, "type") <- attr(V, "type")
  
  # df for t-stats: use fixest df if available; otherwise fall back
  df_r <- tryCatch(fit$fixest$df.residual, error = function(e) NULL)
  if (is.null(df_r)) df_r <- fit$df.residual
  if (is.null(df_r) || !is.finite(df_r)) {
    df_r <- length(fit$residuals) - length(coef(fit))
  }
  
  tval <- coef(fit) / se0
  pval <- 2 * stats::pt(abs(tval), df = df_r, lower.tail = FALSE)
  
  ct <- cbind(
    Estimate = coef(fit),
    `Std. Error` = se0,
    `t value` = tval,
    `Pr(>|t|)` = pval
  )
  
  attr(ct, "type") <- attr(se0, "type")
  
  fit$se <- se0
  fit$coeftable <- ct
  
  fit
}


#' Variance–covariance matrix for \code{ferols} objects
#'
#' @description
#' Computes variance–covariance matrices for estimates generarted by
#' \code{\link{ferols}}. The interface mirrors \code{fixest::vcov.fixest},
#' but currently supports only Huber M-estimation sandwich estimators
#' with either no clustering or one-way clustering.
#'
#' @param object A \code{ferols} object.
#' @param vcov Variance–covariance specification. Either \code{"iid"},
#'   \code{"cluster"}, or a one-sided formula such as \code{~ id} or
#'   \code{cluster ~ id}. Defaults to clustering on the first fixed effect
#'   if present.
#' @param cluster Deprecated. Use \code{vcov} instead.
#' @param se Deprecated. Use \code{vcov} instead.
#' @param ssc Small-sample correction. Not configurable yet and 
#'  must be \code{NULL}.
#' @param ... Additional arguments (currently ignored).
#'
#' @return A variance–covariance matrix with coefficient names as row and
#'   column names. The matrix carries a \code{"type"} attribute describing
#'   the estimator.
#'
#' @details
#' This function is experimental and implements only a subset of the
#' variance estimators available in \code{fixest}. Multi-way clustering,
#' alternative small-sample corrections, and user-supplied covariance
#' estimators are not yet supported.
#'
#' @export
vcov.ferols <- function(
  object, vcov = NULL, cluster = NULL, se = NULL, ssc = NULL, ...
) {
  vcov <- fixest:::oldargs_to_vcov(se, cluster, vcov) 
  
  # --- Argument checks --------------------------------------------------------
  if (!is.null(ssc)) {
    stop("You set 'ssc' but small sample correction is not yet configurable")
  }
  
  if(is.null(vcov)) {vcov <- "iid"}
  
  supported_char <- c("iid", "cluster")
  if (is.character(vcov)) {
    if (!(vcov %in% supported_char)) {
      stop(sprintf(
        "VCOV type '%s' is not supported yet in ferols(). Supported: %s.",
        vcov, paste(supported_char, collapse = ", ")
      ))
    }
    if (vcov == "cluster") {
      if (length(object$fixef_vars) < 1) stop(paste(
        "vcov = 'cluster' but no fixed effects presents. Consider setting",
        "the desired cluster explicitly by providind 'vcov = cluster ~ var'"
      )) else {
        vcov <- as.formula(sprintf("cluster ~ %s", object$fixef_vars[1]))
      }
    }
  }
  
  cl_var <- NULL
  if (inherits(vcov, "formula")) {
    # If there is no lhs, fixest defaults to "cluster"
    if (length(vcov) == 2) {
      vcov <- stats::as.formula(paste0("cluster ~ ", deparse(vcov[[2]])))
    }
    # for now, support only cluster ~ var1 meaning clustering
    if (all.vars(vcov)[1] != "cluster") stop(
      "Currently only 'cluster' is supported as the lhs of a vcov formula"
    ) 
    cl_vars <- fixest:::fml2varnames(vcov, combine_fun = TRUE)
    if (length(cl_vars) < 1 || length(cl_vars) > 1) {
      stop("ferols() currently supports clustering on one variable only.")
    } else cl_var <- cl_vars
  }
  
  please_use_str <- sprintf(
    "Please use instead a formula or vcov='keyword' with keyword in '%s'.",
    paste(supported_char, collapse = "', '")
  )
  
  if (inherits(vcov, "fixest_vcov_request")) {
    stop(paste(
      "ferols() does not support vcov objects.", please_use_str
    ))
  }
  
  if (is.function(vcov) || is.matrix(vcov) || is.list(vcov)) {
    stop(paste(
      "ferols() does not support custom vcov functions/matrices.", 
      please_use_str
    ))
  }
  
  # --- Construct variace/covariance matrix ------------------------------------
  if (is.character(vcov) && identical(vcov, "iid")) {
    obj <- object
    class(obj) <- setdiff(class(obj), "ferols")
    V <- stats::vcov(obj, vcov = "iid", ...)
  } else {
    # estimation sample data
    mf <- fixest::fixest_data(object, sample = "est")
    
    # FE variables (public)
    fe_vars <- object$fixef_vars
    fe_df <- mf[, fe_vars, drop = FALSE]
    
    # cluster ids
    cl <- mf[[cl_var]]
    g <- as.integer(factor(cl))
    G <- length(unique(g))
    if (G < 2) stop("Need at least 2 clusters.")
    
    # slope design matrix (public, stable)
    X <- stats::model.matrix(object, type = "rhs")
    
    # Get parameters from estimation object
    phi <- object$robust$phi[! is.na(object$robust$phi)]
    scale <- object$robust$scale
    k <- object$robust$k
    r <- stats::residuals(object)
    
    # residualize slopes wrt FE under weights = phi
    Xr <- fixest::demean(X, f = fe_df, weights = phi, na.rm = TRUE)
    z <- r / scale
    psi <- huber_psi(z, k)
    
    XtPhiX <- crossprod(Xr * phi, Xr)
    XtPhiX_inv <- solve(XtPhiX)
    
    kdim <- ncol(Xr)
    M <- matrix(0, kdim, kdim)
    for (gg in seq_len(G)) {
      idx <- which(g == gg)
      xg <- Xr[idx, , drop = FALSE]
      pg <- psi[idx]
      v <- crossprod(xg, pg)
      M <- M + v %*% t(v)
    }
    
    V <- scale^2 * XtPhiX_inv %*% M %*% XtPhiX_inv
    
    # CR1-ish finite sample correction
    n <- object$nobs
    temp_obj <- object
    class(temp_obj) <- setdiff(class(temp_obj), "ferols")
    k <- fixest::degrees_freedom(temp_obj, vcov = vcov, type = "k") 
    V <- (G / (G - 1)) * ((n - 1) / (n - k)) * V
    
    colnames(V) <- colnames(X)
    rownames(V) <- colnames(X)
  }
  
  if (is.character(vcov) && tolower(vcov) == "iid") {
    attr(V, "type") <- "IID"
  } else {
    attr(V, "type") <- paste0(
      "Clustered (", paste(cl_vars, sep = " & "), ")"
    )
  }
  
  V
}


#' Summarize a \code{ferols} model
#'
#' @description
#' Produces a coefficient table and inference summary for a model estimated
#' with \code{\link{ferols}}. The method closely follows
#' \code{fixest::summary.fixest} and supports the same \code{vcov} interface,
#' subject to the limitations of \code{ferols}.
#'
#' @param object A \code{ferols} object.
#' @param vcov Optional variance–covariance specification. If omitted,
#'   the covariance matrix stored in the object is used.
#' @param cluster Deprecated. Use \code{vcov} instead.
#' @param se Deprecated. Use \code{vcov} instead.
#' @param ssc Small-sample correction. Not configurable yet and 
#'  must be \code{NULL}.
#' @param ... Additional arguments passed to \code{fixest::summary.fixest}.
#'
#' @return An object of class \code{"summary.fixest"} containing coefficient
#'   estimates, standard errors, test statistics, and p-values. The standard
#'   errors and coefficient table include a \code{"type"} attribute describing
#'   the variance estimator used.
#'
#' @export
summary.ferols <- function(
  object, vcov = NULL, cluster = NULL, se = NULL, ssc = NULL, ...
) {
  # If user supplies vcov/cluster here, use it; otherwise use our vcov(object)
  V <- if (
    !is.null(vcov) || !is.null(cluster) || !is.null(se) || !is.null(ssc)
  ) {
    vcov(object, vcov = vcov, cluster = cluster, se = se, ssc = ssc, ...)
  } else {
    V <- object$cov.scaled
  }
  
  # delegate to fixest summary
  obj <- object
  class(obj) <- setdiff(class(obj), "ferols")
  s <- summary(obj, vcov = V, ssc = ssc, ...)
  
  # Restore our vcoc
  s$cov.scaled <- V
  
  # Set SE attributes
  attr(s$se, "type") <- attr(V, "type")
  attr(s$coeftable, "type") <- attr(V, "type")
  class(s) <- c("ferols", class(s))
  s
}


#' Print method for \code{ferols} objects
#'
#' @description
#' Prints a short model header followed by the coefficient table obtained
#' from \code{\link{summary.ferols}}. The output format is intentionally
#' similar to \code{fixest} model objects.
#'
#' @param x A \code{ferols} object.
#' @param ... Additional arguments passed to \code{print.summary.fixest}.
#'
#' @return Invisibly returns the model object.
#'
#' @export
print.ferols <- function(x, ...) {
  cat("ferols() - Fixed-effects robust IRLS M regression (Huber loss)\n")
  if (!is.null(x$robust)) {
    cat(sprintf(
      paste(
        "  efficiency: %.1f%% | k: %.4f | scale est: %s",
        "| scale: %.4g | iter: %d\n"
      ),
      100 * x$robust$efficiency, x$robust$k, x$robust$scale_est,
      x$robust$scale, x$robust$iter
    ))
  }
  s <- summary.ferols(x)
  class(s) <- setdiff(class(s), "ferols")
  print(s, ...)
  invisible(x)
}