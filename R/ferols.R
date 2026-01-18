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
  out <- rep(1, length(z))
  nz <- (z != 0)
  out[nz] <- pmin(1, k / abs(z[nz]))
  out
}

huber_phi <- function(z, k) as.numeric(abs(z) <= k)

madn <- function(x) {
  stats::median(
    abs(x - stats::median(x, na.rm = TRUE)), na.rm = TRUE
  ) / stats::qnorm(0.75)
}

mmfeqr.brs <- function(y, X, fe = NULL) {
  y <- as.numeric(y)
  X <- as.matrix(X)
  n <- length(y)
  if (nrow(X) != n) stop("y and X must have compatible rows.")
  k <- ncol(X)
  if (k == 0) stop("no indepdendent variables present in model")

  # -------- Step 1: Absorb the FE with more levels and include the other ------
  if (is.null(fe)) {
    Xd = X
    yd = y
  } else {
    if (length(fe) > 1) {
      n_levels <- unlist(lapply(fe, function(x) nlevels(factor(x))))
      idx_absorbed_fe <- which(n_levels == max(n_levels)[1])
      dummy_fe <- fe[, -idx_absorbed_fe, drop = FALSE]
      for(col in 1:ncol(dummy_fe)) {
        X <- cbind(X, stats::model.matrix(~ factor(as.vector(dummy_fe[,col]))))
      }
      fe <- fe[, idx_absorbed_fe, drop = FALSE]
    }
    yd <- fixest::demean(y, fe)
    Xd <- fixest::demean(X, fe) 
  }
  
  # -------- Step 2: location LS on demeaned data (no intercept) ---------------
  # LS: solve (X'X)b = X'y
  b <- qr.coef(qr(Xd), yd)
  b[!is.finite(b)] <- 0
  e <- as.numeric(yd - Xd %*% b) # demeaned residual
  
  # -------- Step 3: scale proxy from "signed residual" regression -------------
  # r_i = 2 * e_i * ( 1{e_i>=0} - mean_w(1{e_i>=0}) )
  r_raw <- 2 * e * (as.numeric(e >= 0) - mean(e >= 0))
  
  # demean r_raw if FE are present
  if (!is.null(fe)) {
    rd <- fixest::demean(r_raw, fe)
  } else {rd <- r_raw}
  bs <- qr.coef(qr(Xd), rd)
  bs[!is.finite(bs)] <- 0
  
  # ------- Step 4: Calculate bt and r for tau = 0.5 ---------------------------
  shat <- as.numeric(X %*% bs)
  u <- e / shat
  q <- as.numeric(stats::quantile(u, probs = 0.5))
  bt <- b + bs*q
  r <- e - q * shat
  
  p <- (2*n - (n - k)) / (2*n)
  s <- stats::quantile(abs(r), probs = p) / stats::qnorm(0.75)
  list(r = r, b = bt[1:ncol(X)], s = s)
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
#'  estimate the scale. Defaults to `"ols"`.
#'    - `"ols"`: Uses `fixest:feols()` with absorbed fixed effects to derive 
#'      the initial residuals.
#'    - `"lad_rq"`: Uses `quantreg:rq.fit()` with absorbed fixed effects 
#'      to estimate a median quantile regression (aka least absolute deviation
#'      (LAD) regression) to derive the initial residuals.
#'    - `"lad_mm"`: Uses an method of moments algorithm for the median quantile 
#'      regression introduced by Machado and Santos Silva 
#'      (2019, https://doi.org/10.1016/j.jeconom.2019.04.009) and extended by
#'      Rios-Avila, Siles, and Canavire-Bacarreza 
#'      (2024, https://docs.iza.org/dp17262.pdf). This should yield similar
#'      scales as the approaches of the Stata packages `robreg` 
#'      (https://github.com/benjann/robreg) and `robtwfe` 
#'      (https://github.com/dveenman/robtwfe). This algorithm is currently
#'      woirk-in-process.
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
#' @export
ferols <- function(
    fml, data, family = "huber",  efficiency = 0.95, 
    scale_est = "ols", scale = NULL, tol = 1e-10, max_iter = 200, 
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
  allowed_scale_ests <- c("ols", "lad_rq", "lad_mm")
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
    fit <- fixest::feols(fml, data = data, ...)
    r <- stats::residuals(fit)
    if (is.null(scale)) scale <- madn(r)
    beta_old <- stats::coef(fit)
  } else {
    mf_main <- stats::model.frame(fml_split[[1]], data)
    y <- as.numeric(stats::model.response(mf_main))
    X <- stats::model.matrix(fml_split[[1]], data = mf_main)
    if (length(fml_split) == 2) { # FE present 
      if ("(Intercept)" %in% colnames(X)) {
        X <- X[, colnames(X) != "(Intercept)", drop = FALSE]
      }
      fe_vars <- all.vars(fml_split[[2]])
      ridx <- as.integer(rownames(mf_main))
      fe_df <- data[ridx, fe_vars, drop = FALSE]
    } else {
      fe_df <- NULL
      fe_vars <- NULL
    }
    if (scale_est == "lad_rq") {
      if (length(fml_split) == 2) { # FE present 
        y_tilde <- fixest::demean(y, f = fe_df, na.rm = TRUE)
        X_tilde <- fixest::demean(X, f = fe_df, na.rm = TRUE)
      } else {
        y_tilde <- y
        X_tilde <- X
      }
      fit <- quantreg::rq.fit(x = X_tilde, y = y_tilde)
      r <- stats::residuals(fit)
      if (is.null(scale)) scale <- madn(r)
      beta_old <- stats::coef(fit)
    } else { # scale_est == "lad_mm"
      brs <- mmfeqr.brs(y, X, fe_df)
      r <- brs$r
      if (is.null(scale)) scale <- brs$s
      beta_old <- brs$b
    }
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
    fit_new <- fixest::feols(fml, data = data, weights = w, ...)
    
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
    
    r <- stats::residuals(fit)
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
    phi <- object$robust$phi
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
    n <- nrow(Xr)
    p <- ncol(Xr)
    V <- (G / (G - 1)) * ((n - 1) / (n - p)) * V
    
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