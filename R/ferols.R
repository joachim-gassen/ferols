# ------------------------------------------------------------------------------
# A fixed effect robust Huber-M estimator
#
# Inspired by and co-developed with the 'robhdfe' packge by David Veenman 
# (https://github.com/dveenman/robhdfe)
#
# Heavily draws from the fixest package and calls feols() as a workhose
# in the IRLS process. UI should mimic the one of fixest as close as possible-
# 
# See LICENSE file for license
# ------------------------------------------------------------------------------


# --- Helper functions ---------------------------------------------------------

# These are not be exported by the package

capture_side_effects <- function(expr, quiet = TRUE) {
  msgs <- character()
  warns <- character()
  
  value <- withCallingHandlers(
    force(expr),
    message = function(m) {
      msgs <<- c(msgs, conditionMessage(m))
      if (quiet) invokeRestart("muffleMessage")
    },
    warning = function(w) {
      warns <<- c(warns, conditionMessage(w))
      if (quiet) invokeRestart("muffleWarning")
    }
  )
  list(value = value, messages = msgs, warnings = warns)
}

relay_side_effects <- function(x) {
  for (m in x$messages) message(m)
  for (w in x$warnings) warning(w, call. = FALSE)
}

huber_k_from_eff <- function(eff_target, tol = 1e-12) {
  eff_fun <- function(k) {
    # Asymptotic efficiency of Huber M-estimator under N(0,1) relative to OLS
    # See Page 27 in Maronna et al. Robust Statistics 
    D <- 2*(k^2 * (1 - pnorm(k)) + pnorm(k) - 0.5 - k*dnorm(k))
    N <- (pnorm(k) - pnorm(-k))^2
    N/D
  }
  stats::uniroot(
    function(k) eff_fun(k) - eff_target, 
    lower = 0.0001, upper = 10, tol = 1e-12
  )$root
}

huber_psi <- function(z, k) ifelse(abs(z) <= k, z, k * sign(z))

huber_w   <- function(z, k) {
  abs_z <- pmax(abs(z), sqrt(.Machine$double.eps))
  pmin(1, k / abs_z)
}

huber_phi <- function(z, k) as.numeric(abs(z) <= k)

madn <- function(x) {
  as.numeric(stats::quantile(
    abs(x - stats::quantile(x, 1/2, na.rm = TRUE, type = 2)), 
    1/2, type = 2, na.rm = TRUE
  )) / stats::qnorm(0.75)
}

# ------------------------------------------------------------------------------
# Simple LAD based on quantreg::rq. Used by default when no fixed effects are
# to be absorbed. Seems to be what Stata's robreg does as well.
# ------------------------------------------------------------------------------
lad_rq <- function(fml, data, adj_rlm = FALSE) {
  fml_split <- fixest:::fml_split(fml, raw = TRUE)
  if (length(fml_split) == 2) {
    stop(paste(
      "'lad_rq' selected for scale estimation but fixed effects are present",
      "in the model's fomula. Use 'scale_est = \"lad_rsc\"' (the default)", 
      "for initial scale/location estimation if fixed effects are to be",  
      "absorbed."
    ))
  }
  
  fit_rq <- quantreg::rq(fml, data = data, tau = 0.5)
  y <- as.numeric(data[, as.character(fml[[2]]), drop = TRUE])
  resid_tau <- y - stats::fitted(fit_rq, na.rm = FALSE)
  
  if (adj_rlm) {
    s <- median(abs(resid_tau))/0.6745
  } else {
    n <- length(fit_rq$y)
    df_initial <- n - length(fit_rq$coefficients)
    pprob <- (2 * n - df_initial) / (2 * n)
    s <- as.numeric(
      stats::quantile(
        abs(resid_tau), probs = pprob, type = 2, names = FALSE, na.rm = TRUE)
    ) / stats::qnorm(0.75)
  }
  
  list(
    b = stats::coef(fit_rq),
    r = resid_tau,
    s = s
  )
}

# ------------------------------------------------------------------------------
# Method of moments quantile estimator as introduced by
# Machado and Santos Silva (2019) 
# https://doi.org/10.1016/j.jeconom.2019.04.009
# Mirrors the approach of Ben Jann in his Stata moremeta code
# https://github.com/benjann/moremata/blob/master/source/mm_aqreg.mata
# ------------------------------------------------------------------------------

lad_mm_ms <- function(
  fml, data, adj_rlm = FALSE, fixef.rm = "perfect_fit", 
  warn = TRUE, notes = TRUE
) {
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
  # User receives the side effect messages of the first stage by default
  fit_loc <- fixest::feols(
    fml_loc, data, fixef.rm = fixef.rm, warn = warn, notes = notes
  )
  yhat_loc <- as.numeric(stats::fitted(fit_loc, na.rm = FALSE))
  y <- as.numeric(data[, as.character(fml[[2]]), drop = TRUE])
  e <- y - yhat_loc
  Ipos <- as.numeric(e >= 0)
  Ibar <- mean(Ipos, na.rm = TRUE)
  data$r_raw <-  2 * e * (Ipos - Ibar)

  # --- Estimate scale ---------------------------------------------------------
  fml_scl <- fml_loc
  fml_scl[[2]] <- substitute(r_raw)
  fit_scl <- fixest::feols(
    fml_scl, data, fixef.rm = fixef.rm,warn = FALSE, notes = FALSE
  )
  denom <- as.numeric(stats::fitted(fit_scl, na.rm = FALSE))
  
  u <- e / denom
  qhat <- as.numeric(stats::quantile(
    u, probs = 0.5, type = 2, names = FALSE, na.rm = TRUE
  ))
  resid_tau <- y - (yhat_loc + qhat * denom)
  
  if (adj_rlm) {
    s <- median(abs(resid_tau))/0.6745
  } else {
    df_initial <- fixest::degrees_freedom(fit_loc, type = "resid") 
    n <- fit_loc$nobs
    pprob <- (2 * n - df_initial) / (2 * n)
    s <- as.numeric(
      stats::quantile(
        abs(resid_tau), probs = pprob, type = 2, names = FALSE, na.rm = TRUE)
    ) / stats::qnorm(0.75)
  }

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
# estimates that are identical within to numerical precision to the ones
# generated by the Machado and Santos Silva (2019) approach.
# ------------------------------------------------------------------------------

lad_mm_rsc <- function(
  fml, data, adj_rlm = FALSE, fixef.rm = "perfect_fit", 
  warn = TRUE, notes = TRUE
) {
  # --- Estimate location ------------------------------------------------------
  # User receives the side effect messages of the first stage by default
  fit_loc <- fixest::feols(
    fml, data, warn = warn, notes = notes, fixef.rm = fixef.rm
  )
  yhat_loc <- as.numeric(stats::fitted(fit_loc, na.rm = FALSE))
  y <- as.numeric(data[, as.character(fml[[2]]), drop = TRUE])
  e <- y - yhat_loc
  Ipos <- as.numeric(e >= 0)
  Ibar <- mean(Ipos, na.rm = TRUE)
  data$r_raw <-  2 * e * (Ipos - Ibar)
  
  # --- Estimate scale ---------------------------------------------------------
  fml_scl <- fml
  fml_scl[[2]] <- substitute(r_raw)
  fit_scl <- fixest::feols(
    fml_scl, data, warn = FALSE, notes = FALSE, fixef.rm = fixef.rm
  )
  denom <- as.numeric(stats::fitted(fit_scl, na.rm = FALSE))
  
  u <- e / denom
  qhat <- as.numeric(stats::quantile(
    u, probs = 0.5, type = 2, names = FALSE, na.rm = TRUE
  ))
  resid_tau <- y - (yhat_loc + qhat * denom)
  
  if (adj_rlm) {
    s <- median(abs(resid_tau))/0.6745
  } else {
    df_initial <- fixest::degrees_freedom(fit_loc, type = "resid") 
    n <- fit_loc$nobs
    pprob <- (2 * n - df_initial) / (2 * n)
    s <- as.numeric(
      stats::quantile(
        abs(resid_tau), probs = pprob, type = 2, names = FALSE, na.rm = TRUE)
    ) / stats::qnorm(0.75)
  }

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
#' iteratively reweighted least squares (IRLS) for Huber M-estimation. 
#' Fixed effects are absorbed using the approach of `fixest::demean`.
#' Standard errors can be clustered via multiple levels using a Huber 
#' (psi/phi) sandwich estimator.
#'
#' @param fml A fixest formula of the form `y ~ x1 + x2 | fe1 + fe2`.
#' @param data A data.frame.
#' @param family Robust loss family. Currently only `"huber"`.
#' @param scale_est Algorithm to estimate the residuals that are used to 
#'  estimate the (inital) location and scale. Defaults to `"lad_rq"` when no 
#'  fixed effects are present. Otherwise, it defaults to `"lad_mm_rsc"`.
#'    - `"ols"`: Uses `fixest:feols()` with absorbed fixed effects to derive 
#'      the initial residuals.
#'    - `"lad_rq"`: Uses absolute residuals from a normal median/LAD regression 
#'       (calling `quantreg::rq()`) to estimate the scale. Can only be used when
#'       no fixed effects are present.   
#'    - `"lad_mm_ms"`: Uses an method of moments algorithm for the median 
#'      quantile regression introduced by Machado and Santos Silva 
#'      (2019, https://doi.org/10.1016/j.jeconom.2019.04.009). 
#'      This should yield similar scale and beta estimates as the approaches of 
#'      the Stata packages `robreg` (https://github.com/benjann/robreg) and 
#'      `robhdfe` (https://github.com/dveenman/robhdfe). 
#'    - `"lad_mm_rsc"`: Method of moments algorithm as extended by
#'      Rios-Avila, Siles, and Canavire-Bacarreza 
#'      (2024, https://docs.iza.org/dp17262.pdf). This algorithm is
#'      faster as it absorbs multi-dimensional fixed effects. Its scale and beta
#'      estimates should be virtually identical (within numerical precision) to 
#'      the ones generated by `"lad_mm_ms"`.
#' @param scale A positive numerical value to override the estimated scale.
#' @param scale_update Should the scale estimate be updated in every IRLS step?
#'  Defaults to `FALSE` like Stata's `robreg`. Set to `TRUE` to get closer to
#'  `MASS::rlm()` behavior.
#' @param efficiency Target normal-efficiency in (0.68, 1). Default is 0.95.
#' @param k Huber Loss function parameter. If `NULL` (the default), it will be
#'   determined by the `efficiency` parameter. If set, it will over-rule the 
#'   `efficiency` parameter.
#' @param tol Convergence tolerance on relative coefficient change.
#' @param max_iter Maximum IRLS iterations.
#' @param adj_rlm Whether or not the algorithm should be adjusted to be closer
#'  ro `MASS::rlm()`. Setting this to `TRUE` implies 
#'  - scale estimates to be calculated by `median(abs(resid_tau))/0.6745` 
#'    instead by a function that factors in the degrees of freedom, and
#'  - the convergence criterion to be based on residuals instead of parameter
#'    changes.
#' @param rm_singletons Whether singletons should be removed when estimating 
#'    inintial coefficients and scale as well as during IRLS steps. Defaults to
#'    `TRUE` (`fixest::feols()`  default). Set to `FALSE` to mimic the current
#'    behavior of Stata's `robtwfe`.
#' @param vcov Variance–covariance specification. Either \code{"iid"}, 
#'  \code{"hetero"}, \code{"cluster"}, \code{"twoway"} or a formula 
#'  such as \code{~ id} or  \code{cluster ~ id + time}. 
#'  Defaults to heteroskedasticity-robust (White) standard errors.
#' @param ssc Small-sample correction. Passed on to 
#'  \code{"fixest::vocov.fixest"}.
#' @param cluster (deprecated) clustering variable name(s) (string).
#'  Note that this argument is deprecated, you should use `vcov` instead.
#' @param se (deprecated) A string ("standard", or "cluster") to indicated how
#'  the standard errors should be calculated. Note that this argument is 
#'  deprecated, you should use `vcov` instead.
#' @param ... Further arguments passed to `fixest::feols`. The `weights` and
#' `fixef.rm` arguments cannot be used as 
#'  - `weight` is set by the IRLS process and regression weights are not 
#'    implemented yet by `ferols()`, and
#'  - `fixef.rm` is controlled by `rm_singletons`.
#'
#' @details
#' `ferols()` implements an IRLS process based on initial location and scale
#' estimates to derive a robust M estimate. During both stages, it uses 
#' `fixest::feols` and its demeaning approach to absorb high-dimensional 
#' fixed effects. This approach generates coefficient and standard error 
#' estimates that are very close to either `MASS::rlm()` or Stata's `robreg  m`.
#' By default, `ferols()` mimics the algorithm of `robreg`. To generate
#' estimates that are identical to `MASS::rlm()` use a call like 
#' `ferols(fml, data, k = 1.345, scale_est = "ols", adj_rlm = TRUE, scale_update = TRUE, tol = 1e-4, vcov = "iid")`.
#' 
#' For the first-step quantile regression that is used to obtain the scale 
#' estimate, `ferols()` implements the MM-QR estimation from 
#' [Machado and Santos Silva(2019)](https://www.sciencedirect.com/science/article/pii/S0304407619300648)
#' to combine quantile estimation with fixed effects, which 
#' [Rios-Avila, Siles, and Canavire-Bacarreza (2024)](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=4944894)
#' show can be extended to the multidimensional fixed effects setting. 
#'
#' `ferols()` drops singleton observations by default. This implementation choice 
#' affects not only the degrees-of-freedom adjustment in the standard error 
#' calculation (see [Correia 2015](https://scorreia.com/research/singletons.pdf)), 
#' but in this case also ensures the scale parameter is correctly estimated with 
#' the MM-QR estimation. Singleton observations produce regression residuals that 
#' are equal to (or very close to) zero, causing the scale parameter to be 
#' understated when singleton observations are retained.
#' 
#' By default, `ferols()` relays the messages and warnings from the initial 
#' location as well as the final IRLS estimation to the user. This behavior can
#' be changed by setting the `warn` and `notes` arguments of `fixest::feols()` 
#' 
#' @return An object of class `c("ferols", fixest")` with an additional list
#'   `robust` containing the following members:
#'   - `family`: The family of the loss function (currently only `"huber"`)
#'   - `scale_est`: The function to generate the initial estimates of location
#'     and scale. See discussion in the argument list above.
#'   - `efficiency`: the desired efficiency under normality, determining the 
#'     `k`tuning parameter.
#'   - `k`: The tuning parameter of the Huber loss function.
#'   - `scale`: The scale estimate for standardizing the residuals. If 
#'      `scale_update = FALSE` (the default), it is the outcome of the initial 
#'      scale estimation. Otherwise, it is the scale estimation resulting from
#'      the final IRLS step.
#'   - `converged`: Whether the IRLS process converged.
#'   - `iter`: The number of IRLS runs it took for convergence or `max_iter` if
#'     the algorithm did not converge.
#'   - `weights`: The weights resulting from the Huber loss weighting of the 
#'     final estimation step. These weights are informative as they identify 
#'     the 'bulk of the data' that remains influential for the robust regression
#'     estimates.
#'
#' @examples
#' df <- generate_panel_data(n_units = 50, n_time = 10, seed = 42)
#' ferols(y ~ x + z | i + t, vcov = ~i, data = df)
#' fixest::feols(y ~ x + z | i + t, vcov = ~i, data = df)
#' 
#' @export
ferols <- function(
    fml, data, family = "huber",  efficiency = 0.95, k = NULL,
    scale_est = NULL, scale = NULL, scale_update = FALSE, 
    tol = 1e-10, max_iter = 200, adj_rlm = FALSE, rm_singletons = TRUE,
    vcov = NULL, ssc = NULL, cluster = NULL, se = NULL, ...
) {
  if (!requireNamespace("fixest", quietly = TRUE)) {
    stop("Package 'fixest' is required to run 'ferols()'.")
  }
  
  # --- Argument checks --------------------------------------------------------
  dots <- list(...)
  if ("warn" %in% names(dots)) warn <- dots$warn else warn <- TRUE
  if ("notes" %in% names(dots)) notes <- dots$notes else notes <- TRUE
  if ("weights" %in% dots) {
    stop("Argument 'weights' is not supported by 'ferols()'.")
  }
  if ("fixef.rm" %in% dots) {
    stop(paste(
      "The exclusion of singletons is controlled by the 'rm_singletons'", 
      "argument."
    ))
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
  if (is.null(scale_est)) {
    if (length(fml_split) == 1) scale_est <- "lad_rq" 
    else scale_est <- "lad_mm_rsc"
  }
  allowed_scale_ests <- c("ols", "lad_mm_ms", "lad_mm_rsc", "lad_rq")
  if (! scale_est %in% allowed_scale_ests) {
    stop(sprintf(
      paste(
        "Scale estimation alogorithm '%s' is not supported.", 
        "Supported algorithms: '%s'"
      ), scale_est, paste(allowed_scale_ests, collapse = "', '")
    ))
  } 
  if (rm_singletons) fixef.rm = "perfect_fit" else fixef.rm = "none"
  
  # --- Initial fit to estimate scale and starting beta ------------------------
  
  
  if (scale_est == "ols")  {
    fit <- fixest::feols(fml, data = data, fixef.rm = fixef.rm, ...)
    r <- stats::residuals(fit, na.rm = FALSE)
    if (is.null(scale)) {
      if (adj_rlm) scale <- median(abs(r), na.rm = TRUE)/0.6745
      else scale <- madn(r)      
    }
    beta_old <- stats::coef(fit)
  } else if (scale_est == "lad_rq") {
    brs <- lad_rq(fml, data, adj_rlm = adj_rlm)
    r <- brs$r
    if (is.null(scale)) scale <- brs$s
    beta_old <- brs$b
  } else if (scale_est == "lad_mm_rsc") {
    brs <- lad_mm_rsc(
      fml, data, fixef.rm = fixef.rm, adj_rlm = adj_rlm, 
      warn = warn, notes = notes
    )
    r <- brs$r
    if (is.null(scale)) scale <- brs$s
    beta_old <- brs$b
  } else if (scale_est == "lad_mm_ms") {
    brs <- lad_mm_ms(
      fml, data, fixef.rm = fixef.rm, adj_rlm = adj_rlm,
      warn = warn, notes = notes
    )
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

  r_old <- r
  if (is.null(k)) k <- huber_k_from_eff(efficiency)
  
  # --- IRLS loop --------------------------------------------------------------
  converged <- FALSE
  iter_done <- 0L
  
  for (it in seq_len(max_iter)) {
    iter_done <- it
    z <- r_old / scale
    w <- huber_w(z, k)
    res_feols <- capture_side_effects(fixest::feols(
      fml, data = data, weights = w, fixef.rm = fixef.rm, ...
    ))
    fit <- res_feols$value 
    beta_new <- stats::coef(fit)
    if (anyNA(beta_new)) stop(
      "IRLS produced NA coefficients; check collinearity/data issues."
    )
    r_new <- stats::residuals(fit, na.rm = FALSE)
    
    # Convergence check
    if (adj_rlm) diff <- sqrt(
      sum((r_old - r_new)^2, na.rm = TRUE) / 
        max(1e-20, sum(r_old^2, na.rm = TRUE))
    )
    else diff <- max(abs(beta_new - beta_old) / pmax(1, abs(beta_old)))
    if (diff <= tol) {
      converged <- TRUE
      break
    }
    
    if (scale_update) {
      if (adj_rlm) scale <- median(abs(r_new), na.rm = TRUE)/0.6745
      else scale <- madn(r_new)
    }

    beta_old <- beta_new
    r_old <- r_new
  }
  
  relay_side_effects(res_feols)
  if (!converged) warning(
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
    weights = w
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
  attr(se0, "vcov_type") <- attr(V, "vcov_type")
  
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
  
  attr(ct, "vcov_type") <- attr(se0, "vcov_type")
  
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
#'  \code{"hetero"}, \code{"cluster"}, \code{"twoway"} or a formula 
#'  such as \code{~ id} or  \code{cluster ~ id + time}. 
#'  Defaults to heteroskedasticity-robust (White) standard errors.
#' @param cluster Deprecated. Use \code{vcov} instead.
#' @param se Deprecated. Use \code{vcov} instead.
#' @param ssc Small-sample correction. Passed on to 
#'  \code{"fixest::vocov.fixest"}.
#' @param ... Additional arguments (passed on to \code{"fixest::vocov.fixest"}).
#'
#' @return A variance–covariance matrix with coefficient names as row and
#'   column names. The matrix carries a \code{"type"} attribute describing
#'   the estimator.
#'
#' @details
#' This function is experimental and implements only a subset of the
#' variance estimators available in \code{fixest}. Several panel robust 
#' estimators and user-supplied functions and covariance matrices are not 
#' supported.
#'
#' @export
vcov.ferols <- function(
    object, vcov = NULL, cluster = NULL, se = NULL, ssc = NULL, ...
) {
  vcov <- fixest:::oldargs_to_vcov(se, cluster, vcov) 
  if(is.null(vcov)) {vcov <- "hetero"}
  supported_char <- c("iid", "twoway", "cluster", "hetero")
  if (is.character(vcov)) {
    vcov <- tolower(vcov)
    if (!(vcov %in% supported_char)) {
      stop(sprintf(
        "VCOV type '%s' is not supported yet in ferols(). Supported: %s.",
        vcov, paste(supported_char, collapse = ", ")
      ))
    }
    if (tolower(vcov) %in% c("cluster", "twoway")) {
      if (length(object$fixef_vars) < 1) stop(paste(
        "vcov indicates clustering but no fixed effects presents.", 
        "Consider setting the desired cluster explicitly by providing",
        "'vcov = cluster ~ var'"
      )) else {
        if (vcov == "cluster") {
          vcov <- as.formula(sprintf("cluster ~ %s", object$fixef_vars[1]))
        } else {
          vcov <- as.formula(
            sprintf("cluster ~ %s", paste(object$fixef_vars, collapse = " + "))
          )
        }
      }
    }
  } 
  if (inherits(vcov, "formula")) {
    # If there is no lhs, fixest defaults to "cluster"
    if (length(vcov) == 2) {
      vcov <- stats::as.formula(paste0("cluster ~ ", deparse(vcov[[2]])))
    }
    if (all.vars(vcov)[1] != "cluster") stop(
      "Only 'cluster' is supported as the lhs of a vcov formula"
    ) 
    cl_vars <- fixest:::fml2varnames(vcov, combine_fun = TRUE)
  } 
  if (! is.character(vcov) & ! inherits(vcov, "formula")) {
    please_use_str <- sprintf( 
      "Please use instead a formula or vcov='keyword' with keyword in '%s'.", 
      paste(supported_char, collapse = "', '") 
    ) 
    if (inherits(vcov, "fixest_vcov_request")) { 
      stop(paste( "ferols() does not support vcov objects.", please_use_str )) 
    } 
    if (is.function(vcov) || is.matrix(vcov) || is.list(vcov)) { 
      stop(
        paste( 
          "ferols() does not support custom vcov functions/matrices.", 
          please_use_str 
        )) 
    }
    stop(paste("Unknown vcov type.",  please_use_str))
  }
  
  mf <- fixest::fixest_data(object, sample = "estimation")
  
  fe_df <- if (length(object$fixef_vars)) {
    mf[, object$fixef_vars, drop = FALSE] 
  } else NULL
  
  X  <- stats::model.matrix(object, type = "rhs")
  r  <- stats::residuals(object)
  k  <- object$robust$k
  sc <- object$robust$scale
  z   <- r / sc
  psi <- huber_psi(z, k) 
  phi <- huber_phi(z, k)
  
  if (vcov == "iid") {
    # We are on our own. Calculating MASS::rlm() model-type standard errors
    Xr <- if (!is.null(fe_df)) {
      fixest::demean(X, f = fe_df, na.rm = TRUE) 
    } else X
    
    w   <- huber_w(z, k)
    n <- nrow(Xr)
    p <- fixest::degrees_freedom(object, type = "k")
    rdf <- n - p
      
    # summary.rlm does: S = sum((wresid*w)^2)/rdf, w = psi(wresid/s)
    # Here take wresid ≈ r, and w = huber_w(r/s, k)
    S <- sum((r * w)^2) / rdf
    mn <- mean(phi)
    kappa <- 1 + p * stats::var(phi) / (n * mn^2)
    stddev <- sqrt(S) * (kappa / mn)
    
    cov_unscaled <- solve(crossprod(Xr))
    V <- (stddev^2) * cov_unscaled
    colnames(V) <- colnames(Xr)
    rownames(V) <- colnames(Xr)
    attr(V, "vcov_type") <- paste0("IID model-based")
  } else {
    Xr <- if (!is.null(fe_df)) {
      fixest::demean(X, f = fe_df, weights = phi, na.rm = TRUE) 
    } else X
    S <- Xr * psi
    
    # Inject bread and scale and then relate to vcov.fixest()
    # bread: scale * (Xr' Phi Xr)^(-1)
    XtPhiX_inv <- solve(crossprod(Xr * as.numeric(phi), Xr))
    bread_rob  <- sc * XtPhiX_inv
    
    tmp <- object
    tmp$scores  <- S
    tmp$sigma2  <- 1
    tmp$cov.iid <- bread_rob
    class(tmp) <- setdiff(class(tmp), "ferols")
    V <- vcov(tmp, vcov = vcov, ssc = ssc, ...)
  
    # Set attributes so that downstream functions report appropriate SE type
    if (is.character(vcov)) {
      attr(V, "vcov_type") <- "Heteroskedasticity-robust"
    } else {
      attr(V, "vcov_type") <- paste0(
        "Clustered (", paste(cl_vars, collapse = " & "), ")"
      )
    }
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

  if (! object$robust$converged) {
    warning(
      "IRLS convergence not achieved (increase max_iter or relax tol)."
    )
  }
  
  # delegate to fixest summary
  obj <- object
  class(obj) <- setdiff(class(obj), "ferols")
  s <- summary(obj, vcov = V, ssc = ssc, ...)
  
  # Restore our vcoc
  s$cov.scaled <- V
  
  # Set SE attributes
  attr(s$se, "vcov_type") <- attr(V, "vcov_type")
  attr(s$coeftable, "vcov_type") <- attr(V, "vcov_type")
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
  if (is.null(x$summary) || ! x$summary) s <- summary.ferols(x)
  else s <- x
  class(s) <- setdiff(class(s), "ferols")
  print(s, ...)
  invisible(x)
}