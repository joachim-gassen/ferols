# ------------------------------------------------------------------------------
# Code to generate a data frame containing simulated panel data for testing
# ferols() or for comparing properties of estimator
#
# Creates observations with a unit and time dimension and errors
# that are auto-correlated within units. A subset of the observations
# receives a multiplied error so that the errors are mixed normal.
# 
# See LICENSE file for license
# ------------------------------------------------------------------------------

#' Generate Panel Data
#' 
#' @description 
#' Generates a data frame containing simulated panel data for testing the
#' \code{ferols()} function. The data contains observations with unit and 
#' time dimensions and two-way clustered composite errors drawn from a mixed 
#' distribution.
#' @param n_units Integer. Number of cross-sectional units. Default is 1000.
#' @param n_time Integer. Number of time periods. Default is 10.
#' @param e_mult Numeric. Multiplier to set the kurtosis fo the composite error. 
#' `e_share` of the composite error will be multiplied by `e_mult`. 
#' Default is 3, resulting in mild kurtosis of the composite error.
#' @param e_share Numeric. Share of observations affected by the error 
#' multiplier. Default is 0.1.
#' @param share_missing Numeric. Percentage of missing data to introduce into 
#' the panel. Default is 0 (no missing data, balanced panel). Needs to be 
#' between 0 and smaller than 1.
#' @param mean_u Numeric. Mean of the unit-level random effects. Default is 0.
#' @param sigma_u Numeric. Standard deviation of the unit-level random effects.
#' Default is 1/3.
#' @param rho_u Numeric. Unit-specific autocorrelation coefficient for the 
#' error term. Default is 0.5.
#' @param mean_t Numeric. Mean of the time-level random effects. Default is 0.
#' @param sigma_t Numeric. Standard deviation of the time-level random effects.
#' Default is 1/3.
#' @param rho_t Numeric. Time-specific autocorrelation coefficient for the 
#' error term. Default is 0.5.
#' @param sigma_e Numeric. Standard deviation of the idiosyncratic noise 
#' component. Default is 1 - sigma_u - sigma_t.
#' @param beta Numeric. Coefficient for the independent variable x. Default is 1.
#' @param mu_x Numeric. Mean of the independent variable x. Default is 0.
#' @param sigma_x Numeric. Standard deviation of the independent variable x.
#' Default is 1.
#' @param n_z Integer. Number of control variables z to include. Default is 1.
#' @param beta_z Numeric vector. Coefficients for the control variable(s) z. 
#' Default is 0 (no effect). Needs to be of length `n_z`.
#' @param mu_z Numeric vector. Mean of the control variable(s) z. Default is 0. 
#' Needs to be of length `n_z`.
#' @param sigma_z Numeric vector. Standard deviation of the control 
#' variable(s) z. Default is 1. Needs to be of length `n_z`.
#' @param seed Integer. Random seed for reproducibility.
#' 
#' @return A data frame containing the simulated panel data.
#' 
#' @examples
#' df <- generate_panel_data(n_units = 50, n_time = 10, seed = 42)
#' head(df)
#' ferols(y ~ x + z | i + t, vcov = ~i, data = df)
#' fixest::feols(y ~ x + z | i + t, vcov = ~i, data = df)
#' 
#' @export
generate_panel_data <- function(
  n_units = 1000, n_time = 10, e_mult = 3, e_share = 0.1, share_missing = 0,
  mean_u = 0, sigma_u = 1/3, rho_u = 0.5,
  mean_t = 0, sigma_t = 1/3, rho_t = 0.5,
  sigma_e = 1 - sigma_u - sigma_t,
  beta = 1, mu_x = 0, sigma_x = 1,
  n_z = 1, beta_z = 0, mu_z = 0, sigma_z = 1,
  seed = NULL
) {

  if (length(mu_z) != n_z) {
    stop("Length of mu_z must be equal to n_z. Each control variable needs a mean.")
  }
  if (length(sigma_z) != n_z) {
    stop("Length of sigma_z must be equal to n_z. Each control variable needs a standard deviation.")
  }
  if (length(beta_z) != n_z) {
    stop("Length of beta_z must be equal to n_z. Each control variable needs a coefficient.")
  }
  if (share_missing < 0 | share_missing >= 1) {
    stop("share_missing must be between 0 and less than 1.")
  }

  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  unit_var_names <- paste0(
    "i_", sprintf(
      paste0("%0", ceiling(log10(n_units + 0.1)), "d"), 
      seq_len(n_units)
    )
  )
  time_var_names <- paste0(
    "t_", sprintf(
      paste0("%0", ceiling(log10(n_time + 0.1)), "d"), 
      seq_len(n_time)
    )
  )
  panel <- expand.grid(t = time_var_names, i = unit_var_names)
  
  panel$x <- rnorm(nrow(panel), mean = mu_x, sd = sigma_x)
  
  z_var_names <- if(n_z == 1) "z" else paste0(
    "z", sprintf(paste0("%0", ceiling(log10(n_z + 0.1)), "d"), seq_len(n_z))
  )
  for (i in z_var_names) {
    panel[i] <- rnorm(nrow(panel), mean = mu_z, sd = sigma_z)
  }
  
  u_fe <- rnorm(n_units, mean = mean_u, sd = sigma_u)
  t_fe <- rnorm(n_time, mean = mean_t, sd = sigma_t)
  
  panel$e <- NA
  for (i in seq_len(n_units)) {
    e_unit <- numeric(n_time)
    e_unit[1] <- rnorm(1, mean = 0, sd = sigma_e)
    for (t in seq(2, n_time)) {
      e_unit[t] <- rho_u * e_unit[t - 1] + 
        rnorm(1, mean = 0, sd = sqrt(1 - rho_u^2) * sigma_e)
    }
    panel$e[panel$i == unit_var_names[i]] <- e_unit
  }
  
  # --- Petersen-style but two-way clustered error: e_it = u_i + r_t + eta_it --
  
  # firm component (AR(1) across t within i)
  for (i in seq_len(n_units)) {
    u <- numeric(n_time)
    u[1] <- rnorm(1, mean = 0, sd = sigma_u)
    for (t in seq(2, n_time)) {
      u[t] <- rho_u * u[t - 1] + 
        rnorm(1, mean = 0, sd = sqrt(1 - rho_u^2) * sigma_u)
    }
    panel$u[panel$i == unit_var_names[i]] <- u
  }

  # time component (AR(1) across i within t)
  for (p in seq_len(n_time)) {
    r <- numeric(n_units)
    r[1] <- rnorm(1, mean = 0, sd = sigma_t)
    for (t in seq(2, n_units)) {
      r[t] <- rho_t * r[t - 1] + 
        rnorm(1, mean = 0, sd = sqrt(1 - rho_t^2) * sigma_t)
    }
    panel$r[panel$t == time_var_names[p]] <- r
  }
  
  # idiosyncratic component 
  panel$eta <- rnorm(n_time*n_units, 0, sigma_e)
  
  # composite error 
  panel$e <- panel$u + panel$r + panel$eta
  
  
  # Add outliers
  n_mixed <- floor(e_share*nrow(panel))
  mixed_indices <- sample(seq_len(nrow(panel)), n_mixed)
  panel$outlier <- FALSE
  panel$outlier[mixed_indices] <- TRUE
  panel$e[mixed_indices] <- panel$e[mixed_indices] * e_mult 

  panel$y <- u_fe[match(panel$i, unit_var_names)] + 
    t_fe[match(panel$t, time_var_names)] + 
    beta * panel$x + 
    as.vector(as.matrix(panel[, z_var_names, drop = FALSE]) %*% beta_z) +
    panel$e 

  if (share_missing > 0) {
    n_missing <- floor(share_missing * nrow(panel))
    missing_indices <- sample(seq_len(nrow(panel)), n_missing)
    panel$y[missing_indices] <- NA
  }

  panel <- panel[, c("i", "t", "y", "x", z_var_names, "e", "outlier")]
  return(panel)
}