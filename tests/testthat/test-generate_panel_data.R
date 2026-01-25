# tests/testthat/test-generate_panel_data.R

test_that("input validation checks fire", {
  expect_error(
    generate_panel_data(n_units = 5, n_time = 2, n_z = 2, mu_z = 0, sigma_z = c(1, 1), seed = 1),
    regexp = "Length of mu_z must be equal to n_z",
    info = "mu_z length must match n_z."
  )

  expect_error(
    generate_panel_data(n_units = 5, n_time = 2, n_z = 2, mu_z = c(0, 0), sigma_z = 1, seed = 1),
    regexp = "Length of sigma_z must be equal to n_z",
    info = "sigma_z length must match n_z."
  )

  expect_error(
    generate_panel_data(n_units = 5, n_time = 2, share_missing = -0.01, seed = 1),
    regexp = "share_missing must be between 0 and less than 1",
    info = "shate_missing must be in [0, 1)."
  )

  expect_error(
    generate_panel_data(n_units = 5, n_time = 2, share_missing = 1, seed = 1),
    regexp = "share_missing must be between 0 and less than 1",
    info = "share_missing must be in [0, 1)."
  )
})

test_that("returns a data.frame with expected shape and columns", {
  n_units <- 12
  n_time  <- 5
  n_z     <- 2

  d <- generate_panel_data(
    n_units = n_units,
    n_time  = n_time,
    n_z     = n_z,
    mu_z    = rep(0, n_z),
    sigma_z = rep(1, n_z),
    beta_z  = rep(0, n_z),
    seed    = 123
  )

  expect_s3_class(d, "data.frame")
  expect_equal(nrow(d), n_units * n_time)

  z_names <- paste0("z", sprintf(paste0("%0", ceiling(log10(n_z + 0.1)), "d"), seq_len(n_z)))
  expect_named(d, c("i", "t", "y", "x", z_names, "e", "outlier"))

  expect_true(is.factor(d$i))
  expect_true(is.factor(d$t))
  expect_type(d$y, "double")
  expect_type(d$x, "double")
  expect_type(d$e, "double")
  expect_type(d$outlier, "logical")
})

test_that("unit/time IDs are correctly zero-padded", {
  n_units <- 12
  n_time  <- 5
  d <- generate_panel_data(
    n_units = n_units, n_time = n_time,
    n_z = 1, mu_z = 0, sigma_z = 1, beta_z = 0,
    seed = 1
  )

  unit_width <- ceiling(log10(n_units + 0.1))
  time_width <- ceiling(log10(n_time + 0.1))

  expected_units <- paste0("i_", sprintf(paste0("%0", unit_width, "d"), seq_len(n_units)))
  expected_times <- paste0("t_", sprintf(paste0("%0", time_width, "d"), seq_len(n_time)))

  expect_setequal(unique(d$i), expected_units)
  expect_setequal(unique(d$t), expected_times)
})

test_that("seed makes output reproducible", {
  args <- list(
    n_units = 10, n_time = 4,
    n_z = 2, mu_z = c(0, 0), sigma_z = c(1, 1), beta_z = c(0.1, -0.2),
    seed = 999
  )

  d1 <- do.call(generate_panel_data, args)
  d2 <- do.call(generate_panel_data, args)

  expect_identical(d1, d2)
})

test_that("outlier flag count matches e_share", {
  n_units <- 20
  n_time  <- 6
  e_share <- 0.15

  d <- generate_panel_data(
    n_units = n_units, n_time = n_time,
    e_share = e_share,
    n_z = 1, mu_z = 0, sigma_z = 1, beta_z = 0,
    seed = 42
  )

  expected_n <- floor(e_share * nrow(d))
  expect_equal(sum(d$outlier), expected_n)
})

test_that("pct_missing introduces missing y at expected rate", {
  n_units <- 30
  n_time  <- 4
  share_missing <- 0.2

  d <- generate_panel_data(
    n_units = n_units, n_time = n_time,
    share_missing = share_missing,
    n_z = 1, mu_z = 0, sigma_z = 1, beta_z = 0,
    seed = 7
  )

  expected_n <- floor(share_missing * nrow(d))
  expect_equal(sum(is.na(d$y)), expected_n)
})

test_that("y equals FE + beta*x + Z*beta_z + e when sigma_u/sigma_t are zero", {
  # Turn off unit/time FE to make y check simpler and deterministic conditional on RNG.
  n_units <- 8
  n_time  <- 3
  beta    <- 1.7
  beta_z  <- c(0.5, -0.25)
  n_z     <- length(beta_z)

  d <- generate_panel_data(
    n_units = n_units, n_time = n_time,
    mean_u = 0, sigma_u = 0,
    mean_t = 0, sigma_t = 0,
    beta = beta,
    n_z = n_z, mu_z = rep(0, n_z), sigma_z = rep(1, n_z), beta_z = beta_z,
    seed = 101
  )

  z_names <- paste0("z", sprintf(paste0("%0", ceiling(log10(n_z + 0.1)), "d"), seq_len(n_z)))
  z_mat <- as.matrix(d[, z_names, drop = FALSE])

  y_expected <- beta * d$x + as.numeric(z_mat %*% beta_z) + d$e
  expect_equal(d$y, y_expected, tolerance = 1e-12)
})
