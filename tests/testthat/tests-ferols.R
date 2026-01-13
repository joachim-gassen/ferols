# Load test data shipped in inst/extdata
get_test_df <- function() {
  path <- system.file("extdata", "ferols_test_data.rds", package = "ferols")
  if (path == "") stop("ferols_test_data.rds not found in package extdata")
  df <- readRDS(path)
}

test_df <- get_test_df()
assign("test_df", test_df, envir = .GlobalEnv)
on.exit(rm("test_df", envir = .GlobalEnv), add = TRUE)

set.seed(1)
n_firm <- 50
n_period <- 20
firm <- rep(seq_len(n_firm), each = n_period)
period <- rep(seq_len(n_period), times = n_firm)

x <- rnorm(n_firm * n_period)
x2 <- rnorm(n_firm * n_period)
fe_firm <- rnorm(n_firm)[firm]
fe_period <- rnorm(n_period)[period]
y <- 1 * x + 0 * x2 + fe_firm + fe_period + rnorm(n_firm * n_period, sd = 0.5)

sim_df <- data.frame(
  y = y, x = x, x2 = x2, firm = firm, period = period, group = firm
)
assign("sim_df", sim_df, envir = .GlobalEnv)
on.exit(rm("sim_df", envir = .GlobalEnv), add = TRUE)


test_that("ferols runs and returns a ferols/fixest object", {
  fit <- ferols(y_cont ~ x, data = test_df, vcov = "iid")
  
  expect_s3_class(fit, "ferols")
  expect_true(inherits(fit, "fixest"))
  expect_true(is.matrix(fit$cov.scaled))
  expect_equal(ncol(fit$cov.scaled), length(coef(fit)))
  
  s <- summary(fit)
  expect_true(is.numeric(s$se))
  expect_equal(attr(s$se, "type"), "IID")
})


test_that("default vcov clusters on first fixed effect", {
  fit <- ferols(y_cont ~ x | i + t, data = test_df)
  
  expect_true(is.matrix(fit$cov.scaled))
  type <- attr(fit$cov.scaled, "type")
  expect_true(is.character(type))
  expect_true(grepl("Clustered", type))
  expect_true(grepl("Clustered (i)", type, fixed = TRUE))
})


test_that("vcov specifications ~i and cluster = 'i' are accepted", {
  f1 <- ferols(y_cont ~ x | i + t, data = test_df, vcov = ~ i)
  f2 <- ferols(y_cont ~ x | i + t, data = test_df, cluster = "i")
  
  se1 <- summary(f1)$se
  se2 <- summary(f2)$se
  expect_equal(unname(se1), unname(se2), tolerance = 1e-10)
  expect_true(grepl("Clustered", attr(se1, "type")))
})


test_that("efficiency 0.7 runs", {
  fit <- ferols(y_cont ~ x | i + t, data = test_df, vcov = ~ i, efficiency = 0.7)
  expect_s3_class(fit, "ferols")
})


test_that("unsupported features error cleanly", {
  expect_error(ferols(y_cont ~ 1 | i + t | p ~ x, data = test_df))
  expect_error(ferols(y_cont ~ x | i + t, data = test_df, vcov = ~ i + t))
  expect_error(ferols(y_cont ~ x | i + t, data = test_df, vcov = "hetero"))
  expect_error(ferols(y_cont ~ x | i + t, data = test_df, vcov = matrix(1, 1, 1)))
  expect_error(ferols(y_cont ~ x | i + t, data = test_df, ssc = fixest::ssc(fixef.K = "none")))
  expect_error(ferols(y_cont ~ x | i + t, data = test_df, efficiency = 0.6))
  expect_error(ferols(y_cont ~ x | i + t, data = test_df, cluster = c("i", "t")))
})


test_that("efficiency=1 gives coefficients close to feols", {
  fit_fe <- fixest::feols(y_cont ~ x | i + t, data = test_df, vcov = "iid")
  fit_fr <- ferols(y_cont ~ x | i + t, data = test_df, vcov = "iid", efficiency = 1)
  
  expect_equal(unname(coef(fit_fr)), unname(coef(fit_fe)), tolerance = 1e-8)
})


test_that("simple twfe simulation returns sensible estimate", {
  fit <- ferols(y ~ x + x2 | firm + period, vcov = ~ group, data = sim_df)
  
  expect_true(abs(coef(fit)["x"] - 1) < 0.1)
})

test_that("data.save = TRUE allows vcov after data removal", {
  assign("temp_test_df", get_test_df(), envir = .GlobalEnv)

  fit <- ferols(y_cont  ~ x | i + t, data = temp_test_df, data.save = TRUE)
  
  # sanity check: data was saved by fixest
  expect_true(!is.null(fit$data))
  
  # check whether data is identical to test_id
  expect_equal(fit$data, test_df, ignore_attr = TRUE)
  
  # remove original data
  rm("temp_test_df", envir = .GlobalEnv)

  # vcov should still work because data.save = TRUE
  V <- vcov(fit)
  expect_true(is.matrix(V))
  expect_true(all(dim(V) > 0))
})

