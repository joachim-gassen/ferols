# Load test data shipped in inst/extdata
cap <- function(expr) paste(capture.output(expr), collapse = "\n")

test_df <- generate_panel_data(seed = 42)
assign("test_df", test_df, envir = .GlobalEnv)
on.exit(rm("test_df", envir = .GlobalEnv), add = TRUE)

test_df_missing <- generate_panel_data(seed = 42, share_missing = 0.3)
assign("test_df_missing", test_df_missing, envir = .GlobalEnv)
on.exit(rm("test_df_missing", envir = .GlobalEnv), add = TRUE)


test_that("ferols runs and returns a ferols/fixest object", {
  fit <- ferols(y ~ x, data = test_df, vcov = "iid")
  
  expect_s3_class(fit, "ferols")
  expect_true(inherits(fit, "fixest"))
  expect_true(is.matrix(fit$cov.scaled))
  expect_equal(ncol(fit$cov.scaled), length(coef(fit)))
  
  s <- summary(fit)
  expect_true(is.numeric(s$se))
  expect_equal(attr(s$se, "type"), "IID")
})


test_that("default vcov is IID and default scale_est is 'lad_mm_rsc'", {
  expect_no_warning(fit <- ferols(y ~ x | i + t, data = test_df))
  expect_true(is.matrix(fit$cov.scaled))
  type <- attr(fit$cov.scaled, "type")
  expect_true(is.character(type))
  expect_true(grepl("IID", type))
  expect_true(grepl("Standard-errors: IID", cap(fit), fixed = TRUE))
  expect_true(grepl("Standard-errors: IID", cap(summary(fit)), fixed = TRUE))
  out <- cap(fit)
  expect_true(grepl("scale est: lad_mm_rsc", out, fixed = TRUE))
})


test_that("vcov specifications ~i and cluster = 'i' are accepted", {
  expect_no_warning(
    f1 <- ferols(y ~ x | i + t, data = test_df, vcov = ~ i)
  )
  expect_no_warning(
    f2 <- ferols(y ~ x | i + t, data = test_df, cluster = "i")
  )
  
  expect_no_warning(se1 <- summary(f1)$se)
  expect_no_warning(se2 <- summary(f2)$se)
  expect_equal(unname(se1), unname(se2), tolerance = 1e-10)
  expect_true(grepl("Clustered", attr(se1, "type")))
  expect_true(grepl("Standard-errors: Clustered (i)", cap(f1), fixed = TRUE))
  expect_true(grepl("Standard-errors: Clustered (i)", cap(summary(f1)), fixed = TRUE))
})

test_that("unsupported features error cleanly", {
  expect_error(ferols(y ~ 1 | i + t | p ~ x, data = test_df))
  expect_error(ferols(y ~ 1 | i + t, data = test_df, weights = runif(nrow(test_df))))
  expect_error(ferols(y ~ x | i + t, data = test_df, scale_est = "unknown"))
  expect_error(ferols(y ~ x | i + t, data = test_df, vcov = ~ i + t))
  expect_error(ferols(y ~ x | i + t, data = test_df, vcov = "hetero"))
  expect_error(ferols(y ~ x | i + t, data = test_df, vcov = matrix(1, 1, 1)))
  expect_error(ferols(y ~ x | i + t, data = test_df, ssc = fixest::ssc(fixef.K = "none")))
  expect_error(ferols(y ~ x | i + t, data = test_df, efficiency = 0.6))
  expect_error(ferols(y ~ x | i + t, data = test_df, cluster = c("i", "t")))
})

test_that("efficiency 0.7 runs", {
  fit <- ferols(y ~ x | i + t, data = test_df, vcov = ~ i, efficiency = 0.7)
  expect_s3_class(fit, "ferols")
})

test_that("scale_est='ols' yields an estimate", {
  expect_no_warning(
    fit <- ferols(y ~ x, data = test_df, scale_est = "ols")
  )
  expect_s3_class(fit, "ferols")
  expect_no_warning(
    fit <- ferols(y ~ x | i + t, data = test_df, scale_est = "ols")
  )
  expect_s3_class(fit, "ferols")
})

test_that("scale_est='lad_mm_ms' yields an estimate", {
  expect_no_warning(
    fit <- ferols(y ~ x, data = test_df, scale_est = "lad_mm_ms")
  )
  expect_s3_class(fit, "ferols")
  expect_no_warning(
    fit <- ferols(y ~ x | i + t, data = test_df, scale_est = "lad_mm_ms")
  )
  expect_s3_class(fit, "ferols")
})

test_that("scale_est='lad_mm_rsc' yields an estimate", {
  expect_no_warning(
    fit <- ferols(y ~ x, data = test_df, scale_est = "lad_mm_rsc")
  )
  expect_s3_class(fit, "ferols")
  expect_no_warning(
    fit <- ferols(y ~ x | i + t, data = test_df, scale_est = "lad_mm_rsc")
  )
  expect_s3_class(fit, "ferols")
})

test_that(
  "'lad_mm_rsc' and 'lad_mm_ms' produce virtually identical estimates", {
    fit_rsc <- ferols(y ~ x, data = test_df, scale_est = "lad_mm_rsc")
    fit_ms <- ferols(y ~ x, data = test_df, scale_est = "lad_mm_ms")
    expect_equal(fit_rsc$robust$scale, fit_ms$robust$scale)
    expect_equal(fit_rsc$coefficients, fit_ms$coefficients)
    expect_equal(fit_rsc$se, fit_ms$se)
  }
)

test_that("efficiency=1 gives coefficients close to feols", {
  fit_fe <- fixest::feols(y ~ x | i + t, data = test_df, vcov = "iid")
  fit_fr <- ferols(y ~ x | i + t, data = test_df, vcov = "iid", efficiency = 1)
  
  expect_equal(unname(coef(fit_fr)), unname(coef(fit_fe)), tolerance = 1e-8)
})

test_that("print.ferols always prints the model info line (eff/k/scale/iter)", {
  out1 <- cap(ferols(y ~ x | i + t, vcov = ~ i, data = test_df))
  fit <- ferols(y ~ x | i + t, vcov = ~ i, data = test_df)
  out2 <- cap(fit)
  out3 <- cap(summary(fit))
  
  for (o in c(out1, out2, out3)) {
    expect_true(grepl(
      "ferols() - Fixed-effects robust IRLS M regression (Huber loss)", 
      o, fixed=TRUE
    ))
    expect_true(grepl("efficiency:", o, fixed = TRUE))
    expect_true(grepl("\\bk:\\s*[-+]?[0-9]*\\.?[0-9]+", o, perl = TRUE))
    expect_true(grepl("scale est: (ols|lad)", o))
    expect_true(grepl(
      "\\bscale:\\s*[-+]?[0-9]*\\.?[0-9]+([eE][-+]?[0-9]+)?", o, perl = TRUE
    ))
    expect_true(grepl("\\biter:\\s*[0-9]+", o, perl = TRUE))
  }
})

test_that("data.save = TRUE allows vcov after data removal", {
  assign("temp_test_df", generate_panel_data(seed = 42), envir = .GlobalEnv)

  fit <- ferols(y  ~ x | i + t, data = temp_test_df, data.save = TRUE)
  
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

test_that("scale parameter is obeyed.", {
  out <- cap(
    ferols(y ~ x | i + t, vcov = ~ i, scale = 0.987, data = test_df)
  )
  expect_true(grepl("scale: 0.987", out, fixed = TRUE))
})

test_that("estimation passes on data with missing values", {
  fit <- ferols(y ~ x + z | i + t, vcov = ~i, data = test_df_missing)
  expect_s3_class(fit, "ferols")
  fit <- ferols(
    y ~ x + z | i + t, vcov = ~i, scale_est = "lad_mm_ms", 
    data = test_df_missing
  )
  expect_s3_class(fit, "ferols")
  fit <- ferols(
    y ~ x + z | i + t, vcov = ~i, scale_est = "ols", data = test_df_missing
  )
  expect_s3_class(fit, "ferols")
})

test_that(
  paste(
    "'lad_mm_rsc' with clustering returns scale, parameter estimates and SEs", 
    "that are virtually identical to those of Stata's robreg m" 
  ), {
    robreg_scale <- 0.668803236125195
    robreg_x <- 1.0029391
    robreg_se <- sqrt(0.00004447)
    fit <- ferols(y ~ x + z | i + t, data = test_df, vcov = ~i)
    expect_equal(
      fit$robust$scale, robreg_scale, tolerance = 1e-4, ignore_attr = TRUE
    )
    expect_equal(
      fit$coefficients[1], robreg_x, tolerance = 1e-4, ignore_attr = TRUE
    )
    expect_equal(
      fit$se[1], robreg_se, tolerance = 1e-4, ignore_attr = TRUE
    )

    robreg_scale <- 0.6909039341817896
    robreg_x <- 1.0033825
    robreg_se <- sqrt(0.00007045)
    fitm <- ferols(y ~ x + z | i + t, data = test_df_missing, vcov = ~i)
    expect_equal(
      fitm$robust$scale, robreg_scale, tolerance = 1e-4, ignore_attr = TRUE
    )
    expect_equal(
      fitm$coefficients[1], robreg_x, tolerance = 1e-4, ignore_attr = TRUE
    )
    expect_equal(
      fitm$se[1], robreg_se, tolerance = 1e-4, ignore_attr = TRUE
    )
  }
)

