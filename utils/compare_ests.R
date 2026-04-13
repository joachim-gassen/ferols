library(ferols)

RUNS <- 500
VARS <- c("run_time", "est", "se", "scale")
STATA_PATH <- "/Applications/Stata/StataSE.app/Contents/MacOS/stata-se"
FORCE_RERUN <- FALSE
set.seed(42)

do_stata_file <- function(do_file, spath = STATA_PATH, log = FALSE) {
  if(.Platform$OS.type != "unix") stop(
    "This only works on unix type OSes (inlcuding MacOS)."
  )
  if (log) system(paste(spath, "-b", do_file))
  else system(sprintf("%s < %s > /dev/null 2>&1", spath, do_file))
}

run_rlm <- function(fml, df, ...) {
  time_in <- Sys.time()
  fit <- MASS::rlm(fml, df, ...)
  time_out <- Sys.time()
  list(
    run_time = as.vector(difftime(time_out, time_in, units = "secs")),
    est = coef(fit)[["x"]],
    se = summary(fit)$coef["x", 2],
    scale = fit$s
  )
}

run_ferols <- function(fml, df, ...) {
  time_in <- Sys.time()
  fit <- ferols(fml, df, ...)
  time_out <- Sys.time()
  list(
    run_time = as.vector(difftime(time_out, time_in, units = "secs")),
    est = coef(fit)[["x"]],
    se = fit$se[["x"]],
    scale = fit$robust$scale
  )
}

run_stata <- function(df, do_file) {
  haven::write_dta(df, "temp/df.dta")
  do_stata_file(do_file)
  as.list(read.delim("temp/stata_out.tsv"))
}

run_ests <- function(
    rlm_fml = NULL, ferols_fml = rlm_fml, 
    robreg_do_file = NULL, robhdfe_do_file = NULL,
    rlm_parms = list(), ferols_parms = list(), data_parms = list()
  ) {
  if (
    is.null(rlm_fml) && is.null(ferols_fml) && 
    is.null(robhdfe_do_file) && is.null(robreg_do_file)
  ) {
    stop("At least one estimation method needs to be specified.")
  }
  if (!"seed" %in% names(data_parms)) {
    data_parms <- c(data_parms, list(seed = sample(seq_len(1e9), 1)))
  }
  df <- do.call(generate_panel_data, data_parms)
  df$g <- rep(sprintf("g_%02d", 1:50), each = nrow(df)/50)
  force(df)
  rl <- list(data_parms = data_parms)
  if (!is.null(rlm_fml)) {
    rl_rlm <- do.call(run_rlm, c(rlm_fml, list(df), rlm_parms))
    rl <- c(rl, rlm = list(c(rlm_parms, rl_rlm)))
  } 
  if (!is.null(robreg_do_file)) {
    rl_robreg <- run_stata(df, robreg_do_file)
    rl <- c(rl, robreg = list(c(list(do_file = robreg_do_file), rl_robreg)))
  } 
  if (!is.null(robhdfe_do_file)) {
    rl_robhdfe <- run_stata(df, robhdfe_do_file)
    rl <- c(rl, robhdfe = list(c(list(do_file = robhdfe_do_file), rl_robhdfe)))
  } 
  if (!is.null(ferols_fml)) {
    rl_ferols <- do.call(run_ferols, c(ferols_fml, list(df), ferols_parms))
    rl <- c(rl, ferols = list(c(ferols_parms, rl_ferols)))
  } 
  
  rl
}

comp_model <- function(model) {
  if (model == "rlm") {
    out <- pbapply::pbreplicate(
      RUNS, run_ests(
        rlm_fml = y ~ x + z + i + t,
        ferols_fml = y ~ x + z | i + t,
        ferols_parms = list(
          k = 1.345, scale_est = "ols", adj_rlm = TRUE, 
          scale_update = TRUE, tol = 1e-4, vcov = "iid", 
          data.save = TRUE
        )
      ), simplify = FALSE  
    )
  } else if (model == "robreg") {
    out <- pbapply::pbreplicate(
      RUNS, run_ests(
        data_parms = list(),
        ferols_fml = y ~ x + z | i + t, 
        robreg_do_file = "utils/run_robreg.do",
        ferols_parms = list(vcov = ~i, scale_est = "lad_mm_ms", data.save = TRUE)
      ), simplify = FALSE  
    )
  } else if (model == "robhdfe") {
    out <- pbapply::pbreplicate(
      RUNS, run_ests(
        data_parms = list(),
        ferols_fml = y ~ x + z | i + t, 
        robhdfe_do_file = "utils/run_robhdfe.do",
        ferols_parms = list(vcov = ~i, data.save = TRUE)
      ), simplify = FALSE  
    )
  } else if (model == "robreg_robhdfe") {
    out <- pbapply::pbreplicate(
      RUNS, run_ests(
        data_parms = list(),
        robreg_do_file = "utils/run_robreg.do",
        robhdfe_do_file = "utils/run_robhdfe.do"
      ), simplify = FALSE  
    )
  } else stop("Unknown model comparison")
  out
}

create_comparison_runs_df <- function(run_list) {
  stats_to_df <- function(l, vars = VARS) {
    rl <- do.call(rbind, lapply(l, function(l) if (all(vars %in% names(l))) unlist(l[vars])))
    df <- cbind(model = rownames(rl), as.data.frame(rl))
    rownames(df) <- NULL
    df
  }
  extract_df <- function(l) {
    do.call(
      rbind, 
      lapply(seq_along(l), function(x) cbind(run = x, stats_to_df(l[[x]])))
    )
  }
  extract_df(run_list)
}

create_comparison_runs_rdiff <- function(run_df) {
  m <- unique(run_df$model)
  rdiff_vars <- VARS[VARS != "run_time"]
  run_wide <- reshape(run_df, direction = "wide", idvar="run", timevar="model")
  within(data.frame(run = run_wide$run), {
    for (v in rdiff_vars)
      assign(
        paste0("rdiff_", v),
        (run_wide[[paste0(v, ".", m[2])]] - run_wide[[paste0(v, ".", m[1])]])/
          run_wide[[paste0(v, ".", m[1])]]
      )
    v <- NULL
  })
}

plot_rdiffs <- function(rdiff_df) {
  rdiff_vars <- VARS[VARS != "run_time"]
  rdiff_cols <- paste0("rdiff_", rdiff_vars)
  xr <- range(unlist(rdiff_df[rdiff_cols]), na.rm = TRUE)
  op <- par(
    mfrow = c(length(rdiff_cols), 1), mar = c(3,4,2,1)
  )                
  for (nm in rdiff_cols) {
    hist(rdiff_df[[nm]], xlim = xr, main = nm, xlab = "", ylab = "")
  }
  mtext("Relative Difference", side = 1, outer = TRUE, line = 1)
  par(op)
}

create_comparison_stats <- function(rdiff_df) {
  rdiff_vars <- VARS[VARS != "run_time"]
  rdiff_cols <- paste0("rdiff_", rdiff_vars)
  data.frame(
    var = rdiff_vars,
    runs = as.numeric(lapply(rdiff_df[rdiff_cols], function(x) sum(!is.na(x)))),
    rdiff_min = as.numeric(lapply(rdiff_df[rdiff_cols], min)),
    rdiff_max = as.numeric(lapply(rdiff_df[rdiff_cols], max)),
    rdiff_amean = as.numeric(lapply(rdiff_df[rdiff_cols], function(x) mean(abs(x))))
  )
}

comp_stats <-data.frame()

for (m in c("rlm", "robreg", "robhdfe", "robreg_robhdfe")) {
  fname <- sprintf("temp/comp_%s_%d.rds", m, RUNS)
  if (file.exists(fname) & ! FORCE_RERUN) {
    message(sprintf("Data file '%s' exists. Skipping.", fname))
    run_list <- readRDS(fname)
  } else {
    message(sprintf("Running Benchmark '%s' (n = %d)", m, RUNS))
    run_list <- comp_model(m)
    message(sprintf("Saving benchmark data to '%s'.", fname))
    saveRDS(run_list, fname)
  }
  run_df <- create_comparison_runs_df(run_list)
  rdiff_df <- create_comparison_runs_rdiff(run_df)
  plot_rdiffs(rdiff_df)
  comp_stats <- rbind(
    comp_stats, cbind(model = m, create_comparison_stats(rdiff_df))
  )
}

saveRDS(comp_stats, "utils/comp_stats.rds")

