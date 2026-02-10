library(ferols)

RUNS <- 100
VARS <- c("run_time", "est", "se", "scale")
STATA_PATH <- "/Applications/Stata/StataSE.app/Contents/MacOS/stata-se"
BMARK_RLM <- FALSE
BMARK_ROBTWFE <- FALSE
BMARK_ROBREG <- TRUE

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
    run_time = as.vector(time_out - time_in),
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
    run_time = as.vector(time_out - time_in),
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
    robreg_do_file = NULL, robtwfe_do_file = NULL,
    rlm_parms = list(), ferols_parms = list(), data_parms = list()
  ) {
  if (
    is.null(rlm_fml) && is.null(ferols_fml) && 
    is.null(robtwfe_do_file) && is.null(robreg_do_file)
  ) {
    stop("At least one estimation method needs to be speic")
  }
  if (!"seed" %in% names(data_parms)) {
    data_parms <- c(data_parms, list(seed = sample(seq_len(1e9), 1)))
  }
  df <- do.call(generate_panel_data, data_parms)
  force(df)
  rl <- list(data_parms = data_parms)
  if (!is.null(rlm_fml)) {
    rl_rlm <- do.call(run_rlm, c(rlm_fml, list(df), rlm_parms))
    rl <- c(rl, rlm = list(c(rlm_parms, rl_rlm)))
  } 
  if (!is.null(ferols_fml)) {
    rl_ferols <- do.call(run_ferols, c(ferols_fml, list(df), ferols_parms))
    rl <- c(rl, ferols = list(c(ferols_parms, rl_ferols)))
  } 
  if (!is.null(robreg_do_file)) {
    rl_robreg <- run_stata(df, robreg_do_file)
    rl <- c(rl, robreg = list(c(list(do_file = robreg_do_file), rl_robreg)))
  } 
  if (!is.null(robtwfe_do_file)) {
    rl_robtwfe <- run_stata(df, robtwfe_do_file)
    rl <- c(rl, robtwfe = list(c(list(do_file = robtwfe_do_file), rl_robtwfe)))
  } 
  
  rl
}

out <- list()

if (BMARK_RLM) {
  out <- c(out, pbapply::pbreplicate(
    RUNS, run_ests(
      rlm_fml = y ~ x + z + i + t,
      ferols_fml = y ~ x + z | i + t,
      ferols_parms = list(
        k = 1.345, scale_est = "ols", adj_rlm = TRUE, 
        scale_update = TRUE, tol = 1e-4, vcov = "iid", 
        data.save = TRUE
      )
    ), simplify = FALSE  
  ))
}

if (BMARK_ROBTWFE) {
  out <- c(out, pbapply::pbreplicate(
    RUNS, run_ests(
      ferols_fml = y ~ x + z | i + t, 
      robtwfe_do_file = "utils/run_robtwfe.do",
      ferols_parms = list(vcov = ~i, data.save = TRUE)
    ), simplify = FALSE  
  ))
}

if (BMARK_ROBREG) {
  out <- c(out, pbapply::pbreplicate(
    RUNS, run_ests(
      ferols_fml = y ~ x + z | i + t, 
      robreg_do_file = "utils/run_robreg.do",
      ferols_parms = list(vcov = ~i, k = 1.3449975, data.save = TRUE)
    ), simplify = FALSE  
  ))
}


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
out_df <- extract_df(out)

m <- sort(unique(out_df$model))
w <- reshape(out_df, direction = "wide", idvar="run", timevar="model")
delta_df <- within(data.frame(run = w$run), {
  for (v in VARS)
    assign(
      paste0("d", v),
      w[[paste0(v, ".", m[2])]] - w[[paste0(v, ".", m[1])]]
    )
  v <- NULL
})

diff_cols <- paste0("d", VARS[VARS != "run_time"])
xr <- range(unlist(delta_df[diff_cols]), na.rm = TRUE)
op <- par(
  mfrow = c(length(diff_cols), 1), mar = c(3,4,2,1)
)                
for (nm in diff_cols) {
  hist(delta_df[[nm]], xlim = xr, main = nm, xlab = "", ylab = "")
}
mtext("Difference", side = 1, outer = TRUE, line = 1)
par(op)
