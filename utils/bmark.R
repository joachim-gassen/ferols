# Util code to benchmark the processing times of various R/Stata packages
# that generate M regression with Huber loss.

# coef, se, and scale estimates from MASS::rlm() will differ from
# ferols() and Stata packages as rlm() uses a different algorithm.
# Check utils/compare_ests.R to see how ferols() can be parameterized
# to return rlm() type estimates.

# This requires the tidyverse as my base R data wrangling skills are 
# rusty...

suppressPackageStartupMessages(suppressMessages({
  library(ferols)
  library(dplyr)
  library(tidyr)
}))

RUNS <- 10
NT <- seq(100, 1000, 100)
NI <- seq(100, 1000, 100)

STATA_PATH <- "/Applications/Stata/StataSE.app/Contents/MacOS/stata-se"
BMARK_FILE <- "utils/bmark_est_res.rds"

# Running this benchmark takes >48h so only do this if really necessary
FORCE_RERUN <- FALSE

set.seed(42)

if (file.exists(BMARK_FILE) & !FORCE_RERUN) {
  res <- readRDS(BMARK_FILE)
} else {
  do_stata_file <- function(do_file, spath = STATA_PATH, log = FALSE) {
    if(.Platform$OS.type != "unix") stop(
      "This only works on unix type OSes (inlcuding MacOS)."
    )
    if (log) system(paste(spath, "-b", do_file))
    else system(sprintf("%s < %s > /dev/null 2>&1", spath, do_file))
  }
  
  run_stata <- function(do_file) {
    haven::write_dta(sim_dta, "temp/df.dta")
    do_stata_file(do_file)
    as.list(read.delim("temp/stata_out.tsv"))
  }
  
  run_ferols <- function(fml, ...) {
    time_in <- Sys.time()
    fit <- ferols(fml, sim_dta, ...)
    time_out <- Sys.time()
    list(
      model = "ferols",
      run_time = as.vector(difftime(time_out, time_in, units = "secs")),
      est = coef(fit)[["x"]],
      se = fit$se[["x"]],
      scale = fit$robust$scale
    )
  }
  
  run_rlm <- function(fml, ...) {
    tryCatch({
      time_in <- Sys.time()
      fit <- MASS::rlm(fml, sim_dta, ...)
      time_out <- Sys.time()
      list(
        model = "rlm",
        run_time = as.vector(difftime(time_out, time_in, units = "secs")),
        est = coef(fit)[["x"]],
        se = summary(fit)$coef["x", 2],
        scale = fit$s
      )
    }, error = function(e) {
      list(
        model = "rlm",
        run_time = NA_real_,
        est = NA_real_,
        se = NA_real_,
        scale = NA_real_
      )
    })
  }
  
  run_robreg <- function() {
    c(list(model = "robreg"), run_stata("utils/run_robreg.do"))
  }
  
  run_robhdfe <- function() {
    c(list(model = "robhdfe"), run_stata("utils/run_robhdfe.do"))
  }
  
  run_ests <- function(data_parms = list()) {
    if (!"seed" %in% names(data_parms)) {
      data_parms <- c(data_parms, list(seed = sample(seq_len(1e9), 1)))
    }
    sim_dta <<- do.call(generate_panel_data, data_parms)
    bind_rows(
      run_ferols(y ~ x + z | i + t, vcov = ~i),
      run_rlm(y ~ x + z + i + t),
      run_robreg(),
      run_robhdfe()
    )
  }
  
  run_iteration <- function(n_units, n_time) {
    message(
      sprintf("Running iteration n_units:%d, n_time:%d...", n_units, n_time)
    )
    rv <- bind_rows(
      pbapply::pbreplicate(
        RUNS, run_ests(
          data_parms = list(n_units = n_units, n_time = n_time)
        ), simplify = FALSE  
      )
    )
    bind_cols(
      n_units = n_units, n_time = n_time, 
      run = rep(seq_len(RUNS), each = nrow(rv)/RUNS), rv
    )
  }
  
  res <- tibble()
  for (i in 1:length(NI)) {
    res <- bind_rows(res, run_iteration(NI[i] , NT[i]))
  }
  saveRDS(res, BMARK_FILE)
}

bmark_table <- res %>%
  mutate(run_time = ifelse(!is.na(est), run_time, NA)) %>%
  group_by(model, n_units, n_time) %>%
  summarise(
    runs = n(),
    mn_run_time = mean(run_time),
    .groups = "drop"
  ) %>% 
  pivot_wider(
    id_cols = c(n_units, n_time, runs), names_from = model,
    values_from = mn_run_time
  )
