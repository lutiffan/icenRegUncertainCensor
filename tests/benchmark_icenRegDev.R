#!/usr/bin/env Rscript
# Runtime benchmarks for icenRegDev fork (standalone process).

pkg_root <- local({
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg)) {
    normalizePath(file.path(dirname(sub("^--file=", "", file_arg)), ".."))
  } else if (file.exists("DESCRIPTION")) {
    normalizePath(".")
  } else {
    normalizePath("..")
  }
})

suppressPackageStartupMessages({
  if (!requireNamespace("survival", quietly = TRUE)) {
    stop("survival is not installed")
  }
  library(survival)
  if (requireNamespace("icenRegDev", quietly = TRUE)) {
    library(icenRegDev)
  } else if (requireNamespace("pkgload", quietly = TRUE) &&
             file.exists(file.path(pkg_root, "DESCRIPTION"))) {
    pkgload::load_all(pkg_root, quiet = TRUE)
  } else {
    stop("icenRegDev is not installed")
  }
})

source(file.path(pkg_root, "tests", "benchmark_common.R"))

bench <- benchmark_setup("icenRegDev")

results <- rbind(
  time_fit(
    icenRegDev::ic_sp(bench$form, data = bench$sim_data, lcProb = NULL,
                      bs_samples = 0, controls = bench$ctrl),
    "icenRegDev: lcProb NULL, bs_samples = 0"
  ),
  time_fit(
    icenRegDev::ic_sp(bench$mix$form_mix, data = bench$mix$dat_mix, lcProb = bench$mix$lc_mix,
                      bs_samples = 0, controls = bench$ctrl),
    "icenRegDev: lcProb mixture, bs_samples = 0"
  ),
  time_fit(
    icenRegDev::ic_sp(bench$form, data = bench$sim_data, lcProb = NULL,
                      bs_samples = 100, useMCores = FALSE, controls = bench$ctrl),
    "icenRegDev: lcProb NULL, bs_samples = 100",
    reps = 1L
  )
)

results$package <- "icenRegDev"
out <- benchmark_output_path("icenRegDev")
saveRDS(results, out)
cat("Wrote", out, "\n")
print(results, row.names = FALSE)
