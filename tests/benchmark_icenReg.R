#!/usr/bin/env Rscript
# Runtime benchmarks for CRAN icenReg (standalone process).

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
  if (!requireNamespace("icenReg", quietly = TRUE)) {
    stop("icenReg is not installed")
  }
  library(survival)
  library(icenReg)
})

source(file.path(pkg_root, "tests", "benchmark_common.R"))

bench <- benchmark_setup("icenReg")

results <- rbind(
  time_fit(
    ic_sp(bench$form, data = bench$sim_data, bs_samples = 0, controls = bench$ctrl),
    "icenReg: bs_samples = 0"
  ),
  time_fit(
    ic_sp(bench$form, data = bench$sim_data, bs_samples = 100, useMCores = FALSE, controls = bench$ctrl),
    "icenReg: bs_samples = 100",
    reps = 1L
  )
)

results$package <- "icenReg"
out <- benchmark_output_path("icenReg")
saveRDS(results, out)
cat("Wrote", out, "\n")
print(results, row.names = FALSE)
