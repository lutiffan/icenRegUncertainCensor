if (!requireNamespace("testthat", quietly = TRUE)) {
  stop("testthat required")
}

run_benchmark_script <- function(script_path) {
  rscript <- file.path(R.home("bin"), "Rscript")
  out <- system2(rscript, script_path, stdout = TRUE, stderr = TRUE)
  status <- attr(out, "status")
  if (is.null(status)) status <- 0L
  list(status = status, output = out)
}

test_that("runtime comparison: icenReg vs icenRegDev (separate processes)", {
  skip_if_not_installed("icenReg")
  skip_if_not(requireNamespace("icenRegDev", quietly = TRUE))

  pkg_root <- if (file.exists("DESCRIPTION")) "." else ".."
  bench_icenReg <- normalizePath(file.path(pkg_root, "tests", "benchmark_icenReg.R"))
  bench_icenRegDev <- normalizePath(file.path(pkg_root, "tests", "benchmark_icenRegDev.R"))

  for (f in c(bench_icenReg, bench_icenRegDev)) {
    skip_if_not(file.exists(f), paste("missing", f))
  }

  out_cran_rds <- file.path(pkg_root, "tests", "benchmark_icenReg.rds")
  out_dev_rds <- file.path(pkg_root, "tests", "benchmark_icenRegDev.rds")
  on.exit(unlink(c(out_cran_rds, out_dev_rds)), add = TRUE)

  cran_run <- run_benchmark_script(bench_icenReg)
  expect_equal(cran_run$status, 0L)

  dev_run <- run_benchmark_script(bench_icenRegDev)
  expect_equal(dev_run$status, 0L)

  cran <- readRDS(out_cran_rds)
  dev <- readRDS(out_dev_rds)
  combined <- rbind(cran, dev)

  baseline <- cran$median_sec[cran$scenario == "icenReg: bs_samples = 0"]
  baseline_floor <- max(baseline, 0.001)
  combined$ratio_vs_icenReg_bs0 <- combined$median_sec / baseline_floor

  cat("\n=== Runtime comparison (seconds) ===\n")
  print(combined[, c("package", "scenario", "median_sec", "min_sec", "max_sec", "ratio_vs_icenReg_bs0")],
        row.names = FALSE)
  cat("\nRatios use max(icenReg bs_samples = 0, 0.001s) as baseline.\n\n")

  dev_null <- dev$median_sec[dev$scenario == "icenRegDev: lcProb NULL, bs_samples = 0"]
  expect_lt(dev_null, 5 * baseline_floor)

  bs_cran <- cran$median_sec[cran$scenario == "icenReg: bs_samples = 100"]
  bs_dev <- dev$median_sec[dev$scenario == "icenRegDev: lcProb NULL, bs_samples = 100"]
  expect_lt(bs_dev, 5 * bs_cran)
})
