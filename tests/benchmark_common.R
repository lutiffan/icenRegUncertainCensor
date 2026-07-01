# Shared data and controls for icenReg / icenRegDev runtime benchmarks.
# Sourced by benchmark_icenReg.R and benchmark_icenRegDev.R (separate processes).

benchmark_pkg_root <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg)) {
    return(normalizePath(file.path(dirname(sub("^--file=", "", file_arg)), "..")))
  }
  if (file.exists("DESCRIPTION")) {
    return(normalizePath("."))
  }
  if (file.exists("../DESCRIPTION")) {
    return(normalizePath(".."))
  }
  stop("Could not locate package root")
}

time_fit <- function(expr, label, reps = 3L) {
  times <- numeric(reps)
  for (r in seq_len(reps)) {
    gc()
    times[r] <- system.time(force(expr))["elapsed"]
  }
  data.frame(
    scenario = label,
    median_sec = median(times),
    min_sec = min(times),
    max_sec = max(times),
    stringsAsFactors = FALSE
  )
}

benchmark_setup <- function(pkg = c("icenReg", "icenRegDev")) {
  pkg <- match.arg(pkg)
  simIC_weib <- get("simIC_weib", envir = asNamespace(pkg))
  makeCtrls_icsp <- get("makeCtrls_icsp", envir = asNamespace(pkg))
  set.seed(2026)
  list(
    sim_data = simIC_weib(n = 100, inspections = 5, inspectLength = 1),
    form = survival::Surv(l, u, type = "interval2") ~ x1 + x2,
    ctrl = makeCtrls_icsp(maxIter = 500),
    mix = local({
      n_mix <- 100
      x1 <- rnorm(n_mix)
      t_exact <- exp(0.2 * x1 + rnorm(n_mix, sd = 0.4))
      y_mix <- cbind(t_exact, t_exact)
      y_mix[sample(n_mix, 30), 2] <- Inf
      dat_mix <- data.frame(l = y_mix[, 1], u = y_mix[, 2], x1 = x1)
      lc_mix <- rep(0, n_mix)
      exact_idx <- which(y_mix[, 1] == y_mix[, 2])
      lc_mix[exact_idx] <- 0.5
      list(
        dat_mix = dat_mix,
        lc_mix = lc_mix,
        form_mix = survival::Surv(l, u, type = "interval2") ~ x1
      )
    })
  )
}

benchmark_output_path <- function(package) {
  file.path(benchmark_pkg_root(), "tests", paste0("benchmark_", package, ".rds"))
}
