#!/usr/bin/env Rscript
# Illustrates Turnbull baseline from ic_sp via plot() and getSCurves().
#
# Observations: (1,3], (3,5], (5,7]  with default B = c(0, 1)
# adjustIntervals() separates shared endpoints -> three Turnbull intervals.
# Each observation maps to one interval (l_ind == r_ind); S(t) reaches 0 at t = 7.

library(icenRegDev)

lower <- c(1, 3, 5)
upper <- c(3, 5, 7)
obs_labels <- c("Obs 1: (1, 3]", "Obs 2: (3, 5]", "Obs 3: (5, 7]")
dat <- data.frame(l = lower, u = upper)

# Default B = c(0, 1): open left, closed right
fit <- ic_sp(Surv(l, u, type = "interval2") ~ 1, data = dat,
             bs_samples = 0, controls = makeCtrls_icsp(maxIter = 500))

dat_adj <- icenRegDev:::adjustIntervals(c(0, 1), cbind(lower, upper))
mi <- icenRegDev:::findMaximalIntersections(dat_adj[, 1], dat_adj[, 2])
cat("Adjusted intervals (internal):\n")
print(dat_adj)
cat("\nTurnbull intervals (mi_l, mi_r):\n")
print(rbind(mi_l = mi$mi_l, mi_r = mi$mi_r))
cat("\nPer-observation indices:\n")
print(data.frame(obs = 1:3, l_ind = mi$l_inds, r_ind = mi$r_inds))
cat("\nBaseline masses p_hat:\n")
print(fit$p_hat)

curves <- getSCurves(fit)
cat("\ngetSCurves() baseline S at Turnbull cutpoints:\n")
print(curves$S_curves$baseline)

outfile <- "examples/turnbull_two_intervals.png"

png(outfile, width = 900, height = 650, res = 120)
par(mar = c(5, 5, 4, 8))

plot(fit, xlab = "time", ylab = expression(S(t)),
     main = "ic_sp baseline (default B = c(0, 1): open left, closed right)",
     col = "steelblue", lwd = 2)

tb <- t(fit$T_bull_Intervals)
ph <- fit$p_hat
s_at_tb <- 1 - c(0, cumsum(ph))  # S just before each drop

# Turnbull interval boundaries
abline(v = c(tb[, 1], tb[nrow(tb), 2]), col = "firebrick", lty = 2, lwd = 1.5)
for (j in seq_len(nrow(tb))) {
  points(tb[j, 2], s_at_tb[j + 1], pch = 19, col = "firebrick", cex = 1.2)
  text(mean(c(tb[j, 1], tb[j, 2])), -0.07,
       labels = sprintf("TB %d: [%g, %g]\np=%.2f", j - 1L, tb[j, 1], tb[j, 2], ph[j]),
       col = "firebrick", cex = 0.78, xpd = NA)
}

# Observation intervals (horizontal bars)
y_obs <- c(0.88, 0.78, 0.68)
for (i in seq_along(obs_labels)) {
  segments(lower[i], y_obs[i], upper[i], y_obs[i], lwd = 4,
           col = adjustcolor("gray30", 0.75))
  points(c(lower[i], upper[i]), rep(y_obs[i], 2), pch = 124, cex = 1.2, col = "gray30")
  text(7.4, y_obs[i], obs_labels[i], adj = 0, cex = 0.82, col = "gray30")
}

# Note: S reaches 0 at the end of the last Turnbull interval
segments(7, 0.02, 7, s_at_tb[length(s_at_tb)], col = "darkgreen", lwd = 2, lty = 3)
text(7.15, 0.35, "S(t) = 0\nat t = 7", col = "darkgreen", cex = 0.85, adj = 0)

legend("topright",
       legend = c("ic_sp baseline S(t)", "Turnbull interval boundary", "Censoring interval"),
       col = c("steelblue", "firebrick", "gray30"),
       lty = c(1, 2, 1), lwd = c(2, 1.5, 4), bty = "n", cex = 0.9)

dev.off()
message("Saved: ", normalizePath(outfile, mustWork = FALSE))
