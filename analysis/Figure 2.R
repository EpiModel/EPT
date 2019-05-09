## EPT analyses
## Coverage by partnership type on incidence
## Figure 2 ---------------------------------------------------------
## Line plot: Incidence and SA interval/HR Interval
rm(list = ls())
library("EpiModelHIV")
library("dplyr")
library("ggplot2")
library("viridis")
library("gridExtra")


#Baseline: 8000,
# All partners
# Main partners: 8023-8032, cas cov:8039:8048, inst cov: 8055:8064

tiff(filename = "analysis/Figure 2.tiff", height = 8, width = 11, units = "in", res = 250)

par(mfrow = c(2, 2), mar = c(3, 3, 2, 1.2), oma = c(0, 0, 2, 0), mgp = c(2,1,0))
sims <- c(8000, 8002, 8004, 8006, 8008, 8010)
pal <- viridis::viridis(n = length(sims), option = "D")
for (i in seq_along(sims)) {
  fn <- list.files("data", pattern = as.character(sims[i]), full.names = TRUE)
  load(fn)
  #sim <- mutate_epi(sim, ir100.gcct = ir100.ct + ir100.gc)
  plot(sim, y = "ir100.gcct", add = i > 1,
       mean.col = pal[i], qnts.col = pal[i], qnts.alpha = 0.3, qnts = 0.0,
       main = "Provision to All Partners",
       xlab = "Week", ylab = "IR per 100 PYAR", ylim = c(0, 12))
}
legend("bottomleft", legend = c("0% coverage", "20% coverage", "40% coverage",
                                "60% coverage" , "80% coverage", "100% coverage"),
       col = pal, lwd = 3, cex = 0.85, bty = "n")
sims <- c(8000, 8024, 8026, 8028, 8030, 8032)
pal <- viridis::viridis(n = length(sims), option = "D")
for (i in seq_along(sims)) {
  fn <- list.files("data", pattern = as.character(sims[i]), full.names = TRUE)
  load(fn)
  #sim <- mutate_epi(sim, ir100.gcct = ir100.ct + ir100.gc)
  plot(sim, y = "ir100.gcct", add = i > 1,
       mean.col = pal[i], qnts.col = pal[i], qnts.alpha = 0.3, qnts = 0.0,
       main = "Provision to Main Partners",
       xlab = "Week", ylab = "IR per 100 PYAR", ylim = c(0, 12))
}
legend("bottomleft", legend = c("0% coverage", "20% coverage", "40% coverage",
                                "60% coverage" , "80% coverage", "100% coverage"),
       col = pal, lwd = 3, cex = 0.85, bty = "n")

sims <- c(8000, 8040, 8042, 8044, 8046, 8048)
pal <- viridis::viridis(n = length(sims), option = "D")
for (i in seq_along(sims)) {
  fn <- list.files("data", pattern = as.character(sims[i]), full.names = TRUE)
  load(fn)
  #sim <- mutate_epi(sim, ir100.gcct = ir100.ct + ir100.gc)
  plot(sim, y = "ir100.gcct", add = i > 1,
       mean.col = pal[i], qnts.col = pal[i], qnts.alpha = 0.3, qnts = 0.0,
       main = "Provision to Casual Partners",
       xlab = "Week", ylab = "IR per 100 PYAR", ylim = c(0, 12))
}
legend("bottomleft", legend = c("0% coverage", "20% coverage", "40% coverage",
                                "60% coverage" , "80% coverage", "100% coverage"),
       col = pal, lwd = 3, cex = 0.85, bty = "n")

sims <- c(8000, 8056, 8058, 8060, 8062, 8064)
pal <- viridis::viridis(n = length(sims), option = "D")
for (i in seq_along(sims)) {
  fn <- list.files("data", pattern = as.character(sims[i]), full.names = TRUE)
  load(fn)
  #sim <- mutate_epi(sim, ir100.gcct = ir100.ct + ir100.gc)
  plot(sim, y = "ir100.gcct", add = i > 1,
       mean.col = pal[i], qnts.col = pal[i], qnts.alpha = 0.3, qnts = 0.0,
       main = "Provision to One-Time Partners",
       xlab = "Week", ylab = "IR per 100 PYAR", ylim = c(0, 12))
}
legend("bottomleft", legend = c("0% coverage", "20% coverage", "40% coverage",
                                "60% coverage" , "80% coverage", "100% coverage"),
       col = pal, lwd = 3, cex = 0.85, bty = "n")
# mtext("Figure 2: Incidence Rates (IR) by EPT Coverage with Provision Limited by Type of Partnership",
#       outer = TRUE, cex = 1.5)

dev.off()
