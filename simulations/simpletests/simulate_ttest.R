# This work is licensed under a Creative Commons Attribution 4.0 International License.
# https://creativecommons.org/licenses/by/4.0/
# Author: Roland Schäfer, Freie Universität Berlin, 2017
# roland.schaefer@fu-berlin.de

nsim <- 10000
n    <- 100
mean <- 0
sd   <- 1

# Data structure for results.
sims <- rep(NA, nsim)

# Simulations.
for (i in 1:nsim) {
  a <- rnorm(n, mean = mean, sd = sd)
  b <- rnorm(n, mean = mean, sd = sd)
  p <- t.test(a,b)$p.value
  sims[i] <- p
}

par(mfrow=c(1,1))
plot(sims, pch=19, cex=0.1, col = "darkgreen")
plot(ecdf(sims))
hist(sims, breaks = 100)
