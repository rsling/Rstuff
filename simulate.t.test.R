
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
plot(sims)
hist(sims, breaks = 100)
