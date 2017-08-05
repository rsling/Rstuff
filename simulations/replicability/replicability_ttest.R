# This work is licensed under a Creative Commons Attribution 4.0 International License.
# https://creativecommons.org/licenses/by/4.0/
# Author: Roland Schäfer, Freie Universität Berlin, 2017
# roland.schaefer@fu-berlin.de

rm(list = ls())
set.seed(3675)

# General settings.
nsim  <- 1000
n     <- 20
alpha <- 0.05

# For false null simulation.
mu1   <- 0.5
mu2   <- 1.2
sd    <- 0.8

true.null.results  <- matrix(rep(NA, nsim * 2), nrow = nsim, byrow = T)
false.null.results <- matrix(rep(NA, nsim * 2), nrow = nsim, byrow = T)

colnames(true.null.results)  <- c("t.statistic", "p.value")
colnames(false.null.results) <- c("t.statistic", "p.value")


# TRUE NULL

for (i in 1:nsim) {
  .s1                    <- rnorm(n, mean = mu1, sd = sd)
  .s2                    <- rnorm(n, mean = mu1, sd = sd)
  .test                  <- t.test(.s1, .s2, conf.level = alpha)
  true.null.results[i,]  <- c(.test$statistic, .test$p.value)
}
true.null.reject.idx         <- which(true.null.results[1:nrow(true.null.results)-1, "p.value"] <= alpha)
true.null.reject.count       <- length(true.null.reject.idx)
true.null.reject.replicate   <- length(which(true.null.results[true.null.reject.idx+1, "p.value"] <= alpha))
true.null.reject.nreplicate  <- true.null.reject.count - true.null.reject.replicate

# FALSE NULL

for (i in 1:nsim) {
  .s1                    <- rnorm(n, mean = mu1, sd = sd)
  .s2                    <- rnorm(n, mean = mu2, sd = sd)
  .test                  <- t.test(.s1, .s2, conf.level = alpha)
  false.null.results[i,] <- c(.test$statistic, .test$p.value)
}
false.null.reject.idx         <- which(false.null.results[1:nrow(false.null.results)-1, "p.value"] <= alpha)
false.null.reject.count       <- length(false.null.reject.idx)
false.null.reject.replicate   <- length(which(false.null.results[false.null.reject.idx+1, "p.value"] <= alpha))
false.null.reject.nreplicate  <- false.null.reject.count - false.null.reject.replicate

cat("#############################################\n")
cat("### Simulation of 'replicability': t-test ###\n")
cat("#############################################\n\n")

cat("\nNOTE: Replication success/failure rates depend on power!\n")
cat("      Power analysis for situation with false null is given below.\n")

cat("\nNOTE: For true null, both samples are taken from N(0,1).\n")

cat("\nSettings for this simulation:\n\n")
cat("nsim  =", nsim, "\n")
cat("n     =", n, "\n")
cat("alpha =", alpha, "\n")
cat("mu1   =", mu1, "\n")
cat("mu2   =", mu2, "\n")
cat("sd    =", sd, "\n")


cat("\n\n### TRUE NULL\n\n")

cat("Rejections of TRUE null:   ", true.null.reject.count, " (", round(100*true.null.reject.count/(nsim-1), 2), "%)\n", sep="")
cat("Successful replications:   ", true.null.reject.replicate, " (", round(100*true.null.reject.replicate/(true.null.reject.count), 2), "%)\n", sep="")
cat("Unsuccessful replications: ", true.null.reject.nreplicate, " (", round(100*true.null.reject.nreplicate/(true.null.reject.count), 2), "%)\n", sep="")

cat("\n\n### FALSE NULL\n\n")
print(power.t.test(n = n, delta = mu2-mu1, sd = sd, sig.level = alpha))

cat("Rejections of FALSE null:  ", false.null.reject.count, " (", round(100*false.null.reject.count/(nsim-1), 2), "%)\n", sep="")
cat("Successful replications:   ", false.null.reject.replicate, " (", round(100*false.null.reject.replicate/(false.null.reject.count), 2), "%)\n", sep="")
cat("Unsuccessful replications: ", false.null.reject.nreplicate, " (", round(100*false.null.reject.nreplicate/(false.null.reject.count), 2), "%)\n", sep="")
