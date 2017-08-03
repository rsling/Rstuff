# This work is licensed under a Creative Commons Attribution 4.0 International License.
# https://creativecommons.org/licenses/by/4.0/
# Author: Roland Schäfer, Freie Universität Berlin, 2017
# roland.schaefer@fu-berlin.de

require(lme4)
require(boot)
library(mvtnorm)
require(MuMIn)
require(fmsb)

rm(list = ls())
set.seed(2707)

source("simulate_glmm_varint_fun.R")

use.saved   <- F
nsim        <-  10
J           <-  10
I           <-  10
beta1       <-   0.8
beta2       <-   1
alpha0      <-  -0.5
sigma_a     <-   0.6
do.r2       <- T

colfunc    <- colorRampPalette(c("gold", "darkblue"))
lwd        <- 3
lwd.small  <- 2
lwd.null   <- 1.5
lty.null   <- 5


if (use.saved) {
  load("simulate_glmm_varint.RData")
} else {
  
  # Matrices for results.
  glmm.raneffs.alpha           <- as.data.frame(matrix(rep(NA, J * nsim), nrow = nsim, byrow = T))
  glmm.fixeffs                 <- as.data.frame(matrix(rep(NA, 3 * nsim), nrow = nsim, byrow = T))
  glmm.p                       <- as.data.frame(matrix(rep(NA, 3 * nsim), nrow = nsim, byrow = T))
  glm.coefs                    <- as.data.frame(matrix(rep(NA, 3 * nsim), nrow = nsim, byrow = T))
  glm.p                        <- as.data.frame(matrix(rep(NA, 3 * nsim), nrow = nsim, byrow = T))
  glm.f.coefs                  <- as.data.frame(matrix(rep(NA, (J+3-1) * nsim), nrow = nsim, byrow = T))
  glm.f.p                      <- as.data.frame(matrix(rep(NA, (J+3-1) * nsim), nrow = nsim, byrow = T))
  r.squared                    <- as.data.frame(matrix(rep(NA, 4 * nsim), nrow = nsim, byrow = T))
  sigmas                       <- as.data.frame(matrix(rep(NA, 1 * nsim), nrow = nsim, byrow = T))
  colnames(glmm.raneffs.alpha) <- paste0("group", 1:J)
  colnames(glmm.fixeffs)       <- c('alpha0', 'beta1', 'beta2')
  colnames(glmm.p)             <- c('alpha0', 'beta1', 'beta2')
  colnames(glm.coefs)          <- c('alpha0', 'beta1', 'beta2')
  colnames(glm.p)              <- c('alpha0', 'beta1', 'beta2')
  colnames(glm.f.coefs)        <- c('alpha0', 'beta1', 'beta2', paste0("group", 2:J))
  colnames(glm.f.p)            <- c('alpha0', 'beta1', 'beta2', paste0("group", 2:J))
  colnames(r.squared)          <- c('glm', 'glm.f', 'marginal', 'conditional')
  colnames(sigmas)             <- "sigma"

  for (i in 1:nsim) {
    cat("Simulation run", i, "...\n")
    
    # In all except the first run, we re-use the alphas to get comparable results.
    if (i == 1) {
      .run <- sim.glmm.varint(J = J, I = I,
                              beta1 = beta1, beta2 = beta2, alpha0 = alpha0,
                              sigma_a = sigma_a)
    } else {
      .run <- sim.glmm.varint(J = J, I = I,
                              beta1 = beta1, beta2 = beta2, alpha0 = alpha0,
                              sigma_a = sigma_a,
                              raneffs = .run$raneffs)
    }

    # Get results.
    glmm.raneffs.alpha[i,]      <- ranef(.run$glmm)$group[,1]
    glmm.fixeffs[i,]            <- fixef(.run$glmm)
    glmm.p[i,]                  <- coef(summary(.run$glmm))[,4]
    glm.coefs[i,]               <- coef(.run$glm)
    glm.p[i,]                   <- coef(summary(.run$glm))[,4]
    glm.f.coefs[i,]             <- coef(.run$glm.f)
    glm.f.p[i,]                 <- coef(summary(.run$glm.f))[,4]
    if (do.r2) {
      r.squared[i,1]   <- NagelkerkeR2(.run$glm)
      r.squared[i,2]   <- NagelkerkeR2(.run$glm.f)
      r.squared[i,3:4] <- suppressMessages(r.squaredGLMM(.run$glmm))
    }
    sigmas[i,]                  <- as.data.frame(VarCorr(.run$glmm))[,"sdcor"]
  }
}

# Save the alphas as actually used.
true.raneffs <- .run$raneffs

par(mfrow=c(1,1))
plot(density(glmm.fixeffs$alpha0),
     xlim = c(min(c(alpha0, beta1, beta2, as.matrix(glmm.fixeffs))),
              max(c(alpha0, beta1, beta2, as.matrix(glmm.fixeffs)))),
     ylim = c(min(c(density(glmm.fixeffs$alpha0)$y, density(glmm.fixeffs$beta1)$y, density(glmm.fixeffs$beta2)$y)),
              max(c(density(glmm.fixeffs$alpha0)$y, density(glmm.fixeffs$beta2)$y, density(glmm.fixeffs$beta2)$y)*1.2) ),
     col = "darkorange", lwd = lwd,
     main = "Estimates of fixed effects in GLMM",
     xlab = "Estimates")
lines(density(glmm.fixeffs$beta1),
      col = "darkgreen", lwd = lwd)
lines(density(glmm.fixeffs$beta2),
      col = "darkred", lwd = lwd)
abline(v = alpha0, col = "darkorange", lwd = lwd, lty = 3)
abline(v = beta1, col = "darkgreen", lwd = lwd, lty = 3)
abline(v = beta2, col = "darkred", lwd = lwd, lty = 3)
legend("top",
       legend = c("alpha0", "beta1", "beta2"),
       col = c("darkorange", "darkgreen", "darkred"),
       lwd = lwd)

plot(density(sigmas$sigma),
     xlim = c( min(sigmas$sigma),
               max(sigmas$sigma)),
     ylim = c( min(density(sigmas$sigma)$y),
               max(density(sigmas$sigma)$y)),
     col = "darkorange", lwd = lwd,
     main = "Variance estimates for random effect",
     xlab = "Estimates")
abline(v = sigma_a, col = "darkorange", lwd = lwd, lty = 3)
legend("topleft",
       legend = "sigma(alpha)",
       col = "darkorange",
       lwd = lwd)


par(mfrow=c(2,4))
alphas.sample.plot <- sort(sample(1:nrow(true.raneffs), size = 8, replace = F))
for (i in alphas.sample.plot) {
  plot(density(glmm.raneffs.alpha[,i]),
       xlim = c( min(0, c(true.raneffs[i,"alpha"], as.matrix(glmm.raneffs.alpha[,i]))),
                 max(c(0, true.raneffs[i,"alpha"], as.matrix(glmm.raneffs.alpha[,i])))),
       xlab = "Predicted alpha", main = paste0("J_", i),
       lwd = lwd.small, col = colfunc(nrow(true.raneffs))[i])
  abline(v = true.raneffs[i,"alpha"], lwd = lwd.small, col = colfunc(nrow(true.raneffs))[i], lty = 3)
  abline(v = 0, lwd = lwd.null, col = "gray", lty = lty.null)
}
par(mfrow=c(1,1))


if (do.r2) {
  par(mfrow=c(1,1))
  plot(density(r.squared[,3]),
       xlim = c( min(c(r.squared[,3], r.squared[,4])),
                 max(c(r.squared[,3], r.squared[,4]))),
       ylim = c( min(c(density(r.squared[,3])$y, density(r.squared[,4])$y)),
                 max(c(density(r.squared[,3])$y, density(r.squared[,4])$y))),
       col = "darkorange", lwd = lwd,
       main = "Estimates of R-squared in GLMM",
       xlab = "Estimates")
  lines(density(r.squared[,4]), col = "darkgreen", lwd = lwd)
  legend("topright",
         legend = c("marginal R-squared", "conditional R-squared"),
         col = c("darkorange", "darkgreen"),
         lwd = lwd)
}

cat("\n\n Difference between marginal and conditional\n R-squared in GLMM (95% interval)\n")
print(quantile(ecdf(r.squared[,4]-r.squared[,3]), probs = c(0.025, 0.975)))
cat("\n")


if (!use.saved) save.image(file = "simulate_glmm_varints.RData")
