# This work is licensed under a Creative Commons Attribution 4.0 International License.
# https://creativecommons.org/licenses/by/4.0/
# Author: Roland Schäfer, Freie Universität Berlin, 2017
# roland.schaefer@fu-berlin.de

require(lme4)
require(boot)
library(mvtnorm)
require(MuMIn)

rm(list = ls())
set.seed(1507)

source("simulate_glmm_varintslope2lp_fun.R")

nsim       <-  10
J          <-  10
I          <-  10
beta1      <-   1
beta2      <-   0.8
alpha0     <-  -0.5
gamma_a    <-   1.2
gamma_b    <-  -0.6
sigma_a    <-   0.6
sigma_b    <-   0.4
rho        <-  -0.5

colfunc    <- colorRampPalette(c("gold", "darkblue"))
lwd        <- 3
lwd.small  <- 2
lwd.null   <- 1.5
lty.null   <- 5
do.r2      <- T

  
# Matrices for results: GLMM.
glmm.raneffs.alpha           <- as.data.frame(matrix(rep(NA, J * nsim), nrow = nsim, byrow = T))
glmm.raneffs.beta            <- as.data.frame(matrix(rep(NA, J * nsim), nrow = nsim, byrow = T))
glmm.fixeffs                 <- as.data.frame(matrix(rep(NA, 5 * nsim), nrow = nsim, byrow = T))
glmm.p                       <- as.data.frame(matrix(rep(NA, 5 * nsim), nrow = nsim, byrow = T))
glm.coefs                    <- as.data.frame(matrix(rep(NA, 4 * nsim), nrow = nsim, byrow = T))
glm.p                        <- as.data.frame(matrix(rep(NA, 4 * nsim), nrow = nsim, byrow = T))
r.squared                    <- as.data.frame(matrix(rep(NA, 4 * nsim), nrow = nsim, byrow = T))
Sigmas                       <- as.data.frame(matrix(rep(NA, 3 * nsim), nrow = nsim, byrow = T))
colnames(glmm.raneffs.alpha) <- paste0("group", 1:J)
colnames(glmm.raneffs.beta)  <- paste0("group", 1:J)
colnames(glmm.fixeffs)       <- c('alpha0', 'beta1', 'beta2', "gamma_a", "gamma_b")
colnames(glmm.p)             <- c('alpha0', 'beta1', 'beta2', "gamma_a", "gamma_b")
colnames(glm.coefs)           <- c('alpha0', 'beta1', 'beta2', "gamma_a")
colnames(glm.p)               <- c('alpha0', 'beta1', 'beta2', "gamma_a")
colnames(r.squared)          <- c('glm', 'glm.f', 'marginal', 'conditional')
colnames(Sigmas)             <- c('sigma_a', 'sigma_b', 'rho')


for (i in 1:nsim) {
  cat("Simulation run", i, "...\n")
  
  # In all except the first run, we re-use the alphas to get comparable results.
  if (i == 1) {
    .run <- sim.glmm.varintslope2lp(J = J, I = I,
                                    beta1 = beta1, beta2 = beta2, alpha0 = alpha0,
                                    gamma_a = gamma_a, gamma_b = gamma_b,
                                    sigma_a = sigma_b, sigma_b = sigma_b, rho = rho)
  } else {
    .run <- sim.glmm.varintslope2lp(J = J, I = I,
                                    beta1 = beta1, beta2 = beta2, alpha0 = alpha0,
                                    gamma_a = gamma_a, gamma_b = gamma_b,
                                    sigma_a = sigma_b, sigma_b = sigma_b, rho = rho,
                                    raneffs = .run$raneffs)
  }
  
  # Get normal GLMM results.
  glmm.raneffs.alpha[i,]   <- ranef(.run$glmm)$group[,1]
  glmm.raneffs.beta[i,]    <- ranef(.run$glmm)$group[,2] 
  glmm.fixeffs[i,]         <- fixef(.run$glmm)
  glmm.p[i,]               <- coef(summary(.run$glmm))[,4]
  glm.coefs[i,]            <- coef(.run$glm)
  glm.p[i,]                <- coef(summary(.run$glm))[,4]
  if (do.r2) {
    r.squared[i,1]         <- NagelkerkeR2(.run$glm)
    r.squared[i,]          <- suppressMessages(r.squaredGLMM(.run$glmm))
  }

  if (is.nan(as.data.frame(VarCorr(.run$glmm))[3,"sdcor"]))
    warning('Covariance is NaN!')
  else
    Sigmas[i,]                <- as.data.frame(VarCorr(.run$glmm))[,"sdcor"]
}

# Save the alphas as actually used.
true.raneffs <- .run$raneffs

# NaN were turned to NA in loop, remove and report.
n.Sigmas <- nrow(Sigmas) 
Sigmas   <- Sigmas[complete.cases(Sigmas),]
m.Sigmas <- n.Sigmas-nrow(Sigmas)
cat("\nFailed runs (NaN in variances):", m.Sigmas)



# OUTPUT

cat("\n\n ### Sample GLMM output\n")
print(summary(.run$glmm))
cat("\n\n ### Sample GLM output (ignore raneffs)\n")
print(summary(.run$glm))
cat("\n\n ### Sample GLM output (raneffs as fixeffs)\n")
print(summary(.run$glm))
cat("\n\n")

plot.fixeffs(glmm.fixeffs, c('alpha0', 'beta1', 'beta2', "gamma_a", "gamma_b"), c(alpha0, beta1, beta2, gamma_a, gamma_b),
             c("darkorange", "darkgreen", "darkred", "darkblue", "darkbrown"), lwd = lwd)
plot.raneff.variance(Sigmas, c("sigma_a", "sigma_b", "rho"), c(sigma_a, sigma_b, rho),
                     c("darkorange", "darkgreen", "darkred"), lwd = lwd)
plot.raneffs(true.raneffs, glmm.raneffs.alpha, "alpha", sample.size = 8, mfrow = c(2,4), lwd = lwd)
plot.raneffs(true.raneffs, glmm.raneffs.alpha, "beta", sample.size = 8, mfrow = c(2,4), lwd = lwd)
print.raneff.variance(Sigmas, c(sigma_a, sigma_b, rho))
print.fixeff.comp(glmm.p, glm.p)
print.fixeff.p.comp(glmm.p, glm.p)

if (do.r2) {
  plot.r2(r.squared, c("darkorange", "darkgreen"), lwd = lwd)
  print.r2.comp(r.squared)
}

save.image(file = "simulate_glmm_varintslope2lp.RData")
