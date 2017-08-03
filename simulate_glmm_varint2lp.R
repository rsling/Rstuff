require(lme4)
require(boot)
library(mvtnorm)
require(MuMIn)

rm(list = ls())
set.seed(2707)

source("simulate_glmm_varint2lp_fun.R")
source("utils.R")

nsim       <-  10
J          <-  10
I          <-  10
beta1      <-   1
beta2      <-   0.8
alpha0     <-  -0.5
gamma      <-   1.2
sigma      <-   0.6

colfunc    <- colorRampPalette(c("gold", "darkblue"))
lwd        <- 3
lwd.small  <- 2
lwd.null   <- 1.5
lty.null   <- 5
do.r2      <- T
  
# Matrices for results: GLMM.
glmm.raneffs.alpha            <- as.data.frame(matrix(rep(NA, J * nsim), nrow = nsim, byrow = T))
glmm.fixeffs                  <- as.data.frame(matrix(rep(NA, 4 * nsim), nrow = nsim, byrow = T))
glmm.p                        <- as.data.frame(matrix(rep(NA, 4 * nsim), nrow = nsim, byrow = T))
glm.coefs                     <- as.data.frame(matrix(rep(NA, 4 * nsim), nrow = nsim, byrow = T))
glm.p                         <- as.data.frame(matrix(rep(NA, 4 * nsim), nrow = nsim, byrow = T))
glm.f.coefs                   <- as.data.frame(matrix(rep(NA, (2*(J-1)+4) * nsim), nrow = nsim, byrow = T))
glm.f.p                       <- as.data.frame(matrix(rep(NA, (2*(J-1)+4) * nsim), nrow = nsim, byrow = T))
r.squared                     <- as.data.frame(matrix(rep(NA, 4 * nsim), nrow = nsim, byrow = T))
sigmas                        <- as.data.frame(matrix(rep(NA, 1 * nsim), nrow = nsim, byrow = T))
colnames(glmm.raneffs.alpha)  <- paste0("group", 1:J)
colnames(glmm.fixeffs)        <- c('alpha0', 'beta1', 'beta2', "gamma")
colnames(glmm.p)              <- c('alpha0', 'beta1', 'beta2', "gamma")
colnames(glm.coefs)           <- c('alpha0', 'beta1', 'beta2', "gamma")
colnames(glm.p)               <- c('alpha0', 'beta1', 'beta2', "gamma")
colnames(glm.f.coefs)         <- c('alpha0', 'beta1', 'beta2', "gamma", paste0("group", char.seq(2, J)), paste0("x_gamma:group", char.seq(2, J)))
colnames(glm.f.p)             <- c('alpha0', 'beta1', 'beta2', "gamma", paste0("group", char.seq(2, J)), paste0("x_gamma:group", char.seq(2, J)))
colnames(r.squared)           <- c('glm', 'glm.f', 'marginal', 'conditional')
colnames(sigmas)              <- c('sigma')


for (i in 1:nsim) {
  cat("Simulation run", i, "...\n")
  
  # In all except the first run, we re-use the alphas to get comparable results.
  if (i == 1) {
    .run <- sim.glmm.varint2lp(J = J, I = I,
                               beta1 = beta1, beta2 = beta2, alpha0 = alpha0,
                               gamma = gamma, sigma = sigma)
  } else {
    .run <- sim.glmm.varint2lp(J = J, I = I,
                               beta1 = beta1, beta2 = beta2, alpha0 = alpha0,
                               gamma = gamma, sigma = sigma,
                               raneffs = .run$raneffs)
  }
  
  # Get normal GLMM results.
  glmm.raneffs.alpha[i,] <- ranef(.run$glmm)$group[,1]
  glmm.fixeffs[i,]       <- fixef(.run$glmm)
  glmm.p[i,]             <- coef(summary(.run$glmm))[,4]
  glm.coefs[i,]          <- coef(.run$glm)
  glm.p[i,]              <- coef(summary(.run$glm))[,4]
  glm.f.coefs[i,]        <- coef(.run$glm.f)
  
  # A lot of NAs here, so we need to be more clever.
  .names <- c("alpha0", "beta1", "beta2", "gamma", names(coef(summary(.run$glm.f))[,4])[-c(1:4)])
  glm.f.p[i,.names]       <- coef(summary(.run$glm.f))[,4]

  if (do.r2) {
    r.squared[i,1]        <- NagelkerkeR2(.run$glm)
    r.squared[i,2]        <- NagelkerkeR2(.run$glm.f)
    r.squared[i,3:4]      <- suppressMessages(r.squaredGLMM(.run$glmm))
  }
  sigmas[i,]              <- as.data.frame(VarCorr(.run$glmm))[,"sdcor"]
  
}

# Save the alphas as actually used.
true.raneffs <- .run$raneffs


# OUTPUT

plot.fixeffs(glmm.fixeffs, c("alpha0", "beta1", "beta2", "gamma"), c(alpha0, beta1, beta2, gamma),
             c("darkorange", "darkgreen", "darkred", "darkblue"), lwd = lwd)
plot.raneff.variance(sigmas, "sigma", sigma, "darkorange", lwd = lwd)
plot.raneffs(true.raneffs, glmm.raneffs.alpha, "alpha", sample.size = 8, mfrow = c(2,4), lwd = lwd)
print.raneff.variance(sigmas, sigma)
print.fixeff.comp(glmm.p, glm.p)
print.fixeff.p.comp(glmm.p, glm.p)

if (do.r2) {
  plot.r2(r.squared, c("darkorange", "darkgreen"), lwd = lwd)
  print.r2.comp(r.squared)
}

save.image(file = "simulate_glmm_varint2lp.RData")
