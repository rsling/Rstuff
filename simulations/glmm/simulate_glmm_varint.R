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
source("utils.R")

fileprefix  <- NULL # "./output/var.int"
nsim        <-  10
J           <-  10
I           <-  10
beta1       <-   0.8
beta2       <-   1
alpha0      <-  -0.5
sigma       <-   0.6
do.r2       <- T

colfunc    <- colorRampPalette(c("gold", "darkblue"))
lwd        <- 3
lwd.small  <- 2
lwd.null   <- 1.5
lty.null   <- 5
lty        <- c(1:4,6:10)

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
                            sigma = sigma)
  } else {
    .run <- sim.glmm.varint(J = J, I = I,
                            beta1 = beta1, beta2 = beta2, alpha0 = alpha0,
                            sigma = sigma,
                            raneffs = .run$raneffs)
  }

  # Get results.
  glmm.raneffs.alpha[i,] <- ranef(.run$glmm)$group[,1]
  glmm.fixeffs[i,]       <- fixef(.run$glmm)
  glmm.p[i,]             <- coef(summary(.run$glmm))[,4]
  glm.coefs[i,]          <- coef(.run$glm)
  glm.p[i,]              <- coef(summary(.run$glm))[,4]
  glm.f.coefs[i,]        <- coef(.run$glm.f)
  glm.f.p[i,]            <- coef(summary(.run$glm.f))[,4]
  if (do.r2) {
    r.squared[i,1]       <- NagelkerkeR2(.run$glm)$R2
    r.squared[i,2]       <- NagelkerkeR2(.run$glm.f)$R2
    r.squared[i,3:4]     <- suppressMessages(r.squaredGLMM(.run$glmm))
  }
  sigmas[i,]             <- as.data.frame(VarCorr(.run$glmm))[,"sdcor"]
}

# Save the alphas as actually used.
true.raneffs <- .run$raneffs


# OUTPUT

if (!is.null(fileprefix)) sink(paste0(fileprefix, '.txt'))
dump.parameters()

cat("\n\n ### Sample GLMM output\n")
print(summary(.run$glmm))
cat("\n\n ### Sample GLM output (ignore raneffs)\n")
print(summary(.run$glm))
cat("\n\n ### Sample GLM output (raneffs as fixeffs)\n")
print(summary(.run$glm))
cat("\n\n")

plot.fixeffs(glmm.fixeffs, c("alpha0", "beta1", "beta2"), c(alpha0, beta1, beta2),
             c("darkorange", "darkgreen", "darkred"), lwd = lwd, lty = lty,
             fileprefix = fileprefix)
plot.raneff.variance(sigmas, "sigma", sigma, "darkorange", lwd = lwd, lty = lty,
                     fileprefix = fileprefix)
plot.raneffs(true.raneffs, glmm.raneffs.alpha, "alpha", sample.size = 8, mfrow = c(2,4),
             lwd = lwd, lty.null = lty.null, colfunc = colfunc,
             fileprefix = fileprefix)

print.raneff.variance(sigmas, sigma)
print.fixeff.comp(glmm.p, glm.p)
print.fixeff.p.comp(glmm.p, glm.p)

if (do.r2) {
  plot.r2(r.squared, c("darkorange", "darkgreen"), lwd = lwd, lty = lty,
          fileprefix = fileprefix)
  print.r2.comp(r.squared)
}

if (!is.null(fileprefix)) sink()

save.image(file = "simulate_glmm_varints.RData")
