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
set.seed(9817)

source("simulate_glmm_varint2lp_fun.R")
source("utils.R")

fileprefix  <- "./output/var.int.2level.j=10.i=10"
nsim       <-  1000
J          <-  10
I          <-  10
beta1      <-   0.8
beta2      <-  -1.3
alpha0     <-  -0.5
gamma      <-   1.2
sigma      <-   0.6
do.r2      <- T

colfunc    <- colorRampPalette(c("gold", "darkblue"))
lwd        <- 3
lwd.small  <- 2
lwd.null   <- 1.5
lty.null   <- 5
lty        <- c(1,2,4,6:10)

  
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
colnames(glmm.raneffs.alpha)  <- paste0("group", char.seq(1, J))
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
    r.squared[i,1]        <- NagelkerkeR2(.run$glm)$R2
    r.squared[i,2]        <- NagelkerkeR2(.run$glm.f)$R2
    r.squared[i,3:4]      <- suppressMessages(r.squaredGLMM(.run$glmm))
  }
  sigmas[i,]              <- as.data.frame(VarCorr(.run$glmm))[,"sdcor"]
  
}

# Save the alphas as actually used.
true.raneffs <- .run$raneffs


# OUTPUT

if (!is.null(fileprefix)) sink(paste0(fileprefix, '.txt'))
cat("\n\n ##################################\n")
cat(" Parameters used in this simulation\n")
cat(" ##################################\n")
if (!is.null(nsim)) cat("\nnsim =", nsim)
if (!is.null(J)) cat("\nJ =", J)
if (!is.null(I)) cat("\nI =", I)
if (!is.null(beta1)) cat("\nbeta1 =", beta1)
if (!is.null(beta2)) cat("\nbeta2 =", beta2)
if (!is.null(alpha0)) cat("\nalpha0 =", alpha0)
if (!is.null(gamma)) cat("\ngamma =", gamma)
if (!is.null(sigma)) cat("\nsigma =", sigma)
cat(date(), "\n\n")

cat("\n\n ### Sample GLMM output\n")
print(summary(.run$glmm))
cat("\n\n ### Sample GLM output (ignore raneffs)\n")
print(summary(.run$glm))
cat("\n\n ### Sample GLM output (raneffs as fixeffs)\n")
print(summary(.run$glm))
cat("\n\n")

plot.fixeffs(glmm.fixeffs, c("alpha0", "beta1", "beta2", "gamma"), c(alpha0, beta1, beta2, gamma),
             c("darkgray", "darkgreen", "darkred", "darkblue"), lwd = lwd, lty = lty,
             fileprefix = fileprefix)
plot.raneff.variance(sigmas, "sigma", sigma, "darkgray", lwd = lwd, lty = lty,
                     fileprefix = fileprefix)
plot.raneffs(true.raneffs, glmm.raneffs.alpha, "alpha", sample.size = 8, mfrow = c(2,4),
             lwd = lwd, lty.null = lty.null, colfunc = colfunc,
             fileprefix = fileprefix)


if (!is.null(fileprefix)) this.fileprefix <- paste0(fileprefix, '_estimates') else this.fileprefix <- NULL
plot.fixeff.comparison(glmm.fixeffs, glm.coefs,
                       l.col = c("gray", "black"),
                       p.col = c("darkgreen", "darkblue"),
                       pch   = c(15, 16),
                       main = "Comparison of fixed effects estimates",
                       fileprefix = this.fileprefix)
if (!is.null(fileprefix)) this.fileprefix <- paste0(fileprefix, '_pvalues') else this.fileprefix <- NULL
plot.fixeff.comparison(glmm.p, glm.p, 
                       l.col = c("gray", "black"),
                       p.col = c("darkgreen", "darkblue"),
                       pch   = c(15, 16),
                       main = "Comparison of p-values for fixed effects",
                       fileprefix = this.fileprefix)


print.raneff.variance(sigmas, sigma)
print.fixeff.comp(glmm.fixeffs, glm.coefs)
print.fixeff.p.comp(glmm.p, glm.p)

if (do.r2) {
  plot.r2(r.squared, c("darkgray", "darkgreen"), lwd = lwd, lty = lty)
  print.r2.comp(r.squared)
}

cat("\n\n ### DUMP OF TRUE RANDOM EFFECTS AND PREDICTIONS \n\n")
alpha.analysis <- dump.raneffs(true.raneffs, glmm.raneffs.alpha)

cat("\n\n ### DUMP OF WARNINGS \n\n")
print(warnings())

if (!is.null(fileprefix)) sink()

if (!is.null(fileprefix)) save.image(file = paste0(fileprefix, ".RData"))
