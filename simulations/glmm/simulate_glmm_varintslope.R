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
set.seed(9298)

source("simulate_glmm_varintslope_fun.R")
source("utils.R")

fileprefix  <- "./output/var.int.slope.j=50.i=50"
nsim       <-  1000
J          <-  50
I          <-  50
beta1      <-   0.8
beta2      <-  -1.3
alpha0     <-  -0.5
sigma_a    <-   0.6
sigma_b    <-   0.4
rho        <-  -0.5
do.raneff  <- T
do.fixeff  <- F
do.r2      <- T

colfunc    <- colorRampPalette(c("gold", "darkblue"))
lwd        <- 3
lwd.small  <- 2
lwd.null   <- 1.5
lty.null   <- 5
lty        <- c(1,2,4,6:10)
  
# Matrices for results.
glmm.raneffs.alpha           <- as.data.frame(matrix(rep(NA, J * nsim), nrow = nsim, byrow = T))
glmm.raneffs.beta            <- as.data.frame(matrix(rep(NA, J * nsim), nrow = nsim, byrow = T))
glmm.fixeffs                 <- as.data.frame(matrix(rep(NA, 3 * nsim), nrow = nsim, byrow = T))
glmm.p                       <- as.data.frame(matrix(rep(NA, 3 * nsim), nrow = nsim, byrow = T))
glm.coefs                    <- as.data.frame(matrix(rep(NA, 3 * nsim), nrow = nsim, byrow = T))
glm.p                        <- as.data.frame(matrix(rep(NA, 3 * nsim), nrow = nsim, byrow = T))
glm.f.coefs                  <- as.data.frame(matrix(rep(NA, (J+3-1) * nsim), nrow = nsim, byrow = T))
glm.f.p                      <- as.data.frame(matrix(rep(NA, (J+3-1) * nsim), nrow = nsim, byrow = T))
r.squared                    <- as.data.frame(matrix(rep(NA, 4 * nsim), nrow = nsim, byrow = T))
Sigmas                       <- as.data.frame(matrix(rep(NA, 3 * nsim), nrow = nsim, byrow = T))
colnames(glmm.raneffs.alpha) <- paste0("group", char.seq(1, J))
colnames(glmm.raneffs.beta)  <- paste0("group", 1:J)
colnames(glmm.fixeffs)       <- c('alpha0', 'beta1', 'beta2')
colnames(glmm.p)             <- c('alpha0', 'beta1', 'beta2')
colnames(glm.coefs)          <- c('alpha0', 'beta1', 'beta2')
colnames(glm.p)              <- c('alpha0', 'beta1', 'beta2')
colnames(glm.f.coefs)        <- c('alpha0', 'beta1', 'beta2', paste0("group", 2:J))
colnames(glm.f.p)            <- c('alpha0', 'beta1', 'beta2', paste0("group", 2:J))
colnames(r.squared)          <- c('marginal', 'conditional')
colnames(Sigmas)             <- c('sigma_a', 'sigma_b', 'rho')


for (i in 1:nsim) {
  cat("Simulation run", i, "...\n")
  
  # In all except the first run, we re-use the alphas to get comparable results.
  if (i == 1) {
    .run <- sim.glmm.varintslope(J = J, I = I,
                                 beta1 = beta1, beta2 = beta2, alpha0 = alpha0,
                                 sigma_a = sigma_b, sigma_b = sigma_b, rho = rho)
  } else {
    .run <- sim.glmm.varintslope(J = J, I = I,
                                 beta1 = beta1, beta2 = beta2, alpha0 = alpha0,
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
  glm.f.coefs[i,]          <- coef(.run$glm.f)
  glm.f.p[i,]              <- coef(summary(.run$glm.f))[,4]

  if (do.r2) {
    r.squared[i,1]    <- NagelkerkeR2(.run$glm)$R2
    r.squared[i,2]    <- NagelkerkeR2(.run$glm.f)$R2
    r.squared[i,3:4]  <- suppressMessages(r.squaredGLMM(.run$glmm))
  }
  if (is.nan(as.data.frame(VarCorr(.run$glmm))[3,"sdcor"]))
    warning('Covariance is NaN!')
  else
    Sigmas[i,]                <- as.data.frame(VarCorr(.run$glmm))[,"sdcor"]

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
if (!is.null(sigma_a)) cat("\nsigma_a =", sigma_a)
if (!is.null(sigma_b)) cat("\nsigma_b =", sigma_b)
if (!is.null(rho)) cat("\nrho =", rho)
cat(date(), "\n\n")

# NaN were turned to NA in loop, remove.
n.Sigmas <- nrow(Sigmas) 
Sigmas   <- Sigmas[complete.cases(Sigmas),]
m.Sigmas <- n.Sigmas-nrow(Sigmas)
cat("\nFailed runs (NaN in variances):", m.Sigmas)

cat("\n\n ### Sample GLMM output\n")
print(summary(.run$glmm))
cat("\n\n ### Sample GLM output (ignore raneffs)\n")
print(summary(.run$glm))
cat("\n\n ### Sample GLM output (raneffs as fixeffs)\n")
print(summary(.run$glm))
cat("\n\n")

plot.fixeffs(glmm.fixeffs, c("alpha0", "beta1", "beta2"), c(alpha0, beta1, beta2),
             c("darkgray", "darkgreen", "darkred"), lwd = lwd, lty = lty,
             fileprefix = fileprefix)
plot.raneff.variance(Sigmas, c("sigma_a", "sigma_b", "rho"), c(sigma_a, sigma_b, rho),
                     c("darkgray", "darkgreen", "darkred"), lwd = lwd, lty = lty,
                     fileprefix = fileprefix)
plot.raneffs(true.raneffs, glmm.raneffs.alpha, "alpha", sample.size = 8, mfrow = c(2,4),
             lwd = lwd, lty.null = lty.null, colfunc = colfunc,
             fileprefix = fileprefix)


if (!is.null(fileprefix)) this.fileprefix <- paste0(fileprefix, '_estimates') else this.fileprefix <- NULL
plot.fixeff.comparison(glmm.fixeffs, glm.coefs, glm.f.coefs,
                       l.col = c("gray", "black"),
                       p.col = c("darkgreen", "darkblue", "darkred"),
                       pch   = c(15, 16, 18),
                       main = "Comparison of fixed effects estimates",
                       fileprefix = this.fileprefix)
if (!is.null(fileprefix)) this.fileprefix <- paste0(fileprefix, '_pvalues') else this.fileprefix <- NULL
plot.fixeff.comparison(glmm.p, glm.p, glm.f.p,
                       l.col = c("gray", "black"),
                       p.col = c("darkgreen", "darkblue", "darkred"),
                       pch   = c(15, 16, 18),
                       main = "Comparison of p-values for fixed effects",
                       fileprefix = this.fileprefix)


print.raneff.variance(Sigmas, c(sigma_a, sigma_b, rho))
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
