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
setwd("~/Workingcopies/Rstuff/simulations/glmm")
set.seed(2708)

source("simulate_glmm_varint_fun.R")
source("utils.R")

fileprefix  <- "./output/var.int_j=5.i=20"
nsim        <-  1000
J           <-  5
I           <-  20
beta1       <-   0.8
beta2       <-  -1.3
alpha0      <-   0 # -0.5
sigma       <-   1.5
do.r2       <- T

colfunc    <- colorRampPalette(c("darkgreen", "darkblue"))
lwd        <- 3
lwd.small  <- 2
lwd.null   <- 1.5
lty.null   <- 5
lty        <- c(1,2,4,6:10)

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
colnames(glmm.raneffs.alpha) <- paste0("group", char.seq(1, J))
colnames(glmm.fixeffs)       <- c('alpha0', 'beta1', 'beta2')
colnames(glmm.p)             <- c('alpha0', 'beta1', 'beta2')
colnames(glm.coefs)          <- c('alpha0', 'beta1', 'beta2')
colnames(glm.p)              <- c('alpha0', 'beta1', 'beta2')
colnames(glm.f.coefs)        <- c('alpha0', 'beta1', 'beta2', paste0("group", char.seq(2, J)))
colnames(glm.f.p)            <- c('alpha0', 'beta1', 'beta2', paste0("group", char.seq(2, J)))
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


cat("\n\n ##################################\n")
cat(" Parameters used in this simulation\n")
cat(" ##################################\n")
if (!is.null(nsim)) cat("\nnsim =", nsim)
if (!is.null(J)) cat("\nJ =", J)
if (!is.null(I)) cat("\nI =", I)
if (!is.null(beta1)) cat("\nbeta1 =", beta1)
if (!is.null(beta2)) cat("\nbeta2 =", beta2)
if (!is.null(alpha0)) cat("\nalpha0 =", alpha0)
if (!is.null(sigma)) cat("\nsigma =", sigma)
cat(date(), "\n\n")

cat("\n\n ### Sample GLMM output\n\n")
print(summary(.run$glmm))
cat("\n\n ### Sample GLM output (ignore raneffs)\n")
print(summary(.run$glm))
cat("\n\n ### Sample GLM output (raneffs as fixeffs)\n")
print(summary(.run$glm))
cat("\n\n")

plot.fixeffs(glmm.fixeffs, c("alpha0", "beta1", "beta2"), c(alpha0, beta1, beta2),
             c("darkgray", "darkgreen", "darkred"), lwd = lwd, lty = lty,
             fileprefix = fileprefix)
plot.raneff.variance(sigmas, "sigma", sigma, "darkgray", lwd = lwd, lty = lty,
                     fileprefix = fileprefix)


if (!is.null(fileprefix)) this.fileprefix <- paste0(fileprefix, '_random') else this.fileprefix <- NULL
plot.effs(true.raneffs[-c(1),], glmm.raneffs.alpha[,-c(1)]+glmm.fixeffs[,1], "alpha", sample.size = 4, nrow = 2,
          lwd = lwd, lty.null = lty.null, colfunc = colfunc,
          xlim.perc = c(0, 1),
          plot.0 <- F,
          fileprefix = this.fileprefix)

if (!is.null(fileprefix)) this.fileprefix <- paste0(fileprefix, '_fixed') else this.fileprefix <- NULL
plot.effs(true.raneffs[-c(1),], glm.f.coefs[,-c(1:3)]+glm.f.coefs[,1], "alpha", sample.size = 4, nrow = 2,
          lwd = lwd, lty.null = lty.null, colfunc = colfunc,
          xlim.perc = c(0.5, 1),
          plot.0 <- F,
          fileprefix = this.fileprefix)


if (!is.null(fileprefix)) this.fileprefix <- paste0(fileprefix, '_estimates') else this.fileprefix <- NULL
plot.fixeff.comparison(glmm.fixeffs, glm.coefs, glm.f.coefs,
                       l.col = c("gray", "black"),
                       p.col = c("darkgreen", "darkblue", "darkred"),
                       pch   = c(15, 16, 18),
                       main = "Comparison of fixed effects estimates",
                       fileprefix = this.fileprefix)

plot.fixeff.p.comp(glmm.p, glm.p, c("beta1", "beta2"), log = T)

if (!is.null(fileprefix)) this.fileprefix <- paste0(fileprefix, '_pvalues') else this.fileprefix <- NULL
plot.fixeff.comparison(glmm.p, glm.p, glm.f.p,
                       l.col = c("gray", "black"),
                       p.col = c("darkgreen", "darkblue", "darkred"),
                       pch   = c(15, 16, 18),
                       main = "Comparison of p-values for fixed effects",
                       fileprefix = this.fileprefix)


if (!is.null(fileprefix)) pdf(paste0(fileprefix, "_fixeff.p-comparison.pdf"))
par(mfrow=c(1,2))
plot.fixeff.p.comp(glmm.p, glm.f.p, c("beta1"), 
                   lims = c(0, 1),
                   col = c("darkgreen"),
                   pch = 18, lty = 2,
                   xlab="GLMM", ylab="GLM (ignoring grouping factor)",
                   main="beta1",
                   lines = "lowess")
plot.fixeff.p.comp(glmm.p, glm.f.p, c("beta2"), 
                   lims = c(0, 0.3),
                   col = c("darkblue"),
                   pch = 18, lty = 2,
                   xlab="GLMM", ylab="GLM (ignorign grouping factor)",
                   main="beta2",
                   lines = "lowess")
if (!is.null(fileprefix)) dev.off()
par(mfrow=c(1,1))



# par(mfrow=c(1,2))
# plot.fixeff.p.comp(glm.p, glm.f.p, c("beta1"), 
#                    lims = c(0, 1),
#                    col = c("darkgreen"),
#                    pch = 18, lty = 2,
#                    xlab="GLM (ignoring grouping factor)", ylab="GLM (grouping factor as fixed)",
#                    main="beta1",
#                    lines = "lowess")
# plot.fixeff.p.comp(glm.p, glm.f.p, c("beta2"), 
#                    lims = c(0, 0.3),
#                    col = c("darkblue"),
#                    pch = 18, lty = 2,
#                    xlab="GLM (ignoring grouping factor)", ylab="GLM (grouping factor as fixed)",
#                    main="beta2",
#                    lines = "lowess")
# par(mfrow=c(1,1))


print.raneff.variance(sigmas, sigma)
print.fixeff.comp(glmm.fixeffs, glm.coefs, glm.f.coefs[,1:3])
print.fixeff.p.comp(glmm.p, glm.p, glm.f.p[,1:3])

if (do.r2) {
  plot.r2(r.squared, c("darkgray", "darkgreen"), lwd = lwd, lty = lty,
          fileprefix = fileprefix)
  print.r2.comp(r.squared)
}

cat("\n\n ### DUMP OF TRUE RANDOM EFFECTS AND PREDICTIONS \n\n")
alpha.analysis <- dump.raneffs(true.raneffs, glmm.raneffs.alpha)

cat("\n\n ### DUMP OF WARNINGS \n\n")
print(warnings())

if (!is.null(fileprefix)) sink()

if (!is.null(fileprefix)) save.image(file = paste0(fileprefix, ".RData"))
