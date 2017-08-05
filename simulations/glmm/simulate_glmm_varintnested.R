# This work is licensed under a Creative Commons Attribution 4.0 International License.
# https://creativecommons.org/licenses/by/4.0/
# Author: Roland Schäfer, Freie Universität Berlin, 2017
# roland.schaefer@fu-berlin.de

require(lme4)
require(boot)
require(MuMIn)
require(fmsb)

rm(list = ls())
set.seed(1131)

source("simulate_glmm_varintnested_fun.R")
source("utils.R")

fileprefix  <- "./output/var.int.nested.j0=10.j1=10.i=20"
nsim        <- 1000
J0          <- 10
J1          <- 10
I           <- 20
beta1       <-  0.8
beta2       <- -1.3
alpha       <- -0.5
sigma0      <-  0.8
sigma1      <-  0.4
nested      <- T    # This creates nested data, nothing in GLMM spec changes.
do.r2       <- T

colfunc     <- colorRampPalette(c("gold", "darkblue"))
lwd         <- 3
lwd.small   <- 3
lwd.null    <- 1.5
lty.null    <- 5
lty         <- c(1,2,4,6:10)

# Matrices for results.
glmm.raneffs.group0           <- as.data.frame(matrix(rep(NA, J0 * nsim),      nrow = nsim, byrow = T))
glmm.raneffs.group1           <- as.data.frame(matrix(rep(NA, J1 * J0 * nsim), nrow = nsim, byrow = T))
glmm.fixeffs                  <- as.data.frame(matrix(rep(NA, 3 * nsim),       nrow = nsim, byrow = T))
glmm.p                        <- as.data.frame(matrix(rep(NA, 3 * nsim),       nrow = nsim, byrow = T))
glm.coefs                     <- as.data.frame(matrix(rep(NA, 3 * nsim),       nrow = nsim, byrow = T))
glm.p                         <- as.data.frame(matrix(rep(NA, 3 * nsim),       nrow = nsim, byrow = T))
r.squared                     <- as.data.frame(matrix(rep(NA, 4 * nsim),       nrow = nsim, byrow = T))
raneff.var                    <- as.data.frame(matrix(rep(NA, 2 * nsim),       nrow = nsim, byrow = T))
colnames(glmm.raneffs.group0) <- paste0("group", 1:J0)
colnames(glmm.raneffs.group1) <- paste0("group", 1:(J1*J0))
colnames(glmm.fixeffs)        <- c('alpha0', 'beta1', 'beta2')
colnames(glmm.p)              <- c('alpha0', 'beta1', 'beta2')
colnames(glm.coefs)           <- c('alpha0', 'beta1', 'beta2')
colnames(glm.p)               <- c('alpha0', 'beta1', 'beta2')
colnames(r.squared)           <- c('glm', 'glm.f', 'marginal', 'conditional')
colnames(raneff.var)          <- c('sigma0', 'sigma1')


# RUN SIMULATIONS

for (i in 1:nsim) {
  cat("Simulation run", i, "...\n")
  
  # In all except the first run, we re-use the alphas to get comparable results.
  if (i == 1) {
    .run <- sim.glmm.varintnested(nested = nested,
                                  J0 = J0, J1 = J1, I = I,
                                  beta1 = beta1, beta2 = beta2,
                                  alpha = alpha,
                                  sigma0 = sigma0, sigma1 = sigma1)
  } else {
    .run <- sim.glmm.varintnested(nested = nested,
                                  J0 = J0, J1 = J1, I = I,
                                  beta1 = beta1, beta2 = beta2,
                                  alpha = alpha,
                                  sigma0 = sigma0, sigma1 = sigma1,
                                  alphas0 = .run$alphas0, alphas1 = .run$alphas1)
  }

  # Checker whether groups are correctly aligned in data frames.
  if ( !all( unique(as.character(.run[["alphas0"]]$group0)) %in% rownames(ranef(.run$glmm)$group0))  )
    stop('Level names in pre-specified list alpha0 is not a subset of those in ranef(glmm)$group0.')

  if ( !all( unique(as.character(.run[["alphas1"]]$group1)) %in% rownames(ranef(.run$glmm)$group1))  )
    stop('Level names in pre-specified list alpha1 is not a subset of those in ranef(glmm)$group1.')
  
  # Get results.
  glmm.raneffs.group0[i, ]  <- unlist(ranef(.run$glmm)$group0)
  glmm.raneffs.group1[i, ]  <- unlist(ranef(.run$glmm)$group1)
  glmm.fixeffs[i,]          <- fixef(.run$glmm)
  glmm.p[i,]                <- coef(summary(.run$glmm))[,4]
  glm.coefs[i,]             <- coef(.run$glm)
  glm.p[i,]                 <- coef(summary(.run$glm))[,4]
  if (do.r2) {
    r.squared[i,1]          <- NagelkerkeR2(.run$glm)$R2
    r.squared[i,3:4]        <- suppressMessages(r.squaredGLMM(.run$glmm))
  }
  raneff.var[i,1]           <- as.numeric(VarCorr(.run$glmm)$group0)
  raneff.var[i,2]           <- as.numeric(VarCorr(.run$glmm)$group1)

}
# Save the alphas as actually used.
alphas0             <- .run[["alphas0"]]
alphas1             <- .run[["alphas1"]]


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

plot.fixeffs(glmm.fixeffs, c("alpha0", "beta1", "beta2"), c(alpha, beta1, beta2),
             c("darkgray", "darkgreen", "darkred"), lwd = lwd, lty = lty,
             fileprefix = fileprefix)
plot.raneff.variance(raneff.var, c("sigma0", "sigma1"), c(sigma0, sigma1),
                     c("darkgray", "darkgreen"), lwd = lwd, lty = lty,
                     fileprefix = fileprefix)
plot.raneffs(alphas0, glmm.raneffs.group0, "group0", sample.size = 8, mfrow = c(2,4), 
             lwd = lwd, lty.null = lty.null, colfunc = colfunc,
             fileprefix = fileprefix)
plot.raneffs(alphas1, glmm.raneffs.group1, "group1", sample.size = 8, mfrow = c(2,4),
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


print.raneff.variance(raneff.var, c(sigma0, sigma1))
print.fixeff.comp(glmm.fixeffs, glm.coefs)
print.fixeff.p.comp(glmm.p, glm.p)

if (do.r2) {
  plot.r2(r.squared, c("darkgray", "darkgreen"), lwd = lwd, lty = lty)
  print.r2.comp(r.squared)
}

if (!is.null(fileprefix)) sink()

cat("\n\n ### DUMP OF WARNINGS \n\n")
print(warnings())

save.image(file = "simulate_glmm_varintnested.RData")
