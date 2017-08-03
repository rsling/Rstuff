# This work is licensed under a Creative Commons Attribution 4.0 International License.
# https://creativecommons.org/licenses/by/4.0/
# Author: Roland Schäfer, Freie Universität Berlin, 2017
# roland.schaefer@fu-berlin.de

require(lme4)
require(boot)
require(MuMIn)
require(fmsb)

rm(list = ls())
set.seed(1255)

source("simulate_glmm_varintnested_fun.R")

use.saved       <- F
nsim            <- 10
J0              <- 10
J1              <- 10
I               <- 10
beta1           <-  0.8
beta2           <- -1.1
alpha           <- -0.5
sigma0          <-  0.8
sigma1          <-  0.4
nested          <- T    # This creates nested data, nothing in GLMM spec changes.
do.r2           <- T

colfunc         <- colorRampPalette(c("gold", "darkblue"))
lwd             <- 3
lwd.small       <- 3
lwd.null        <- 1.5
lty.null        <- 5


if (use.saved) {
  load("simulate_glmm_varintnested.RData")
} else {
  
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
      r.squared[i,3:4]        <- suppressMessages(r.squaredGLMM(.run$glmm))
    }
    raneff.var[i,1]           <- as.numeric(VarCorr(.run$glmm)$group0)
    raneff.var[i,2]           <- as.numeric(VarCorr(.run$glmm)$group1)
  
  }
  # Save the alphas as actually used.
  alphas0             <- .run[["alphas0"]]
  alphas1             <- .run[["alphas1"]]
}


# PLOTS

par(mfrow=c(1,1))
plot(density(glmm.fixeffs$alpha0),
     xlim = c(min(c(alpha, beta1, beta2, as.matrix(glmm.fixeffs))),
              max(c(alpha, beta1, beta2, as.matrix(glmm.fixeffs)))),
     ylim = c(min(c(density(glmm.fixeffs$alpha0)$y, density(glmm.fixeffs$beta1)$y, density(glmm.fixeffs$beta2)$y)),
              max(c(density(glmm.fixeffs$alpha0)$y, density(glmm.fixeffs$beta1)$y, density(glmm.fixeffs$beta2)$y)) ),
     col = "darkorange", lwd = lwd,
     main = "Estimates of fixed effects in GLMM",
     xlab = "Estimates")
lines(density(glmm.fixeffs$beta1),
      col = "darkgreen", lwd = lwd)
lines(density(glmm.fixeffs$beta2),
      col = "darkred", lwd = lwd)
abline(v = alpha, col = "darkorange", lwd = lwd, lty = 3)
abline(v = beta1, col = "darkgreen", lwd = lwd, lty = 3)
abline(v = beta2, col = "darkred", lwd = lwd, lty = 3)
legend("top",
       legend = c("beta2", "alpha", "beta1"),
       col = c("darkred", "darkorange", "darkgreen"),
       lwd = lwd)



plot(density(raneff.var$sigma0),
     xlim = c(min(as.matrix(raneff.var))*0.75, max(as.matrix(raneff.var))*1.5),
     ylim = c( min(c(density(raneff.var$sigma0)$y, density(raneff.var$sigma1)$y)),
               max(c(density(raneff.var$sigma0)$y, density(raneff.var$sigma1)$y)) ),
     col = "darkorange", lwd = lwd,
     main = "Variance estimates for random effects",
     xlab = "Estimates")
lines(density(raneff.var$sigma1), col = "darkgreen", lwd = lwd)
abline(v = sigma0^2, col = "darkorange", lwd = lwd, lty = 3)
abline(v = sigma1^2, col = "darkgreen", lwd = lwd, lty = 3)
legend("topright",
       legend = c("sigma0 (outer groups)", "sigma1 (inner groups)"),
       col = c("darkorange", "darkgreen"),
       lwd = lwd)




par(mfrow=c(2,4))
alphas0.sample.plot <- sort(sample(1:nrow(unique(alphas0)), size = 8, replace = F))
for (i in alphas0.sample.plot) {
  true <- unique(alphas0)[order(unique(alphas0$group0)),][i,2]
  plot(density(glmm.raneffs.group0[,i]),
       xlim = c( min(0, c(true, as.matrix(glmm.raneffs.group0[,i]))),
                 max(c(0, true, as.matrix(glmm.raneffs.group0[,i])))),
       xlab = "Predicted alpha", main = paste0("J0_", i),
       lwd = lwd.small, col = colfunc(nrow(unique(alphas0)))[i])
  abline(v = true, lwd = lwd.small, col = colfunc(nrow(unique(alphas0)))[i], lty = 3)
  abline(v = 0, lwd = lwd.null, col = "gray", lty = lty.null)
}
par(mfrow=c(1,1))



par(mfrow=c(2,4))
alphas1.sample.plot <- sort(sample(1:nrow(unique(alphas1)), size = 8, replace = F))
for (i in 1:nrow(unique(alphas1))) {
  if (i %in% alphas1.sample.plot) {
    true <- unique(alphas1)[order(unique(alphas1$group1)),][i,2]
    plot(density(glmm.raneffs.group1[,i]),
         xlim = c( min(c(0, true, as.matrix(glmm.raneffs.group1[,i]))),
                   max(c(0, true, as.matrix(glmm.raneffs.group1[,i])))),
         xlab = "Predicted alpha", main = paste0("J1_", unique(alphas1)[i,1]),
       lwd = lwd.small, col = colfunc(8)[ match(i, alphas1.sample.plot)  ])
    abline(v = true, lwd = lwd.small, col = colfunc(8)[ match(i, alphas1.sample.plot) ], lty = 3)
    abline(v = 0, lwd = lwd.null, col = "gray", lty = lty.null)
  }
}
par(mfrow=c(1,1))


# Plot R2 comparisons.
if (do.r2) plot.r2(r.squared)

# Print R2 comparison.
print.r2.comp(r.squared)

# Print comaprison between R2 measures.
print.fixeff.comp(glmm.p, glm.p)

if (!use.saved) save.image(file = "simulate_glmm_varintnested.RData")
