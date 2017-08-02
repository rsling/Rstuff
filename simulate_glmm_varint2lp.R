require(lme4)
require(boot)
library(mvtnorm)
require(MuMIn)

rm(list = ls())
set.seed(2707)

source("simulate_glmm_varint2lp_fun.R")

use.saved  <- T
nsim       <-  10
J          <-  20
I          <-  20
beta0      <-   1
beta1      <-   0.8
alpha0     <-  -0.5
gamma_a    <-   1.2
gamma_b    <-  -0.6
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

if (use.saved) {
  load("simulate_glmm_varint2lp.RData")
} else {
  
  # Matrices for results: GLMM.
  raneffs.alpha            <- as.data.frame(matrix(rep(NA, J * nsim), nrow = nsim, byrow = T))
  raneffs.beta             <- as.data.frame(matrix(rep(NA, J * nsim), nrow = nsim, byrow = T))
  fixefs                   <- as.data.frame(matrix(rep(NA, 5 * nsim), nrow = nsim, byrow = T))
  r.squared                <- as.data.frame(matrix(rep(NA, 2 * nsim), nrow = nsim, byrow = T))
  Sigmas                   <- as.data.frame(matrix(rep(NA, 3 * nsim), nrow = nsim, byrow = T))
  colnames(raneffs.alpha)  <- paste0("group", 1:J)
  colnames(raneffs.beta)   <- paste0("group", 1:J)
  colnames(fixefs)         <- c('alpha0', 'beta0', 'beta1', "gamma_a", "gamma_b")
  colnames(r.squared)      <- c('marginal', 'conditional')
  colnames(Sigmas)         <- c('sigma_a', 'sigma_b', 'covar_ab')
  
  
  for (i in 1:nsim) {
    cat("Simulation run", i, "...\n")
    
    # In all except the first run, we re-use the alphas to get comparable results.
    if (i == 1) {
      .run <- sim.glmm.varint2lp(J = J, I = I,
                                    beta1 = beta1, alpha0 = alpha0,
                                    gamma_a = gamma_a, gamma_b = gamma_b,
                                    sigma_a = sigma_b, sigma_b = sigma_b, rho = rho,
                                    do.raneff = do.raneff, do.fixeff = do.fixeff)
    } else {
      .run <- sim.glmm.varint2lp(J = J, I = I,
                                    beta1 = beta1, alpha0 = alpha0,
                                    gamma_a = gamma_a, gamma_b = gamma_b,
                                    sigma_a = sigma_b, sigma_b = sigma_b, rho = rho,
                                    do.raneff = do.raneff, do.fixeff = do.fixeff,
                                    raneffs = .run$raneffs)
    }
    
    # Get normal GLMM results.
    raneffs.alpha[i,]         <- ranef(.run$glmm)[["group"]][,1]
    raneffs.beta[i,]          <- ranef(.run$glmm)[["group"]][,2] 
    fixefs[i,]                <- fixef(.run[["glmm"]])
    if (do.r2) r.squared[i,]  <- r.squaredGLMM(.run[["glmm"]])
    if (is.nan(as.data.frame(VarCorr(.run[["glmm"]]))[3,"sdcor"]))
      warning('Covariance is NaN!')
    else
      Sigmas[i,]                <- as.data.frame(VarCorr(.run[["glmm"]]))[,"sdcor"]
    
  }

  # Save the alphas as actually used.
  true.raneffs <- .run$raneffs
  
  # NaN were turned to NA in loop, remove.
  n.Sigmas <- nrow(Sigmas) 
  Sigmas   <- Sigmas[complete.cases(Sigmas),]
  m.Sigmas <- n.Sigmas-nrow(Sigmas)
}


cat("\nFailed runs (NaN in variances):", m.Sigmas)


par(mfrow=c(1,2))
plot(density(fixefs$alpha0),
     xlim = c(min(c(alpha0, beta0, beta1, as.matrix(fixefs)[,1:3])), max(c(alpha0, beta0, beta1, as.matrix(fixefs)[,1:3]))),
     ylim = c(min(c(density(fixefs$alpha0)$y, density(fixefs$beta0)$y, density(fixefs$beta1)$y)),
              max(c(density(fixefs$alpha0)$y, density(fixefs$beta0)$y, density(fixefs$beta1)$y)*1.2) ),
     col = "darkorange", lwd = lwd,
     main = "Estimates of fixed effects in GLMM\n(first level)",
     xlab = "Estimates")
lines(density(fixefs$beta0),
      col = "darkgreen", lwd = lwd)
lines(density(fixefs$beta1),
      col = "darkred", lwd = lwd)
abline(v = alpha0, col = "darkorange", lwd = lwd, lty = 3)
abline(v = beta0, col = "darkgreen", lwd = lwd, lty = 3)
abline(v = beta1, col = "darkred", lwd = lwd, lty = 3)
legend("top",
       legend = c("alpha0", "beta0", "beta1"),
       col = c("darkorange", "darkgreen", "darkred"),
       lwd = lwd)

plot(density(fixefs$gamma_a),
     xlim = c(min(c(gamma_a, gamma_b, as.matrix(fixefs)[,4:5])), max(c(gamma_a, gamma_b, as.matrix(fixefs)[,4:5]))),
     ylim = c(min(c(density(fixefs$gamma_a)$y, density(fixefs$gamma_b)$y)),
              max(c(density(fixefs$gamma_a)$y, density(fixefs$gamma_b)$y)*1.2)),
     col = "darkorange", lwd = lwd,
     main = "Estimates of fixed effects in GLMM\n(2nd level)",
     xlab = "Estimates")
lines(density(fixefs$gamma_b),
      col = "darkgreen", lwd = lwd)
abline(v = gamma_a, col = "darkorange", lwd = lwd, lty = 3)
abline(v = gamma_b, col = "darkgreen", lwd = lwd, lty = 3)
legend("top",
       legend = c("gamma_a", "gamma_b"),
       col = c("darkorange", "darkgreen"),
       lwd = lwd)
par(mfrow=c(1,1))



plot(density(Sigmas$sigma_a),
     xlim = c( min( c( Sigmas$sigma_a , Sigmas$sigma_b, Sigmas$covar_ab)), max( c( Sigmas$sigma_a , Sigmas$sigma_b, Sigmas$covar_ab))),
     ylim = c( min( c( density(Sigmas$sigma_a)$y, density(Sigmas$sigma_b)$y, density(Sigmas$covar_ab)$y)),
               max( c( density(Sigmas$sigma_a)$y, density(Sigmas$sigma_b)$y, density(Sigmas$covar_ab)$y))),
     col = "darkorange", lwd = lwd,
     main = "Variance estimates for random effects",
     xlab = "Estimates")
lines(density(Sigmas$sigma_b), col = "darkgreen", lwd = lwd)
lines(density(Sigmas$covar_ab), col = "darkred", lwd = lwd)
abline(v = sigma_a^2, col = "darkorange", lwd = lwd, lty = 3)
abline(v = sigma_b^2, col = "darkgreen", lwd = lwd, lty = 3)
abline(v = rho, col = "darkred", lwd = lwd, lty = 3)
legend("topleft",
       legend = c("sigma(alpha)", "sigma(beta)", "rho"),
       col = c("darkorange", "darkgreen", "darkred"),
       lwd = lwd)


par(mfrow=c(2,4))
alphas.sample.plot <- sort(sample(1:nrow(true.raneffs), size = 8, replace = F))
for (i in alphas.sample.plot) {
  plot(density(raneffs.alpha[,i]),
       xlim = c( min(0, c(true.raneffs[i,"alpha"], as.matrix(raneffs.alpha[,i]))), max(c(0, true.raneffs[i,"alpha"], as.matrix(raneffs.alpha[,i])))),
       xlab = "Predicted alpha", main = paste0("J_", i),
       lwd = lwd.small, col = colfunc(nrow(true.raneffs))[i])
  abline(v = true.raneffs[i,"alpha"], lwd = lwd.small, col = colfunc(nrow(true.raneffs))[i], lty = 3)
  abline(v = 0, lwd = lwd.null, col = "gray", lty = lty.null)
}
par(mfrow=c(1,1))


par(mfrow=c(2,4))
for (i in alphas.sample.plot) {
  plot(density(raneffs.beta[,i]),
       xlim = c( min(0, c(true.raneffs[i,"beta"], as.matrix(raneffs.beta[,i]))), max(c(0, true.raneffs[i,"beta"], as.matrix(raneffs.beta[,i])))),
       xlab = "Predicted beta", main = paste0("J_", i),
       lwd = lwd.small, col = colfunc(nrow(true.raneffs))[i])
  abline(v = true.raneffs[i,"beta"], lwd = lwd.small, col = colfunc(nrow(true.raneffs))[i], lty = 3)
  abline(v = 0, lwd = lwd.null, col = "gray", lty = lty.null)
}
par(mfrow=c(1,1))


if (do.r2) {
  par(mfrow=c(1,1))
  plot(density(r.squared[,1]),
       xlim = c( min(c(r.squared[,1], r.squared[,2])), max(c(r.squared[,1], r.squared[,2]))),
       ylim = c( min(c(density(r.squared[,1])$y, density(r.squared[,2])$y)), max(c(density(r.squared[,1])$y, density(r.squared[,2])$y))),
       col = "darkorange", lwd = lwd,
       main = "Estimates of R-squared",
       xlab = "Estimates")
  lines(density(r.squared[,2]), col = "darkgreen", lwd = lwd)
  legend("topright",
         legend = c("marginal R-squared", "conditional R-squared"),
         col = c("darkorange", "darkgreen"),
         lwd = lwd)
}

cat("\n\n Difference between marginal and conditional R-squared (95% interval)\n")
print(quantile(ecdf(r.squared[,2]-r.squared[,1]), probs = c(0.025, 0.975)))
cat("\n")

if (!use.saved) save.image(file = "simulate_glmm_varint2lp.RData")
