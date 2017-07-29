require(lme4)
require(boot)
require(MuMIn)
require(fmsb)

rm(list = ls())
set.seed(3979)

source("simulate_glmm_varintnested_fun.R")

nsim <- 20

# We use fixed random effects.
J0              <- 20
J1              <- 20
I               <- 20
thin            <-  0  # Remove on third of observations, creating unevenly represented groups. NOT WORKING!!!
beta1           <-  0.8
beta2           <- -1.1
alpha           <- -0.5
sigma0          <-  0.8
sigma1          <-  0.4

nested          <- T

do.raneff       <- T
do.raneff.slash <- F
do.fixeff       <- F
do.r2           <- T

colfunc         <- colorRampPalette(c("gold", "darkblue"))
lwd             <- 3
lwd.small       <- 3
lwd.null        <- 1.5
lty.null        <- 5

# Matrices for results: GLMM.
raneffs.group0           <- as.data.frame(matrix(rep(NA, J0 * nsim),      nrow = nsim, byrow = T))
raneffs.group1           <- as.data.frame(matrix(rep(NA, J1 * J0 * nsim), nrow = nsim, byrow = T))
fixefs                   <- as.data.frame(matrix(rep(NA, 3 * nsim),       nrow = nsim, byrow = T))
r.squared                <- as.data.frame(matrix(rep(NA, 2 * nsim),       nrow = nsim, byrow = T))
raneff.var               <- as.data.frame(matrix(rep(NA, 2 * nsim),       nrow = nsim, byrow = T))
colnames(raneffs.group0) <- paste0("group", 1:J0)
colnames(raneffs.group1) <- paste0("group", 1:(J1*J0))
colnames(fixefs)         <- c('intcpt', 'beta1', 'beta2')
colnames(r.squared)      <- c('marginal', 'conditional')
colnames(raneff.var)     <- c('sigma0', 'sigma1')


# Matrices for results: GLMM with slash.
if (do.raneff.slash) {
  slash.raneffs.group0           <- as.data.frame(matrix(rep(NA, J0 * nsim),      nrow = nsim, byrow = T))
  slash.raneffs.group1           <- as.data.frame(matrix(rep(NA, J1 * J0 * nsim), nrow = nsim, byrow = T))
  slash.fixefs                   <- as.data.frame(matrix(rep(NA, 3 * nsim),       nrow = nsim, byrow = T))
  slash.r.squared                <- as.data.frame(matrix(rep(NA, 2 * nsim),       nrow = nsim, byrow = T))
  slash.raneff.var               <- as.data.frame(matrix(rep(NA, 2 * nsim),       nrow = nsim, byrow = T))
  colnames(slash.raneffs.group0) <- paste0("group", 1:J0)
  colnames(slash.raneffs.group1) <- paste0("group", 1:(J1*J0))
  colnames(slash.fixefs)         <- c('intcpt', 'beta1', 'beta2')
  colnames(slash.r.squared)      <- c('marginal', 'conditional')
  colnames(slash.raneff.var)     <- c('sigma0', 'sigma1')
}


# RUN SIMULATIONS


for (i in 1:nsim) {
  cat("Simulation run", i, "...\n")
  
  # In all except the first run, we re-use the alphas to get comparable results.
  if (i == 1) {
    .run <- sim.glmm.varintnested(nested = nested,
                                  J0 = J0, J1 = J1, I = I, thin = thin,
                                  beta1 = beta1, beta2 = beta2,
                                  alpha = alpha,
                                  sigma0 = sigma0, sigma1 = sigma1,
                                  do.raneff = do.raneff, do.raneff.slash = do.raneff.slash, do.fixeff = do.fixeff)
  } else {
    .run <- sim.glmm.varintnested(nested = nested,
                                  J0 = J0, J1 = J1, I = I, thin = thin,
                                  beta1 = beta1, beta2 = beta2,
                                  alpha = alpha,
                                  sigma0 = sigma0, sigma1 = sigma1,
                                  alphas0 = .run[["alphas0"]], alphas1 = .run[["alphas1"]],
                                  do.raneff = do.raneff, do.raneff.slash = do.raneff.slash, do.fixeff = do.fixeff)
  }
  
  # Checker whether groups are correctly aligned in data frames.
  if ( !all( unique(as.character(.run[["alphas0"]]$group0)) %in% rownames(ranef(.run[["glmm"]])$group0))  )
    stop('Level names in pre-specified list alpha0 is not a subset of those in ranef(glmm)$group0.')

  if ( !all( unique(as.character(.run[["alphas1"]]$group1)) %in% rownames(ranef(.run[["glmm"]])$group1))  )
    stop('Level names in pre-specified list alpha1 is not a subset of those in ranef(glmm)$group1.')
  

  # Get normal GLMM results.
  raneffs.group0[i, ]       <- unlist(ranef(.run[["glmm"]])$group0)
  raneffs.group1[i, ]       <- unlist(ranef(.run[["glmm"]])$group1)
  fixefs[i,]                <- fixef(.run[["glmm"]])
  if (do.r2) r.squared[i,]  <- r.squaredGLMM(.run[["glmm"]])
  raneff.var[i,1]           <- as.numeric(VarCorr(.run[["glmm"]])$group0)
  raneff.var[i,2]           <- as.numeric(VarCorr(.run[["glmm"]])$group1)

  # Get GLMM with slash syntax results.  
  if (do.raneff.slash) {
    slash.raneffs.group0[i, ]      <- unlist(ranef(.run[["glmm.slash"]])$group0)
    slash.raneffs.group1[i, ]      <- unlist(ranef(.run[["glmm.slash"]])$group1)
    slash.fixefs[i,]               <- fixef(.run[["glmm.slash"]])
    if (do.r2) slash.r.squared[i,] <- r.squaredGLMM(.run[["glmm.slash"]])
    slash.raneff.var[i,1]          <- as.numeric(VarCorr(.run[["glmm.slash"]])$group0)
    slash.raneff.var[i,2]          <- as.numeric(VarCorr(.run[["glmm.slash"]])$group1)
  }
}

# Save the alphas as actually used.
alphas0             <- .run[["alphas0"]]
alphas1             <- .run[["alphas1"]]


# PLOTS

par(mfrow=c(1,1))
plot(density(fixefs$intcpt),
     xlim = c(min(c(alpha, beta1, beta2, as.matrix(fixefs))), max(c(alpha, beta1, beta2, as.matrix(fixefs)))),
     ylim = c(min(c(density(fixefs$intcpt)$y, density(fixefs$beta1)$y, density(fixefs$beta2)$y)),
              max(c(density(fixefs$intcpt)$y, density(fixefs$beta1)$y, density(fixefs$beta2)$y)) ),
     col = "darkorange", lwd = lwd,
     main = "Estimates of fixed effects in GLMM",
     xlab = "Estimates")
lines(density(fixefs$beta1),
      col = "darkgreen", lwd = lwd)
lines(density(fixefs$beta2),
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
     ylim = c( min(c(density(raneff.var$sigma0)$y, density(raneff.var$sigma1)$y)), max(c(density(raneff.var$sigma0)$y, density(raneff.var$sigma1)$y)) ),
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
  plot(density(raneffs.group0[,i]),
       xlim = c( min(0, c(true, as.matrix(raneffs.group0[,i]))), max(c(0, true, as.matrix(raneffs.group0[,i])))),
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
    plot(density(raneffs.group1[,i]),
         xlim = c( min(c(0, true, as.matrix(raneffs.group1[,i]))), max(c(0, true, as.matrix(raneffs.group1[,i])))),
         xlab = "Predicted alpha", main = paste0("J1_", unique(alphas1)[i,1]),
       lwd = lwd.small, col = colfunc(8)[ match(i, alphas1.sample.plot)  ])
    abline(v = true, lwd = lwd.small, col = colfunc(8)[ match(i, alphas1.sample.plot) ], lty = 3)
    abline(v = 0, lwd = lwd.null, col = "gray", lty = lty.null)
  }
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

save.image(file = "simulate_glmm_varintnested.RData")
