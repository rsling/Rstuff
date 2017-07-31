char.seq <- function(start, end, by = 1, pad = 4, pad.char = "0") {
  formatC(seq(start, end, by), width = pad, format = "d", flag = pad.char)
}

sim.glmm.varint <- function(
  J                =  50,            # Number of groups.
  I                =  50,            # Number of obs. per group.
  alpha0           = -0.5,           # Overall intercept.
  beta1            =  0.8,           # Fixed effect coefficient for binary predictor.
  beta2            =  1,             # Fixed effect coefficient for continuous predictor.
  
  sigma_a          =  0.5,           # Varying intercept SD.

  raneffs          = NULL,           # Specify this if you want to use constant ranefs across sims.
                                     # sigma_a, sigma_b, rho, Sigma are IGNORED if it is specified.

  do.raneff        = T,              # Whether random effects model should be run.
  do.fixeff        = F               # Whether fixed effects model should be run, ignoring random effect structure.
) {
  
  # Total number of observations.
  N =  J * I
  
  # Create group labels in a full observations-size data frame.
  groups <- data.frame(group = as.factor(sort(rep(char.seq(1, J), I))))

  # Make random effects if not passed in call.
  if (is.null(raneffs)) {
    
    # Var-covar-matrix first.
    raneffs           <- data.frame( char.seq(1, J), rnorm(J, mean = 0, sd = sigma_a))
    colnames(raneffs) <- c("group", "alpha")
  }

  observations <- merge(groups, raneffs)
  
  # Put together the data frame with everything in it.
  observations <- cbind(observations,
                        data.frame(
                          i      = rep(c(1:I), J),
                          x1     = rbinom(N, 1, prob = 0.5),
                          x2     = rnorm(N, mean = 0, sd = 1)
                          )
  )
  
  # Generate the data using the actual model.
  observations <- within(observations,
                         y <- rbinom(N, 1, prob = inv.logit(alpha0 + alpha + beta1 * x1 + beta2 * x2))
  )
  
  # Calculate random effects model.
  raneff.glmer <- NULL
  if (do.raneff) {
    raneff.glmer    <-  glmer(y ~ factor(x1) + x2 + (1 | group),
                              data = observations,
                              family=binomial(link=logit))
  }
  
  # Fixed effects model.
  fixeff.glm <- NULL
  if (do.fixeff) {
    fixeff.glm    <-  glm(y ~ factor(x1) + x2, data = observations,
                          family=binomial(link=logit))
  }

  # Return results.
  list(
    raneffs      = raneffs,
    observations = observations,
    glmm         = raneff.glmer,
    glm          = fixeff.glm
  )
}

.test.simulate.glmm.varint <- sim.glmm.varint()
