# This work is licensed under a Creative Commons Attribution 4.0 International License.
# https://creativecommons.org/licenses/by/4.0/
# Author: Roland Schäfer, Freie Universität Berlin, 2017
# roland.schaefer@fu-berlin.de

require(lme4)
require(boot)
source("utils.R")


sim.glmm.varint <- function(
  J                =  50,            # Number of groups.
  I                =  50,            # Number of obs. per group.
  alpha0           =   0,            # Overall intercept.
  beta1            =  0.8,           # Fixed effect coefficient for binary predictor.
  beta2            =  1,             # Fixed effect coefficient for continuous predictor.
  sigma            =  0.5,           # Varying intercept SD.
  raneffs          = NULL            # Specify this if you want to use constant ranefs across sims.
                                     # sigma is IGNORED if it is specified.
) {
  
  # Total number of observations.
  N =  J * I
  
  # Create group labels in a full observations-size data frame.
  groups <- data.frame(group = as.factor(sort(rep(char.seq(1, J), I))))

  # Make random effects if not passed in call.
  if (is.null(raneffs)) {
    
    # Var-covar-matrix first.
    raneffs           <- data.frame( char.seq(1, J), rnorm(J, mean = 0, sd = sigma))
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
  raneff.glmer    <-  glmer(y ~ factor(x1) + x2 + (1 | group),
                            data = observations,
                            family=binomial(link=logit))

  # Fixed effects model.
  fixeff.glm    <-  glm(y ~ factor(x1) + x2, data = observations,
                        family=binomial(link=logit))

  # Fixed effects model.
  fixeff.glm.f    <-  glm(y ~ factor(x1) + x2 + group, data = observations,
                        family=binomial(link=logit))

  # Return results.
  list(
    raneffs      = raneffs,
    observations = observations,
    glmm         = raneff.glmer,
    glm          = fixeff.glm,
    glm.f        = fixeff.glm.f
  )
}

# .test.simulate.glmm.varint <- sim.glmm.varint()
