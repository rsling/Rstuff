# This work is licensed under a Creative Commons Attribution 4.0 International License.
# https://creativecommons.org/licenses/by/4.0/
# Author: Roland Schäfer, Freie Universität Berlin, 2017
# roland.schaefer@fu-berlin.de

# NOTE:
# Fixed effects model including raneffs as fixeffs was DISABLED
# because it won't converge, and R session hangs.

require(lme4)
require(boot)
source("utils.R")

sim.glmm.varintnested <- function(
  nested           = T,              # Nested or non nested data (group1 within group0).
  J0               =  10,            # Number of groups level 1.
  J1               =  20,            # Number of groups level 2 (per level group0).
  I                =  50,            # Number of obs. per group.
  beta1            =  0.8,           # Fixed effect coefficient 1 (binary factor).
  beta2            = -1.1,           # Fixed effect coefficient 2 (numeric).
  alpha            = -0.5,           # Overall intercept.
  sigma0           =  0.2,           # Intercept SD level 1.
  sigma1           =  0.5,           # Intercept SD level 2.
  alphas0          = NULL,           # Specify if you want to use constant ranefs across sims. J0, sigma0 are ignored.
  alphas1          = NULL            # Specify if you want to use constant ranefs across sims. J1, sigma1 are ignored.
) {
  
  # Total number of observations.
  N =  J0 * J1 * I
  
  # Create group labels.
  groups0 <- data.frame(group0 = as.factor(sort(rep(char.seq(1, J0), I * J1))))
  if (nested) {
    groups1 <- as.factor( paste0( sort(rep(char.seq(1, J0), I * J1)) , "_" , rep(sort(rep(char.seq(1, J1), I)), J0)))
  } else {
    groups1 <- as.factor( rep(char.seq(1, J1), I))
  }
  groups1 <- data.frame(group1 = groups1)

  # Random effects vary around 0 with SD sigma0 and sigma1.
  # The number of required intercepts differs between nested and non-nested!
  if (is.null(alphas0)) {
    alphas0 <- data.frame(
      group0 = levels(groups0$group0),
      alpha0 = rnorm(length(levels(groups0$group0)), mean = 0, sd = sigma0)
    )
    alphas0 <- merge(groups0, alphas0)
  }

  if (is.null(alphas1)) {
    alphas1 <- data.frame(
      group1 = levels(groups1$group1),
      alpha1 = rnorm(length(levels(groups1$group1)), mean = 0, sd = sigma1)
    )
    alphas1 <- merge(groups1, alphas1)
  }

  # Put together the data frame with everything in it.
  observations <- data.frame(
    group0 = alphas0$group0,
    group1 = alphas1$group1,
    i      = rep(c(1:I), J0*J1),
    x1     = rbinom(N, 1, prob = 0.5),
    x2     = rnorm(N, mean = 0, sd = 1),
    alpha0 = alphas0$alpha0,
    alpha1 = alphas1$alpha1
  )

  # Generate the data using the actual model.
  observations <- within(observations,
                         y <- rbinom(N, 1, prob = inv.logit(alpha + alpha0 + alpha1 + beta1 * x1 + beta2 * x2))
  )
  
  # Calculate random effects model.
  raneff.glmer    <-  glmer(y ~ factor(x1) + x2 + (1 | group0) + (1 | group1),
                            data = observations,
                            family=binomial(link=logit))

  # Fixed effects model, ignoring raneff.
  fixeff.glm    <-  glm(y ~ factor(x1) + x2, data = observations,
                        family=binomial(link=logit))

  # Return results.
  list(
    alphas0      = alphas0,
    alphas1      = alphas1,
    observations = observations,
    glmm         = raneff.glmer,
    glm          = fixeff.glm
  )
}

# .test.simulate.glmm.varintnested <- sim.glmm.varintnested()