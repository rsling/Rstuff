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
  nested           = T,              # Nested or non nested data (group2 within group1).
  J1               =  20,            # Number of groups level 1.
  J2               =  20,            # Number of groups level 2 (per level group1).
  I                =  20,            # Number of obs. per group.
  beta1            =  0.8,           # Fixed effect coefficient 1 (binary factor).
  beta2            = -1.3,           # Fixed effect coefficient 2 (numeric).
  alpha0           =    0,           # Overall intercept.
  sigma1           =  0.2,           # Intercept SD level 1.
  sigma2           =  0.5,           # Intercept SD level 2.
  alphas1          = NULL,           # Specify if you want to use constant ranefs across sims. J1, sigma2 are ignored.
  alphas2          = NULL            # Specify if you want to use constant ranefs across sims. J2, sigma2 are ignored.
) {
  
  # Total number of observations.
  N =  J1 * J2 * I
  
  # Create group labels.
  groups1 <- data.frame(group1 = as.factor(sort(rep(char.seq(1, J1), I * J2))))
  if (nested) {
    groups2 <- as.factor( paste0( sort(rep(char.seq(1, J1), I * J2)) , "_" , rep(sort(rep(char.seq(1, J2), I)), J1)))
  } else {
    groups2 <- as.factor( rep(char.seq(1, J2), I))
  }
  groups2 <- data.frame(group2 = groups2)

  # Random effects vary around 0 with SD sigma2 and sigma2.
  # The number of required intercepts differs between nested and non-nested!
  if (is.null(alphas1)) {
    alphas1 <- data.frame(
      group1 = levels(groups1$group1),
      alpha1 = rnorm(length(levels(groups1$group1)), mean = 0, sd = sigma2)
    )
    alphas1 <- merge(groups1, alphas1)
  }

  if (is.null(alphas2)) {
    alphas2 <- data.frame(
      group2 = levels(groups2$group2),
      alpha2 = rnorm(length(levels(groups2$group2)), mean = 0, sd = sigma2)
    )
    alphas2 <- merge(groups2, alphas2)
  }

  # Put together the data frame with everything in it.
  observations <- data.frame(
    group1 = alphas1$group1,
    group2 = alphas2$group2,
    i      = rep(c(1:I), J1*J2),
    x1     = rbinom(N, 1, prob = 0.5),
    x2     = rnorm(N, mean = 0, sd = 1),
    alpha1 = alphas1$alpha1,
    alpha2 = alphas2$alpha2
  )

  # Generate the data using the actual model.
  observations <- within(observations,
                         y <- rbinom(N, 1, prob = inv.logit(alpha0 + alpha1 + alpha2 + beta1 * x1 + beta2 * x2))
  )
  
  # Calculate random effects model.
  raneff.glmer    <-  glmer(y ~ factor(x1) + x2 + (1 | group1) + (1 | group2),
                            data = observations,
                            family=binomial(link=logit))

  # # Calculate random effects model.
  # raneff.glmer.slash <-  glmer(y ~ factor(x1) + x2 + (1 | group1) + (1 | group2) + (1 | group1:group2),
  #                           data = observations,
  #                           family=binomial(link=logit))
  
  # Fixed effects model, ignoring raneff.
  fixeff.glm    <-  glm(y ~ factor(x1) + x2, data = observations,
                        family=binomial(link=logit))

  # Return results.
  list(
    alphas1      = alphas1,
    alphas2      = alphas2,
    observations = observations,
    glmm         = raneff.glmer,
#    glmm.slash   = raneff.glmer.slash,
    glm          = fixeff.glm
  )
}

.test.simulate.glmm.varintnested <- sim.glmm.varintnested()
