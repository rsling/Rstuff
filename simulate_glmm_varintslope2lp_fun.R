# This work is licensed under a Creative Commons Attribution 4.0 International License.
# https://creativecommons.org/licenses/by/4.0/
# Author: Roland Schäfer, Freie Universität Berlin, 2017
# roland.schaefer@fu-berlin.de

require(lme4)
require(boot)
library(mvtnorm)
source("utils.R")

sim.glmm.varintslope2lp <- function(
  J                =  50,            # Number of groups.
  I                =  50,            # Number of obs. per group.
  alpha0           = -0.5,           # Overall intercept.
  beta1            =  1,             # Overall slope offset for continuous predictor.
  beta2            =  0.8,           # Fixed effect coefficient 1 (binary factor).
  gamma_a          =  2,             # Second-level intercept.
  gamma_b          = -0.6,           # Second-level coefficient.
  sigma_a          =  0.5,           # Intercept SD.
  sigma_b          =  0.2,           # Slope SD — varying slope is for numeric regressor.
  rho              =  0.4,           # Intcpt-slope covariance.
  raneffs          = NULL            # Specify this if you want to use constant ranefs across sims.
                                     # sigma_a, sigma_b, rho, Sigma are IGNORED if it is specified.
) {

  # Total number of observations.
  N =  J * I

  # Create group labels in a full observations-size data frame.
  groups <- data.frame(group = as.factor(sort(rep(char.seq(1, J), I))))

  # Make random effects if not passed in call.
  Sigma <- NULL
  if (is.null(raneffs)) {
  
    # Var-covar-matrix first.
    Sigma             <- matrix(c(sigma_a^2, rep(sigma_a * sigma_b * rho, 2), sigma_b^2),
                                nrow = 2, byrow = TRUE)
    raneffs           <- data.frame(char.seq(1, J),
                                    rmvnorm(J, mean = c(0, 0), sigma = Sigma),
                                    rnorm(J))
    colnames(raneffs) <- c("group", "alpha", "beta", "x_gamma")

    # Calculate second-level model.
    raneffs <- within(raneffs, alpha_modelled <- alpha + gamma_a * x_gamma)
    raneffs <- within(raneffs,  beta_modelled <- beta  + gamma_b * x_gamma)
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
                         y <- rbinom(N, 1,
                                     prob = inv.logit(alpha0 + alpha_modelled + beta2 * x1 + (beta1 + beta_modelled) * x2))
  )
  
  # Calculate random effects model.
  raneff.glmer    <-  glmer(y ~ factor(x1) + x2 + x_gamma + x2 : x_gamma + (1 + x2 | group),
                            data = observations,
                            family=binomial(link=logit))

  # Fixed effects model, ignoring raneff.
  fixeff.glm    <-  glm(y ~ factor(x1) + x2 + x_gamma, data = observations,
                        family=binomial(link=logit))

  # # Fixed effects model, including raneffs as fixeffs.
  # fixeff.glm.f <-  glm(y ~ factor(x1) + x2 + x_gamma + group + x2 : x_gamma : group, data = observations,
  #                       family=binomial(link=logit))

  # Return results.
  list(
    raneffs      = raneffs,
    Sigma        = Sigma,
    observations = observations,
    glmm         = raneff.glmer,
    glm          = fixeff.glm
#    glm.f        = fixeff.glm.f
  )
}

.test.simulate.glmm.varintslope2lp <- sim.glmm.varintslope2lp()
