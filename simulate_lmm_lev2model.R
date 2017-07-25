require(lme4)

# Inspired from (but heavily modified):
# http://anythingbutrbitrary.blogspot.de/2012/10/hierarchical-linear-models-and-lmer.html

rm(list = ls())
set.seed(6854)   # date +%s%N | md5sum | tr -dC '[^0-9]' | cut -c1-4

# Units = Groups 
#                counter:    i
#                intercept:  alpha_i
#                slope:      beta_i
#                true value: a_i
#
# Measurements  counter:     j
#
# Fixeff measurements:       x_ij.
#
# Responses:                 y_ij.

###### GROUPS (or UNITS) ######

N <- 30  #  Number of groups.
unit.df = data.frame(unit = c(1:N), a = rnorm(N))

# Unit slope and intercept. Linearly related to a_i.
unit.df <-  within(unit.df, {
  E.alpha.given.a <-  1 - 0.15 * a
  E.beta.given.a <-  3 + 0.3 * a
})


# Based on our choices of qq and ss above, slopes vary more than intercepts.
# Also we've specified a correlation of .9.9 (i.e., 0.09/√[0.04∗0.25] or r*q*s),
# so when the intercept's random effect is above its mean (of zero),
# the same is likely true for the slope. 
#

library(mvtnorm)
q = 0.2  # Variance in intercepts.
r = 0.9  # Covariance.
s = 0.5  # Variance in slopes.
cov.matrix <-  matrix(c(q^2, r * q * s, r * q * s, s^2), nrow = 2, byrow = TRUE)
random.effects <-  rmvnorm(N, mean = c(0, 0), sigma = cov.matrix)

# Add random variance (random effects) to idealised intercepts and slopes.
unit.df$alpha = unit.df$E.alpha.given.a + random.effects[, 1]
unit.df$beta = unit.df$E.beta.given.a + random.effects[, 2]



###### MEASUREMENTS (or INDIVIDUALS) ######


J = 100     # Number of obs. per GROUP
M = J * N  # Total number of observations.

# Within each group, x values are perfectly the same.
x.grid = seq(-4, 4, by = 8/J)[0:J]
within.unit.df <-  data.frame(unit = sort(rep(c(1:N), J)), j = rep(c(1:J), N), x =rep(x.grid, N))

# Create full data frame with y values (with some standard normal epsilon).
flat.df = merge(unit.df, within.unit.df)
flat.df <-  within(flat.df, y <-  alpha + x * beta + 0.75 * rnorm(n = M))

# Reduce to information we would have in actual experiment.
simple.df <-  flat.df[, c("unit", "a", "x", "y")]

# WHAT DOES THIS MEAN EXACTLY?
# For the purpose of comparison, we'll keep track of the Akaike information criterion (AIC),
# a general-purpose criterion for model comparison. Smaller is better, and there is good
# reason to believe we will find a model with a smaller AIC. For while αi and βi vary,
# they vary not about fixed means, but rather about conditional means given ai.
#
# => OK, they vary around the "idealised" group means introduced at the outset. I think
#   this is what he means.

# Calculate random effects model.
raneff.lmer <-  lmer(y ~ x + a + x : a + (1 + x | unit), data = simple.df)
print(summary(raneff.lmer))

# Compare with simple LM, ignoring random effects.
fixeff.lm <-  lm(y ~ x , data = simple.df)
print(summary(fixeff.lm))
