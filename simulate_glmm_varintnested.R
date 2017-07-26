require(lme4)
require(boot)
require(MuMIn)
require(fmsb)

rm(list = ls())
set.seed(3979)

save.to.disk <- T

J0     <-  10            #  Number of groups level 1.
J1     <-  20            #  Number of groups level 2.
I      <-  50            # Number of obs. per group.
N      <-  J0 * J1 * I   # Total number of observations.
beta1  <-  0.8           # Fixed effect coefficient 1 (binary factor).
beta2  <- -1.1           # Fixed effect coefficient 2 (numeric).

alpha  <- -0.5     # Overall intercept.

sigma0 <-  0.2     # Intercept SD level 1.
sigma1 <-  0.5     # Intercept SD level 2.


# Random effects vary around 0 with SD sigma0 and sigma1.
alphas0 <- rnorm(J0, mean = 0, sd = sigma0)
alphas1 <- rnorm(J1, mean = 0, sd = sigma1)

# Make fake obervational data.
observations <- data.frame(
  group0 = sort(rep(c(1:J0), I * J1)),
  group1 = paste0( sort(rep(c(1:J0), I * J1)) , "_" , rep(sort(rep(c(1:J1), I)), J0)),
  i      = rep(c(1:I), J0*J1),
  x1     = rbinom(N, 1, prob = 0.5),
  x2     = rnorm(N, mean = 0, sd = 1),
  alpha0 = alphas0[sort(rep(1:J0, I* J1 ))],
  alpha1 = rep(alphas1[sort(rep(1:J1, I))], J0)
)


# Generate the data using the actual model.
observations <- within(observations,
                       y <- rbinom(N, 1, prob = inv.logit(alpha + alpha0 + alpha1 + beta1 * x1 + beta2 * x2))
)

if (save.to.disk) sink("simulate_glmm_varintnested.txt")

cat("\n\n Random effects model (two nested effects).\n\n")

# Calculate random effects model.
raneff.glmer <-  glmer(y ~ factor(x1) + x2 + (1 | group0) + (1 | group1),
                       data = observations,
                       family=binomial(link=logit))
print(summary(raneff.glmer))
print(r.squaredGLMM(raneff.glmer))
plot(density(observations$y))

cat("\n\n Random effects model with explicit / notation.\n\n")

raneff.glmer.slash <-  glmer(y ~ factor(x1) + x2 + (1 | group0 / group1),
                       data = observations,
                       family=binomial(link=logit))
print(summary(raneff.glmer.slash))
print(r.squaredGLMM(raneff.glmer.slash))

cat("\n\n Fixed effects-only model.\n\n")

# Comparison to model without raneff.
fixeff.glm <-  glm(y ~ factor(x1) + x2, data = observations,
                   family=binomial(link=logit))
print(summary(fixeff.glm))
print(NagelkerkeR2(fixeff.glm))

if (save.to.disk) sink()
