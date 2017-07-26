require(lme4)
require(boot)

rm(list = ls())
#set.seed(5086)

J      <-  20      #  Number of groups.
I      <-  50      # Number of obs. per GROUP
N      <-  J * I   # Total number of observations.
q      <-  0.2     # Intercept SD.
beta1  <-  0.5     # Fixed effect coefficient.
alpha0 <- -0.5     # Common intercept.


# Random effects vary around 0 with sd q.
alphas <- data.frame(
  group = c(1:N),
  alpha = rnorm(N, mean = 0, sd = q)
)


observations <-  data.frame(
  group = sort(rep(c(1:I), J)),
  i     = rep(c(1:I), J),
  x     = rbinom(N, 1, prob = 0.5)
)

observations <- merge(alphas, observations)  # natural join using "unit"
observations <- within(observations,  y <- rbinom(N, 1, prob = inv.logit(alpha0 + alpha + beta1 * x)) )


# Calculate random effects model.
raneff.glmer <-  glmer(y ~ factor(x) + (1 | group), data = observations,
                       family=binomial(link=logit))
print(summary(raneff.glmer))

# Comparison to model without raneff.
fixeff.glm <-  glm(y ~ factor(x), data = observations,
                   family=binomial(link=logit))
print(summary(fixeff.glm))

par(mfrow=c(2,2))
plot(fixeff.glm)
par(mfrow=c(1,1))

plot(raneff.glmer)


