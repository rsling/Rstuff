require(lme4)
require(boot)
library(mvtnorm)
require(MuMIn)

rm(list = ls())
#set.seed(2707)

J        <-  20      #  Number of groups.
I        <-  50      # Number of obs. per GROUP
N        <-  J * I   # Total number of observations.
sigma_a  <-  0.2     # Intercept SD.
sigma_b  <-  0.5     # Slope SD.
rho      <-  0.3     # Intcp/slope covariance.
alpha0   <- -0.5     # Common intercept.



# Random effects vary around 0 with sd q.
covar  <-  matrix(  c(sigma_a^2, rep(sigma_a * sigma_b * rho, 2), sigma_b^2),
                    nrow = 2, byrow = TRUE)

alphasbetas           <- data.frame( c(1:J), rmvnorm(J, mean = c(0, 0), sigma = covar) )
colnames(alphasbetas) <- c("group", "alpha", "beta")

observations <-  data.frame(
  group = sort(rep(c(1:J), I)),
  i     = rep(c(1:I), J),
  x     = rbinom(N, 1, prob = 0.5)
)

observations <- merge(alphasbetas, observations)  # Natural join using "group".

# Generate the data using the actual model.
observations <- within(observations,
                  y <- rbinom(N, 1, prob = inv.logit(alpha0 + alpha + beta * x))
                )


# Calculate random effects model.
raneff.glmer <-  glmer(y ~ x + (1 + x | group), data = observations,
                       family=binomial(link=logit))
print(summary(raneff.glmer))
print(r.squaredGLMM(raneff.glmer))


print(sqrt(covar))

# # Comparison to model without raneff.
# fixeff.glm <-  glm(y ~ factor(x), data = observations,
#                    family=binomial(link=logit))
# print(summary(fixeff.glm))
# 
# par(mfrow=c(2,2))
# plot(fixeff.glm)
# par(mfrow=c(1,1))
# 
# plot(raneff.glmer)


