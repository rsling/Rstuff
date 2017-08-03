
char.seq <- function(start, end, by = 1, pad = 4, pad.char = "0") {
  formatC(seq(start, end, by), width = pad, format = "d", flag = pad.char)
}

print.fixeff.comp <- function(glmm.p, glm.p, glm.f.p = NULL) {
  cat("\n\n Comparison of p-values of 'fixed' effects\n")
  cat("\n ### Summaries\n")
  cat("\n GLMM\n")
  print(apply(glmm.p, 2, function(r) {summary(ecdf(r))}))
  cat("\n GLM (no random)\n")
  print(apply(glm.p, 2, function(r) {summary(ecdf(r))}))
  if (!is.null(glm.f.p)) {
    cat("\n GLM (random as fixed)\n")
    print(apply(glm.f.p, 2, function(r) {summary(ecdf(r))}))
  }
  cat("\n ### 95% intervals\n")
  cat("\n GLMM\n")
  print(apply(glmm.p, 2, function(r) {quantile(ecdf(r), probs = c(0.025, 0.975))}))
  cat("\n GLM (no random)\n")
  print(apply(glm.p, 2, function(r) {quantile(ecdf(r), probs = c(0.025, 0.975))}))
  if (!is.null(glm.f.p)) {
    cat("\n GLM (random as fixed)\n")
    print(apply(glm.f.p, 2, function(r) {quantile(ecdf(r), probs = c(0.025, 0.975))}))
  }
}


print.r2.comp <- function(r.squared) {
  cat("\n\n Difference between marginal and conditional R-squared in GLMM (95% interval)\n")
  print(quantile(ecdf(r.squared[,4]-r.squared[,3]), probs = c(0.025, 0.975)))
}



plot.r2 <- function(r.squared) {
  par(mfrow=c(1,1))
  plot(density(r.squared[,3]),
       xlim = c( min(c(r.squared[,3], r.squared[,4])),
                 max(c(r.squared[,3], r.squared[,4]))),
       ylim = c( min(c(density(r.squared[,3])$y, density(r.squared[,4])$y)),
                 max(c(density(r.squared[,3])$y, density(r.squared[,4])$y))),
       col = "darkorange", lwd = lwd,
       main = "Estimates of R-squared in GLMM",
       xlab = "Estimates")
  lines(density(r.squared[,4]), col = "darkgreen", lwd = lwd)
  legend("top",
         legend = c("marginal R-squared", "conditional R-squared"),
         col = c("darkorange", "darkgreen"),
         lwd = lwd)
}