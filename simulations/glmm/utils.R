# This work is licensed under a Creative Commons Attribution 4.0 International License.
# https://creativecommons.org/licenses/by/4.0/
# Author: Roland Schäfer, Freie Universität Berlin, 2017
# roland.schaefer@fu-berlin.de


char.seq <- function(start, end, by = 1, pad = 4, pad.char = "0") {
  formatC(seq(start, end, by), width = pad, format = "d", flag = pad.char)
}


print.raneff.variance <- function(raneff.var, true.variance) {
  cat("\n\n Distribution of variance estimates for 'random' effects\n")
  cat("\n ### Summary\n")
  print(apply(raneff.var, 2, function(r) {summary(ecdf(r))}))
  cat("\n ### 95% intervals\n")
  print(apply(raneff.var, 2, function(r) {quantile(ecdf(r), probs = c(0.025, 0.975))}))
  cat("\n True values:\n")
  print(true.variance)
}


print.fixeff.comp <- function(glmm.fixeffs, glm.fixeffs, glm.f.fixeffs = NULL) {
  cat("\n\n Comparison of coefficient estimates of 'fixed' effects\n")
  cat("\n ### Summaries\n")
  cat("\n GLMM\n")
  print(apply(glmm.fixeffs, 2, function(r) {summary(ecdf(r))}))
  cat("\n GLM (no random)\n")
  print(apply(glm.fixeffs, 2, function(r) {summary(ecdf(r))}))
  if (!is.null(glm.f.fixeffs)) {
    cat("\n GLM (random as fixed)\n")
    print(apply(glm.f.p, 2, function(r) {summary(ecdf(r))}))
  }
  cat("\n ### 95% intervals\n")
  cat("\n GLMM\n")
  print(apply(glmm.fixeffs, 2, function(r) {quantile(ecdf(r), probs = c(0.025, 0.975))}))
  cat("\n GLM (no random)\n")
  print(apply(glm.fixeffs, 2, function(r) {quantile(ecdf(r), probs = c(0.025, 0.975))}))
  if (!is.null(glm.f.fixeffs)) {
    cat("\n GLM (random as fixed)\n")
    print(apply(glm.f.p, 2, function(r) {quantile(ecdf(r), probs = c(0.025, 0.975))}))
  }
}

print.fixeff.p.comp <- function(glmm.p, glm.p, glm.f.p = NULL) {
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



plot.r2 <- function(r.squared, cols, lwd) {
  xdens <- c(density(r.squared[,3])$x, density(r.squared[,4])$x)
  ydens <- c(density(r.squared[,3])$y, density(r.squared[,4])$y)
  par(mfrow=c(1,1))
  plot(density(r.squared[,3]),
       xlim = c(min(xdens), max(xdens)),
       ylim = c(min(ydens), max(ydens)),
       col = cols[1], lwd = lwd,
       main = "Estimates of R-squared",
       xlab = "Estimates")
  lines(density(r.squared[,4]), col = cols[2], lwd = lwd)
  legend("topright",
         legend = c("marginal R-squared", "conditional R-squared"),
         col = cols,
         lwd = lwd)
}


plot.raneffs <- function(alphas, raneffs, column_name, sample.size, mfrow, lwd) {
  par(mfrow=mfrow)
  alphas.sample.plot <- sort(sample(1:nrow(unique(alphas)), size = sample.size, replace = F))
  for (i in 1:nrow(unique(alphas))) {
    if (i %in% alphas.sample.plot) {
      true <- unique(alphas)[order(unique(alphas[[column_name]])),][i,2]
      plot(density(raneffs[,i]),
           xlim = c( min(c(0, true, as.matrix(raneffs[,i]))),
                     max(c(0, true, as.matrix(raneffs[,i])))),
           xlab = "Predicted alpha", main = paste0("J1_", unique(alphas)[i,1]),
           lwd = lwd.small, col = colfunc(8)[ match(i, alphas.sample.plot)  ])
      abline(v = true, lwd = lwd.small, col = colfunc(8)[ match(i, alphas.sample.plot) ], lty = 3)
      abline(v = 0, lwd = lwd.null, col = "gray", lty = lty.null)
    }
  }
  par(mfrow=c(1,1))
}


plot.raneff.variance <- function(raneff.var, column_names, true_sigmas, cols, lwd) {
  xlims <- true_sigmas
  ylims <- c()
  for (cn in column_names) {
    .dens <- density(raneff.var[[cn]])
    xlims <- c(xlims, .dens$x)
    ylims <- c(ylims, .dens$y)
  }
  par(mfrow=c(1,1))
  plot(0, type="n",
       xlim = c(min(xlims), max(xlims)*1.1),
       ylim = c(min(ylims), max(ylims)*1.2),
       main = "Variance estimates for random effects",
       xlab = "Estimates",
       ylab = "Density")

  for (i in 1:length(column_names)) {
    lines(density(raneff.var[[column_names[i]]]), col = cols[i], lwd = lwd)
    abline(v = true_sigmas[i], col = cols[i], lwd = lwd, lty = 3)
  }
  legend("topright",
         legend = column_names,
         col = cols,
         lwd = lwd)
}


plot.fixeffs <- function(fixeffs, column_names, true_coefs, cols, lwd) {
  xlims <- true_coefs
  ylims <- c()
  for (cn in column_names) {
    .dens <- density(fixeffs[[cn]])
    xlims <- c(xlims, .dens$x)
    ylims <- c(ylims, .dens$y)
  }
  par(mfrow=c(1,1))
  plot(0, type='n',
       xlim = c(min(xlims), max(xlims)*1.1),
       ylim = c(min(ylims), max(ylims)*1.2),
       main = "Estimates of fixed effects",
       xlab = "Estimates",
       ylab = "Density")
  for (i in 1:length(column_names)) {
    lines(density(fixeffs[[column_names[i]]]),
          col = cols[i], lwd = lwd)
    abline(v = true_coefs[i], col = cols[i], lwd = lwd, lty = 3)
  }
  legend("topright",
         legend = column_names,
         col = cols,
         lwd = lwd)
}
