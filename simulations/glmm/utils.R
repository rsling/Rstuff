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



plot.r2 <- function(r.squared,
                    cols, lwd, lty,
                    fileprefix = NULL) {
  xdens <- c(density(r.squared[,3])$x, density(r.squared[,4])$x)
  ydens <- c(density(r.squared[,3])$y, density(r.squared[,4])$y)
  if (!is.null(fileprefix)) pdf(paste0(fileprefix, "_r.squared_glmm.pdf"))
  par(mfrow=c(1,1))
  plot(density(r.squared[,3]),
       xlim = c(min(xdens), max(xdens)),
       ylim = c(min(ydens), max(ydens)),
       col = cols[1], lwd = lwd, lty = lty[1],
       main = "Nakagawa & Schielzeth's R-squared (GLMM)",
       xlab = "Estimates")
  lines(density(r.squared[,4]), col = cols[2], lwd = lwd, lty = lty[2])
  legend("topright",
         legend = c("marginal R-squared", "conditional R-squared"),
         col = cols,
         lwd = lwd,
         lty = lty)
  if (!is.null(fileprefix)) dev.off()

  if (!is.na(r.squared[1,2])) {
    xdens <- c(density(r.squared[,1])$x, density(r.squared[,2])$x)
    ydens <- c(density(r.squared[,1])$y, density(r.squared[,2])$y)    
  } else {
    xdens <- density(r.squared[,1])$x
    ydens <- density(r.squared[,1])$y
  }
  if (!is.null(fileprefix)) pdf(paste0(fileprefix, "_r.squared_glm.pdf"))
  plot(density(r.squared[,1]),
       xlim = c(min(xdens), max(xdens)),
       ylim = c(min(ydens), max(ydens)),
       col = cols[1], lwd = lwd, lty = lty[1],
       main = "Nagelkerke's R-squared (GLM)",
       xlab = "Estimates")
  .legend <- c("ignore random")
  if (!is.na(r.squared[1,2])) {
    lines(density(r.squared[,2]), col = cols[2], lwd = lwd, lty = lty[2])
    .legend <- c(.legend, "random as fixed")
  }
  legend("topright",
         legend = .legend,
         col = cols,
         lwd = lwd,
         lty = lty)
  if (!is.null(fileprefix)) dev.off()
}


plot.raneffs <- function(alphas, raneffs, column.name, sample.size, mfrow,
                         lwd, lty.null, colfunc,
                         use.first.instead.of.random = T,
                         xlab = "Predicted alpha",
                         fileprefix = NULL) {

  # Check whether nos. of true and simulated raneffs are equal.
  if (!nrow(alphas) == ncol(raneffs)) stop(paste("Numbers of true and simulated random effects not equal: ", nrow(alphas), ncol(raneffs)))

  # Randomly select which effects to plot.
  if (use.first.instead.of.random)
    .selection <- 1:sample.size
  else
    .selection <- sort(sample(x = 1:nrow(alphas), replace = F, size = sample.size))

  if (!is.null(fileprefix)) pdf(paste0(fileprefix, "_raneffs.pdf"))
  par(mfrow=mfrow)

  for (i in .selection) {
    .true <- alphas[i, column.name]
    .dens <- density(raneffs[,i])
    plot(.dens,
         xlim = c( min(c(0, .true, .dens$x)),
                   max(c(0, .true, .dens$x))),
         xlab = xlab, main = paste0("J1_", colnames(raneffs)[i]),
         lwd = lwd.small, col = colfunc(8)[ match(i, .selection)  ])
    abline(v = .true, lwd = lwd.small, col = colfunc(8)[ match(i, .selection) ])
    abline(v = 0, lwd = lwd.null, col = "gray", lty = lty.null)
  }
  par(mfrow=c(1,1))

  if (!is.null(fileprefix)) dev.off()
}

plot.raneff.variance <- function(raneff.var, column_names, true_sigmas,
                                 cols, lwd, lty,
                                 fileprefix = NULL) {
  xlims <- true_sigmas
  ylims <- c()
  for (cn in column_names) {
    .dens <- density(raneff.var[[cn]])
    xlims <- c(xlims, .dens$x)
    ylims <- c(ylims, .dens$y)
  }
  if (!is.null(fileprefix)) pdf(paste0(fileprefix, "_raneff.variance.pdf"))
  par(mfrow=c(1,1))
  plot(0, type="n",
       xlim = c(min(xlims), max(xlims)*1.1),
       ylim = c(min(ylims), max(ylims)*1.2),
       main = "Variance estimates for random effects",
       xlab = "Estimates",
       ylab = "Density")

  for (i in 1:length(column_names)) {
    lines(density(raneff.var[[column_names[i]]]), col = cols[i], lwd = lwd, lty = lty[i])
    abline(v = true_sigmas[i], col = cols[i], lwd = lwd, lty = lty[i])
  }
  legend("topright",
         legend = column_names,
         col = cols,
         lwd = lwd,
         lty = lty)
  if (!is.null(fileprefix)) dev.off()
}


plot.fixeffs <- function(fixeffs, column_names, true_coefs,
                         cols, lwd, lty = lty,
                         fileprefix = NULL) {
  xlims <- true_coefs
  ylims <- c()
  for (cn in column_names) {
    .dens <- density(fixeffs[[cn]])
    xlims <- c(xlims, .dens$x)
    ylims <- c(ylims, .dens$y)
  }
  if (!is.null(fileprefix)) pdf(paste0(fileprefix, "_fixeffs.pdf"))
  par(mfrow=c(1,1))
  plot(0, type='n',
       xlim = c(min(xlims), max(xlims)*1.1),
       ylim = c(min(ylims), max(ylims)*1.2),
       main = "Estimates of fixed effects",
       xlab = "Estimates",
       ylab = "Density")
  for (i in 1:length(column_names)) {
    lines(density(fixeffs[[column_names[i]]]),
          col = cols[i], lwd = lwd, lty = lty[i])
    abline(v = true_coefs[i], col = cols[i], lwd = lwd, lty = lty[i])
  }
  legend("topright",
         legend = column_names,
         col = cols,
         lwd = lwd,
         lty = lty)
  if (!is.null(fileprefix)) dev.off()
}



plot.fixeff.comparison <- function(glmm.fixeffs, glm.coefs, glm.f.coefs = NULL,
                                   l.col, p.col, pch, main, fileprefix = NULL) {
  if (!is.null(glm.f.coefs)) {
    .effects <- intersect(intersect(colnames(glmm.fixeffs), colnames(glm.coefs)), colnames(glm.f.coefs))
    .xlim    <- cbind(as.matrix(glmm.fixeffs[, .effects]), as.matrix(glm.coefs[, .effects]), as.matrix(glm.f.coefs[, .effects]))
    .labels  <- c("GLMM", "GLM (ignore radnom)", "GLM (random as fixed)")
    .models  <- list(glmm = glmm.fixeffs, glm = glm.coefs, glm.f = glm.f.coefs)
  } else {
    .effects <- intersect(colnames(glmm.fixeffs), colnames(glm.coefs))
    .xlim    <- cbind(as.matrix(glmm.fixeffs[, .effects]), as.matrix(glm.coefs[, .effects]))
    .labels  <- c("GLMM", "GLM (ign.)")
    .models  <- list(glmm = glmm.fixeffs, glm = glm.coefs)
  }
  .step <- 100%/%(length(.effects)+1)
  .mstep <- .step%/%5
  .stops <- seq(.step, .step*length(.effects), .step)
  
  if (!is.null(fileprefix)) pdf(paste0(fileprefix, "_fixeff.comparison.pdf"))
  par(mfrow=c(1,1))
  plot(0, type = "n",
       xlim = c( quantile(ecdf(.xlim), 0.005), quantile(ecdf(.xlim), 0.99)),
       ylim = c(0,100),
       yaxt = "n", ann = F
  )
  title(main = main,
        xlab = ""
  )
  axis(side = 2, at = .stops, labels = rev(.effects), tick = t, cex.axis = 0.9)
  
  for (.m in 1:length(.models)) {
    .ypos <- -0.5*((length(.effects)-1)*.mstep)+((.m-1)*.mstep)
    for (i in 1:length(.effects)) {
      .eff  <- .effects[i]
      .ecdf <- ecdf(.models[[.m]][[.eff]])
      .y    <- rev(.stops)[i]+.ypos
#      cat(names(.models)[[.m]], .eff, .y, "\n", sep = "  |  ")
      lines(c(quantile(.ecdf, 0.025), quantile(.ecdf, 0.975)), c(.y, .y),
            lwd = 2, col = l.col[2])
      lines(c(quantile(.ecdf, 0.05), quantile(.ecdf, 0.95)), c(.y, .y),
            lwd = 5, col = l.col[1])
      points(quantile(.ecdf, 0.5), .y, pch = pch[.m], cex = 2, col = p.col[.m])
    }
  }
  legend("topright", legend = rev(.labels),
         col = rev(p.col), pch = rev(pch))
  if (!is.null(fileprefix)) dev.off()
}


dump.raneffs <- function(true.raneffs, estimate, print = T) {
  .eff.name <- colnames(true.raneffs)[2]
  .est <- t(apply(estimate, 2, function(m) { quantile(m, probs = c(0.025, 0.5, 0.975))}))
  .est <- as.data.frame(cbind(true.raneffs[,2], .est))
  colnames(.est) <- c("true", paste0(rep("predicted_", 3), colnames(.est)[2:4]))
  if (print) print(.est)
  .est
}
