# This work is licensed under a Creative Commons Attribution 4.0 International License.
# https://creativecommons.org/licenses/by/4.0/
# Author: Roland Schäfer, Freie Universität Berlin, 2017
# roland.schaefer@fu-berlin.de

rm(list = ls())

setwd("~/Workingcopies/Rstuff/simulations/glmm")

outdir  <- "./output/"

scripts <- c("simulate_glmm_varint",
#             "simulate_glmm_varintnested",
             "simulate_glmm_varintslope",
             "simulate_glmm_varint2lp",
             "simulate_glmm_varintslope2lp")

nsim    <-  1000
Js      <- c(10, 20, 50)
Is      <- c(10, 20, 30)

for (J in Js) {
  for (I in Is) {

    # For the nested case, there are two Js.
    J0  <- J
    J1  <- J

    for (script in scripts) {
      fileprefix  <- paste0(outdir, script, "_j=", I, "_i=", J)
      cat("\nRunning: ", fileprefix, "...\n")
      source(paste0(script, ".R"))
    }
  }
}
