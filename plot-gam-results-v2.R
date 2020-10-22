# Plot GAM to estimate shark trends 
#CJ Brown
# 2018-07-02

rm(list=ls()) # clear history

# Load packages
library(mgcv)
library(dplyr)
library(ggplot2)
library(PlotTools)
library(purrr)
library(tidyr)

load("results/model-runs-gamm-02-07-2019-allruns.rda") 

dfpreds <- expand.grid(choosek = c("free", "one-third"), sample = 1:15)
dfpreds$predchange1 <- dfpreds$predchange2 <- dfpreds$location <- dfpreds$numregs <- NA

icount <- 0
for (i in 1:length(modfits)){
  for (j in 1:2){
    icount <- icount + 1
    dfpreds$predchange1[icount] <- modfits[[i]][[j]]$pred_change[1]
    dfpreds$predchange2[icount] <- modfits[[i]][[j]]$pred_change[2]
    dfpreds$numregs[icount] <- datsamp$numregs[i]
    dfpreds$location[icount] <- as.character(datsamp$location[i])
  }
}

ggplot(dfpreds, aes(x = numregs, y = predchange2, color = location)) + 
  geom_point() + 
  facet_grid(.~choosek)
