# GAM to estimate shark trends 
#CJ Brown
# 2019-07-02

rm(list=ls()) # clear history

# Load packages
library(mgcv)
library(dplyr)
library(ggplot2)
library(PlotTools)
library(purrr)
library(tidyr)

#Wood: "In practice k-1 (or k) sets the upper limit on the degrees of freedom associated with an s smooth"
# ---------------
# parameters
# ---------------
# years to calfculate lincombs over
yrbase1 <- 1970 #oldest year that is in at least one group
yrbase2 <- 1984 #= 3 gens ago = 2017- 10.9*3
theta_fixed <- 2.5 #GAMM requires you fix theta in the negbin. 

load( "results/model-predictions.rda")
load("results/files-to-fit-gam.rda")

yearnorm_max <- max(xin$datmod1[[1]]$yearnorm)
newdat <- expand.grid(gear = "drum",
                      yearnorm = c(xin$midyear[[1]], yearnorm_max),
                      region.num = NA, Leffort = 0)



#
# Fit local and global GAMMs
#
#Takes about 30 seconds
etime <- Sys.time()
modfits <- NULL
    for (i in 1:nsamps){ #loop through each sample
        xtemp <- xin #each xtemp has one sample with 6 different priors, for GAM we can ignore priors
        # and just fit the model once. 
        xtemp$datmod1[[1]] <- xin$datmod1[[1]] %>%
            filter(region.num %in% samples[[i]])
        xtemp$datextra <- list(
              filter(datextra, yearnorm %in% unique(xtemp$datmod1[[1]]$yearnorm))
        )
        nyears_temp <- length(unique(xtemp$datmod1[[1]]$yearnorm))
        k_upper <- c(-1, # very large value so that mgcv chooses df
           round((nyears_temp/3)) + 1) # recommended value, 1/3 length of ts (see Knape 2016)
        #run mixed effects model, or fixed effects if only one region. 
        xtemp$datmod1[[1]]$yearnorm2 <- xtemp$datmod1[[1]]$yearnorm
        xtemp$datmod1[[1]]$regfact <- factor(paste(xtemp$datmod1[[1]]$region.num, 
                                                   xtemp$datmod1[[1]]$gear))
        
        modout <- NULL
        for (ik in 1:2){
            if (length(unique(xtemp$datmod1[[1]]$region.num))>1){
              
              mtemp <- gamm(abundance ~ s(yearnorm, k = k_upper[ik])+ 
                              gear + offset(Leffort), 
                            correlation = corAR1(form =~ yearnorm2 | regfact), 
                            family = negbin(theta = theta_fixed),
                            data = xtemp$datmod1[[1]])

              preds <- predict(mtemp$gam, newdata = newdat)
            } else {
              mtemp <- gam(abundance ~ s(yearnorm, k = k_upper[ik]) + 
                             gear + offset(Leffort), 
                           # correlation = corAR1(form =~ yearnorm2 | regfact), 
                           family = negbin(theta = theta_fixed),
                           data = xtemp$datmod1[[1]])
              #predict(mtemp, newdata = xtemp$datextra[[1]], type = "response")
              # preds <- predict(mtemp$gam, newdata = newdat)
              preds <- predict(mtemp, newdata = newdat)
            }
          #Calculate predicted change 
          
          pred_change <- c(-100*(1-exp(preds[3] - preds[1])),
              -100*(1-exp(preds[3] - preds[2])))
          mtemp$pred_change <- pred_change
          modout <- c(modout, list(mtemp))
        }
        
        modfits <- c(modfits, list(modout))
        print(paste(i,"/",nsamps))
        print(proc.time()/60)
    }

print(Sys.time()-etime)
save(list = ls(), file = "results/model-runs-gamm-02-07-2019-allruns.rda") 

# ----------------- #
# Process results to save a df 
# ----------------- #

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

save(dfpreds, file = "results/gam-results-v2.rda")
