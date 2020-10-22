# Simulations of the effect of time-series non-linearity, LH traits and prior choice on model fits
# CJ Brown 2018-03-26
#v2 has no trends, just random walks

rm(list = ls())
library(INLA)
library(dplyr)
library(ggplot2)

source("sharks-functions.R")
#
# Evaluation functions 
#

#posterior prob of new data
evalfun1 <- function(t, obs, marg){
  log(approx(marg[[t]], xout = obs[t])$y)
}


#
# Params 
#
tmax <- 30
nsims <- 20
sigma_r_true <- log(c(1.02,1.34, 1.66))
nspp <- length(sigma_r_true)
U_prior <- c(log(1.01), log(1.34), log(1.66))
npriors <- length(U_prior)
alpha <- 0.025 #set mean expected r as upper limit.

n_priors <- length(U_prior)

init <- 10
rdiff <- -0.025 #annual trend
theta_disp <- 2
pshock <- 0.1
pmagnitude <- 0.1

#
# Simulate time-series 
#

set.seed(3838)
datout <- NULL
x1 <- data.frame(t = 1:tmax, RW = NA, cpue_nb = NA, 
                 cpue_pois = NA, ispp = NA, isim = NA, cpue_nb2 = NA,
                 cpue_pois2 = NA)

for (ispp in 1:nspp){
  icount <- 1
  for (isim in 1:nsims){
    set.seed(icount+1000+1)
    errors <- rnorm(tmax, sd=1)
    icount <- icount+1
    xtemp <- x1
    xtemp$ispp <- ispp;  xtemp$isim <- isim
    xtemp$RW[1] <- init
    for (t in 1:(tmax-1)){
      xtemp$RW[t+1] <- xtemp$RW[t]*exp(errors[t] * sigma_r_true[ispp]+rdiff)
    }
    set.seed(icount + 3838)
    xtemp$cpue_nb <- rnbinom(tmax, size = theta_disp, mu = xtemp$RW)
    xtemp$cpue_nb2 <- rnbinom(tmax, size = theta_disp, mu = xtemp$RW)
    set.seed(icount + 3838)
    xtemp$cpue_pois <- rpois(tmax, lambda = xtemp$RW)
    xtemp$cpue_pois2 <- rpois(tmax, lambda = xtemp$RW)
    xtemp$cpue_shock <- xtemp$cpue_pois * 
      ifelse(runif(tmax, 0,1)<pshock, pmagnitude, 1)
    datout <- c(datout, list(xtemp))
    
  }
}

dout <- do.call("rbind",datout)
# 
# ggplot(dout, aes(x = t, y = RW, group = isim)) +
#   geom_line(color = grey(0.5, 0.5)) +
#  facet_grid(ispp~., scales = "free") +
# theme_bw()

#
# Fit models  and evaluate
#

dffits <- dout %>% select(ispp, isim) %>% distinct()
dffits <- rbind(dffits, dffits, dffits)
dffits$iprior <- rep(1:3, each = length(datout))
dffits$ifit <- rep(1:length(datout), npriors)
dffits$eval_nb <- dffits$eval_pois <- dffits$eval_shock <- dffits$eval_nb_mean <- 
  dffits$eval_pois_mean <- dffits$eval_shock_mean <- dffits$po_negbin <- dffits$po_pois <- 
 dffits$po_shock <- NA

nfits <- nrow(dffits)

#Loop over TS and fit models
etime <- Sys.time()
for (ifit in 1:nfits){
#for (ifit in 1:60){
  #select dataframe
  dattemp <- datout[[dffits$ifit[ifit]]]
  
  #specify models
  form_negbin <- cpue_nb ~ 1 +
    f(t, model = "rw1", scale.model = FALSE,constr=TRUE,
      hyper = list(theta = list(prior="pc.prec", 
                  param=c(U_prior[dffits$iprior[ifit]], alpha))))
  form_pois <- cpue_pois ~ 1 +
    f(t, model = "rw1", scale.model = FALSE,constr=TRUE,
      hyper = list(theta = list(prior="pc.prec", 
                                param=c(U_prior[dffits$iprior[ifit]], alpha))))
  form_shock <- cpue_shock ~ 1 +
    f(t, model = "rw1", scale.model = FALSE,constr=TRUE,
      hyper = list(theta = list(prior="pc.prec", 
                                param=c(U_prior[dffits$iprior[ifit]], alpha))))
  
  #Fit models to negbin and pois distributed data 
  mnegbin <- inla(form_negbin, data = dattemp, family = "nbinomial",
                   control.predictor = list(compute = TRUE, link = 1),
                  control.compute = list(return.marginals = TRUE, po= TRUE, waic = TRUE))
  
  mpois <- inla(form_pois, data = dattemp, family = "nbinomial",
             control.predictor = list(compute = TRUE, link = 1),
            control.compute = list(return.marginals = TRUE, po= TRUE, waic = TRUE))
  
  mshock <- inla(form_shock, data = dattemp, family = "nbinomial",
                 control.predictor = list(compute = TRUE, link = 1),
                 control.compute = list(return.marginals = TRUE, po= TRUE, waic = TRUE))
  
  #
  #Evaluate fits
  #
  
  #in sample: po: 
  dffits$po_negbin[ifit] <-tryCatch(sum(log(mnegbin$po$po)), 
                                    error = function(e) NA)
  dffits$po_pois[ifit] <- tryCatch(sum(log(mpois$po$po)), 
                                   error = function(e) NA)
  dffits$po_shock[ifit] <- tryCatch(sum(log(mshock$po$po)), 
                                   error = function(e) NA)
  
  #Compare mean to marginal predictions
  dffits$eval_nb_mean[ifit] <- tryCatch(
    sum(purrr::map_dbl(1:tmax, evalfun1, dattemp$RW, mnegbin$marginals.fitted.values)), 
    error = function(e) NA)
  
  dffits$eval_pois_mean[ifit] <- tryCatch(
    sum(purrr::map_dbl(1:tmax, evalfun1, dattemp$RW, mpois$marginals.fitted.values)), 
    error = function(e) NA)
  
  dffits$eval_shock_mean[ifit] <- tryCatch(
    sum(purrr::map_dbl(1:tmax, evalfun1, dattemp$RW, mshock$marginals.fitted.values)), 
    error = function(e) NA)
  
  #Likelihood of new data 
  dffits$eval_nb[ifit] <- tryCatch(sum(dnbinom(dattemp$cpue_nb2, 
              size = mnegbin$summary.hyperpar$`0.5quant`[1], 
              mu = mnegbin$summary.fitted.values$`0.5quant`, log = TRUE)), 
              error = function(e) NA)
  dffits$eval_pois[ifit] <- tryCatch(sum(dnbinom(dattemp$cpue_pois2, 
              size = mpois$summary.hyperpar$`0.5quant`[1], 
              mu = mpois$summary.fitted.values$`0.5quant`, log = TRUE)),
              error = function(e) NA)
  
  rm(mnegbin, mpois, mshock)
  print(Sys.time() -etime)
  print(ifit)

}
save(dffits, dout, file = "results/simulation-study-results-v4.rda")
#
# Plot 
#
# ggplot(dffits, aes(x = iprior, y = po_pois)) + 
  

