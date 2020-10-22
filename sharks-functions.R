#Function for sharks trends ms
# CJ Brown 2018-03-25
# ---------------
# Functions
# ---------------

#Interpolates posterior distribution from one model fit (sampled data) to full 
#model fit then integrates to calculate KL divergence.
calc_kl <- function(xfull, pfull, xsamp, psamp){
    fx <- approxfun(xsamp, psamp, yleft = min(pfull), yright = min(pfull))
    dkl <- pfull * log(pfull/fx(xfull))
    fdkl <- approxfun(xfull, dkl)
    dkl <- integrate(fdkl, min(xfull), max(xfull), stop.on.error = FALSE)$value
    return(dkl)
    }

#Creates predictions data frame and extracts predictions from model object.
create_datpred <- function(model, dat, ipred1, modnam){
    datpred <- dat
    datpred$abundance <- model$summary.fitted.values[ipred1, 4]
    datpred$CI975 <- model$summary.fitted.values[ipred1, 5]
    datpred$CI025 <- model$summary.fitted.values[ipred1, 3]
    datpred$modnam <- modnam
    return(datpred)
}

#Core function, fits the random walk models given a dataset and prior selection.
fitmod <- function(datmod1, datextra, ipred1, modnam = "m1", prior = "pc.prec", param, midyear = c(0, 0.5)){

    #Set lincombs.
    ntimes <- length(unique(datmod1$yearnorm))
    yearsnorm <- rep(NA, ntimes)
    yearsnorm_lc1 <- yearsnorm
    yearsnorm_lc1[which(datextra$yearnorm==midyear[1])] <- -1
    yearsnorm_lc1[ntimes] <- 1
    yearsnorm_lc2 <- yearsnorm
    yearsnorm_lc2[which(datextra$yearnorm==midyear[2])] <- -1
    yearsnorm_lc2[ntimes] <- 1
    nregions <- length(na.omit(unique(datmod1$region.num)))

    if (grepl("loggamma", modnam)){
        rw1prior <- list(theta = list(prior="loggamma", param=param))
    } else {
        rw1prior <- list(theta = list(prior="pc.prec", param=param))
    }
    if (nregions >1){
        form1 <- abundance ~ 1 +
            f(yearnorm, model = "rw1", scale.model = FALSE, constr = TRUE,
            hyper = rw1prior) +
            f(region.num, model = "iid") +
            gear + offset(Leffort)
        lc1 <- inla.make.lincomb(yearnorm = yearsnorm_lc1, region.num = 0, gearnet=0)
        names(lc1) <- "popn_decline"
        lc2 <- inla.make.lincomb(yearnorm = yearsnorm_lc2, region.num = 0, gearnet=0)
        names(lc2) <- "popn_decline2"
    } else {
        form1 <- abundance ~ 1 +
            f(yearnorm, model = "rw1", scale.model = FALSE,constr = TRUE,
            hyper = rw1prior) +
            gear + offset(Leffort)
        lc1 <- inla.make.lincomb(yearnorm = yearsnorm_lc1, gearnet=0)
        names(lc1) <- "popn_decline"
        lc2 <- inla.make.lincomb(yearnorm = yearsnorm_lc2, gearnet=0)
        names(lc2) <- "popn_decline2"
    }

    m1 <- try(
        inla(form1, data = datmod1, family = "nbinomial",
         lincomb = c(lc1, lc2),
       control.predictor = list(compute = TRUE, link = 1),
    control.compute =list(waic = FALSE, cpo = TRUE),
    control.inla = list(lincomb.derived.only=TRUE),
    control.family = list(prior="pc.mgamma", param = lambda_negbin))
    , TRUE)

    #Plots for debugging
    # plot(datmod1$yearnorm[datmod1$gear == "drum"], datmod1$abundance[datmod1$gear == "drum"])
    # points(datmod1$yearnorm[datmod1$gear == "drum"], m1$summary.fitted.values$mean[datmod1$gear == "drum"], pch = 16)
    # lines(datextra$yearnorm,m1$summary.fitted$mean[ipred1]*exp(3.5))

    #check lincomb
    # m1$summary.lincomb.derived$mean
     # -m1$summary.linear.predictor$mean[ipred1[datextra$yearnorm==midyear[1]]] +
     # m1$summary.linear.predictor$mean[ipred1[length(ipred1)]]
     # -m1$summary.linear.predictor$mean[ipred1[datextra$yearnorm==midyear[2]]] +
     # m1$summary.linear.predictor$mean[ipred1[length(ipred1)]]

     if (class(m1) != 'try-error'){
        datpred1 <- create_datpred(m1, datextra, ipred1, modnam)
        # hyperpar <- inla.hyperpar(m1)
    } else {
        datpred1 <- datextra
        # hyperpar <- NULL
    }

    mout <- list(model = m1, datpred = datpred1, hyperparam = NULL)
    print(modnam)
    return(mout)
}

#Core function, fits the random walk models given a dataset and prior selection.
#fits with AR1 on region/gear
fitmod_AR1 <- function(datmod1, datextra, ipred1, modnam = "m1", prior = "pc.prec", param, midyear = c(0, 0.5)){
  
  #Set lincombs.
  ntimes <- length(unique(datmod1$yearnorm))
  yearsnorm <- rep(NA, ntimes)
  yearsnorm_lc1 <- yearsnorm
  yearsnorm_lc1[which(datextra$yearnorm==midyear[1])] <- -1
  yearsnorm_lc1[ntimes] <- 1
  yearsnorm_lc2 <- yearsnorm
  yearsnorm_lc2[which(datextra$yearnorm==midyear[2])] <- -1
  yearsnorm_lc2[ntimes] <- 1
  nregions <- length(na.omit(unique(datmod1$region.num)))
  
  ar1prior_region_1 <- ar1prior_region_1 <- list(
    theta1 = list(prior="pc.prec", param = c(1, 0.01)), 
    theta2 = list(prior = "pc.cor0", param = c(0.5, 0.5)))
  
  if (grepl("loggamma", modnam)){
    rw1prior <- list(theta = list(prior="loggamma", param=param))
  } else {
    rw1prior <- list(theta = list(prior="pc.prec", param=param))
  }
  if (nregions >1){
    form1 <- abundance ~ 1 +
      f(yearnorm, model = "rw1", scale.model = FALSE, constr = TRUE,
        hyper = rw1prior) +
      f(yearnorm2, model = "ar1", replicate = region_gear, constr = FALSE, 
        hyper = ar1prior_region_1) +
      # f(region.num, model = "iid") +
      gear + offset(Leffort)
    # lc1 <- inla.make.lincomb(yearnorm = yearsnorm_lc1, yearnorm2 = 0, region.num = 0, gearnet=0)
    lc1 <- inla.make.lincomb(yearnorm = yearsnorm_lc1, yearnorm2 = 0, gearnet=0)
    names(lc1) <- "popn_decline"
    # lc2 <- inla.make.lincomb(yearnorm = yearsnorm_lc2, yearnorm2 = 0, region.num = 0, gearnet=0)
    lc2 <- inla.make.lincomb(yearnorm = yearsnorm_lc2, yearnorm2 = 0, gearnet=0)
    names(lc2) <- "popn_decline2"
  } else {
    form1 <- abundance ~ 1 +
      f(yearnorm, model = "rw1", scale.model = FALSE,constr = TRUE,
        hyper = rw1prior) +
      f(yearnorm2, model = "ar1", replicate = region_gear, constr = FALSE, 
        hyper = ar1prior_region_1) +
      gear + offset(Leffort)
    lc1 <- inla.make.lincomb(yearnorm = yearsnorm_lc1, yearnorm2 = 0, gearnet=0)
    names(lc1) <- "popn_decline"
    lc2 <- inla.make.lincomb(yearnorm = yearsnorm_lc2, yearnorm2 = 0, gearnet=0)
    names(lc2) <- "popn_decline2"
  }
  
  m1 <- try(
    inla(form1, data = datmod1, 
         # family = "nbinomial",
         family = "poisson", 
         lincomb = c(lc1, lc2),
         control.predictor = list(compute = TRUE, link = 1),
         control.compute =list(waic = FALSE, cpo = TRUE),
         control.inla = list(lincomb.derived.only=TRUE)
         # control.family = list(prior="pc.mgamma", param = lambda_negbin)
         )
    , TRUE)
  
  # plot(datmod1$yearnorm[datmod1$gear == "drum"], datmod1$abundance[datmod1$gear == "drum"])
  # points(datmod1$yearnorm[datmod1$gear == "drum"], m1$summary.fitted.values$mean[datmod1$gear == "drum"], pch = 16)
  # lines(datextra$yearnorm,m1$summary.fitted$mean[ipred1]*exp(3.5))
  
  #check lincomb
  # m1$summary.lincomb.derived$mean
  # -m1$summary.linear.predictor$mean[ipred1[datextra$yearnorm==midyear[1]]] +
  # m1$summary.linear.predictor$mean[ipred1[length(ipred1)]]
  # -m1$summary.linear.predictor$mean[ipred1[datextra$yearnorm==midyear[2]]] +
  # m1$summary.linear.predictor$mean[ipred1[length(ipred1)]]
  
  if (class(m1) != 'try-error'){
    datpred1 <- create_datpred(m1, datextra, ipred1, modnam)
    # hyperpar <- inla.hyperpar(m1)
  } else {
    datpred1 <- datextra
    # hyperpar <- NULL
  }
  
  mout <- list(model = m1, datpred = datpred1, hyperparam = NULL)
  print(modnam)
  return(mout)
}

#Utilty function to ploy hyperparam distribution (for RW precision)
plot.hyperpar <- function(hyperpar, U, alpha){
    x <- hyperpar$marginals.hyperpar[[2]][,1]
    pr1 <- hyperpar$marginals.hyperpar[[2]][,2]
    pr2 <- pc.prior(x, U, alpha)
    ylim <- c(0, max(c(pr1, pr2)))
    xlim <- c(0, max(c(pr1[pr1>0.05], pr2[pr2>0.05])))
    plot(x, pr1, xlim = xlim, ylim = ylim,type = 'l')
    lines(x, pr2, col = "red")
}


#Fits without gear effect. Core function, fits the random walk models given a dataset and prior selection.
fitmod_drum <- function(datmod1, datextra, ipred1, modnam = "m1", prior = "pc.prec", param, midyear = c(0, 0.5)){
  
  #Set lincombs.
  ntimes <- length(unique(datmod1$yearnorm))
  yearsnorm <- rep(NA, ntimes)
  yearsnorm_lc1 <- yearsnorm
  yearsnorm_lc1[which(datextra$yearnorm==midyear[1])] <- -1
  yearsnorm_lc1[ntimes] <- 1
  yearsnorm_lc2 <- yearsnorm
  yearsnorm_lc2[which(datextra$yearnorm==midyear[2])] <- -1
  yearsnorm_lc2[ntimes] <- 1
  nregions <- length(na.omit(unique(datmod1$region.num)))
  
  if (grepl("loggamma", modnam)){
    rw1prior <- list(theta = list(prior="loggamma", param=param))
  } else {
    rw1prior <- list(theta = list(prior="pc.prec", param=param))
  }
  if (nregions >1){
    form1 <- abundance ~ 1 +
      f(yearnorm, model = "rw1", scale.model = FALSE, constr = TRUE,
        hyper = rw1prior) +
      f(region.num, model = "iid") +
      offset(Leffort)
    lc1 <- inla.make.lincomb(yearnorm = yearsnorm_lc1, region.num = 0)
    names(lc1) <- "popn_decline"
    lc2 <- inla.make.lincomb(yearnorm = yearsnorm_lc2, region.num = 0)
    names(lc2) <- "popn_decline2"
  } else {
    form1 <- abundance ~ 1 +
      f(yearnorm, model = "rw1", scale.model = FALSE, constr = TRUE,
        hyper = rw1prior) +
       offset(Leffort)
    lc1 <- inla.make.lincomb(yearnorm = yearsnorm_lc1)
    names(lc1) <- "popn_decline"
    lc2 <- inla.make.lincomb(yearnorm = yearsnorm_lc2)
    names(lc2) <- "popn_decline2"
  }
  
  m1 <- try(
    inla(form1, data = datmod1, family = "nbinomial",
         lincomb = c(lc1, lc2),
         control.predictor = list(compute = TRUE, link = 1),
         control.compute =list(waic = FALSE, cpo = TRUE),
         control.inla = list(lincomb.derived.only=TRUE),
         control.family = list(prior="pc.mgamma", param = lambda_negbin))
    , TRUE)
  
  if (class(m1) != 'try-error'){
    datpred1 <- create_datpred(m1, datextra, ipred1, modnam)
    # hyperpar <- inla.hyperpar(m1)
  } else {
    datpred1 <- datextra
    # hyperpar <- NULL
  }
  
  mout <- list(model = m1, datpred = datpred1, hyperparam = NULL)
  print(modnam)
  return(mout)
}
