# Look at impact of prior choice on analysis of shark trends
#CJ Brown
# 2019-06-09
# Same as v4, but using the AR1 model for residual region errors

rm(list=ls()) # clear history

# Load packages
library(INLA)
library(dplyr)
library(ggplot2)
library(PlotTools)
library(purrr)
library(tidyr)
#inla.doc("loggamma")
#inla.list.models("prior")

# ---------------
# parameters
# ---------------
# years to calfculate lincombs over
yrbase1 <- 1970 #oldest year that is in at least one group
yrbase2 <- 1984 #= 3 gens ago = 2017- 10.9*3
lambda_negbin <- 7 #prior for pc.mgamma prior for negative binomial.

# ---------------
# Data
# ---------------

#  x <- seq(0.01, 1, length.out = 100)
#  plot(1/x, pgamma(x, 1, 0.1))
# lines(1/x, pgamma(x, 1, 0.001))

source("sharks-functions.R")

load("databund.rda")
nrow(databund)
nrow(databund)
years <- seq(min(databund$year), max(databund$year), by = 1)
yearsnorm <- years - min(years)+1

databund$yearnorm <- with(databund, year - min(years)+ 1)
groups <- levels(databund$group)
ngrps <- length(groups)
databund$regsite <- paste0(databund$region, databund$site)
databund$regsite.num <- as.numeric(factor(databund$regsite))
regionord <- data.frame(region = unique(databund$region),
    region.num = c(6, 1, 4, 5, 11, 3, 10, 7, 8, 9, 2)) %>%
    arrange(region.num)
databund <- databund %>% left_join(regionord)
 dattiger <- filter(databund, group == "tiger")

#Get rid of site data.
datuse <- dattiger %>% group_by(gear, year, yearnorm, region, region.num) %>%
     summarize(abundance = sum(abundance), effort = sum(effort)) %>%
    ungroup() %>%
    mutate(yearnorm2 = yearnorm,
      Leffort = log(effort), cpue = abundance/effort, 
           region_gear = as.numeric(factor(paste(region, gear))))

#Years available at all regions
datuse %>% group_by(region, region.num) %>% summarize(miny = min(year), maxy = max(yearnorm)) %>%
    arrange(region.num)
datuse %>% group_by(region, gear) %>% summarize(n = n()) %>% data.frame()

 ggplot(datuse, aes(x = year, y = cpue, color = region)) +
 geom_line()+ facet_grid(gear~., scales = "free")

#dataframes for matching once fitting is complete
datyrs <- datuse %>% select(year, yearnorm) %>% distinct()
#Prediction data frames
datextra <- datuse %>% select(gear, yearnorm, yearnorm2, region.num, region_gear) %>%
mutate(gear = "drum", region.num = NA, region_gear = NA, yearnorm2 = NA) %>%
 distinct() %>%
    mutate(abundance = NA, Leffort = 0)
dattemp <- datuse[,names(datextra)]
datmod1 <- rbind(dattemp, datextra)
#indices for prediction dfs.
ipred1 <- (nrow(dattemp)+1) : nrow(datmod1)

# ---------------
# Specify samples
# ---------------
samples <- list(
    11, 1, 3,
    c(8,9,11), c(1:3), c(3,4,6), c(1,5,11),
    c(6:11), c(1:6), c(3:8), c(1,3,5,7,9,11),
    c(3:11), c(1:9), c(2:10), c(1:11)
    # c(6:9), c(1:4), c(3:6), c(1,4,6,9),
    # c(5:9), c(1:5), c(3:7), c(1,3,5,7,9)
)
datsamp <- data.frame(numregs = map_dbl(samples, length),
    location = c("South", "North", "Middle",
        rep(c("South", "North", "Middle", "Extremes"),2),
        "South", "North", "Middle", "All")
)
nsamps <- length(samples)

# ---------------
# Models
# ---------------
mu <- 1.27 #finite growth rate our estimate = exp(0.24)
U <- log(mu)#*sqrt(2) #= sd of intrinsic rate of increase
Ufast <- log(1.66)#*sqrt(2) #= sd of intrinsic rate of increase
Uslow <- log(1.02)#*sqrt(2) #= sd of intrinsic rate of increase
alpha <- 0.025 #set mean expected r as upper limit.

#priors
xin <- list(
    datmod1 = list(datmod1),
    datextra = list(datextra),
    ipred1 = list(ipred1),
    midyear = list(c(datyrs$yearnorm[datyrs$year == yrbase1], datyrs$yearnorm[datyrs$year == yrbase2])),
    #modnam = c(list("Life-history"), list("Slow"), list("Fast"), list("loggamma-wide"), list("loggamma-default")),
    #param = c(list(c(U, alpha)), list(c(Uslow,alpha)), list(c(Ufast,alpha)), list(c(1,0.1)), list(c(1, 1E-05)))
    modnam = c(list("Life-history"), list("Slow"), list("Fast"), list("loggamma-default")),
    param = c(list(c(U, alpha)), list(c(Uslow,alpha)), list(c(Ufast,alpha)), list(c(1, 1E-05)))
)

# save(xin, samples, datsamp, nsamps, datextra, file ="results/files-to-fit-gam.rda")

# ---------------
# Plot priors
# ---------------
#redo with inla.pc.prior code
x <- seq(0.01, 1000, length.out = 50000)
sdx <- 1/sqrt(x)
datpriors <- data.frame(sdx=rep(sdx, 3),
    pr = c(inla.pc.dprec(x, xin$param[[1]][1], xin$param[[1]][2]),
    inla.pc.dprec(x, xin$param[[2]][1], xin$param[[2]][2]),
    inla.pc.dprec(x, xin$param[[3]][1], xin$param[[3]][2])),
    modnam = c(rep(xin$modnam[[1]], 50000),
            rep(xin$modnam[[2]], 50000),
            rep(xin$modnam[[3]], 50000)
    ))

# ggplot(datpriors, aes(x = sdx, y = pr, color = modnam)) +
#     geom_line() +
#     theme_bw()

#
# Fit local and global models
#
#Takes about 25 mins for 15 models on my 2013 3 GhZ dual core macbook pro. Single region models are slowest (not sure why that is). 38 minutes for 90 models.
#6.5 minutes on my 2018 Dell with R3.6.0! 
etime <- Sys.time()
modfits <- NULL
    for (i in 1:nsamps){
        xtemp <- xin
        xtemp$datmod1[[1]] <- xin$datmod1[[1]] %>%
            filter(region.num %in% samples[[i]])
        datextra2 <- filter(datextra, yearnorm %in% unique(xtemp$datmod1[[1]]$yearnorm))
        xtemp$datmod1[[1]] <- rbind(xtemp$datmod1[[1]], datextra2)
        xtemp$ipred1[[1]] <- (nrow(xtemp$datmod1[[1]])-nrow(datextra2)+1):
            nrow(xtemp$datmod1[[1]])
        xtemp$datextra <- list(datextra2)
        mtemp <- pmap(xtemp, fitmod_AR1)
        datpredtemp <- do.call("rbind",map(mtemp, function(x) x$datpred))
        modout <- list(mout = mtemp, datpred = datpredtemp)
        modfits <- c(modfits, list(modout))
        print(paste(i,"/",nsamps))
        print(proc.time()/60)
    }
print(Sys.time()-etime)
save(list = ls(), file = "sharks-trends-inla-results/model-runs-19-06-2019-allruns.rda")
# ~500 mb

  load("/model-runs-19-06-2019-allruns.rda")

#
# Save predictions
#

datpred <- NULL
icount <- 0
for (imod in 1:nrow(datsamp)){
  icount=icount+1
  for (iprior in 1:length(xin$modnam)){

      dpred <- modfits[[imod]]$mout[[iprior]]$datpred
        n <- nrow(dpred)
        n2 <- nrow(modfits[[imod]]$mout[[iprior]]$model$summary.linear.predictor)
        dpred$Q05 <- modfits[[imod]]$mout[[iprior]]$model$summary.linear.predictor$`0.5quant`[(n2-n+1):n2]
        dpred$Q975 <- modfits[[imod]]$mout[[iprior]]$model$summary.linear.predictor$`0.975quant`[(n2-n+1):n2]
        dpred$Q025 <- modfits[[imod]]$mout[[iprior]]$model$summary.linear.predictor$`0.025quant`[(n2-n+1):n2]
        
        xrw1 <- modfits[[imod]]$mout[[iprior]]$model$marginals.hyperpar$`Precision for yearnorm`[,1]
        prw1 <- modfits[[imod]]$mout[[iprior]]$model$marginals.hyperpar$`Precision for yearnorm`[,2]
        
        # Calculate KL divergence. Hard to get right, because
        # where to set max bound depends on the length of the prior's tail. 
        # nprec <- 100000
        # prec <- seq(0.05, 100000, length.out = nprec)
        # if (iprior < 4){
        #   pr <- inla.pc.dprec(prec, xin$param[[iprior]][1],  xin$param[[iprior]][2])
        # } else{
        #   pr <- dgamma(log(prec), xin$param[[iprior]][1],  xin$param[[iprior]][2])
        #   }
        # dpred$kl_rw1 <- calc_kl(prec, pr, xrw1, prw1)
        # 
        
        dpred$prior <- xin$modnam[[iprior]]
        dpred$numregs <- datsamp$numregs[imod]
        dpred$location <- datsamp$location[[imod]]
        datpred <- c(datpred, list(dpred))
    }
}
datpred <- do.call("rbind", datpred)
save(datpred, xin, datsamp, file = "results/model-predictions-19-06-2019.rda")

# ---------------
# Data processing
# ---------------
# precision on lincomb

lincomb_metrics <- function(thismod, fullmod){
    # thismod <- modfits[[2]]$mout[[1]]$model
    # fullmod <- modfits[[15]]$mout[[1]]$model
    thismarg1 <- thismod$marginals.lincomb.derived[[1]]
    thismarg2 <- thismod$marginals.lincomb.derived[[2]]
    fullmarg1 <- fullmod$marginals.lincomb.derived[[1]]
    fullmarg2 <- fullmod$marginals.lincomb.derived[[2]]
    ipred1 <- is.na(thismod$cpo$failure)
    ipred2 <- is.na(fullmod$cpo$failure)
    kl_lc1 <- calc_kl(fullmarg1[,1], fullmarg1[,2], thismarg1[,1], thismarg1[,2])
    kl_lc2 <- calc_kl(fullmarg2[,1], fullmarg2[,2], thismarg2[,1], thismarg2[,2])

    perc_change <- -100*(1-exp(thismod$summary.lincomb.derived[,c(4,5,6)]))
    perc_change$kld <- c(kl_lc1, kl_lc2)
    perc_change$varz <- var(diff(thismod$summary.fitted.values$`0.5quant`[ipred1]))
    perc_change$rmse <- sqrt(sum((thismod$summary.fitted.values$`0.5quant`[ipred1] -
        fullmod$summary.fitted.values$`0.5quant`[ipred2])^2/sum(ipred1)))
    perc_change$p_30 <- round(pnorm(log(1-0.3),
                              thismod$summary.lincomb.derived$mean,
                              thismod$summary.lincomb.derived$sd, lower.tail = FALSE),3)
    perc_change$p_50 <- round(pnorm(log(1-0.5),
                              thismod$summary.lincomb.derived$mean,
                              thismod$summary.lincomb.derived$sd, lower.tail = FALSE), 3)
    perc_change$p_80 <- round(pnorm(log(1-0.8),
                              thismod$summary.lincomb.derived$mean,
                              thismod$summary.lincomb.derived$sd, lower.tail = FALSE), 3)
    perc_change <- tibble::rownames_to_column(perc_change, var = "lc")
    return(perc_change)
}

combine_metrics <- function(modlist, fullmod, dats, modnam){
    m <- modlist %>% transpose() %>% pluck("model")
    dout <- map(m, ~lincomb_metrics(.x, fullmod)) %>%
         do.call("rbind",.)
    dout$samprow <- rep(1:nrow(dats), each = 2)
    dats$samprow <- 1:nrow(dats)
    dout2 <- dout %>%
        gather(stat, value,-lc, -samprow) %>%
        unite(lc_stat, lc, stat, sep = ".") %>%
        spread(lc_stat, value) %>%
        inner_join(dats)
    dout2$modnam <- modnam
    return(dout2)
}

#All this shennanigans is just to extract some statse from each model
#and match it to corrrect scenario.
full_fits <- modfits %>% transpose() %>% pluck("mout") %>% transpose()
all_region_model <- full_fits[[1]][[length(samples)]]$model

full_stats <- map2(full_fits,
        unlist(xin$modnam), ~combine_metrics(.x, all_region_model, datsamp, .y)) %>%
     do.call("rbind",.)

save(full_stats, file = "results/stats_model_fits-19-06-2019.rda")
# ---------------
# Plots
# ---------------


full_stats$modnam <- factor(full_stats$modnam)
beststats <- filter(full_stats, numregs ==11 & modnam =="Life-history")

dodgewdth <- 1

dev.new()
ggplot(full_stats, aes(x = numregs, y = popn_decline2.0.5quant, color = location)) +
    facet_wrap(~modnam) +
  geom_rect(aes(xmin = 0, xmax = 11, ymin = beststats$popn_decline2.0.025quant, 
                ymax= beststats$popn_decline2.0.975quant), color = NA, 
            fill = "grey90") +
    geom_point(position=position_dodge(width=dodgewdth), size = 2) +
    geom_linerange(aes(ymin = popn_decline2.0.025quant, ymax = popn_decline2.0.975quant), position=position_dodge(width=dodgewdth)) +
    theme_classic() +
  geom_line(position=position_dodge(width=dodgewdth), alpha = 0.7)

    # geom_hline(yintercept = beststats$popn_decline2.0.025quant, linetype = 2, 
    #            color = "grey30") +
    # geom_hline(yintercept = beststats$popn_decline2.0.975quant, linetype = 2, 
    #            color = "grey30") +
    # geom_hline(yintercept = beststats$popn_decline2.0.5quant, 
    #            color = "grey30")


dev.new()
ggplot(full_stats, aes(x = numregs, y = 1-popn_decline.p_50, color = location)) +
  facet_wrap(~modnam) +
  geom_rect(aes(xmin = 0, xmax = 12, ymin = 0.75, ymax = 1), color = NA,
            fill = "grey90")+
  geom_rect(aes(xmin = 0, xmax = 12, ymin = 0.95, ymax = 1), color = NA,
            fill = "grey80")+
  geom_point(position=position_dodge(width=dodgewdth), size = 3) +
  geom_line(position=position_dodge(width=dodgewdth)) + 
  theme_bw() +
  ylab("Probability of detecting decline \n
       of 50% or greater") + 
  theme_classic()


dev.new()
ggplot(full_stats, aes(x = numregs, y = popn_decline.0.5quant, color = location)) +
    facet_wrap(~modnam) +
    geom_point(position=position_dodge(width=dodgewdth)) +
    geom_linerange(aes(ymin = popn_decline.0.025quant, ymax = popn_decline.0.975quant), position=position_dodge(width=dodgewdth)) +
    theme_bw() +
    ylim(-100, 30) +
    geom_hline(yintercept = beststats$popn_decline.0.025quant, linetype = 2) +
    geom_hline(yintercept = beststats$popn_decline.0.975quant, linetype = 2) +
    geom_hline(yintercept = beststats$popn_decline.0.5quant)


dev.new()
ggplot(full_stats, aes(x = numregs, y = popn_decline2.kld, color = location)) +
    facet_wrap(~modnam) +
    geom_point(size = 2, position=position_dodge(width=dodgewdth))+
    theme_bw()

dev.new()
ggplot(full_stats, aes(x = numregs, y = popn_decline.kld, color = location)) +
    facet_wrap(~modnam) +
    geom_point(size = 2, position=position_dodge(width=dodgewdth))+
    theme_bw()

dev.new()
ggplot(full_stats, aes(x = numregs, y = varz, color = location)) +
    facet_wrap(~modnam) +
    geom_point(size = 2, position=position_dodge(width=dodgewdth))+
    theme_bw()

dev.new()
ggplot(full_stats, aes(x = numregs, y = rmse, color = location)) +
    facet_wrap(~modnam) +
    geom_point(size = 2, position=position_dodge(width=dodgewdth))+
    theme_bw()


full_stats_arrange <- full_stats %>% arrange(popn_decline2.kld)

full_stats_arrange[1:10,]

# ---------------
# TS plots
# ---------------

#str(modfits[[1]],1)
# str(modfits[[1]]$mout[[1]],1)
xin$modnam
imod <- 15
datsamp
firstyr <- min(datuse$year)-1
dev.new()
par(mfrow = c(2,3))
for (iprior in 1:6){
    dpred <- modfits[[imod]]$mout[[iprior]]$datpred
    n <- nrow(dpred)
    n2 <- nrow(modfits[[imod]]$mout[[iprior]]$model$summary.fitted.values)
    dpred$Q05 <- modfits[[imod]]$mout[[iprior]]$model$summary.fitted.values$`0.5quant`[(n2-n+1):n2]
    print(1000*var(diff(modfits[[imod]]$mout[[iprior]]$model$summary.fitted.values$`0.5quant`[(n2-n+1):n2])))
    plot(dpred$yearnorm+firstyr, dpred$Q05, main = xin$modnam[[iprior]], type = 'l')
    abline(v = yrbase1)
    abline(v = yrbase2, lty = 2)
}
