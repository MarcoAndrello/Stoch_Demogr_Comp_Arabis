# Build and analyse IPM

rm(list=ls())

library(popbio)
library(IPMpack)
library(tidyverse)
library(lme4)
library(gridExtra)
library(fields)
library(tictoc)
library(abind)

# Replication parameters
iparallel <- 1
set.seed(iparallel)
startclim <- 1
endclim <- 2
nboot <- 200
nrepl <- (endclim-startclim + 1) * nboot

# Various definitions and source functions
v.site <- c("BRU", "CHA", "VIL", "LAU", "GAL", "PIC")
v.year <- c(2008:2013)
source("01 - Functions_IPM.R")

# Load non-climatic environmental data
load(paste0(getwd(),"/../09 - PCA Soil and vegetation/PCA_soil_veg.RData"))
soilveg$SiteQuad <- row.names(soilveg)
soilveg <-
    as_tibble(soilveg) %>%
    add_column(Site=substr(.[["SiteQuad"]],1,3)) %>%
    add_column(Quad=paste0("Q",substr(.[["SiteQuad"]],5,5))) %>%
    # Reorder sites
    mutate(Site = factor(Site,levels=c("BRU","CHA","VIL","LAU","GAL","PIC"))) %>%
    mutate(SiteQuad=paste0(Site, "-", Quad)) %>%
    arrange(Site) %>%
    # Keep only the PCA Axes
    select(Site, Quad, SiteQuad, Axis1, Axis2)
soilveg %>% group_by(Site) %>% summarise(Axis1=mean(Axis1),Axis2=mean(Axis2)) %>% mutate(Site=as.character(Site)) -> soilveg_site

# Definition of matrix structure
min.size <- 0.5
max.size <- 50.5
n <- 50
b <- min.size + c(0:n)*(max.size-min.size)/n
y <- 0.5*(b[1:n]+b[2:(n+1)])
h <- y[2]-y[1]

# Representative stages (sizes) to keep track of
repr.stage <- c(1,5,10,15)

# Final results
r <-
Lexp.mean <-
Lexp.sd <- array(NA,c(length(v.site),nrepl))
w <- array(NA,c(length(y),
                length(v.site),
                length(v.year),
                nrepl))
ws <- array(NA,c(length(y),
                 length(v.site),
                 nrepl))
elast <- array(NA,c(length(v.site),
              16, # Elasticity to mean and variability of 8 vital rates
              nrepl))
dimnames(elast) <- list(site=v.site,
                       stat=c("E.surv.mu",
                              "E.surv.si",
                              "E.retro.growth.mu",
                              "E.retro.growth.si",
                              "E.stasis.mu",
                              "E.stasis.si",
                              "E.progr.growth.mu",
                              "E.progr.growth.si",
                              "E.fec.mu",
                              "E.fec.si",
                              "E.F0.mu",
                              "E.F0.si",
                              "E.F1.mu",
                              "E.F1.si",
                              "E.F2.mu",
                              "E.F2.si"),
                       repl=c(1:nrepl)
)

lambda.det <-
  mean.surv <-
  mean.progr <-
  mean.retro <-
  mean.stasis <-
  mean.fec <-
  mean.pfruit <-
  mean.nfruits.net <-
  mean.nfruits <-
  mean.recr.size <- array(NA,c(length(v.site),
                             length(v.year),
                             nrepl))
mat.site <- list(matrix(0,n,n), # 2008
                 matrix(0,n,n), # 2009
                 matrix(0,n,n), # 2010
                 matrix(0,n,n), # 2011
                 matrix(0,n,n), # 2012
                 matrix(0,n,n)) # 2013
mat.G <-
mat.P <-
mat.F <-
mat.K <-  list(mat.site, # BRU
              mat.site, # CHA
              mat.site, # VIL
              mat.site, # LAU
              mat.site, # GAL
              mat.site) # PIC
surv <-
retro <-
stasis <-
progr <-
pfruit <-
nfruits <- array(NA,
                     c(length(repr.stage),
                       length(v.site),
                       length(v.year),
                       nrepl)
                )

cnet <- array(NA,c(7, # Number of life-cycle components (survival, retro.growth, stasis, prog.growth, repro, fec, recr)
                   4, # Number of type of SLTRE contributions (mean, elast, cv, corr)
                   length(v.site),
                   nrepl))



i.repl <- 1
i.clim <- 1
for (i.clim in startclim : endclim) { # Loop on clim
  # Load climatic data
  load(paste0(getwd(),"/../01 - iButton data and imputation/01 - Resampled data/02 - Monthly/data.monthly.",i.clim,".RData"))
  monthly %>%
    ungroup %>%
    rename(Site=Pop) %>%
    group_by(Site,Year) %>%
    summarise(T.mean=mean(AvgMeanTemp),T.range=mean(AvgAmp)) %>%
      extend.env.dataset %>%
      left_join(soilveg_site,by="Site") %>%
      mutate(T.mean.sq = T.mean^2,
             T.mean_f.sq = T.mean_f^2,
             T.mean_m.sq = T.mean_m^2) ->
      env


  i.boot <- 1
  for (i.boot in 1 : nboot) {
    cat("Clim",i.clim,"; Boot",i.boot,"of",nboot,"\n"); flush.console()
    # Load models
    load(paste0(getwd(),"/../../noBackup/bootstraps/Models_Clim",i.clim,"_Boot",i.boot,".RData")) 
    
    matrices.G <-
    matrices.P <-
    matrices.F <-
    matrices.K <- list(mat.site,   # BRU
                         mat.site, # CHA
                         mat.site, # VIL
                         mat.site, # LAU
                         mat.site, # GAL
                         mat.site) # PIC
    
    # The following variables are needed to calculate stochastic elasticities AND SLTRE with F0, F1 and F2
    array.all <- array(NA, c(n,length(v.year)) )
    pfruit.all      <-
    log.nfruits.all <-
    nsd             <-  # Needed also to store the year-specific newb.size.distr to average life-expectancy
      list(array.all, # BRU
           array.all, # CHA
           array.all, # VIL
           array.all, # LAU
           array.all, # GAL
           array.all) # PIC

    i.site <- 1
    # Here loop on site
    for (i.site in 1 : length(v.site)) {
      site <- v.site[i.site]
      
      i.year <- 1
      for (i.year in 1 : length(v.year) ) {
        year <- v.year[i.year]
        # Generate new data
        env %>% filter(Site==v.site[i.site],Year==v.year[i.year]) -> envsy
        data.frame(Nb_Tot.x.or = y,
                   Nb_Tot.x.sq.or = y^2,
                   # logNbtot.x.or = y,
                   T.mean.or = pull(envsy,T.mean),
                   T.mean_f.or = pull(envsy,T.mean_f),
                   T.mean_m.or = pull(envsy,T.mean_m),
                   T.mean.sq.or = pull(envsy,T.mean.sq),
                   T.mean_f.sq.or = pull(envsy,T.mean_f.sq),
                   T.mean_m.sq.or = pull(envsy,T.mean_m.sq),
                   T.range.or = pull(envsy,T.range),
                   T.range_f.or = pull(envsy,T.range_f),
                   T.range_m.or = pull(envsy,T.range_m),
                   Axis1.or = pull(envsy,Axis1),
                   Axis2.or = pull(envsy,Axis2),
                   Fate = 0,
                   # logNbtot.y = 0,
                   Nb_Tot.y = 0,
                   Fru = 0,
                   logTotalSil = 0) %>%
          as_tibble ->
          newdat
        # Build kernels
        nsd[[i.site]][,i.year] <- h * newb.size.distr(newdat,site,year)
        nsd_year <- nsd[[i.site]][,i.year]
        res.build.IPM <- build.IPM(newdat,site,year,nsd[[i.site]][,i.year])
        # Matrices
        matrices.G[[i.site]][[i.year]] <- res.build.IPM$G
        matrices.P[[i.site]][[i.year]] <- res.build.IPM$P
        matrices.F[[i.site]][[i.year]] <- res.build.IPM$F
        matrices.K[[i.site]][[i.year]] <- res.build.IPM$K
        mat.G[[i.site]][[i.year]] <- mat.G[[i.site]][[i.year]] + matrices.G[[i.site]][[i.year]] # Store it to have a mean matrix to show
        mat.P[[i.site]][[i.year]] <- mat.P[[i.site]][[i.year]] + matrices.P[[i.site]][[i.year]]
        mat.F[[i.site]][[i.year]] <- mat.F[[i.site]][[i.year]] + matrices.F[[i.site]][[i.year]]
        mat.K[[i.site]][[i.year]] <- mat.K[[i.site]][[i.year]] + matrices.K[[i.site]][[i.year]]
        # Stage-specific vital rates for all stages (used for elasticities) (used also for SLTRE with F0, F1 and F2)
        pfruit.all[[i.site]][,i.year] <- res.build.IPM$pfruit
        log.nfruits.all[[i.site]][,i.year] <- res.build.IPM$log.nfruits
        # Stage-specific vital rates for 4 representative stages (stored)
        surv[,i.site,i.year,i.repl] <- colSums(res.build.IPM$P)[repr.stage]
        for (k in 1 : length(repr.stage)) {
            retro[k,i.site,i.year,i.repl] <- sum(res.build.IPM$G[(1:(repr.stage[k]-1)),repr.stage[k]])
            stasis[k,i.site,i.year,i.repl] <- res.build.IPM$G[repr.stage[k],repr.stage[k]]
            progr[k,i.site,i.year,i.repl] <- sum(res.build.IPM$G[(repr.stage[k]+1):n,repr.stage[k]])
        }
        if (any(repr.stage == 1)) retro[which(repr.stage == 1),i.site,i.year,i.repl] <- 0 # For the stage 1, no retrogression is possible, but the sum above puts a value
        pfruit[,i.site,i.year,i.repl] <- pfruit.all[[i.site]][repr.stage,i.year]
        nfruits[,i.site,i.year,i.repl] <- 10^log.nfruits.all[[i.site]][repr.stage,i.year]
        # Life-cycle components
        res.ctvr <- calc.life.cycle.comp(matrices.K[[i.site]][[i.year]],
                                         matrices.G[[i.site]][[i.year]],
                                         matrices.F[[i.site]][[i.year]],
                                         matrices.P[[i.site]][[i.year]],
                                         res.build.IPM$pfruit,
                                         res.build.IPM$log.nfruits)
        lambda.det[i.site,i.year,i.repl] <- res.ctvr$lambda.det
        w[,i.site,i.year,i.repl] <- res.ctvr$w
        mean.surv[i.site,i.year,i.repl] <- res.ctvr$mean.surv
        mean.progr[i.site,i.year,i.repl] <- res.ctvr$mean.progr
        mean.retro[i.site,i.year,i.repl] <- res.ctvr$mean.retro
        mean.stasis[i.site,i.year,i.repl] <- res.ctvr$mean.stasis
        mean.fec[i.site,i.year,i.repl] <- res.ctvr$mean.fec
        mean.pfruit[i.site,i.year,i.repl] <- res.ctvr$mean.pfruit
        mean.nfruits[i.site,i.year,i.repl] <- res.ctvr$mean.nfruits
        mean.nfruits.net[i.site,i.year,i.repl] <- res.ctvr$mean.nfruits.net
        mean.recr.size[i.site,i.year,i.repl] <- weighted.mean(c(1:n),nsd_year)  # Average of size (1:n) weighted by the size distribution of recruits (nsd_year)

      } # End loop on years
      # Stochastic growth rate
      r[i.site,i.repl] <- stochGrowthRateSampleList(nRunIn=1000, tMax=11000, listIPMmatrix=matrices.K[[i.site]])
      # Stochastic life-expectancy
      sle <- stoch.life.expectancy(matrices.P[[i.site]])
      nsd_mean <- rowMeans(nsd[[i.site]])
      Lexp.mean[i.site,i.repl] <- weighted.mean(sle$mean, nsd_mean)
      Lexp.sd[i.site,i.repl] <- weighted.mean(sle$sd, nsd_mean)
      # Stochastic elasticities and stochastic stage distribution
      res.elast <- stoch.sens.Marco(matrices.K[[i.site]],
                                    matrices.G[[i.site]],
                                    matrices.P[[i.site]],
                                    matrices.F[[i.site]],
                                    pfruit.all[[i.site]],
                                    log.nfruits.all[[i.site]],
                                    nsd[[i.site]],
                                    100)
      elast[i.site,,i.repl] <- res.elast$elast
      ws[,i.site,i.repl] <- res.elast$w


    } # End loop on sites
    # SLTRE
    cnet[,,,i.repl] <- calc.SLTRE(matrices.G, matrices.P, pfruit.all, log.nfruits.all, nsd)

    i.repl <- i.repl + 1
  } # End loop on boot
} # End loop on clim

dimnames(cnet) <- list(par=c("S","G-","G=","G+","F0","F1","F2"),
                       stat=c("mu","e","cc","rho"),
                       site = v.site,
                       repl = c(1:nrepl))

save(cnet,
     file=paste0("Results_IPM_cnet_parallel",iparallel,".RData"))


save(nrepl, y,
     lambda.det, w,
     surv, progr, retro, stasis, pfruit, nfruits,
     mean.surv, mean.progr, mean.retro, mean.stasis, mean.fec, mean.pfruit, mean.nfruits, mean.nfruits.net, mean.recr.size,
     r, ws, Lexp.mean, Lexp.sd, elast,
     cnet,
     file=paste0("Results_IPM_parallel",iparallel,".RData"))

save(mat.G, mat.P, mat.F, mat.K,
     file=paste0("Results_IPM_matrices_parallel",iparallel,".RData"))


# ## Compose results
# rm(list=ls())
# library(abind)
# load("Results_IPM_parallel1.RData")
# load("Results_IPM_matrices_parallel1.RData")
# 
# # Vital rates
# f.nfruits <- nfruits
# f.pfruit <- pfruit
# f.progr <- progr
# f.retro <- retro
# f.stasis <- stasis
# f.surv <- surv
# 
# # Life-cycle components
# f.mean.fec    <- mean.fec
# f.mean.nfruits.net <- mean.nfruits.net
# f.mean.nfruits <- mean.nfruits
# f.mean.pfruit <- mean.pfruit
# f.mean.progr <- mean.progr
# f.mean.retro <- mean.retro
# f.mean.stasis <- mean.stasis
# f.mean.recr.size <- mean.recr.size
# f.mean.surv    <- mean.surv
# 
# # Stage distribution
# f.w <- w
# f.ws <- ws
# 
# # Matrices
# f.mat.G <- mat.G
# f.mat.P <- mat.P
# f.mat.F <- mat.F
# f.mat.K <- mat.K
# 
# # Other variables: life-expectancies, growth rates, SLTRE, elasticities.
# f.Lexp.mean <- Lexp.mean
# f.Lexp.sd <- Lexp.sd
# f.r <- r
# f.cnet <- cnet
# f.elast <- elast
# f.lambda.det   <- lambda.det
# 
# for (i.parallel in 2 : 2){
#   load(paste0("Results_IPM_matrices_parallel",i.parallel,".RData"))
#   load(paste0("Results_IPM_parallel",i.parallel,".RData"))
# 
#   f.nfruits <- abind(f.nfruits, nfruits)
#   f.pfruit <- abind(f.pfruit, pfruit)
#   f.progr <- abind(f.progr, progr)
#   f.retro <- abind(f.retro, retro)
#   f.stasis <- abind(f.stasis, stasis)
#   f.surv <- abind(f.surv, surv)
# 
#   f.Lexp.mean <- cbind(f.Lexp.mean,Lexp.mean)
#   f.Lexp.sd <- cbind(f.Lexp.sd,Lexp.sd)
#   f.r <- cbind(f.r,r)
#   f.cnet <- abind(f.cnet, cnet, along=4)
#   f.elast <- abind(f.elast,elast, along=3)
#   f.lambda.det   <- abind(f.lambda.det, lambda.det)
# 
#   f.mean.nfruits.net    <- abind(f.mean.nfruits.net, mean.nfruits.net)
#   f.mean.nfruits    <- abind(f.mean.nfruits, mean.nfruits)
#   f.mean.pfruit    <- abind(f.mean.pfruit, mean.pfruit)
#   f.mean.recr.size    <- abind(f.mean.recr.size, mean.recr.size)
#   f.mean.surv    <- abind(f.mean.surv, mean.surv)
#   f.mean.progr <- abind(f.mean.progr, mean.progr)
#   f.mean.retro <- abind(f.mean.retro, mean.retro)
#   f.mean.stasis       <- abind(f.mean.stasis, mean.stasis)
#   f.mean.fec    <- abind(f.mean.fec, mean.fec)
# 
#   f.w <- abind(f.w, w)
#   f.ws <- abind(f.ws, ws)
# 
#   for (i.site in 1 : 6){
#     for (i.year in 1 : 6){
#       f.mat.G[[i.site]][[i.year]] <- f.mat.G[[i.site]][[i.year]] + mat.G[[i.site]][[i.year]]
#       f.mat.P[[i.site]][[i.year]] <- f.mat.P[[i.site]][[i.year]] + mat.P[[i.site]][[i.year]]
#       f.mat.F[[i.site]][[i.year]] <- f.mat.F[[i.site]][[i.year]] + mat.F[[i.site]][[i.year]]
#       f.mat.K[[i.site]][[i.year]] <- f.mat.K[[i.site]][[i.year]] + mat.K[[i.site]][[i.year]]
#     }
#   }
# }
# 
# f.nfruits -> nfruits
# f.pfruit -> pfruit
# f.progr -> progr
# f.retro -> retro
# f.stasis -> stasis
# f.surv -> surv
# 
# f.mean.nfruits.net   -> mean.nfruits.net
# f.mean.nfruits -> mean.nfruits
# f.mean.pfruit    -> mean.pfruit
# f.mean.recr.size -> mean.recr.size
# f.mean.surv    -> mean.surv
# f.mean.progr -> mean.progr
# f.mean.retro -> mean.retro
# f.mean.stasis       -> mean.stasis
# f.mean.fec    -> mean.fec
# 
# f.w -> w
# f.ws -> ws
# 
# f.mat.G -> mat.G
# f.mat.P -> mat.P
# f.mat.F -> mat.F
# f.mat.K -> mat.K
# 
# f.Lexp.mean -> Lexp.mean
# f.Lexp.sd -> Lexp.sd
# f.r -> r
# f.cnet -> cnet
# f.elast -> elast
# f.lambda.det   -> lambda.det
# 
# # Divide the "summed" matrices by the number of replicates (2000)
# for (i.site in 1 : 6){
#   for (i.year in 1 : 6){
#     mat.G[[i.site]][[i.year]] <- mat.G[[i.site]][[i.year]] / 2000
#     mat.P[[i.site]][[i.year]] <- mat.P[[i.site]][[i.year]] / 2000
#     mat.F[[i.site]][[i.year]] <- mat.F[[i.site]][[i.year]] / 2000
#     mat.K[[i.site]][[i.year]] <- mat.K[[i.site]][[i.year]] / 2000
#   }
# }
# 
# 
# 
# save(
#      nfruits,
#      pfruit,
#      progr,
#      retro,
#      stasis,
#      surv,
# 
#      mean.nfruits.net,
#      mean.nfruits,
#      mean.pfruit,
#      mean.recr.size,
#      mean.surv,
#      mean.progr,
#      mean.retro,
#      mean.stasis,
#      mean.fec,
# 
#      w,
#      ws,
# 
#      mat.G,
#      mat.P,
#      mat.F,
#      mat.K,
# 
#      y,
#      Lexp.mean,
#      Lexp.sd,
#      r,
#      cnet,
#      elast,
#      lambda.det,
# 
#      file="Results_IPM.RData")



