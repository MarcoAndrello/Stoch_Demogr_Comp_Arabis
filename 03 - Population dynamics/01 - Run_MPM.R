# Build and analyse MPM
# Marco Andrello
# Last modification: 07/01/2020

rm(list=ls())

library(popbio)
library(IPMpack)
library(tidyverse)
library(lme4)
library(gridExtra)
library(fields)
library(abind)

# Replication parameters
iparallel <- 2
set.seed(iparallel)
startclim <- 6
endclim <- 10
nboot <- 200
nrepl <- (endclim-startclim + 1) * nboot

# Various definitions and source functions
v.site <- c("BRU", "CHA", "VIL", "LAU", "GAL", "PIC")
v.year <- c(2008:2013)
source("00.1 - Functions_MPM.R")
source("00.2 - Functions_SLTRE.R")
source("00.3 - Functions_Elasticity_SLTRE.R")

# Load non-climatic environmental data
load(paste0(getwd(),"/../01 - Environmental data/02 - SoilVeg/PCA_soil_veg.RData")) 
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

####################################################################################
# Set up objects to store results
####################################################################################

# Deterministic lambda
lambda.det <- array(NA,c(length(v.site),
  length(v.year),
  nrepl))

# Life-cycle components
mean.surv <-
    mean.retro <-
    mean.stasis <-
    mean.progr <-
    mean.pfruit <-
    mean.nfruits <-
    mean.recr.size <- array(NA,c(length(v.site),
                                 length(v.year),
                                 nrepl))
dimnames(mean.surv) <- 
    dimnames(mean.retro) <- 
    dimnames(mean.stasis) <- 
    dimnames(mean.progr) <- 
    dimnames(mean.pfruit) <- 
    dimnames(mean.nfruits) <- 
    dimnames(mean.recr.size) <-  list(site = v.site,
                                      year = v.year,
                                      repl = c(1:nrepl))
# Stochastic growth rate r = log-lambda (should be a = log-lambda)
r <- array(NA,c(length(v.site),
                nrepl))
dimnames(r) <-  list(site = v.site,
                     repl = c(1:nrepl))
# # Stable stage distribution as average of the iterated population vector in a stochastic framework
ws <- array(NA,c(length(y),
                 length(v.site),
                 nrepl))
dimnames(ws) <- list(stage = c(1:50),
                     site = v.site,
                     repl = c(1:nrepl))
# Elasticity to mean and standard deviation of life-cycle components
elast <- array(NA,c(length(v.site),    
                    14, 
                    nrepl))
dimnames(elast) <- list(site = v.site,
                        stat = c("E.surv.mu",
                                 "E.surv.si",
                                 "E.retro.growth.mu",
                                 "E.retro.growth.si",
                                 "E.stasis.mu",
                                 "E.stasis.si",
                                 "E.progr.growth.mu",
                                 "E.progr.growth.si",
                                 "E.F0.mu",
                                 "E.F0.si",
                                 "E.F1.mu",
                                 "E.F1.si",
                                 "E.F2.mu",
                                 "E.F2.si"),
                        repl = c(1:nrepl)
)
# Object containing the six projection matrices (one for each year) for a single site
mat.site <- list(matrix(0,n,n), # 2008
                 matrix(0,n,n), # 2009
                 matrix(0,n,n), # 2010
                 matrix(0,n,n), # 2011
                 matrix(0,n,n), # 2012
                 matrix(0,n,n)) # 2013
names(mat.site) <- v.year
# Objects containing the different projection matrices (G, P, F and K) for each site and year
mat.G <-
mat.P <-
mat.F <-
mat.K <-  list(mat.site, # BRU
                   mat.site,  # CHA
                   mat.site,  # VIL
                   mat.site,  # LAU
                   mat.site,  # GAL
                   mat.site)  # PIC
names(mat.G) <- 
names(mat.P) <-
names(mat.F) <-
names(mat.K) <- v.site
# The seven life-cycle components but not averaged over all stages but taken for a set of representative stages
repr.stage <- c(1,5,10,15) # # Representative stages (sizes) to keep track of
surv <-
    retro <-
    stasis <-
    progr <-
    pfruit <-
    nfruits <- array(NA,c(length(repr.stage),
                          length(v.site),
                          length(v.year),
                          nrepl)
    )
dimnames(surv) <-
    dimnames(retro) <-
    dimnames(stasis) <-
    dimnames(progr) <-
    dimnames(pfruit) <-
    dimnames(nfruits) <- list(stage = repr.stage,
                              site = v.site,
                              year = v.year,
                              repl = c(1:nrepl))
# SLTRE: net contributions of each demographic paramters (descriptor of a life-cycle component) to differences in log-lambdas between sites
cnet <- array(NA,c(7, # Number of life-cycle components (survival, retro.growth, stasis, prog.growth, repro, fec, recr)
                   4, # Descriptor (mean, elast, cv, corr)
                   length(v.site),
                   nrepl))
dimnames(cnet) <- list(par = c("S","G-","G=","G+","F0","F1","F2"),
                       descr = c("mu","e","cc","rho"),
                       site = v.site,
                       repl = c(1:nrepl))
log.lambda.R <- array(NA,nrepl)

####################################################################################
####################################################################################


i.repl <- 1
i.clim <- 1
for (i.clim in startclim : endclim) { # Loop on clim
  # Load climatic data
  load(paste0(getwd(),"/../01 - Environmental data/01 - Climatic/01 - Resampled data/02 - Monthly/data.monthly.",i.clim,".RData"))
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
    load(paste0(getwd(),"/../ResultsGLMM/Models_Clim",i.clim,"_Boot",i.boot,".RData")) 
    
    matrices.G <-
    matrices.P <-
    matrices.F <-
    matrices.K <- list(mat.site,   # BRU
                         mat.site, # CHA
                         mat.site, # VIL
                         mat.site, # LAU
                         mat.site, # GAL
                         mat.site) # PIC
    names(matrices.G) <- 
    names(matrices.P) <-
    names(matrices.F) <-
    names(matrices.K) <- v.site
    
    # The following variables are needed to calculate stochastic elasticities AND SLTRE with F0, F1 and F2
    array.all <- array(NA, c(n,length(v.year)) )
    dimnames(array.all) <- list(stage = c(1:n),
                                year = v.year)
    pfruit.all      <-
    log.nfruits.all <-
    nsd             <-  # Needed also to store the year-specific newb.size.distr to average life-expectancy. Also for the SLTRE?
      list(array.all, # BRU
           array.all, # CHA
           array.all, # VIL
           array.all, # LAU
           array.all, # GAL
           array.all) # PIC
    names(pfruit.all)      <-
    names(log.nfruits.all) <- 
    names(nsd)             <- 
        v.site

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
                   Nb_Tot.y = 0,
                   Fru = 0,
                   logTotalSil = 0) %>%
          as_tibble ->
          newdat
        # Calculate the nuwborn size distribution
        nsd[[i.site]][,i.year] <- h * newb.size.distr(newdat,site,year)
        # nsd_year <- nsd[[i.site]][,i.year]
        # Build matrices and calculate life-cycle components
        res.build.matrices <- build.matrices(newdat,site,year,nsd[[i.site]][,i.year])
        # Store matrices of this replicate
        matrices.G[[i.site]][[i.year]] <- res.build.matrices$G
        matrices.P[[i.site]][[i.year]] <- res.build.matrices$P
        matrices.F[[i.site]][[i.year]] <- res.build.matrices$F
        matrices.K[[i.site]][[i.year]] <- res.build.matrices$K
        # Add matrices to the cumulative matrix for all replicates (to have a mean matrix to show)
        mat.G[[i.site]][[i.year]] <- mat.G[[i.site]][[i.year]] + res.build.matrices$G 
        mat.P[[i.site]][[i.year]] <- mat.P[[i.site]][[i.year]] + res.build.matrices$P
        mat.F[[i.site]][[i.year]] <- mat.F[[i.site]][[i.year]] + res.build.matrices$F
        mat.K[[i.site]][[i.year]] <- mat.K[[i.site]][[i.year]] + res.build.matrices$K
        # Lambda and w - needed ?
        lambda.det[i.site,i.year,i.repl] <- res.build.matrices$lambda.det
        # w[,i.site,i.year,i.repl] <- res.build.matrices$w
        # Life-cycle components
        mean.surv[i.site,i.year,i.repl] <- res.build.matrices$mean.surv
        mean.retro[i.site,i.year,i.repl] <- res.build.matrices$mean.retro
        mean.stasis[i.site,i.year,i.repl] <- res.build.matrices$mean.stasis
        mean.progr[i.site,i.year,i.repl] <- res.build.matrices$mean.progr
        mean.pfruit[i.site,i.year,i.repl] <- res.build.matrices$mean.pfruit
        mean.nfruits[i.site,i.year,i.repl] <- res.build.matrices$mean.nfruits
        mean.recr.size[i.site,i.year,i.repl] <- weighted.mean(c(1:n),nsd[[i.site]][,i.year])  # Average of size (1:n) weighted by the size distribution of recruits (nsd_year)
        
        # Stage-specific vital rates for all stages (used for elasticities) (used also for SLTRE with F0, F1 and F2)
        pfruit.all[[i.site]][,i.year] <- res.build.matrices$pfruit
        log.nfruits.all[[i.site]][,i.year] <- res.build.matrices$log.nfruits
        # Stage-specific vital rates for 4 representative stages (stored) - A cosa servono ? A fare la Fig. 1 ? 
        surv[,i.site,i.year,i.repl] <- colSums(res.build.matrices$P)[repr.stage]
        for (k in 1 : length(repr.stage)) {
            retro[k,i.site,i.year,i.repl] <- sum(res.build.matrices$G[(1:(repr.stage[k]-1)),repr.stage[k]])
            stasis[k,i.site,i.year,i.repl] <- res.build.matrices$G[repr.stage[k],repr.stage[k]]
            progr[k,i.site,i.year,i.repl] <- sum(res.build.matrices$G[(repr.stage[k]+1):n,repr.stage[k]])
        }
        if (any(repr.stage == 1)) retro[which(repr.stage == 1),i.site,i.year,i.repl] <- 0 # For the stage 1, no retrogression is possible, but the sum above puts a value
        pfruit[,i.site,i.year,i.repl] <- pfruit.all[[i.site]][repr.stage,i.year]
        nfruits[,i.site,i.year,i.repl] <- 10^log.nfruits.all[[i.site]][repr.stage,i.year]
      }
      # End loop on years
      # Stochastic growth rate
      r[i.site,i.repl] <- stochGrowthRateSampleList(nRunIn=1000, tMax=11000, listIPMmatrix=matrices.K[[i.site]])
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
      rm(res.build.matrices,res.elast)
    }
    # End loop on sites
    # SLTRE
    res.calc.SLTRE       <- calc.SLTRE(matrices.G, matrices.P, pfruit.all, log.nfruits.all, nsd)
    cnet[,,,i.repl]      <- res.calc.SLTRE$cnet
    log.lambda.R[i.repl] <- res.calc.SLTRE$log.lambda.R
    rm(res.calc.SLTRE)
    i.repl <- i.repl + 1
  } # End loop on boot
} # End loop on clim


save(nrepl, y,
     surv, progr, retro, stasis, pfruit, nfruits,
     mean.surv, mean.progr, mean.retro, mean.stasis, mean.pfruit, mean.nfruits, mean.recr.size,
     r, ws, elast, lambda.det,
     cnet, log.lambda.R,
     file=paste0("Results_MPM_parallel",iparallel,".RData"))

save(mat.G, mat.P, mat.F, mat.K,
      file=paste0("Results_MPM_matrices_parallel",iparallel,".RData"))





