# Build and analyse MPM for Seed Bank model
# Marco Andrello
# Last modification: 21/01/2020

rm(list=ls())

library(popbio)
library(IPMpack)
library(tidyverse)
library(lme4)
library(gridExtra)
library(fields)
library(abind)

# Replication parameters
elev <- c(930,1480,1980,2090,2590,2930)
#germ <- rep(1,6)                 # Germination rate: fixed (0.1, 0.4, 0.7)
germ <- 1-elev/3000              # or site specific   (a vector)
iparallel <- 1
set.seed(iparallel)
startclim <- 1
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

# Stochastic growth rate r = log-lambda (should be a = log-lambda)
r <- array(NA,c(length(v.site),
                nrepl))
dimnames(r) <-  list(site = v.site,
                     repl = c(1:nrepl))

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
# SEED BANK MODEL: Objects containing the different projection matrices
P1 <- S1 <- G1 <- F1 <- K1 <- list()
# SLTRE: net contributions of each demographic paramters (descriptor of a life-cycle component) to differences in log-lambdas between sites
# SEED BANK MODEL: We have an 8th life-cycle component: germination
cnet <- array(NA,c(8, # Number of life-cycle components (survival, retro.growth, stasis, prog.growth, repro, repr.out, recr, germ)
                   4, # Descriptor (mean, elast, cv, corr)
                   length(v.site),
                   nrepl))
dimnames(cnet) <- list(par = c("S","G-","G=","G+","F0","F1","F2","gamma"),
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
    load(paste0("F:/noBackup/Arabis/ResultsGLMM/Models_Clim",i.clim,"_Boot",i.boot,".RData"))
    
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
      # SEED BANK: Modified matrices
      P1[[i.site]] <- S1[[i.site]] <- G1[[i.site]] <- F1[[i.site]] <- K1[[i.site]] <- list()
      
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
        # Calculate the newborn size distribution
        nsd[[i.site]][,i.year] <- h * newb.size.distr(newdat,site,year)
        # Build matrices and calculate life-cycle components
        res.build.matrices <- build.matrices(newdat,site,year,nsd[[i.site]][,i.year])
        # Store matrices of this replicate
        matrices.G[[i.site]][[i.year]] <- res.build.matrices$G
        matrices.P[[i.site]][[i.year]] <- res.build.matrices$P
        matrices.F[[i.site]][[i.year]] <- res.build.matrices$F
        # matrices.K[[i.site]][[i.year]] <- res.build.matrices$K
     
        # Stage-specific vital rates for all stages (used for elasticities and for SLTRE with F0, F1 and F2)
        pfruit.all[[i.site]][,i.year] <- res.build.matrices$pfruit
        log.nfruits.all[[i.site]][,i.year] <- res.build.matrices$log.nfruits
      
        # ........................................................Modify matrices to add a seed bank stage: only to calculate log_lambda, not for the SLTRE
        surv <- colSums(matrices.P[[i.site]][[i.year]])
        # Modify G
        g_sb <- c((1-germ[i.site]),germ[i.site]*nsd[[i.site]][,i.year])
        G1[[i.site]][[i.year]] <- rbind(rep(0,n),
                                        matrices.G[[i.site]][[i.year]])
        G1[[i.site]][[i.year]] <- cbind(g_sb,G1[[i.site]][[i.year]])
        # Modify F
        fj <- colSums(matrices.F[[i.site]][[i.year]]) / (germ[i.site]*(2-germ[i.site]))
        F1[[i.site]][[i.year]] <- nsd[[i.site]][,i.year] %*% t(fj) * germ[i.site]
        F1[[i.site]][[i.year]] <- rbind(fj * (1-germ[i.site]),
                                        F1[[i.site]][[i.year]])
        F1[[i.site]][[i.year]] <- cbind(rep(0,(n+1)),
                                        F1[[i.site]][[i.year]])
        # Modify S and K
        S1[[i.site]][[i.year]] <- matrix(c(1,surv),(n+1),(n+1),byrow=T)
        P1[[i.site]][[i.year]] <- S1[[i.site]][[i.year]]*G1[[i.site]][[i.year]]
        K1[[i.site]][[i.year]] <- P1[[i.site]][[i.year]] + F1[[i.site]][[i.year]]
      }
      # End loop on years
      # Stochastic growth rate
      # r[i.site,i.repl] <- stochGrowthRateSampleList(nRunIn=1000, tMax=11000, listIPMmatrix=matrices.K[[i.site]])
      r[i.site,i.repl] <- stochGrowthRateSampleList(nRunIn=1000, tMax=11000, listIPMmatrix=K1[[i.site]])
      rm(res.build.matrices)
    }
    # End loop on sites
    
    # SLTRE
    # We input the SLTRE funciton the untransformed matrices
    res.calc.SLTRE       <- calc.SLTRE_SeedBank(matrices.G, matrices.P, pfruit.all, log.nfruits.all, nsd, germ=germ)
    cnet[,,,i.repl]      <- res.calc.SLTRE$cnet
    log.lambda.R[i.repl] <- res.calc.SLTRE$log.lambda.R
    rm(res.calc.SLTRE)
    i.repl <- i.repl + 1
  } # End loop on boot
} # End loop on clim


suffix_name_file <- if (var(germ) == 0) mean(germ) else "VarElev"
save(nrepl, y,
     r,
     cnet, log.lambda.R,
     file=paste0("Results_MPM_SeedBank_germ_",suffix_name_file,"_.RData"))






