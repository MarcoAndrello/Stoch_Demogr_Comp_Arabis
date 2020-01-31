# Compose results MPM for SensSurv
# Marco Andrello
# Last modification: 13/01/2020

rm(list=ls())

nparallel <- 5 # Number of parallel runs

library(abind)

load("Results_MPM_SensSurv_parallel1.RData")
load("Results_MPM_matrices_SensSurv_parallel1.RData")

# Vital rates
f.nfruits <- nfruits
f.pfruit <- pfruit
f.progr <- progr
f.retro <- retro
f.stasis <- stasis
f.surv <- surv

# Life-cycle components
f.mean.nfruits <- mean.nfruits
f.mean.pfruit <- mean.pfruit
f.mean.progr <- mean.progr
f.mean.recr.size <- mean.recr.size
f.mean.retro <- mean.retro
f.mean.stasis <- mean.stasis
f.mean.surv    <- mean.surv

# Stage distribution
f.ws <- ws

# Other variables: growth rates, SLTRE, elasticities.
f.r <- r
f.cnet <- cnet
f.elast <- elast
f.lambda.det   <- lambda.det
f.log.lambda.R <- log.lambda.R

# Matrices
f.mat.G <- mat.G
f.mat.P <- mat.P
f.mat.F <- mat.F
f.mat.K <- mat.K

for (i.parallel in 2 : nparallel){
  load(paste0("Results_MPM_matrices_SensSurv_parallel",i.parallel,".RData"))
  load(paste0("Results_MPM_SensSurv_parallel",i.parallel,".RData"))

  # Vital rates
  f.nfruits <- abind(f.nfruits, nfruits)
  f.pfruit <- abind(f.pfruit, pfruit)
  f.progr <- abind(f.progr, progr)
  f.retro <- abind(f.retro, retro)
  f.stasis <- abind(f.stasis, stasis)
  f.surv <- abind(f.surv, surv)
  
  # Life-cycle components
  f.mean.nfruits    <- abind(f.mean.nfruits, mean.nfruits)
  f.mean.pfruit    <- abind(f.mean.pfruit, mean.pfruit)
  f.mean.recr.size    <- abind(f.mean.recr.size, mean.recr.size)
  f.mean.surv    <- abind(f.mean.surv, mean.surv)
  f.mean.progr <- abind(f.mean.progr, mean.progr)
  f.mean.retro <- abind(f.mean.retro, mean.retro)
  f.mean.stasis       <- abind(f.mean.stasis, mean.stasis)
  
  # Stage distribution
  f.ws <- abind(f.ws, ws)

  # Other variables: growth rates, SLTRE, elasticities.
  f.r <- cbind(f.r,r)
  f.cnet <- abind(f.cnet, cnet, along=4)
  f.elast <- abind(f.elast,elast, along=3)
  f.lambda.det   <- abind(f.lambda.det, lambda.det)
  f.log.lambda.R   <- c(f.log.lambda.R, log.lambda.R)
  
  # Matrices
  for (i.site in 1 : 6){
    for (i.year in 1 : 6){
      f.mat.G[[i.site]][[i.year]] <- f.mat.G[[i.site]][[i.year]] + mat.G[[i.site]][[i.year]]
      f.mat.P[[i.site]][[i.year]] <- f.mat.P[[i.site]][[i.year]] + mat.P[[i.site]][[i.year]]
      f.mat.F[[i.site]][[i.year]] <- f.mat.F[[i.site]][[i.year]] + mat.F[[i.site]][[i.year]]
      f.mat.K[[i.site]][[i.year]] <- f.mat.K[[i.site]][[i.year]] + mat.K[[i.site]][[i.year]]
    }
  }
}


# Vital rates
f.nfruits -> nfruits
f.pfruit -> pfruit
f.progr -> progr
f.retro -> retro
f.stasis -> stasis
f.surv -> surv

# Life-cycle components
f.mean.nfruits -> mean.nfruits
f.mean.pfruit    -> mean.pfruit
f.mean.recr.size -> mean.recr.size
f.mean.surv    -> mean.surv
f.mean.progr -> mean.progr
f.mean.retro -> mean.retro
f.mean.stasis       -> mean.stasis

# Stage distribution
f.ws -> ws

# Matrices
f.mat.G -> mat.G
f.mat.P -> mat.P
f.mat.F -> mat.F
f.mat.K -> mat.K

# Other variables: growth rates, SLTRE, elasticities.
f.r -> r
f.cnet -> cnet
f.elast -> elast
f.lambda.det   -> lambda.det
f.log.lambda.R   -> log.lambda.R

# Divide the "summed" matrices by the number of replicates (2000)
for (i.site in 1 : 6){
  for (i.year in 1 : 6){
    mat.G[[i.site]][[i.year]] <- mat.G[[i.site]][[i.year]] / 2000
    mat.P[[i.site]][[i.year]] <- mat.P[[i.site]][[i.year]] / 2000
    mat.F[[i.site]][[i.year]] <- mat.F[[i.site]][[i.year]] / 2000
    mat.K[[i.site]][[i.year]] <- mat.K[[i.site]][[i.year]] / 2000
  }
}



save(
    nfruits,
    pfruit,
    progr,
    retro,
    stasis,
    surv,
    
    mean.nfruits,
    mean.pfruit,
    mean.recr.size,
    mean.surv,
    mean.progr,
    mean.retro,
    mean.stasis,
    
    ws,
    
    mat.G,
    mat.P,
    mat.F,
    mat.K,
    
    y,
    r,
    cnet,
    elast,
    lambda.det,
    log.lambda.R,
    
    file="Results_MPM_SensSurv.RData"
)



