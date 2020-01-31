# Compose results MPM for SensRef
# Marco Andrello
# Last modification: 13/01/2020

rm(list=ls())
scen <- 5
nparallel <- 2 # Number of parallel runs
library(abind)
load(paste0("Results_MPM_SensRef_scen_",scen,"_parallel1.RData"))
f.cnet <- cnet
f.log.lambda.R <- log.lambda.R
for (i.parallel in 2 : nparallel){
  load(paste0("Results_MPM_SensRef_scen_",scen,"_parallel",i.parallel,".RData"))
  f.cnet <- abind(f.cnet, cnet, along=4)
  f.log.lambda.R   <- c(f.log.lambda.R, log.lambda.R)
}
f.cnet -> cnet
f.log.lambda.R   -> log.lambda.R

save(
    y,
    cnet,
    log.lambda.R,
    file=paste0("Results_MPM_SensRef_",scen,".RData")
)
