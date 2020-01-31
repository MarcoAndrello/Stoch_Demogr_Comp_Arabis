# Calculate average coefficients
# Marco Andrello
# Last modification: 12/01/2020

rm(list=ls())

library(tidyverse)
library(lme4)
library(MuMIn)

# Controls replicates
i.parallel <- 1
startclim <- 1
endclim <- 2
nboot <- 200
nrepl <- (endclim-startclim+1)*nboot

sample_model <- function(object,weights) {
  n.models <- length(object)
  i.model <- sample(n.models,1,replace=T,prob=weights)
  object[i.model][[1]]
}

get_null_model <- function(object) {
  object[1][[1]]
}

names_var <- c("surv","growth","fruiting","fruits","recruit")

# Create the object to store the values of coefficients
coef_fixed  <- list()
sd_rand     <- list()
theta       <- list()
r_squared   <- list()
load("G:/noBackup/Arabis/ResultsGLMM/Models_Clim1_Boot1.RData")
for (i.var in 1 : length(names_var)){
  var <- names_var[i.var]
  d.var <- paste0("d.",var)
  # Prepare array for coefficients of fixed effects
  get(d.var)[[length(get(d.var))]] %>%
    fixef %>%
    names -> 
    names_predictors
  tabcoef <- array(NA,c(length(names_predictors),nrepl))
  dimnames(tabcoef) <- list(predictors=names_predictors,repl=c(1:nrepl))
  coef_fixed[[i.var]] <- tabcoef
  # Prepare array for sd of random effects
  tabsd <- array(NA,c(length(as.data.frame(VarCorr(get(d.var)[[1]]))$grp),
                      nrepl))
  # dimnames(tabsd) <- list(group = as.data.frame(VarCorr(get(d.var)[[1]]))$grp,
  #                         repl=c(1:nrepl))
  dimnames(tabsd) <- list(group = paste(as.data.frame(VarCorr(get(d.var)[[1]]))[,1],
                                        as.data.frame(VarCorr(get(d.var)[[1]]))[,2],
                                        as.data.frame(VarCorr(get(d.var)[[1]]))[,3],
                                        sep="_"),
                          repl=c(1:nrepl))
  sd_rand[[i.var]] <- tabsd
  # Prepare array for theta
  if (i.var == 2 | i.var == 5) { # Add theta
      theta[[i.var]] <- rep(NA,nrepl)
      }
  # Prepare array for R squared
  r_squared[[i.var]] <- array(NA,c(2,
                                   nrepl))
  dimnames(r_squared[[i.var]]) <- list(R2 = c("R2m","R2c"),
                                       repl=c(1:nrepl) )
}
names(coef_fixed) <- names_var
names(sd_rand) <- names_var
names(r_squared) <- names_var

# Loop on climate replicates and bootstrap replicates, for all variables
i.repl <- 1
for (i.clim in startclim : endclim){
  for(i.boot in 1 : nboot) {
    # Load models
    cat("Clim",i.clim,"; Boot",i.boot,"of",nboot,"\n"); flush.console()
    load(paste0("G:/noBackup/Arabis/ResultsGLMM/Models_Clim",i.clim,"_Boot",i.boot,".RData"))
    # Loop on variables
    for (i.var in 1 : length(names_var)) {
      var <- names_var[i.var]
      d.var <- paste0("d.",var)
      m.weight <- paste0("m.weight.",var)
      # Sample model
      mod <- sample_model(get(d.var),get(m.weight))
      # Read values of fixed effects
      fe <- fixef(mod)
      id <- match( names(fe),row.names(coef_fixed[[i.var]]) )
      coef_fixed[[i.var]][id,i.repl] <- fe
      # Read values of random effects
      re <- as.data.frame(VarCorr(mod))$sdcor
      sd_rand[[i.var]][,i.repl] <- re
      # Read values of theta
      if (i.var == 2 | i.var == 5) {
          theta[[i.var]][i.repl] <- getME(mod,"glmer.nb.theta")
          }
      # Calculate R squared
      mod.null <- get_null_model(get(d.var))
      if (i.var == 4) r.squaredGLMM(mod,null=mod.null) -> r_squared[[i.var]][,i.repl]
      if (i.var == 1 | i.var == 3) r.squaredGLMM(mod,null=mod.null)["theoretical",] -> r_squared[[i.var]][,i.repl] 
      if (i.var == 2 | i.var == 5) r.squaredGLMM(mod,null=mod.null)["trigamma",] -> r_squared[[i.var]][,i.repl] 
    }
    i.repl <- i.repl + 1
  } # End loop boot
} # End loop clim

save(coef_fixed, sd_rand, theta, r_squared, file=paste0("Coefficients_",i.parallel,".RData"))


