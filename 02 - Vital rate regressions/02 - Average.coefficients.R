# Calculate average coefficients

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
load("C:/Users/mandrell/Desktop/BOOTSTRAPS/Models_Clim1_Boot1.RData")
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
  dimnames(tabsd) <- list(group = as.data.frame(VarCorr(get(d.var)[[1]]))$grp,
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
    load(paste0("C:/Users/mandrell/Desktop/BOOTSTRAPS/Models_Clim",i.clim,"_Boot",i.boot,".RData"))
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



# # Compose results
# rm(list=ls())
# nparallel <- 5
# names_var <- c("surv","growth","fruiting","fruits","recruit")
# 
# load("Coefficients_1.RData")
# f.coef_fixed <- coef_fixed
# f.sd_rand <- sd_rand
# f.theta <- theta
# f.r_squared <- r_squared
# for (i.parallel in 2 : nparallel){
#   for (i.var in 1 : 5) {
#     load(paste0("Coefficients_",i.parallel,".RData"))
#     f.coef_fixed[[i.var]] <- cbind(f.coef_fixed[[i.var]],coef_fixed[[i.var]])
#     f.sd_rand[[i.var]] <- cbind(f.sd_rand[[i.var]],sd_rand[[i.var]])
#     if (i.var == 2 | i.var == 5) f.theta[[i.var]] <- c(f.theta[[i.var]],theta[[i.var]])
#     f.r_squared[[i.var]] <- cbind(f.r_squared[[i.var]],r_squared[[i.var]])
#   }
# }
# 
# f.coef_fixed -> coef_fixed
# f.sd_rand    -> sd_rand
# f.theta      -> theta
# f.r_squared  -> r_squared
# rm(f.coef_fixed, f.sd_rand, f.theta, f.r_squared)
# save(coef_fixed, sd_rand, theta, r_squared, file="Coefficients.RData")
# 
# load("Coefficients.RData")
# nrepl <- dim(sd_rand[[1]])[2]
# 
# ## Fixed factors
# table_res_coef <- array(NA,c(8,length(names_var)*3))
# colnames(table_res_coef) <- rep("",length(names_var)*3)
# colnames(table_res_coef)[(((c(1:5)-1)*3)+1)] <- paste0(names_var,"_Mean")
# colnames(table_res_coef)[(((c(1:5)-1)*3)+2)] <- name.low <-paste0(names_var,"_low")
# colnames(table_res_coef)[(c(1:5)*3)] <- name.high <- paste0(names_var,"_high")
# # Add two lines for recruit
# nrepl <- dim(coef_fixed[[5]])[2] # Recruit
# coef_fixed[[5]] <- rbind(coef_fixed[[5]][1:3,],
#                          rep(NA,nrepl),
#                          rep(NA,nrepl),
#                          coef_fixed[[5]][4:6,])
# # Calculate mean and CI
# for(i.var in 1 : length(names_var)) {
#   table_res_coef[,(((i.var-1)*3)+1)] <- apply(coef_fixed[[i.var]],1,mean,na.rm=T)
#   table_res_coef[,(((i.var-1)*3)+2)] <- apply(coef_fixed[[i.var]],1,quantile,c(0.025),na.rm=T)
#   table_res_coef[,(i.var*3)] <- apply(coef_fixed[[i.var]],1,quantile,c(0.975),na.rm=T)
# }
# # Write csv
# table_res_coef %>% as_tibble %>% round(2) %>% add_column(Predictor=row.names(coef_fixed[[2]])) %>% write.csv("Results_coef_fixed.csv")
# 
# ## Random factors
# table_res_sd <- array(NA,c(3,length(names_var)*3))
# colnames(table_res_sd) <- rep("",length(names_var)*3)
# colnames(table_res_sd)[(((c(1:5)-1)*3)+1)] <- paste0(names_var,"_Mean")
# colnames(table_res_sd)[(((c(1:5)-1)*3)+2)] <- name.low <-paste0(names_var,"_low")
# colnames(table_res_sd)[(c(1:5)*3)] <- name.high <- paste0(names_var,"_high")
# # Add two blank lines for surv and fruiting
# for (i.var in c(1,3)){
#   sd_rand[[i.var]] <- rbind(sd_rand[[i.var]][1:2,],
#                             rep(NA,nrepl))
# }
# # Add theta for growth and recruit
# for (i.var in c(2,5)){
#   sd_rand[[i.var]] <- rbind(sd_rand[[i.var]][1:2,],
#                             theta[[i.var]])
# }
# # Calculate mean and CI
# for(i.var in 1 : length(names_var)) {
#   table_res_sd[,(((i.var-1)*3)+1)] <- apply(sd_rand[[i.var]],1,mean,na.rm=T)
#   table_res_sd[,(((i.var-1)*3)+2)] <- apply(sd_rand[[i.var]],1,quantile,c(0.025),na.rm=T)
#   table_res_sd[,(i.var*3)] <- apply(sd_rand[[i.var]],1,quantile,c(0.975),na.rm=T)
# }
# # Write csv
# table_res_sd %>% as_tibble %>% round(2) %>% add_column(Predictor=row.names(sd_rand[[2]])) %>% write.csv("Results_sd_rand.csv")
# 
# 
# ## R-squared
# table_res_r <- array(NA,c(2,length(names_var)*3))
# colnames(table_res_r) <- rep("",length(names_var)*3)
# colnames(table_res_r)[(((c(1:5)-1)*3)+1)] <- paste0(names_var,"_Mean")
# colnames(table_res_r)[(((c(1:5)-1)*3)+2)] <- name.low <-paste0(names_var,"_low")
# colnames(table_res_r)[(c(1:5)*3)] <- name.high <- paste0(names_var,"_high")
# # Calculate mean and CI
# for(i.var in 1 : length(names_var)) {
#   table_res_r[,(((i.var-1)*3)+1)] <- apply(r_squared[[i.var]],1,mean,na.rm=T)
#   table_res_r[,(((i.var-1)*3)+2)] <- apply(r_squared[[i.var]],1,quantile,c(0.025),na.rm=T)
#   table_res_r[,(i.var*3)] <- apply(r_squared[[i.var]],1,quantile,c(0.975),na.rm=T)
# }
# # Write csv
# table_res_r %>% as_tibble %>% round(2) %>% add_column(r.squared=row.names(r_squared[[2]])) %>% write.csv("Results_r_squared.csv")
# 
# 
# 
# # Correlation between site and Year standard deviation
# i.var <- 1
# plot(sd_rand[[i.var]][1,],sd_rand[[i.var]][2,],xlab="Year within Site",ylab="Site")
# 
# names_vitalr <- c("Survival","Growth","Reproduction","Fecundity","Recruit size")
# par(mfrow=c(3,2))
# for (i.var in 1 : 5) {
# hist(sd_rand[[i.var]][2,],main=names_vitalr[i.var],xlab="Site (SD)")
# }
# 

