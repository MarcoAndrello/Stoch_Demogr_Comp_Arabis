# Analysis


rm(list=ls())

library(tidyverse)
library(lme4)
library(MuMIn)

# Replicate paramters
rep.clim <- 1
start.boot <- 1
n.boot <- 2

# Function to append climatic data of the following yeat to year t
extend.env.dataset <- function(env) {
  env.onlyclim <- ungroup(env) %>% select(-c(Site, Year, Quad, SiteQuad))
  env <-  mutate(env, UniqueID = paste0(SiteQuad,"-",Year))  
  env_dat_follow <- env_dat_mean <- array(NA,dim(env.onlyclim))
  for (i in 1 : dim(env)[1]) {
    current_UniqueID <- env$UniqueID[i]
    current_year <- env$Year[i]
    follow_year <- current_year + 1
    follow_UniqueID <- current_UniqueID
    substr(follow_UniqueID,8,11) <- as.character(follow_year)
    id <- which(env$UniqueID == follow_UniqueID)
    env_dat_follow[i,] <- as.numeric(env.onlyclim[id,])
    env_dat_mean[i,] <- ( as.numeric(env.onlyclim[i,]) + as.numeric(env.onlyclim[id,]) ) / 2
  }
  colnames(env_dat_follow) <- paste0(colnames(env.onlyclim),"_f")
  colnames(env_dat_mean) <- paste0(colnames(env.onlyclim),"_m")
  env <- cbind(as.data.frame(env),env_dat_follow,env_dat_mean)
  env <- as_tibble(env) %>% select(-UniqueID)
  rm(current_UniqueID,follow_UniqueID,current_year,follow_year,env_dat_follow,env_dat_mean)
  return(env)
}

# Function to calculate model weights
calc.m.weight <- function(object){
  id.min <- lapply(object,AICc) %>% unlist %>% which.min()
  AICcBase <- AICc(object[[id.min]])
  deltaAICc <- unlist(lapply(object,AICc)) - AICcBase 
  exp(-0.5 * deltaAICc) / sum(exp(-0.5 * deltaAICc))
}

# Load environmental non-climatic data
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
    # Keep only the PCA Axis
    select(Site, Quad, SiteQuad, Axis1, Axis2)

# Load demographic data
load("Demographic_data.RData")
data <-
  data %>%
  mutate(SiteQuad=paste0(Site, "-", Quad), Year = as.integer(Year)) %>%
  # Join the environmental data to the demographic data
  left_join(soilveg, by="SiteQuad") %>%
  # Change some names and rearrange columns
  mutate(Site = Site.y,
         Quad = Quad.x,
         logNbtot.x = log10(Nb_Tot.x),
         logNbtot.y = log10(Nb_Tot.y),
         logTotalSil = log10(Total_Sil),
         SiteYear = paste0(Site,Year)) %>%
  select(Site, Year, Quad, SiteQuad, SiteYear, ID,
         Nb_Tot.x, logNbtot.x,
         Nb_Tot.y, logNbtot.y,
         Fru, Total_Sil, logTotalSil,
         Fate,
         is.recr, First_stade,
         Axis1, Axis2)

######### here loop on resampled climatic data


# Load resampled climatic data
load(paste0(getwd(),"/../01 - Environmental data/01 - Climatic/01 - Resampled data/02 - Monthly/data.monthly.",rep.clim,".RData"))
monthly <- 
  monthly %>%
  rename(Site=Pop) %>% 
  mutate(SiteQuad = paste0(Site, "-Q", Quad)) %>%
  select(Site, Year, Quad, SiteQuad, AvgMeanTemp, AvgAmp) %>%
  rename(T.mean = AvgMeanTemp, T.range = AvgAmp)

## Set outlier temperature to NA (correct aberrant prediction for VIL 2008)
monthly$T.range[which(monthly$T.range>32)] <- NA

# Extend climatic data
monthly <- extend.env.dataset(monthly)

# Join resampled climatic data to demographic data
data_rep.clim <-
  left_join(data, monthly, by = c("SiteQuad", "Year")) %>%
  rename(Site = Site.x, Quad = Quad.x) %>%
    select(-c(Site.y,Quad.y))


######### here loop on bootstrap

set.seed(rep.clim*5)
i.boot <- 1
for (i.boot in start.boot : n.boot) {
  cat("Boot",i.boot,"of",n.boot,"\n")
  # Boostrap the data
  data_boot <- data_rep.clim[sample(dim(data_rep.clim)[1],replace=T),]
  
  # Here create quadratic variables
  data_boot <- 
    mutate(data_boot, 
         Nb_Tot.x.sq = Nb_Tot.x^2,
         T.mean.sq = T.mean^2,
         T.mean_f.sq = T.mean_f^2,
         T.mean_m.sq = T.mean_m^2)
         
  # Analysis of survival
  # DO NOT Remove case affected by sheep trampling
  # data_boot_surv <- data_boot %>% filter(SiteYear != "LAU2009")
  data_boot_surv <- data_boot
  # Select only necessary columns and filter out NA
  data_boot_surv <-
    data_boot_surv %>%
    select(Site, SiteYear, Nb_Tot.x, Nb_Tot.x.sq, Fate, Axis1, Axis2, T.mean_m, T.mean_m.sq, T.range_m) %>%
    na.omit()
  # Scale variables and save mean and sd
  predictors.surv.mean <- data_boot_surv %>% summarise_at(c("Nb_Tot.x","Nb_Tot.x.sq","T.mean_m","T.mean_m.sq","T.range_m", "Axis1", "Axis2"),mean)
  predictors.surv.sd <-data_boot_surv %>% summarise_at(c("Nb_Tot.x","Nb_Tot.x.sq","T.mean_m","T.mean_m.sq","T.range_m", "Axis1", "Axis2"),sd)
  data_boot_surv_scaled <- data_boot_surv %>% mutate_at(c("Nb_Tot.x","Nb_Tot.x.sq","T.mean_m","T.mean_m.sq","T.range_m", "Axis1", "Axis2"),scale)
  data_boot_surv_scaled$SiteYear <- factor(data_boot_surv_scaled$SiteYear)
  # Fit
  m <- glmer(Fate ~ Nb_Tot.x + Nb_Tot.x.sq + T.mean_m + T.mean_m.sq + T.range_m + Axis1 + Axis2 + (1 | Site) + (1 + Axis2 | SiteYear), family=binomial(),
             data=data_boot_surv_scaled, na.action="na.fail",nAGQ=0, control = glmerControl(calc.derivs = FALSE))
  d.surv <- lapply(dredge(m, evaluate = FALSE), eval)
  m.weight.surv <- calc.m.weight(d.surv)


  save(d.surv,
       m.weight.surv,
       predictors.surv.mean,
       predictors.surv.sd,
       file = paste0(paste0(getwd(),"/../../noBackup/Arabis/ResultsGLMM/SensSurv/Models_Clim",rep.clim,"_Boot",i.boot,".RData")) # The path where the GLMM models will be saved 
  )

  rm(data_boot)
  rm(data_boot_surv,data_boot_surv_scaled)

} #end loop on boot

