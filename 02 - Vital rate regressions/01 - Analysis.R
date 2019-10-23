# Analysis


rm(list=ls())

library(tidyverse)
library(lme4)
library(MuMIn)

# Replicate paramters
rep.clim <- 1
start.boot <- 1
n.boot <- 200

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
load(paste0(getwd(),"/../01 - iButton data and imputation/01 - Resampled data/02 - Monthly/data.monthly.",rep.clim,".RData"))
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
  # Remove case affected by sheep trampling
  data_boot_surv <-
    data_boot %>% filter(SiteYear != "LAU2009")
  # Select only necessary columns and filter out NA
  data_boot_surv <-
    data_boot_surv %>%
    select(Site, SiteYear, Nb_Tot.x, Nb_Tot.x.sq, Fate, Axis1, Axis2, T.mean_m, T.mean_m.sq, T.range_m) %>%
    na.omit()
  # Scale variables and save mean and sd
  predictors.surv.mean <- data_boot_surv %>% summarise_at(c("Nb_Tot.x","Nb_Tot.x.sq","T.mean_m","T.mean_m.sq","T.range_m", "Axis1", "Axis2"),mean)
  predictors.surv.sd <-data_boot_surv %>% summarise_at(c("Nb_Tot.x","Nb_Tot.x.sq","T.mean_m","T.mean_m.sq","T.range_m", "Axis1", "Axis2"),sd)
  data_boot_surv_scaled <- data_boot_surv %>% mutate_at(c("Nb_Tot.x","Nb_Tot.x.sq","T.mean_m","T.mean_m.sq","T.range_m", "Axis1", "Axis2"),scale)
  # Fit
  m <- glmer(Fate ~ Nb_Tot.x + Nb_Tot.x.sq + T.mean_m + T.mean_m.sq + T.range_m + Axis1 + Axis2 + (1 | Site/SiteYear) ,family=binomial(),
             data=data_boot_surv_scaled, na.action="na.fail",nAGQ=0, control = glmerControl(calc.derivs = FALSE))
  d.surv <- lapply(dredge(m, evaluate = FALSE), eval)
  m.weight.surv <- calc.m.weight(d.surv)

  
  # Analysis of growth
  # Select only necessary columns and filter out NA
  data_boot_growth <-
    data_boot %>%
    select(Site, SiteYear, Nb_Tot.x, Nb_Tot.x.sq, Nb_Tot.y, Axis1, Axis2, T.mean_f, T.mean_f.sq, T.range_f) %>%
    na.omit()
  # Scale variables and save mean and sd, subtract 1 from dependent variable
  predictors.growth.mean <- data_boot_growth %>% summarise_at(c("Nb_Tot.x","Nb_Tot.x.sq","T.mean_f","T.mean_f.sq","T.range_f", "Axis1", "Axis2"),mean)
  predictors.growth.sd <-data_boot_growth %>% summarise_at(c("Nb_Tot.x","Nb_Tot.x.sq","T.mean_f","T.mean_f.sq","T.range_f", "Axis1", "Axis2"),sd)
  data_boot_growth_scaled <- data_boot_growth %>% mutate_at(c("Nb_Tot.x","Nb_Tot.x.sq","T.mean_f","T.mean_f.sq","T.range_f", "Axis1", "Axis2"),scale)
  data_boot_growth_scaled$Nb_Tot.y_m1 <- data_boot_growth_scaled$Nb_Tot.y - 1
  # Fit
  m <- glmer.nb(Nb_Tot.y_m1 ~ Nb_Tot.x + Nb_Tot.x.sq + T.mean_f + T.mean_f.sq + T.range_f + Axis1 + Axis2 + (1 | Site/SiteYear),
                data=data_boot_growth_scaled, na.action="na.fail", control = glmerControl(calc.derivs = FALSE), nAGQ=0)
  d.growth <- lapply(dredge(m, evaluate = FALSE), eval)
  m.weight.growth <- calc.m.weight(d.growth)
  
  
  # Analysis of fruiting
  # Select only necessary columns and filter out NA
  data_boot_fruiting <-
    data_boot %>%
    select(Site, SiteYear, Nb_Tot.x, Nb_Tot.x.sq, Fru, Axis1, Axis2, T.mean, T.mean.sq, T.range) %>%
    na.omit()
  # Scale variables and save mean and sd
  predictors.fruiting.mean <- data_boot_fruiting %>% summarise_at(c("Nb_Tot.x", "Nb_Tot.x.sq","T.mean","T.mean.sq","T.range", "Axis1", "Axis2"),mean)
  predictors.fruiting.sd <-data_boot_fruiting %>% summarise_at(c("Nb_Tot.x", "Nb_Tot.x.sq","T.mean","T.mean.sq","T.range", "Axis1", "Axis2"),sd)
  data_boot_fruiting_scaled <- data_boot_fruiting %>% mutate_at(c("Nb_Tot.x", "Nb_Tot.x.sq","T.mean","T.mean.sq","T.range", "Axis1", "Axis2"),scale)
  # Fit
  m <- glmer(Fru ~ Nb_Tot.x + Nb_Tot.x.sq + T.mean + T.mean.sq + T.range + Axis1 + Axis2 + (1 | Site/SiteYear), family=binomial(),
             data=data_boot_fruiting_scaled, na.action="na.fail", nAGQ=0, control = glmerControl(calc.derivs = FALSE))
  d.fruiting <- lapply(dredge(m, evaluate = FALSE), eval)
  m.weight.fruiting <- calc.m.weight(d.fruiting)


  # Analysis of fruits
  # Select only necessary columns and filter out NA
  data_boot_fruits <-
    data_boot %>%
    select(Site, SiteYear, Nb_Tot.x, Nb_Tot.x.sq, Fru, logTotalSil, Axis1, Axis2, T.mean, T.mean.sq, T.range) %>%
    filter(Fru == 1) %>%
    na.omit()
  # Scale variables and save mean and sd
  predictors.fruits.mean <- data_boot_fruits %>% summarise_at(c("Nb_Tot.x", "Nb_Tot.x.sq","T.mean","T.mean.sq","T.range", "Axis1", "Axis2"),mean)
  predictors.fruits.sd <-data_boot_fruits %>% summarise_at(c("Nb_Tot.x", "Nb_Tot.x.sq","T.mean","T.mean.sq","T.range", "Axis1", "Axis2"),sd)
  data_boot_fruits_scaled <- data_boot_fruits %>% mutate_at(c("Nb_Tot.x", "Nb_Tot.x.sq","T.mean","T.mean.sq","T.range", "Axis1", "Axis2"),scale)
  # Fit
  m <- lmer(logTotalSil ~ Nb_Tot.x + Nb_Tot.x.sq + T.mean + T.mean.sq + T.range + Axis1 + Axis2 + (1 | Site/SiteYear),
            data=data_boot_fruits_scaled, na.action="na.fail", REML=F, control = lmerControl(calc.derivs = FALSE))
  d.fruits <- lapply(dredge(m, evaluate = FALSE), eval)
  m.weight.fruits <- calc.m.weight(d.fruits)

  
  # Analysis of recruit size
  # Select only necessary columns and filter out large plants and NA
  data_boot_recruit <-
      data_boot %>%
      select(Site, SiteYear, is.recr, First_stade, Nb_Tot.x, logNbtot.x, Axis1, Axis2, T.mean, T.mean.sq, T.range) %>%
      filter(is.recr == 1) %>%
      filter(First_stade == 0) %>%
      filter(Nb_Tot.x <= 5) %>%
      na.omit()
  # Scale variables and save mean and sd, subtract 1 from dependent variable
  predictors.recruit.mean <- data_boot_recruit %>% summarise_at(c("T.mean","T.mean.sq","T.range", "Axis1", "Axis2"),mean)
  predictors.recruit.sd <-data_boot_recruit %>% summarise_at(c("T.mean","T.mean.sq","T.range", "Axis1", "Axis2"),sd)
  data_boot_recruit_scaled <- data_boot_recruit %>% mutate_at(c("T.mean","T.mean.sq","T.range", "Axis1", "Axis2"),scale)
  data_boot_recruit_scaled$Nb_Tot.x_m1 <- data_boot_recruit_scaled$Nb_Tot.x - 1
  # Fit
  m <- glmer.nb(Nb_Tot.x_m1 ~ T.mean + T.mean.sq + T.range + Axis1 + Axis2 + (1 | Site/SiteYear),
                data=data_boot_recruit_scaled, na.action="na.fail", control = glmerControl(calc.derivs = FALSE), nAGQ=0)
  d.recruit <- lapply(dredge(m, evaluate = FALSE), eval)
  m.weight.recruit <- calc.m.weight(d.recruit)


  save(d.surv,
       d.growth,
       d.fruiting,
       d.fruits,
       d.recruit,
       m.weight.surv,
       m.weight.growth,
       m.weight.fruiting,
       m.weight.fruits,
       m.weight.recruit,
       predictors.surv.mean,
       predictors.growth.mean,
       predictors.fruiting.mean,
       predictors.fruits.mean,
       predictors.recruit.mean,
       predictors.surv.sd,
       predictors.growth.sd,
       predictors.fruiting.sd,
       predictors.fruits.sd,
       predictors.recruit.sd,
       file = paste0(paste0(getwd(),"/../../noBackup/bootstraps/Models_Clim",rep.clim,"_Boot",i.boot,".RData")) # The path where the GLMM models will be saved 
  )

  rm(data_boot)
  rm(data_boot_surv,data_boot_surv_scaled)
  rm(data_boot_growth,data_boot_growth_scaled)
  rm(data_boot_fruiting,data_boot_fruiting_scaled)
  rm(data_boot_fruits,data_boot_fruits_scaled)
  rm(data_boot_recruit,data_boot_recruit_scaled)

} #end loop on boot

