# Analysis


rm(list=ls())

library(tidyverse)
library(lme4)
library(MuMIn)

# Function to append climatic data of the following yeat to year t (same as in the main script, so it can be put in an heklper file)
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

# Load AVERAGE climatic data
load(paste0(getwd(),"/../01 - Environmental data/01 - Climatic/reconstruct_monthly.Rdata"))
monthly <- 
  monthly %>%
  rename(Site=Pop) %>% 
  mutate(SiteQuad = paste0(Site, "-Q", Quad)) %>%
  select(Site, Year, Quad, SiteQuad, AvgDTemp, AvgAmp) %>%
  rename(T.range = AvgAmp) %>%
  rename(T.mean = AvgDTemp)

## Set outlier temperature to NA (correct aberrant prediction for VIL 2008)
monthly$T.range[which(monthly$T.range>32)] <- NA

# Extend climatic data
monthly <- extend.env.dataset(monthly)

# Join AVERAGE climatic data to demographic data
data_rep.clim <-
  left_join(data, monthly, by = c("SiteQuad", "Year")) %>%
  rename(Site = Site.x, Quad = Quad.x) %>%
    select(-c(Site.y,Quad.y))

# Here create quadratic variables
data_rep.clim <- 
    mutate(data_rep.clim, 
           Nb_Tot.x.sq = Nb_Tot.x^2,
           T.mean.sq = T.mean^2,
           T.mean_f.sq = T.mean_f^2,
           T.mean_m.sq = T.mean_m^2)

# Analysis of survival
# Remove case affected by sheep trampling
data_surv <-
    data_rep.clim %>% filter(SiteYear != "LAU2009")
# Select only necessary columns and filter out NA
data_surv <-
    data_surv %>%
    select(Site, SiteYear, SiteQuad, Nb_Tot.x, Nb_Tot.x.sq, Fate, Axis1, Axis2, T.mean_m, T.mean_m.sq, T.range_m) %>%
    na.omit()
# Scale variables
data_surv_scaled <- data_surv %>% mutate_at(c("Nb_Tot.x","Nb_Tot.x.sq","T.mean_m","T.mean_m.sq","T.range_m", "Axis1", "Axis2"),scale)
data_surv_scaled$SiteYear <- factor(data_surv_scaled$SiteYear)
data_surv_scaled$SiteQuad <- factor(data_surv_scaled$SiteQuad)
# Fit
# No random slope
m01 <- glmer(Fate ~ Nb_Tot.x + Nb_Tot.x.sq + T.mean_m + T.mean_m.sq + T.range_m + Axis1 + Axis2 + (1 | Site) + (1 | SiteYear),family=binomial(),
           data=data_surv_scaled, na.action="na.fail",nAGQ=0, control = glmerControl(calc.derivs = FALSE))
m02 <- glmer(Fate ~ Nb_Tot.x + Nb_Tot.x.sq + T.mean_m + T.mean_m.sq + T.range_m + Axis1 + Axis2 + (1 | Site) + (1 | SiteQuad) ,family=binomial(),
            data=data_surv_scaled, na.action="na.fail",nAGQ=0, control = glmerControl(calc.derivs = FALSE))
m03 <- glmer(Fate ~ Nb_Tot.x + Nb_Tot.x.sq + T.mean_m + T.mean_m.sq + T.range_m + Axis1 + Axis2 + (1 | Site) + (1 | SiteYear) + (1 | SiteQuad) ,family=binomial(),
             data=data_surv_scaled, na.action="na.fail",nAGQ=0, control = glmerControl(calc.derivs = FALSE))
AIC(m01,m02,m03)
r.squaredGLMM(m01)
r.squaredGLMM(m02)
r.squaredGLMM(m03) # little improvement, keep only SiteYear
rm(m01,m02,m03)
m0 <- glmer(Fate ~ Nb_Tot.x + Nb_Tot.x.sq + T.mean_m + T.mean_m.sq + T.range_m + Axis1 + Axis2 + (1 | Site) + (1 | SiteYear),family=binomial(),
             data=data_surv_scaled, na.action="na.fail",nAGQ=0, control = glmerControl(calc.derivs = FALSE))

# Random slopes for Sites
ms1 <- glmer(Fate ~ Nb_Tot.x + Nb_Tot.x.sq + T.mean_m + T.mean_m.sq + T.range_m + Axis1 + Axis2 + (1 | SiteYear) + (1 + Nb_Tot.x | Site),family=binomial(),
            data=data_surv_scaled, na.action="na.fail",nAGQ=0, control = glmerControl(calc.derivs = FALSE))
ms2 <- glmer(Fate ~ Nb_Tot.x + Nb_Tot.x.sq + T.mean_m + T.mean_m.sq + T.range_m + Axis1 + Axis2 + (1 | SiteYear) + (1 + Nb_Tot.x.sq | Site),family=binomial(),
             data=data_surv_scaled, na.action="na.fail",nAGQ=0, control = glmerControl(calc.derivs = FALSE))
ms3 <- glmer(Fate ~ Nb_Tot.x + Nb_Tot.x.sq + T.mean_m + T.mean_m.sq + T.range_m + Axis1 + Axis2 + (1 | SiteYear) + (1 + T.mean_m | Site),family=binomial(),
             data=data_surv_scaled, na.action="na.fail",nAGQ=0, control = glmerControl(calc.derivs = FALSE))
ms4 <- glmer(Fate ~ Nb_Tot.x + Nb_Tot.x.sq + T.mean_m + T.mean_m.sq + T.range_m + Axis1 + Axis2 + (1 | SiteYear) + (1 + T.mean_m.sq | Site),family=binomial(),
             data=data_surv_scaled, na.action="na.fail",nAGQ=0, control = glmerControl(calc.derivs = FALSE))
ms5 <- glmer(Fate ~ Nb_Tot.x + Nb_Tot.x.sq + T.mean_m + T.mean_m.sq + T.range_m + Axis1 + Axis2 + (1 | SiteYear) + (1 + T.range_m | Site),family=binomial(),
             data=data_surv_scaled, na.action="na.fail",nAGQ=0, control = glmerControl(calc.derivs = FALSE))
ms6 <- glmer(Fate ~ Nb_Tot.x + Nb_Tot.x.sq + T.mean_m + T.mean_m.sq + T.range_m + Axis1 + Axis2 + (1 | SiteYear) + (1 + Axis1 | Site),family=binomial(),
             data=data_surv_scaled, na.action="na.fail",nAGQ=0, control = glmerControl(calc.derivs = FALSE))
ms7 <- glmer(Fate ~ Nb_Tot.x + Nb_Tot.x.sq + T.mean_m + T.mean_m.sq + T.range_m + Axis1 + Axis2 + (1 | SiteYear) + (1 + Axis2 | Site),family=binomial(),
             data=data_surv_scaled, na.action="na.fail",nAGQ=0, control = glmerControl(calc.derivs = FALSE))

# Random slopes for SiteYear
msy1 <- glmer(Fate ~ Nb_Tot.x + Nb_Tot.x.sq + T.mean_m + T.mean_m.sq + T.range_m + Axis1 + Axis2 + (1 | Site) + (1 + Nb_Tot.x | SiteYear),family=binomial(),
             data=data_surv_scaled, na.action="na.fail",nAGQ=0, control = glmerControl(calc.derivs = FALSE))
msy2 <- glmer(Fate ~ Nb_Tot.x + Nb_Tot.x.sq + T.mean_m + T.mean_m.sq + T.range_m + Axis1 + Axis2 + (1 | Site) + (1 + Nb_Tot.x.sq | SiteYear),family=binomial(),
             data=data_surv_scaled, na.action="na.fail",nAGQ=0, control = glmerControl(calc.derivs = FALSE))
msy3 <- glmer(Fate ~ Nb_Tot.x + Nb_Tot.x.sq + T.mean_m + T.mean_m.sq + T.range_m + Axis1 + Axis2 + (1 | Site) + (1 + T.mean_m | SiteYear),family=binomial(),
             data=data_surv_scaled, na.action="na.fail",nAGQ=0, control = glmerControl(calc.derivs = FALSE))
msy4 <- glmer(Fate ~ Nb_Tot.x + Nb_Tot.x.sq + T.mean_m + T.mean_m.sq + T.range_m + Axis1 + Axis2 + (1 | Site) + (1 + T.mean_m.sq | SiteYear),family=binomial(),
             data=data_surv_scaled, na.action="na.fail",nAGQ=0, control = glmerControl(calc.derivs = FALSE))
msy5 <- glmer(Fate ~ Nb_Tot.x + Nb_Tot.x.sq + T.mean_m + T.mean_m.sq + T.range_m + Axis1 + Axis2 + (1 | Site) + (1 + T.range_m | SiteYear),family=binomial(),
             data=data_surv_scaled, na.action="na.fail",nAGQ=0, control = glmerControl(calc.derivs = FALSE))
msy6 <- glmer(Fate ~ Nb_Tot.x + Nb_Tot.x.sq + T.mean_m + T.mean_m.sq + T.range_m + Axis1 + Axis2 + (1 | Site) + (1 + Axis1 | SiteYear),family=binomial(),
             data=data_surv_scaled, na.action="na.fail",nAGQ=0, control = glmerControl(calc.derivs = FALSE))
msy7 <- glmer(Fate ~ Nb_Tot.x + Nb_Tot.x.sq + T.mean_m + T.mean_m.sq + T.range_m + Axis1 + Axis2 + (1 | Site) + (1 + Axis2 | SiteYear),family=binomial(),
             data=data_surv_scaled, na.action="na.fail",nAGQ=0, control = glmerControl(calc.derivs = FALSE))

( res.surv <- data.frame(Site=AIC(m0,ms1,ms2,ms3,ms4,ms5,ms6,ms7),
                       SiteYear=AIC(m0,msy1,msy2,msy3,msy4,msy5,msy6,msy7)) )

mb1 <- glmer(Fate ~ Nb_Tot.x + Nb_Tot.x.sq + T.mean_m + T.mean_m.sq + T.range_m + Axis1 + Axis2 + (1 | Site) + (1 + Axis1 + Axis2 | SiteYear),family=binomial(),
             data=data_surv_scaled, na.action="na.fail",nAGQ=0, control = glmerControl(calc.derivs = FALSE))
AIC(mb1) # not better than msy7
r.squaredGLMM(m0)
r.squaredGLMM(msy7) # good improvement ! keep msy7 !

# Try to dredge it
d.surv <- lapply(dredge(msy7, evaluate = FALSE), eval)
m.weight.surv <- calc.m.weight(d.surv)

# A note on dredge
# Do I need to remove the randomly varying slope for Axis2 when Axis2 is not included as a fixed effect ?
# After reading this post
# https://stats.stackexchange.com/questions/199377/is-it-reasonable-to-include-a-random-slope-term-in-an-lmer-model-without-the-cor
# I am convinced that i can keep it (note that the post makes the opposite statement, but assuming that the true fixed effect is not zero.
# In our case, the true effect seems to be 0, and model averaging ensures that bad models will not be used anyways).
# So I can dredge it in the usual way.

##

# # Random effect of SiteYear on Axis2
# expand.grid(Axis2 = seq(-1.5,2,0.1),levels(data_surv_scaled$SiteYear)) %>% as_tibble -> newdat 
# newdat %>% mutate(SiteYear = Var2) %>% mutate(Site = substr(SiteYear,1,3), Nb_Tot.x = 0, Nb_Tot.x.sq = 0, T.mean_m = 0, T.mean_m.sq = 0, T.range_m = 0, Axis1 = 0) -> newdat
# newdat$Site <- factor(newdat$Site,levels=levels(data_surv_scaled$Site))
# newdat$y <- predict(msy7,newdat,type="response")
# ggplot(newdat,aes(x=Axis2,y=y,col=SiteYear)) + geom_line()

##
##

# Analysis of growth
# Select only necessary columns and filter out NA
data_growth <-
    data_rep.clim %>%
    select(Site, SiteYear, SiteQuad, Nb_Tot.x, Nb_Tot.x.sq, Nb_Tot.y, Axis1, Axis2, T.mean_f, T.mean_f.sq, T.range_f) %>%
    na.omit()
# Scale variables and subtract 1 from dependent variable
data_growth_scaled <- data_growth %>% mutate_at(c("Nb_Tot.x","Nb_Tot.x.sq","T.mean_f","T.mean_f.sq","T.range_f", "Axis1", "Axis2"),scale)
data_growth_scaled$Nb_Tot.y_m1 <- data_growth_scaled$Nb_Tot.y - 1
data_growth_scaled$SiteYear <- factor(data_growth_scaled$SiteYear)
data_growth_scaled$SiteQuad <- factor(data_growth_scaled$SiteQuad)
# Fit
# No random slope
m01 <- glmer.nb(Nb_Tot.y_m1 ~ Nb_Tot.x + Nb_Tot.x.sq + T.mean_f + T.mean_f.sq + T.range_f + Axis1 + Axis2 + (1 | Site) + (1 | SiteYear),
              data=data_growth_scaled, na.action="na.fail", control = glmerControl(calc.derivs = FALSE), nAGQ=0)
m02 <- glmer.nb(Nb_Tot.y_m1 ~ Nb_Tot.x + Nb_Tot.x.sq + T.mean_f + T.mean_f.sq + T.range_f + Axis1 + Axis2 + (1 | Site) + (1 | SiteQuad),
               data=data_growth_scaled, na.action="na.fail", control = glmerControl(calc.derivs = FALSE), nAGQ=0)
m03 <- glmer.nb(Nb_Tot.y_m1 ~ Nb_Tot.x + Nb_Tot.x.sq + T.mean_f + T.mean_f.sq + T.range_f + Axis1 + Axis2 + (1 | Site) + (1 | SiteYear) + (1 | SiteQuad),
                data=data_growth_scaled, na.action="na.fail", control = glmerControl(calc.derivs = FALSE), nAGQ=0)
AIC(m01,m02,m03)
r.squaredGLMM(m01)
r.squaredGLMM(m02)
r.squaredGLMM(m03) # little improvement, keep only SiteYear
rm(m01,m02,m03) 
m0 <- glmer.nb(Nb_Tot.y_m1 ~ Nb_Tot.x + Nb_Tot.x.sq + T.mean_f + T.mean_f.sq + T.range_f + Axis1 + Axis2 + (1 | Site) + (1 | SiteYear),
                data=data_growth_scaled, na.action="na.fail", control = glmerControl(calc.derivs = FALSE), nAGQ=0)

# Random slopes for Sites
ms1 <- glmer.nb(Nb_Tot.y_m1 ~ Nb_Tot.x + Nb_Tot.x.sq + T.mean_f + T.mean_f.sq + T.range_f + Axis1 + Axis2 + (1 | SiteYear) + (1 + Nb_Tot.x | Site),
              data=data_growth_scaled, na.action="na.fail", control = glmerControl(calc.derivs = FALSE), nAGQ=0)
ms2 <- glmer.nb(Nb_Tot.y_m1 ~ Nb_Tot.x + Nb_Tot.x.sq + T.mean_f + T.mean_f.sq + T.range_f + Axis1 + Axis2 + (1 | SiteYear) + (1 + Nb_Tot.x.sq | Site),
                data=data_growth_scaled, na.action="na.fail", control = glmerControl(calc.derivs = FALSE), nAGQ=0)
ms3 <- glmer.nb(Nb_Tot.y_m1 ~ Nb_Tot.x + Nb_Tot.x.sq + T.mean_f + T.mean_f.sq + T.range_f + Axis1 + Axis2 + (1 | SiteYear) + (1 + T.mean_f | Site),
                data=data_growth_scaled, na.action="na.fail", control = glmerControl(calc.derivs = FALSE), nAGQ=0)
ms4 <- glmer.nb(Nb_Tot.y_m1 ~ Nb_Tot.x + Nb_Tot.x.sq + T.mean_f + T.mean_f.sq + T.range_f + Axis1 + Axis2 + (1 | SiteYear) + (1 + T.mean_f.sq | Site),
                data=data_growth_scaled, na.action="na.fail", control = glmerControl(calc.derivs = FALSE), nAGQ=0)
ms5 <- glmer.nb(Nb_Tot.y_m1 ~ Nb_Tot.x + Nb_Tot.x.sq + T.mean_f + T.mean_f.sq + T.range_f + Axis1 + Axis2 + (1 | SiteYear) + (1 + T.range_f | Site),
                data=data_growth_scaled, na.action="na.fail", control = glmerControl(calc.derivs = FALSE), nAGQ=0)
ms6 <- glmer.nb(Nb_Tot.y_m1 ~ Nb_Tot.x + Nb_Tot.x.sq + T.mean_f + T.mean_f.sq + T.range_f + Axis1 + Axis2 + (1 | SiteYear) + (1 + Axis1 | Site),
                data=data_growth_scaled, na.action="na.fail", control = glmerControl(calc.derivs = FALSE), nAGQ=0)
ms7 <- glmer.nb(Nb_Tot.y_m1 ~ Nb_Tot.x + Nb_Tot.x.sq + T.mean_f + T.mean_f.sq + T.range_f + Axis1 + Axis2 + (1 | SiteYear) + (1 + Axis2 | Site),
                data=data_growth_scaled, na.action="na.fail", control = glmerControl(calc.derivs = FALSE), nAGQ=0)

# Random slopes for SiteYears
msy1 <- glmer.nb(Nb_Tot.y_m1 ~ Nb_Tot.x + Nb_Tot.x.sq + T.mean_f + T.mean_f.sq + T.range_f + Axis1 + Axis2 + (1 | Site) + (1 + Nb_Tot.x | SiteYear),
                data=data_growth_scaled, na.action="na.fail", control = glmerControl(calc.derivs = FALSE), nAGQ=0)
msy2 <- glmer.nb(Nb_Tot.y_m1 ~ Nb_Tot.x + Nb_Tot.x.sq + T.mean_f + T.mean_f.sq + T.range_f + Axis1 + Axis2 + (1 | Site) + (1 + Nb_Tot.x.sq | SiteYear),
                data=data_growth_scaled, na.action="na.fail", control = glmerControl(calc.derivs = FALSE), nAGQ=0)
msy3 <- glmer.nb(Nb_Tot.y_m1 ~ Nb_Tot.x + Nb_Tot.x.sq + T.mean_f + T.mean_f.sq + T.range_f + Axis1 + Axis2 + (1 | Site) + (1 + T.mean_f | SiteYear),
                data=data_growth_scaled, na.action="na.fail", control = glmerControl(calc.derivs = FALSE), nAGQ=0)
msy4 <- glmer.nb(Nb_Tot.y_m1 ~ Nb_Tot.x + Nb_Tot.x.sq + T.mean_f + T.mean_f.sq + T.range_f + Axis1 + Axis2 + (1 | Site) + (1 + T.mean_f.sq | SiteYear),
                data=data_growth_scaled, na.action="na.fail", control = glmerControl(calc.derivs = FALSE), nAGQ=0)
msy5 <- glmer.nb(Nb_Tot.y_m1 ~ Nb_Tot.x + Nb_Tot.x.sq + T.mean_f + T.mean_f.sq + T.range_f + Axis1 + Axis2 + (1 | Site) + (1 + T.range_f | SiteYear),
                data=data_growth_scaled, na.action="na.fail", control = glmerControl(calc.derivs = FALSE), nAGQ=0)
msy6 <- glmer.nb(Nb_Tot.y_m1 ~ Nb_Tot.x + Nb_Tot.x.sq + T.mean_f + T.mean_f.sq + T.range_f + Axis1 + Axis2 + (1 | Site) + (1 + Axis1 | SiteYear),
                data=data_growth_scaled, na.action="na.fail", control = glmerControl(calc.derivs = FALSE), nAGQ=0)
msy7 <- glmer.nb(Nb_Tot.y_m1 ~ Nb_Tot.x + Nb_Tot.x.sq + T.mean_f + T.mean_f.sq + T.range_f + Axis1 + Axis2 + (1 | Site) + (1 + Axis2 | SiteYear),
                data=data_growth_scaled, na.action="na.fail", control = glmerControl(calc.derivs = FALSE), nAGQ=0)

( res.growth <- data.frame(Site=AIC(m0,ms1,ms2,ms3,ms4,ms5,ms6,ms7),
                       SiteYear=AIC(m0,msy1,msy2,msy3,msy4,msy5,msy6,msy7)) )

# msy1 seems the best
mb1 <- glmer.nb(Nb_Tot.y_m1 ~ Nb_Tot.x + Nb_Tot.x.sq + T.mean_f + T.mean_f.sq + T.range_f + Axis1 + Axis2 + (1 + Nb_Tot.x | Site) + (1 + Nb_Tot.x | SiteYear),
                 data=data_growth_scaled, na.action="na.fail", control = glmerControl(calc.derivs = FALSE), nAGQ=0)
AIC(mb1) # not better than msy1
mb2 <- glmer.nb(Nb_Tot.y_m1 ~ Nb_Tot.x + Nb_Tot.x.sq + T.mean_f + T.mean_f.sq + T.range_f + Axis1 + Axis2 + (1 | Site) + (1 + Nb_Tot.x + Nb_Tot.x.sq | SiteYear),
                data=data_growth_scaled, na.action="na.fail", control = glmerControl(calc.derivs = FALSE), nAGQ=0)
AIC(mb2) # not better than msy1
r.squaredGLMM(m0)
r.squaredGLMM(msy1) # little improvement, keep only random intercepts

# # Random effect of SiteYear on Nb_Tot.x
# expand.grid(Nb_Tot.x = seq(-0.7,5.8,0.1),levels(data_growth_scaled$SiteYear)) %>% as_tibble -> newdat 
# newdat %>% mutate(SiteYear = Var2) %>% mutate(Site = substr(SiteYear,1,3), Nb_Tot.x.sq = 0, T.mean_f = 0, T.mean_f.sq = 0, T.range_f = 0, Axis1 = 0, Axis2 = 0) -> newdat
# newdat$Site <- factor(newdat$Site,levels=levels(data_growth_scaled$Site))
# newdat$y <- predict(msy1,newdat)
# ggplot(newdat,aes(x=Nb_Tot.x,y=y,col=SiteYear)) + geom_line()



# Analysis of fruiting
# Select only necessary columns and filter out NA
data_fruiting <-
    data_rep.clim %>%
    select(Site, SiteYear, SiteQuad, Nb_Tot.x, Nb_Tot.x.sq, Fru, Axis1, Axis2, T.mean, T.mean.sq, T.range) %>%
    na.omit()
# Scale variables and save mean and sd
data_fruiting_scaled <- data_fruiting %>% mutate_at(c("Nb_Tot.x", "Nb_Tot.x.sq","T.mean","T.mean.sq","T.range", "Axis1", "Axis2"),scale)
data_fruiting_scaled$SiteYear <- factor(data_fruiting_scaled$SiteYear)
data_fruiting_scaled$SiteQuad <- factor(data_fruiting_scaled$SiteQuad)

# Fit
# No random slope. Testing plot
m01 <- glmer(Fru ~ Nb_Tot.x + Nb_Tot.x.sq + T.mean + T.mean.sq + T.range + Axis1 + Axis2 + (1 | Site) + (1| SiteYear), family=binomial(),
           data=data_fruiting_scaled, na.action="na.fail", nAGQ=0, control = glmerControl(calc.derivs = FALSE))
m02 <- glmer(Fru ~ Nb_Tot.x + Nb_Tot.x.sq + T.mean + T.mean.sq + T.range + Axis1 + Axis2 + (1 | Site) + (1| SiteQuad), family=binomial(),
             data=data_fruiting_scaled, na.action="na.fail", nAGQ=0, control = glmerControl(calc.derivs = FALSE))
m03 <- glmer(Fru ~ Nb_Tot.x + Nb_Tot.x.sq + T.mean + T.mean.sq + T.range + Axis1 + Axis2 + (1 | Site) + (1| SiteYear) + (1| SiteQuad), family=binomial(),
             data=data_fruiting_scaled, na.action="na.fail", nAGQ=0, control = glmerControl(calc.derivs = FALSE))
AIC(m01,m02,m03)
r.squaredGLMM(m01)
r.squaredGLMM(m02)
r.squaredGLMM(m03) # keep also SiteQuad 
rm(m01,m02,m03) 

m0 <- glmer(Fru ~ Nb_Tot.x + Nb_Tot.x.sq + T.mean + T.mean.sq + T.range + Axis1 + Axis2 + (1 | Site) + (1| SiteYear) + (1| SiteQuad), family=binomial(),
             data=data_fruiting_scaled, na.action="na.fail", nAGQ=0, control = glmerControl(calc.derivs = FALSE))

# Random slopes for Sites
ms1 <- glmer(Fru ~ Nb_Tot.x + Nb_Tot.x.sq + T.mean + T.mean.sq + T.range + Axis1 + Axis2 + (1 | SiteYear) + (1 + Nb_Tot.x | Site) + (1| SiteQuad), family=binomial(),
            data=data_fruiting_scaled, na.action="na.fail", nAGQ=0, control = glmerControl(calc.derivs = FALSE))
ms2 <- glmer(Fru ~ Nb_Tot.x + Nb_Tot.x.sq + T.mean + T.mean.sq + T.range + Axis1 + Axis2 + (1 | SiteYear) + (1 + Nb_Tot.x.sq | Site) + (1| SiteQuad), family=binomial(),
             data=data_fruiting_scaled, na.action="na.fail", nAGQ=0, control = glmerControl(calc.derivs = FALSE))
ms3 <- glmer(Fru ~ Nb_Tot.x + Nb_Tot.x.sq + T.mean + T.mean.sq + T.range + Axis1 + Axis2 + (1 | SiteYear) + (1 + T.mean | Site) + (1| SiteQuad), family=binomial(),
             data=data_fruiting_scaled, na.action="na.fail", nAGQ=0, control = glmerControl(calc.derivs = FALSE))
ms4 <- glmer(Fru ~ Nb_Tot.x + Nb_Tot.x.sq + T.mean + T.mean.sq + T.range + Axis1 + Axis2 + (1 | SiteYear) + (1 + T.mean.sq | Site) + (1| SiteQuad), family=binomial(),
             data=data_fruiting_scaled, na.action="na.fail", nAGQ=0, control = glmerControl(calc.derivs = FALSE))
ms5 <- glmer(Fru ~ Nb_Tot.x + Nb_Tot.x.sq + T.mean + T.mean.sq + T.range + Axis1 + Axis2 + (1 | SiteYear) + (1 + T.range | Site) + (1| SiteQuad), family=binomial(),
             data=data_fruiting_scaled, na.action="na.fail", nAGQ=0, control = glmerControl(calc.derivs = FALSE))
ms6 <- glmer(Fru ~ Nb_Tot.x + Nb_Tot.x.sq + T.mean + T.mean.sq + T.range + Axis1 + Axis2 + (1 | SiteYear) + (1 + Axis1 | Site) + (1| SiteQuad), family=binomial(),
             data=data_fruiting_scaled, na.action="na.fail", nAGQ=0, control = glmerControl(calc.derivs = FALSE))
ms7 <- glmer(Fru ~ Nb_Tot.x + Nb_Tot.x.sq + T.mean + T.mean.sq + T.range + Axis1 + Axis2 + (1 | SiteYear) + (1 + Axis2 | Site) + (1| SiteQuad), family=binomial(),
             data=data_fruiting_scaled, na.action="na.fail", nAGQ=0, control = glmerControl(calc.derivs = FALSE))

# Random slopes for SiteYears
msy1 <- glmer(Fru ~ Nb_Tot.x + Nb_Tot.x.sq + T.mean + T.mean.sq + T.range + Axis1 + Axis2 + (1 | Site) + (1 + Nb_Tot.x | SiteYear) + (1| SiteQuad), family=binomial(),
             data=data_fruiting_scaled, na.action="na.fail", nAGQ=0, control = glmerControl(calc.derivs = FALSE))
msy2 <- glmer(Fru ~ Nb_Tot.x + Nb_Tot.x.sq + T.mean + T.mean.sq + T.range + Axis1 + Axis2 + (1 | Site) + (1 + Nb_Tot.x.sq | SiteYear) + (1| SiteQuad), family=binomial(),
             data=data_fruiting_scaled, na.action="na.fail", nAGQ=0, control = glmerControl(calc.derivs = FALSE))
msy3 <- glmer(Fru ~ Nb_Tot.x + Nb_Tot.x.sq + T.mean + T.mean.sq + T.range + Axis1 + Axis2 + (1 | Site) + (1 + T.mean | SiteYear) + (1| SiteQuad), family=binomial(),
             data=data_fruiting_scaled, na.action="na.fail", nAGQ=0, control = glmerControl(calc.derivs = FALSE))
msy4 <- glmer(Fru ~ Nb_Tot.x + Nb_Tot.x.sq + T.mean + T.mean.sq + T.range + Axis1 + Axis2 + (1 | Site) + (1 + T.mean.sq | SiteYear) + (1| SiteQuad), family=binomial(),
             data=data_fruiting_scaled, na.action="na.fail", nAGQ=0, control = glmerControl(calc.derivs = FALSE))
msy5 <- glmer(Fru ~ Nb_Tot.x + Nb_Tot.x.sq + T.mean + T.mean.sq + T.range + Axis1 + Axis2 + (1 | Site) + (1 + T.range | SiteYear) + (1| SiteQuad), family=binomial(),
             data=data_fruiting_scaled, na.action="na.fail", nAGQ=0, control = glmerControl(calc.derivs = FALSE))
msy6 <- glmer(Fru ~ Nb_Tot.x + Nb_Tot.x.sq + T.mean + T.mean.sq + T.range + Axis1 + Axis2 + (1 | Site) + (1 + Axis1 | SiteYear) + (1| SiteQuad), family=binomial(),
             data=data_fruiting_scaled, na.action="na.fail", nAGQ=0, control = glmerControl(calc.derivs = FALSE))
msy7 <- glmer(Fru ~ Nb_Tot.x + Nb_Tot.x.sq + T.mean + T.mean.sq + T.range + Axis1 + Axis2 + (1 | Site) + (1 + Axis2 | SiteYear) + (1| SiteQuad), family=binomial(),
             data=data_fruiting_scaled, na.action="na.fail", nAGQ=0, control = glmerControl(calc.derivs = FALSE))

# Random slopes for SiteQuad
msq1 <- glmer(Fru ~ Nb_Tot.x + Nb_Tot.x.sq + T.mean + T.mean.sq + T.range + Axis1 + Axis2 + (1 | Site) + (1 | SiteYear) + (1 + Nb_Tot.x| SiteQuad), family=binomial(),
              data=data_fruiting_scaled, na.action="na.fail", nAGQ=0, control = glmerControl(calc.derivs = FALSE))
msq2 <- glmer(Fru ~ Nb_Tot.x + Nb_Tot.x.sq + T.mean + T.mean.sq + T.range + Axis1 + Axis2 + (1 | Site) + (1 | SiteYear) + (1 + Nb_Tot.x.sq | SiteQuad), family=binomial(),
              data=data_fruiting_scaled, na.action="na.fail", nAGQ=0, control = glmerControl(calc.derivs = FALSE))
msq3 <- glmer(Fru ~ Nb_Tot.x + Nb_Tot.x.sq + T.mean + T.mean.sq + T.range + Axis1 + Axis2 + (1 | Site) + (1 | SiteYear) + (1 + T.mean | SiteQuad), family=binomial(),
              data=data_fruiting_scaled, na.action="na.fail", nAGQ=0, control = glmerControl(calc.derivs = FALSE))
msq4 <- glmer(Fru ~ Nb_Tot.x + Nb_Tot.x.sq + T.mean + T.mean.sq + T.range + Axis1 + Axis2 + (1 | Site) + (1 | SiteYear) + (1 + T.mean.sq | SiteQuad), family=binomial(),
              data=data_fruiting_scaled, na.action="na.fail", nAGQ=0, control = glmerControl(calc.derivs = FALSE))
msq5 <- glmer(Fru ~ Nb_Tot.x + Nb_Tot.x.sq + T.mean + T.mean.sq + T.range + Axis1 + Axis2 + (1 | Site) + (1 | SiteYear) + (1 + T.range | SiteQuad), family=binomial(),
              data=data_fruiting_scaled, na.action="na.fail", nAGQ=0, control = glmerControl(calc.derivs = FALSE))
msq6 <- glmer(Fru ~ Nb_Tot.x + Nb_Tot.x.sq + T.mean + T.mean.sq + T.range + Axis1 + Axis2 + (1 | Site) + (1 | SiteYear) + (1 + Axis1 | SiteQuad), family=binomial(),
              data=data_fruiting_scaled, na.action="na.fail", nAGQ=0, control = glmerControl(calc.derivs = FALSE))
msq7 <- glmer(Fru ~ Nb_Tot.x + Nb_Tot.x.sq + T.mean + T.mean.sq + T.range + Axis1 + Axis2 + (1 | Site) + (1 | SiteYear) + (1 + Axis2 | SiteQuad), family=binomial(),
              data=data_fruiting_scaled, na.action="na.fail", nAGQ=0, control = glmerControl(calc.derivs = FALSE))


( res.fruiting <- data.frame(Site=AIC(m0,ms1,ms2,ms3,ms4,ms5,ms6,ms7),
                         SiteYear=AIC(m0,msy1,msy2,msy3,msy4,msy5,msy6,msy7),
                         SiteQuad=AIC(m0,msq1,msq2,msq3,msq4,msq5,msq6,msq7)) )
# Best seems msq1
mb1 <- glmer(Fru ~ Nb_Tot.x + Nb_Tot.x.sq + T.mean + T.mean.sq + T.range + Axis1 + Axis2 + (1 | Site) + (1 | SiteYear) + (1 + Nb_Tot.x + Nb_Tot.x.sq | SiteQuad), family=binomial(),
              data=data_fruiting_scaled, na.action="na.fail", nAGQ=0, control = glmerControl(calc.derivs = FALSE))
AIC(mb1) # not better than msq1
mb2 <- glmer(Fru ~ Nb_Tot.x + Nb_Tot.x.sq + T.mean + T.mean.sq + T.range + Axis1 + Axis2 + (1 | Site) + (1 | SiteYear) + (1 + Nb_Tot.x + T.mean.sq | SiteQuad), family=binomial(),
             data=data_fruiting_scaled, na.action="na.fail", nAGQ=0, control = glmerControl(calc.derivs = FALSE))
AIC(mb2) # not better than msq1
mb3 <- glmer(Fru ~ Nb_Tot.x + Nb_Tot.x.sq + T.mean + T.mean.sq + T.range + Axis1 + Axis2 + (1 | Site) + (1 | SiteYear) + (1 + Nb_Tot.x + T.mean | SiteQuad), family=binomial(),
             data=data_fruiting_scaled, na.action="na.fail", nAGQ=0, control = glmerControl(calc.derivs = FALSE))
AIC(mb3) # not better than msq1
mb4 <- glmer(Fru ~ Nb_Tot.x + Nb_Tot.x.sq + T.mean + T.mean.sq + T.range + Axis1 + Axis2 + (1 | Site) + (1 | SiteYear) + (1 + Nb_Tot.x + Axis2 | SiteQuad), family=binomial(),
             data=data_fruiting_scaled, na.action="na.fail", nAGQ=0, control = glmerControl(calc.derivs = FALSE))
AIC(mb4) # not better than msq1
r.squaredGLMM(m0)
r.squaredGLMM(msq1) # good improvement, keep msq1




# Analysis of fruits
# Select only necessary columns and filter out NA
data_fruits <-
    data_rep.clim %>%
    select(Site, SiteYear, SiteQuad, Nb_Tot.x, Nb_Tot.x.sq, Fru, logTotalSil, Axis1, Axis2, T.mean, T.mean.sq, T.range) %>%
    filter(Fru == 1) %>%
    na.omit()
# Scale variables and save mean and sd
data_fruits_scaled <- data_fruits %>% mutate_at(c("Nb_Tot.x", "Nb_Tot.x.sq","T.mean","T.mean.sq","T.range", "Axis1", "Axis2"),scale)
data_fruits_scaled$SiteYear <- factor(data_fruits_scaled$SiteYear)
data_fruits_scaled$SiteQuad <- factor(data_fruits_scaled$SiteQuad)
# Fit
# No random slopes
m01 <- lmer(logTotalSil ~ Nb_Tot.x + Nb_Tot.x.sq + T.mean + T.mean.sq + T.range + Axis1 + Axis2 + (1 | Site) + (1 | SiteYear),
           data=data_fruits_scaled, na.action="na.fail", REML=T, control = lmerControl(calc.derivs = FALSE))
m02 <- lmer(logTotalSil ~ Nb_Tot.x + Nb_Tot.x.sq + T.mean + T.mean.sq + T.range + Axis1 + Axis2 + (1 | Site) + (1 | SiteQuad),
            data=data_fruits_scaled, na.action="na.fail", REML=T, control = lmerControl(calc.derivs = FALSE))
m03 <- lmer(logTotalSil ~ Nb_Tot.x + Nb_Tot.x.sq + T.mean + T.mean.sq + T.range + Axis1 + Axis2 + (1 | Site) + (1 | SiteQuad) + (1 | SiteYear),
          data=data_fruits_scaled, na.action="na.fail", REML=T, control = lmerControl(calc.derivs = FALSE))
AIC(m01,m02,m03)
r.squaredGLMM(m01)
r.squaredGLMM(m02)
r.squaredGLMM(m03) # little improvement, keep only SiteYear
rm(m01,m02,m03)
m0 <- lmer(logTotalSil ~ Nb_Tot.x + Nb_Tot.x.sq + T.mean + T.mean.sq + T.range + Axis1 + Axis2 + (1 | Site) + (1 | SiteYear),
            data=data_fruits_scaled, na.action="na.fail", REML=T, control = lmerControl(calc.derivs = FALSE))

# Random slopes for Sites
ms1 <- lmer(logTotalSil ~ Nb_Tot.x + Nb_Tot.x.sq + T.mean + T.mean.sq + T.range + Axis1 + Axis2 + (1 + Nb_Tot.x | Site)  + (1 | SiteYear),
           data=data_fruits_scaled, na.action="na.fail", REML=T, control = lmerControl(calc.derivs = FALSE))
ms2 <- lmer(logTotalSil ~ Nb_Tot.x + Nb_Tot.x.sq + T.mean + T.mean.sq + T.range + Axis1 + Axis2 + (1 + Nb_Tot.x.sq | Site) + (1 | SiteYear),
            data=data_fruits_scaled, na.action="na.fail", REML=T, control = lmerControl(calc.derivs = FALSE))
ms3 <- lmer(logTotalSil ~ Nb_Tot.x + Nb_Tot.x.sq + T.mean + T.mean.sq + T.range + Axis1 + Axis2 + (1 + T.mean | Site) + (1 | SiteYear),
            data=data_fruits_scaled, na.action="na.fail", REML=T, control = lmerControl(calc.derivs = FALSE))
ms4 <- lmer(logTotalSil ~ Nb_Tot.x + Nb_Tot.x.sq + T.mean + T.mean.sq + T.range + Axis1 + Axis2 + (1 + T.mean.sq | Site) + (1 | SiteYear),
            data=data_fruits_scaled, na.action="na.fail", REML=T, control = lmerControl(calc.derivs = FALSE))
ms5 <- lmer(logTotalSil ~ Nb_Tot.x + Nb_Tot.x.sq + T.mean + T.mean.sq + T.range + Axis1 + Axis2 + (1 + T.range | Site) + (1 | SiteYear),
            data=data_fruits_scaled, na.action="na.fail", REML=T, control = lmerControl(calc.derivs = FALSE))
ms6 <- lmer(logTotalSil ~ Nb_Tot.x + Nb_Tot.x.sq + T.mean + T.mean.sq + T.range + Axis1 + Axis2 + (1 + Axis1 | Site) + (1 | SiteYear),
            data=data_fruits_scaled, na.action="na.fail", REML=T, control = lmerControl(calc.derivs = FALSE))
ms7 <- lmer(logTotalSil ~ Nb_Tot.x + Nb_Tot.x.sq + T.mean + T.mean.sq + T.range + Axis1 + Axis2 + (1 + Axis2 | Site) + (1 | SiteYear),
            data=data_fruits_scaled, na.action="na.fail", REML=T, control = lmerControl(calc.derivs = FALSE))

# Random slopes for SiteYear
msy1 <- lmer(logTotalSil ~ Nb_Tot.x + Nb_Tot.x.sq + T.mean + T.mean.sq + T.range + Axis1 + Axis2 + (1  | Site) + (1 + Nb_Tot.x| SiteYear),
            data=data_fruits_scaled, na.action="na.fail", REML=T, control = lmerControl(calc.derivs = FALSE))
msy2 <- lmer(logTotalSil ~ Nb_Tot.x + Nb_Tot.x.sq + T.mean + T.mean.sq + T.range + Axis1 + Axis2 + (1 | Site) + (1 + Nb_Tot.x.sq | SiteYear),
            data=data_fruits_scaled, na.action="na.fail", REML=T, control = lmerControl(calc.derivs = FALSE))
msy3 <- lmer(logTotalSil ~ Nb_Tot.x + Nb_Tot.x.sq + T.mean + T.mean.sq + T.range + Axis1 + Axis2 + (1 | Site) + (1 + T.mean | SiteYear),
            data=data_fruits_scaled, na.action="na.fail", REML=T, control = lmerControl(calc.derivs = FALSE))
msy4 <- lmer(logTotalSil ~ Nb_Tot.x + Nb_Tot.x.sq + T.mean + T.mean.sq + T.range + Axis1 + Axis2 + (1 | Site) + (1 + T.mean.sq | SiteYear),
            data=data_fruits_scaled, na.action="na.fail", REML=T, control = lmerControl(calc.derivs = FALSE))
msy5 <- lmer(logTotalSil ~ Nb_Tot.x + Nb_Tot.x.sq + T.mean + T.mean.sq + T.range + Axis1 + Axis2 + (1 | Site) + (1 + T.range | SiteYear),
            data=data_fruits_scaled, na.action="na.fail", REML=T, control = lmerControl(calc.derivs = FALSE))
msy6 <- lmer(logTotalSil ~ Nb_Tot.x + Nb_Tot.x.sq + T.mean + T.mean.sq + T.range + Axis1 + Axis2 + (1 | Site) + (1 + Axis1 | SiteYear),
            data=data_fruits_scaled, na.action="na.fail", REML=T, control = lmerControl(calc.derivs = FALSE))
msy7 <- lmer(logTotalSil ~ Nb_Tot.x + Nb_Tot.x.sq + T.mean + T.mean.sq + T.range + Axis1 + Axis2 + (1 | Site) + (1 + Axis2 | SiteYear),
            data=data_fruits_scaled, na.action="na.fail", REML=T, control = lmerControl(calc.derivs = FALSE))

( res.fruits <- data.frame(Site=AIC(m0,ms1,ms2,ms3,ms4,ms5,ms6,ms7),
                           SiteYear=AIC(m0,msy1,msy2,msy3,msy4,msy5,msy6,msy7)) )

# ms1 seems the best 

mb1 <- lmer(logTotalSil ~ Nb_Tot.x + Nb_Tot.x.sq + T.mean + T.mean.sq + T.range + Axis1 + Axis2 + (1 + Nb_Tot.x + Nb_Tot.x.sq | Site) + (1 | SiteYear),
            data=data_fruits_scaled, na.action="na.fail", REML=T, control = lmerControl(calc.derivs = FALSE))
AIC(mb1) # not better than ms1
mb2 <- lmer(logTotalSil ~ Nb_Tot.x + Nb_Tot.x.sq + T.mean + T.mean.sq + T.range + Axis1 + Axis2 + (1 + Nb_Tot.x | Site) + (1 + Nb_Tot.x | SiteYear),
            data=data_fruits_scaled, na.action="na.fail", REML=T, control = lmerControl(calc.derivs = FALSE))
AIC(mb2) # not better than ms1
# So ms1 is the best
r.squaredGLMM(m0)
r.squaredGLMM(ms1) # little improvement, keep m0



# Analysis of recruit size
# Select only necessary columns and filter out large plants and NA
data_recruit <-
    data_rep.clim %>%
    select(Site, SiteYear, SiteQuad, is.recr, First_stade, Nb_Tot.x, logNbtot.x, Axis1, Axis2, T.mean, T.mean.sq, T.range) %>%
    filter(is.recr == 1) %>%
    filter(First_stade == 0) %>%
    filter(Nb_Tot.x <= 5) %>%
    na.omit()
# Scale variables and save mean and sd, subtract 1 from dependent variable
data_recruit_scaled <- data_recruit %>% mutate_at(c("T.mean","T.mean.sq","T.range", "Axis1", "Axis2"),scale)
data_recruit_scaled$Nb_Tot.x_m1 <- data_recruit_scaled$Nb_Tot.x - 1
data_recruit_scaled$SiteYear <- factor(data_recruit_scaled$SiteYear)
data_recruit_scaled$SiteQuad <- factor(data_recruit_scaled$SiteQuad)

# Fit
# No random slopes, testing plot
m01 <- glmer.nb(Nb_Tot.x_m1 ~ T.mean + T.mean.sq + T.range + Axis1 + Axis2 + (1 | Site) + (1 | SiteYear),
              data=data_recruit_scaled, na.action="na.fail", control = glmerControl(calc.derivs = FALSE), nAGQ=0)
m02 <- glmer.nb(Nb_Tot.x_m1 ~ T.mean + T.mean.sq + T.range + Axis1 + Axis2 + (1 | Site) + (1 | SiteQuad),
                data=data_recruit_scaled, na.action="na.fail", control = glmerControl(calc.derivs = FALSE), nAGQ=0)
m03 <- glmer.nb(Nb_Tot.x_m1 ~ T.mean + T.mean.sq + T.range + Axis1 + Axis2 + (1 | Site) + (1 | SiteYear) + (1 | SiteQuad),
                data=data_recruit_scaled, na.action="na.fail", control = glmerControl(calc.derivs = FALSE), nAGQ=0)
AIC(m01,m02,m03)
r.squaredGLMM(m01)
r.squaredGLMM(m02)
r.squaredGLMM(m03) # Keep only SiteYear
rm(m01,m02,m03)

m0 <- glmer.nb(Nb_Tot.x_m1 ~ T.mean + T.mean.sq + T.range + Axis1 + Axis2 + (1 | Site) + (1 | SiteYear),
                data=data_recruit_scaled, na.action="na.fail", control = glmerControl(calc.derivs = FALSE), nAGQ=0)

# Random slope for Site
ms1 <- glmer.nb(Nb_Tot.x_m1 ~ T.mean + T.mean.sq + T.range + Axis1 + Axis2 + (1 + T.mean | Site) + (1 | SiteYear),
               data=data_recruit_scaled, na.action="na.fail", control = glmerControl(calc.derivs = FALSE), nAGQ=0)
ms2 <- glmer.nb(Nb_Tot.x_m1 ~ T.mean + T.mean.sq + T.range + Axis1 + Axis2 + (1 + T.mean.sq | Site) + (1 | SiteYear),
                data=data_recruit_scaled, na.action="na.fail", control = glmerControl(calc.derivs = FALSE), nAGQ=0)
ms3 <- glmer.nb(Nb_Tot.x_m1 ~ T.mean + T.mean.sq + T.range + Axis1 + Axis2 + (1 + T.range | Site) + (1 | SiteYear),
                data=data_recruit_scaled, na.action="na.fail", control = glmerControl(calc.derivs = FALSE), nAGQ=0)
ms4 <- glmer.nb(Nb_Tot.x_m1 ~ T.mean + T.mean.sq + T.range + Axis1 + Axis2 + (1 + Axis1 | Site) + (1 | SiteYear),
                data=data_recruit_scaled, na.action="na.fail", control = glmerControl(calc.derivs = FALSE), nAGQ=0)
ms5 <- glmer.nb(Nb_Tot.x_m1 ~ T.mean + T.mean.sq + T.range + Axis1 + Axis2 + (1 + Axis2 | Site) + (1 | SiteYear),
                data=data_recruit_scaled, na.action="na.fail", control = glmerControl(calc.derivs = FALSE), nAGQ=0)

# Random slope for SiteYear
msy1 <- glmer.nb(Nb_Tot.x_m1 ~ T.mean + T.mean.sq + T.range + Axis1 + Axis2 + (1 | Site) + (1 + T.mean | SiteYear),
                data=data_recruit_scaled, na.action="na.fail", control = glmerControl(calc.derivs = FALSE), nAGQ=0)
msy2 <- glmer.nb(Nb_Tot.x_m1 ~ T.mean + T.mean.sq + T.range + Axis1 + Axis2 + (1 | Site) + (1 + T.mean.sq | SiteYear),
                data=data_recruit_scaled, na.action="na.fail", control = glmerControl(calc.derivs = FALSE), nAGQ=0)
msy3 <- glmer.nb(Nb_Tot.x_m1 ~ T.mean + T.mean.sq + T.range + Axis1 + Axis2 + (1 | Site) + (1 + T.range | SiteYear),
                data=data_recruit_scaled, na.action="na.fail", control = glmerControl(calc.derivs = FALSE), nAGQ=0)
msy4 <- glmer.nb(Nb_Tot.x_m1 ~ T.mean + T.mean.sq + T.range + Axis1 + Axis2 + (1 | Site) + (1 + Axis1 | SiteYear),
                data=data_recruit_scaled, na.action="na.fail", control = glmerControl(calc.derivs = FALSE), nAGQ=0)
msy5 <- glmer.nb(Nb_Tot.x_m1 ~ T.mean + T.mean.sq + T.range + Axis1 + Axis2 + (1 | Site) + (1 + Axis2 | SiteYear),
                data=data_recruit_scaled, na.action="na.fail", control = glmerControl(calc.derivs = FALSE), nAGQ=0)

( res.recruits <- data.frame(Site=AIC(m0,ms1,ms2,ms3,ms4,ms5),
                           SiteYear=AIC(m0,msy1,msy2,msy3,msy4,msy5)) )

 # Keep only random intercept !

save(res.surv, res.growth, res.fruiting, res.fruits, res.recruits, file="Res_TestRandomFactors.RData")
