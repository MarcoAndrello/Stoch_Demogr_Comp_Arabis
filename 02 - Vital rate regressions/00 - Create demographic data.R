rm(list=ls())

library(tidyverse)

load(paste0(getwd(),"/demographic_datasets/Data.growth.global.2017_02_01.RData"))
data.growth <- as_tibble(data)
load(paste0(getwd(),"/demographic_datasets/Data.survival.global.2017_01_31.RData"))
data.surv <- as_tibble(data)
load(paste0(getwd(),"/demographic_datasets/Data.fruiting.global.2017_07_24.RData"))
data.fruiting <- as_tibble(data)
load(paste0(getwd(),"/demographic_datasets/Data.recruits.size.2017_02_06.RData"))
data.recruit <- as_tibble(data)
data.recruit$Year <- factor(data.recruit$Year)
rm(data)

# Checking
anti_join(data.recruit,data.surv)
# All records in recr are also appearing in surv
# Attempting join
data <- 
  left_join(data.surv,data.recruit) %>%
  left_join(data.growth) %>%
  left_join(data.fruiting)


 # Check joins
filter(data,!is.na(data$Nb_Tot)) %>% mutate(comp= Nb_Tot- Nb_Tot.x) %>% summary() # Nb_Tot and Nb_Tot.x are the same
filter(data,is.na(data$Nb_Tot.y)) %>% filter(.$Fate == 1) # These are "missing" individuals
filter(data, Year != 2008) %>% filter(is.na(.$Total_Sil))

data <- 
  data %>% 
  add_column(Quad=substr(.[["ID"]],5,6)) %>%
  mutate(is.recr = as.numeric(!is.na(Nb_Tot))) %>%
  select(Site, Year, Quad, ID, Nb_Tot.x, Fru, Total_Sil, Fate, Nb_Tot.y, is.recr)

data$ID1 <- paste0(data$ID,"-",data$Year)

# Import information on First_stade
dat <- read.csv(paste0(getwd(),"/demographic_datasets/data corrigé2.csv"))
dat %>% as_tibble -> dat

filter(dat,!is.na(First_stade)) %>%
    mutate(ID1 = paste(Site, Quad, Ind, Year, sep="-")) ->
    datfilt

# Missing individuals in data (recruits of 2014)
recr2014 <- 
anti_join(datfilt, data, by="ID1") %>%
    mutate(ID=paste(Site,Quad,Ind,sep="-"), Nb_Tot.x=Nb_Tot, Fru=Total_Sil,
           Fate=NA, Nb_Tot.y=NA, is.recr=1, ID1=paste(ID,Year,sep="-")) %>%
    select(Site, Year, Quad, ID, Nb_Tot.x, Fru, Total_Sil, Fate, Nb_Tot.y, is.recr, ID1) 
recr2014$Fru[recr2014$Fru>0] <- 1
# Add them to data
data <- rbind(data,recr2014)
# Join datafilt (containing First_stade) to data
left_join(data, select(datfilt, ID1, First_stade), by="ID1") -> data2
data2 %>% filter(!is.na(First_stade)) %>% print(n=250)

data <- data2
save(data,file="Demographic_data.RData")
