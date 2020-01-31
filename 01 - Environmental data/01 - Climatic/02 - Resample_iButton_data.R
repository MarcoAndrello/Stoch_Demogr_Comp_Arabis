# Generate resampled iButton data
# Marco Andrello
# 05/06/2019

rm(list=ls())

library(tidyverse)
library(lubridate)
library(hms)
library(magrittr)
library(broom)


load("reconstruct_consensus.RData")

id1 <- which(reconstruct$MeanTemp_Type=="Predicted")
id2 <- which(reconstruct$MinTemp_Type=="Predicted")
id3 <- which(reconstruct$MaxTemp_Type=="Predicted")
id4 <- which(reconstruct$AmpTemp_Type=="Predicted")
id5 <- which(reconstruct$Freeze_Type=="Predicted")

for (rep in 1 : 50) {
  
  cat("Rep",rep,"\n"); flush.console()
  
  sampled <- reconstruct
  sampled$MeanTemp[id1] <- rnorm(length(id1),sampled$MeanTemp_Pred[id1],sampled$MeanTemp_PredSE[id1])
  sampled$MinTemp[id2] <- rnorm(length(id2),sampled$MinTemp_Pred[id2],sampled$MinTemp_PredSE[id2])
  sampled$MaxTemp[id3] <- rnorm(length(id3),sampled$MaxTemp_Pred[id3],sampled$MaxTemp_PredSE[id3])
  sampled$AmpTemp[id4] <- rnorm(length(id4),sampled$AmpTemp_Pred[id4],sampled$AmpTemp_PredSE[id4])
  sampled$Freeze[id5] <- rbinom(length(id5),1,sampled$Prob_Freeze[id5])
  
  save(sampled,file=paste0(getwd(),"/01 - Resampled data/01 - Daily/data.daily.",rep,".RData"))
  
  # Generate the monthly estimates
  monthly <-
    sampled %>%
    group_by(Pop, Quad, Year) %>%
    summarise(AvgMeanTemp  = mean(MeanTemp, na.rm = TRUE),
              AvgMinTemp  = mean(MinTemp, na.rm = TRUE),
              AvgMaxTemp  = mean(MaxTemp, na.rm = TRUE),
              AvgAmp    = mean(AmpTemp, na.rm = TRUE),
              NbFreeze  = sum(Freeze, na.rm = TRUE))
  
  save(monthly,file=paste0(getwd(),"/01 - Resampled data/02 - Monthly/data.monthly.",rep,".RData"))
  
}


