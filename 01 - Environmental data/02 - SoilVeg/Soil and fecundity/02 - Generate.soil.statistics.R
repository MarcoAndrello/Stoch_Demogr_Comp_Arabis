rm(list=ls())
library(MCMCglmm)

dataCN <- read.csv("CN.csv",h=T,sep=",")
dataCN$CNratio <- dataCN$C / dataCN$N

datapH <- read.csv("pH.csv",h=T,sep=",")

all(levels(dataCN$SiteQuad) == levels(datapH$SiteQuad))

soil.stat <- data.frame(SiteQuad=levels(dataCN$SiteQuad))
soil.stat$N_mean <- as.vector(by(dataCN$N,dataCN$SiteQuad,mean))
soil.stat$C_mean <- as.vector(by(dataCN$C,dataCN$SiteQuad,mean))
soil.stat$CNratio_mean <- as.vector(by(dataCN$CNratio,dataCN$SiteQuad,mean))
soil.stat$pH <- as.vector(by(datapH$pH,datapH$SiteQuad,mean))

soil.stat$Site <- factor(substr(soil.stat$SiteQuad,1,3))

save(soil.stat,file="Soil.statistics.RData")

