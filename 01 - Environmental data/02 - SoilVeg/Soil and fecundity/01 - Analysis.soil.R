rm(list=ls())
library(MCMCglmm)

dataCN <- read.csv("CN.csv",h=T,sep=",")

# Nitrogen 
hist(dataCN$N)
summary(dataCN$N)
m1 <- MCMCglmm(N ~ 1,data=dataCN,random= ~site + SiteQuad,
               nitt=101000, thin=50, burnin=1000,verbose=F)
summary(m1)
autocorr.diag(m1$VCV)
plot(m1$VCV)

vSite <- m1$VCV[, "site"]/(m1$VCV[, "SiteQuad"] + m1$VCV[, "site"] + m1$VCV[, "units"])
mean(vSite)
vSiteQuad <- m1$VCV[, "SiteQuad"]/(m1$VCV[, "SiteQuad"] + m1$VCV[, "site"] + m1$VCV[, "units"])
mean(vSiteQuad)
vUnits <- m1$VCV[, "units"]/(m1$VCV[, "SiteQuad"] + m1$VCV[, "site"] + m1$VCV[, "units"])
mean(vUnits)



# Carbon
hist(dataCN$C)
summary(dataCN$C)
m1 <- MCMCglmm(C ~ 1,data=dataCN,random= ~site + SiteQuad,
               nitt=101000, thin=50, burnin=1000,verbose=F)
summary(m1)
autocorr.diag(m1$VCV)
plot(m1$VCV)

vSite <- m1$VCV[, "site"]/(m1$VCV[, "SiteQuad"] + m1$VCV[, "site"] + m1$VCV[, "units"])
mean(vSite)
vSiteQuad <- m1$VCV[, "SiteQuad"]/(m1$VCV[, "SiteQuad"] + m1$VCV[, "site"] + m1$VCV[, "units"])
mean(vSiteQuad)
vUnits <- m1$VCV[, "units"]/(m1$VCV[, "SiteQuad"] + m1$VCV[, "site"] + m1$VCV[, "units"])
mean(vUnits)

cor.test(dataCN$C,dataCN$N,method="sp")

# pH
datapH <- read.csv("pH.csv",h=T,sep=",")
hist(datapH$pH)
summary(datapH$pH)

m1 <- MCMCglmm(pH ~ 1,data=datapH,random= ~Site + SiteQuad,
               nitt=101000, thin=50, burnin=1000,verbose=F)
autocorr.diag(m1$VCV)
plot(m1$VCV)
summary(m1)

vSite <- m1$VCV[, "Site"]/(m1$VCV[, "SiteQuad"] + m1$VCV[, "Site"] + m1$VCV[, "units"])
mean(vSite)
vSiteQuad <- m1$VCV[, "SiteQuad"]/(m1$VCV[, "SiteQuad"] + m1$VCV[, "Site"] + m1$VCV[, "units"])
mean(vSiteQuad)
vUnits <- m1$VCV[, "units"]/(m1$VCV[, "SiteQuad"] + m1$VCV[, "Site"] + m1$VCV[, "units"])
mean(vUnits)



# Boxplots
par(mfrow=c(2,2))
boxplot(C~site,dataCN,main="Soil: carbon content")
boxplot(N~site,dataCN,main="Soil: nitrogen content")
boxplot(C/N~site,dataCN,main="Soil: carbon/nitrogen ratio")
boxplot(pH~site,datapH,ylim=c(4,8),ylab="pH",main="Soil acidity")
