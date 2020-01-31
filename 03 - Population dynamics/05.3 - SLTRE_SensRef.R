# SLTRE for the Seed Bank analysis
# Marco Andrello
# Last modification: 23/01/2020

rm(list=ls())

library(tidyverse)
library(gridExtra)
library(RColorBrewer)
library(abind)

v.scen <- c(1,5)

v.site <- c("BRU", "CHA", "VIL", "LAU", "GAL", "PIC")
v.year <- c(2008:2013)
elev <- c(930,1480,1980,2090,2590,2930)
names(elev) <- v.site 

# Set themes for plotting
theme_set(theme_classic())
theme_update(axis.text.x = element_text(size=7, color="black"),
             axis.text.y = element_text(size=7, color="black"),
             axis.title.x = element_text(size=7, color="black"),
             axis.title.y = element_text(size=7, color="black"),
             panel.spacing.x = unit(0.45, "cm"),
             plot.margin = margin(0,0.5,0,0, unit="cm"),
             strip.text.x = element_text(size = 6.5),
             strip.background.x = element_rect(linetype = "blank"),
             legend.position = "none")

png("Figure_2a_SensRef.png",width=16,height=8,units="cm",res=300)
par(xpd=NA,cex.axis=0.7, cex.main=0.8, mar=c(3,3,2,1),mfrow=c(1,2))
pal <- brewer.pal(4,"Accent")
for (i.scen in 1 : length(v.scen)) {
  load(paste0("Results_MPM_SensRef_scen_",v.scen[i.scen],"_.RData"))
  cnet <- apply(cnet,c(1,2,3),mean)
  C_k <- array(NA,c(4,length(v.site))) # Total effects of each descriptor
  dimnames(C_k) <- list(dimnames(cnet)[[2]],dimnames(cnet)[[3]])
  for (i.descr in 1 : 4) {
    C_k[i.descr,] <- colSums(abs(cnet[,i.descr,]))
  }
  barplot(C_k,col=pal,border=NA,las=1,main=paste0("Ref = ",v.site[v.scen[i.scen]]))
}
dev.off()

png("Figure_2b_SensRef.png",width=16,height=8,units="cm",res=300)
par(xpd=NA,cex.axis=0.7, cex.main=0.8, mar=c(3,3,2,1),mfrow=c(1,2))
pal <- brewer.pal(7,"Set1")
for (i.scen in 1 : length(v.scen)) {
  load(paste0("Results_MPM_SensRef_scen_",v.scen[i.scen],"_.RData"))
  cnet <- apply(cnet,c(1,2,3),mean)
  C_l <- array(NA,c(7,length(v.site))) # Total effects of each life-cycle component
  dimnames(C_l) <- list(dimnames(cnet)[[1]],dimnames(cnet)[[3]])
  for (i.lcc in 1 : 7) {
    C_l[i.lcc,] <- colSums(abs(cnet[i.lcc,,]))
  }
  barplot(C_l,col=pal,border=NA,las=1,main=paste0("Ref = ",v.site[v.scen[i.scen]]))
}
dev.off()


