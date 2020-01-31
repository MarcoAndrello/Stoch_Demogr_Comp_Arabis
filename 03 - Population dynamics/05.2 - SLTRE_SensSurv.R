# SLTRE, SensSurv
# Marco Andrello
# Last modification: 13/01/2020

rm(list=ls())

library(tidyverse)
library(gridExtra)
library(RColorBrewer)
library(abind)

load("Results_MPM_SensSurv.RData")

v.site <- c("BRU", "CHA", "VIL", "LAU", "GAL", "PIC")
v.year <- c(2008:2013)
v.stage <- c(1,5,10,15)
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



# Take the mean of cnet
cnet_2000 <- cnet
cnet <- apply(cnet_2000,c(1,2,3),mean)

# Total effects of each descriptor..................
# By summing the absolute values of its contributions over the seven life-cycle components in each site
C_k <- array(NA,c(4,length(v.site)))
dimnames(C_k) <- list(dimnames(cnet)[[2]],dimnames(cnet)[[3]])
for (i.descr in 1 : 4) {
    C_k[i.descr,] <- colSums(abs(cnet[,i.descr,]))
}

# Calculating the percent effect
C_k.perc <- apply(C_k, 2, function(x) {x/sum(x)*100})
C_k.perc
colSums(C_k.perc[2:4,])

# Plot Figure
png("Figure_2a_SensSurv.png",width=8,height=8,units="cm",res=300)
par(xpd=NA,cex=0.7)
pal <- brewer.pal(4,"Accent")
barplot(C_k,col=pal,border=NA,las=1)
legend(0,-20,c(expression(italic("\u03BC")),expression(italic(e)),
               expression(italic(CV)),expression(rho)),fill=pal,bty="n",ncol=4)
dev.off()


# Total effects of each life-cycle component..............
# By summing the absolute values of its contributions over the four descriptors in each site
C_l <- array(NA,c(7,length(v.site)))
dimnames(C_l) <- list(dimnames(cnet)[[1]],dimnames(cnet)[[3]])
for (i.lcc in 1 : 7) {
    C_l[i.lcc,] <- colSums(abs(cnet[i.lcc,,]))
}

# Calculating the percent effect
C_l.perc <- apply(C_l, 2, function(x) {x/sum(x)*100})

# Plot Figure
png("Figure_2b_SensSurv.png",width=8,height=8,units="cm",res=300)
par(xpd=NA,cex=0.7)
pal <- brewer.pal(7,"Set1")
barplot(C_l,col=pal,border=NA,las=1)
legend(0,-15,c("S","G-","G=","G+","F0","F1","F2"),fill=pal,bty="n",ncol=4)
dev.off()


# Site-specific plots....................................
r.plot <- r %>% apply(1,mean,na.rm=T)
v.site.elev <- c("BRU (930m)", "CHA (1480m)", "VIL (1980m)", "LAU (2090m)", "GAL (2590m)", "PIC (2930m)")
png("Fig_2_BySite_SensSurv.png",width=15,height=9,units="cm",res=300)
pal <- brewer.pal(4, "Accent")
par(mfrow=c(2,3), mar=c(3,3,2,1), cex.axis=0.75)
for (i.site in 1 : length(v.site)) {
  c1 <- c2 <- c3 <- cnet[,,i.site]
  c2[c2>0] <- 0 
  c3[c3<0] <- 0 
  myRange <- c( min(rowSums(c2)), max(rowSums(c3)) )
  myRange <- c(-0.25,0.45)
  barplot(t(c2),col=pal,space=0,border=NA,axisnames=FALSE, ylim=myRange,main=v.site.elev[i.site],las=1)
  barplot(t(c3),col=pal,space=0,border=NA,axisnames=FALSE, add=T, axes=F)
  axis(1, at=c(0:6)+0.5, lab= row.names(cnet), las=1)
  text(5.5,0.4,bquote(log~lambda[s]*''==''*.(round(r.plot[i.site],1))))
}
dev.off()


