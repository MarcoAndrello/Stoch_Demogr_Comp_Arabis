# SLTRE
# Marco Andrello
# Last modification: 25/12/2019

rm(list=ls())

library(tidyverse)
library(gridExtra)
library(RColorBrewer)
library(abind)

load("Results_MPM.RData")

v.site <- c("BRU", "CHA", "VIL", "LAU", "GAL", "PIC")
v.year <- c(2008:2013)
v.stage <- c(1,5,10,15)
elev <- c(930,1480,1980,2090,2590,2930)
names(elev) <- v.site 

# Set themes for plotting
theme_set(theme_classic())
theme_update(axis.text.x = element_text(size=7, color="black"),
             axis.text.y = element_text(size=7, color="black"),
             axis.title.x = element_text(size=7, color="black", margin = margin(0.5,0,0,0, unit="cm")),
             axis.title.y = element_text(size=7, color="black"),
             panel.spacing.x = unit(0.45, "cm"),
             plot.margin = margin(0,0.5,0,0, unit="cm"),
             strip.text.x = element_text(colour="white"),
             strip.background.x = element_rect(linetype = "blank"),
             legend.position = "none")


# Take the mean of cnet
cnet_2000 <- cnet
cnet <- apply(cnet_2000,c(1,2,3),mean)
# Reminder: dimensions of cnet are:
# par = c("S","G-","G=","G+","F0","F1","F2"),
# stat =c("mu","e","cc","rho"),
# site = v.site

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

# # Plot Figure
# png("Figure_2a.png",width=8,height=8,units="cm",res=300)
# par(xpd=NA,cex=0.7, mar=c(3,3,2,1))
# pal <- brewer.pal(4,"Accent")
# barplot(C_k.perc,col=pal,border=NA,las=1)
# dev.off()

# Plot Figure no_perc
png("Figure_2a_noPerc.png",width=8,height=8,units="cm",res=300)
par(xpd=NA,cex=0.7, mar=c(3,3,2,1))
pal <- brewer.pal(4,"Accent")
barplot(C_k,col=pal,border=NA,las=1)
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

# # Plot Figure
# png("Figure_2b.png",width=8,height=8,units="cm",res=300)
# par(xpd=NA,cex=0.7, mar=c(3,3,2,1))
# pal <- brewer.pal(7,"Set1")
# barplot(C_l.perc,col=pal,border=NA,las=1)
# dev.off()

# Plot Figure no Perc
png("Figure_2b_noPerc.png",width=8,height=8,units="cm",res=300)
par(xpd=NA,cex=0.7, mar=c(3,3,2,1))
pal <- brewer.pal(7,"Set1")
barplot(C_l,col=pal,border=NA,las=1)
dev.off()

# Site-specific plots....................................
r.plot <- r %>% apply(1,mean,na.rm=T)
v.site.elev <- c("BRU (930m)", "CHA (1480m)", "VIL (1980m)", "LAU (2090m)", "GAL (2590m)", "PIC (2930m)")
png("Fig_2_BySite.png",width=16,height=9,units="cm",res=300)
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

# Plot legends
png("Fig_2_Legend.png",width=16,height=2.75,units="cm",res=600)
par(mar=c(0,0,0,0))
plot.new()
pal <- brewer.pal(4,"Accent")
legend(0,0.5,c(expression(italic("\u03BC")),expression(italic(e)),
               expression(italic(CV)),expression(rho)),fill=pal,bty="n",ncol=4)
pal <- brewer.pal(7,"Set1")
legend(0.5,0.5,c("S",expression(G^'-'),expression(G^'='),expression(G^'+'),expression(F[0]),expression(F[1]),expression(F[2])),fill=pal,bty="n",ncol=4)
dev.off()

# Regression of SLTRE contributions with elevation.......................

# 1 Combine observed values into a data.frame

# First rearrange parameters in the first dimension
temp <- abind(cnet_2000[,1,,],cnet_2000[,2,,],cnet_2000[,3,,],cnet_2000[,4,,],along=1)
temp.names <- expand.grid(dimnames(cnet_2000)[[1]], dimnames(cnet_2000)[[2]])
dimnames(temp) <- list(Var = paste0(temp.names$Var1,"_",temp.names$Var2),
                       Site = v.site,
                       Repl = c(1:dim(temp)[3]))
i <- 1
data <- 
data.frame(
    Site  = dimnames(temp)[[2]],
    med   = temp[i,,] %>% apply(1,mean,na.rm=T),
    low   = temp[i,,] %>% apply(1,quantile,0.025,na.rm=T),
    high  = temp[i,,] %>% apply(1,quantile,0.975,na.rm=T),
    Rate  = dimnames(temp)[[1]][i]
)

for (i in 2 : 28) {
    data <- rbind(data,
                  data.frame(
                      Site  = dimnames(temp)[[2]],
                      med   = temp[i,,] %>% apply(1,mean,na.rm=T),
                      low   = temp[i,,] %>% apply(1,quantile,0.025,na.rm=T),
                      high  = temp[i,,] %>% apply(1,quantile,0.975,na.rm=T),
                      Rate  = dimnames(temp)[[1]][i]
                  ))
    
    
}

# 2 Predict values from linear regressions

newdat <- data.frame(elev=seq(900,3000,100))
nrepl <- 2000
ypred <- array(NA,c(dim(newdat)[1],nrepl,28))
coef.b <- array(NA,c(nrepl,28))
# For each variable, fit a linear model for each replicate
for (i.var in 1 : 28) {
    cat(i.var,"\n"); flush.console()
    for (i.repl in 1 : nrepl) {
        yobs <- temp[i.var,,i.repl]
        mod <- lm(yobs~elev)
        coef.b[i.repl,i.var] <- coef(mod)[2]
        ypred[,i.repl,i.var] <- predict(mod,newdat)
    }
}


# Calculate mean and 95%CI of the coefficients
data.coef <- data.frame(
    med = as.vector(apply(coef.b,2,mean)),
    low = as.vector(apply(coef.b,2,quantile,0.025)),
    high = as.vector(apply(coef.b,2,quantile,0.975)),
    Rate = row.names(temp)
)
# Calculate mean and 95%CI of the predicted values
data.pred <- data.frame(
    elev = newdat$elev,
    med = as.vector(apply(ypred,c(1,3),mean)),
    low = as.vector(apply(ypred,c(1,3),quantile,0.025)),
    high = as.vector(apply(ypred,c(1,3),quantile,0.975)),
    Rate = rep(row.names(temp),each=dim(newdat)[1])
)

# Significant relationships are those with concordant signs among lower and upper bounds of CI
data.coef <- mutate(data.coef,significant=sign(low)*sign(high)) 
data.coef$linetype <- "solid"
data.coef$linetype[data.coef$significant==-1] <- "dotted"

# Reorder rates 
data.pred$Rate <- as.character(data.pred$Rate)
data.pred$Rate <- factor(data.pred$Rate,
                         levels=c("S_mu", "G-_mu", "G=_mu", "G+_mu", "F0_mu", "F1_mu", "F2_mu",
                                  "S_e", "G-_e", "G=_e", "G+_e", "F0_e", "F1_e", "F2_e",
                                  "S_cc", "G-_cc", "G=_cc", "G+_cc", "F0_cc", "F1_cc", "F2_cc",
                                  "S_rho", "G-_rho", "G=_rho", "G+_rho", "F0_rho", "F1_rho", "F2_rho")
)


save(temp,data,data.pred,data.coef,elev,file="Data_for_Figure_S8.RData")  


# 4 Plot

load("Data_for_Figure_S8.RData")
data$elev <- elev[match(data$Site,names(elev))]

png("Figure_S7_unlabelled.png",width=11,height=20,units="cm",res=600)
ggplot(data,aes(x=elev,y=med, ymin=low, ymax=high)) +
    geom_ribbon(data=data.pred,fill="gray") +
    geom_line(data=data.pred,colour="blue",aes(linetype=Rate),size=1) +
    scale_linetype_manual(values=data.coef$linetype) +
    geom_errorbar(width=0,size=0.5) +
    geom_point(size=1.5) + 
    labs(y="",x="") +
    geom_hline(yintercept=0,linetype="dashed") +
    scale_x_continuous(breaks=c(1000,2000,3000)) +
    labs(y="",x="Elevation (m)") +
    facet_wrap(~Rate,scales="free_y",nrow=7,dir="v")
dev.off()





