# Seed bank model
rm(list=ls())

library(popbio)
library(IPMpack)
library(tidyverse)
library(lme4)
library(gridExtra)
library(fields)
library(tictoc)
library(abind)
library(RColorBrewer)

# Various definitions and source functions
v.site <- c("BRU", "CHA", "VIL", "LAU", "GAL", "PIC")
v.year <- c(2008:2013)
num.sites <- length(v.site)
num.years <- length(v.year)
source("01 - Functions_IPM.R")
source("05.1 - Calc.SLTRE.sb.R")

load("Results_IPM.RData")

k <- 50 # Number of stages

# Constant germination rate across sites

v.g <- seq(0.1,1,0.3)

# Initialize matrices
P <- mat.P
G <- mat.G
F <- mat.F
K <- mat.F
r <- Lexp.mean <- Lexp.sd <- array(NA,c(length(v.g), length(v.site)))
nsd <- array(NA,c(k,num.years))
cnet <- array(NA,c(6, # Number of life-cycle components: seed germination, (survival, retro.growth, stasis, prog.growth, fec)
                   4, # Number of type of SLTRE contributions (mean, elast, cv, corr)
                   length(v.site),
                   length(v.g)))

i.g <- 1
for (i.g in 1 : length(v.g)) {
  g <- v.g[i.g] # Germination rate
  cat(g,"\n"); flush.console()
  P1 <- S1 <- G1 <- F1 <- K1 <- list()
  i.site <- 1
  for (i.site in 1 : num.sites) {
    P1[[i.site]] <- S1[[i.site]] <- G1[[i.site]] <- F1[[i.site]] <- K1[[i.site]] <- list()
    i.year <- 1
    for (i.year in 1 : num.years) {
      
      surv <- colSums(P[[i.site]][[i.year]])
      nsd[,i.year] <- F[[i.site]][[i.year]][,1] / sum(F[[i.site]][[i.year]][,1])
      
      # Modify G
      g_sb <- c((1-g),g*nsd[,i.year])
      G1[[i.site]][[i.year]] <- rbind(rep(0,k),
                                      G[[i.site]][[i.year]])
      G1[[i.site]][[i.year]] <- cbind(g_sb,G1[[i.site]][[i.year]])
      
      # Modify F
      fj <- colSums(F[[i.site]][[i.year]]) / (g*(2-g))
      
      F1[[i.site]][[i.year]] <- nsd[,i.year] %*% t(fj) * g
      F1[[i.site]][[i.year]] <- rbind(fj * (1-g),
                                      F1[[i.site]][[i.year]])
      F1[[i.site]][[i.year]] <- cbind(rep(0,(k+1)),
                                      F1[[i.site]][[i.year]])
      
      # Modify S and K
      S1[[i.site]][[i.year]] <- matrix(c(1,surv),(k+1),(k+1),byrow=T)
      P1[[i.site]][[i.year]] <- S1[[i.site]][[i.year]]*G1[[i.site]][[i.year]]
      K1[[i.site]][[i.year]] <- P1[[i.site]][[i.year]] + F1[[i.site]][[i.year]]
    }
    
    r[i.g, i.site] <- stochGrowthRateSampleList(nRunIn=1000, tMax=11000, listIPMmatrix=K1[[i.site]])
    sle <- stoch.life.expectancy(P1[[i.site]])
    nsd_mean <- rowMeans(nsd)
    nsd_mean_1 <- c((1-g),g*nsd_mean)
    Lexp.mean[i.g, i.site] <- weighted.mean(sle$mean, nsd_mean_1)
    Lexp.sd[i.g, i.site] <- weighted.mean(sle$sd, nsd_mean_1)
    # Stochastic elasticities and stochastic stage distribution
    # res.elast <- stoch.sens.Marco(K1[[i.site]],
    #                               G1[[i.site]],
    #                               P1[[i.site]],
    #                               F1[[i.site]],
    #                               pfruit,
    #                               log.nfruits,
    #                               nsd,
    #                               100)
    # elast[i.site,,i.repl] <- res.elast$elast
    # ws[,i.site,i.repl] <- res.elast$w
  }
  

  cnet[,,,i.g] <- calc.SLTRE.sb(K1,
                                G1,
                                P1,
                                F1)
}

# Examine contributions of SLTRE
cnet_all <- cnet
dimnames(cnet_all) <- list(par=c("S","G-","G=","G+","F","germ"),
                       stat=c("mu","e","cc","rho"),
                       site = v.site,
                       g = v.g)

# For mathematical annotation:
#https://lukemiller.org/index.php/2017/05/r-plotmath-functions-combined-with-variable-values/

pal <- brewer.pal(4, "Accent")
for (i.g in 1 : 4) {
  cnet <- cnet_all[,,,i.g]
  png(paste0(getwd(),"/Figures/Fig_SLTRE_BySite_sb_",v.g[i.g],".png"),width=13,height=9,units="cm",res=300)
  par(mfrow=c(2,3),mar=c(3,3,2,1))
  for (i.site in 1 : num.sites) {
    c1 <- c2 <- c3 <- cnet[,,i.site]
    c2[c2>0] <- 0 
    c3[c3<0] <- 0 
    #myRange <- c( min(rowSums(c2)), max(rowSums(c3)) )
    myRange <- c(-0.2,0.45)
    barplot(t(c2),col=pal,space=0,border=NA,axisnames=FALSE, ylim=myRange,main=v.site[i.site])
    barplot(t(c3),col=pal,space=0,border=NA,axisnames=FALSE, add=T)
    axis(1, at=c(0:5)+0.5, lab= row.names(cnet), las=2, cex.axis=1)
    text(4,0.4,bquote(log~lambda*''==''*.(round(r[i.g,i.site],2))))
    text(4,0.3,bquote(italic(L)[exp]*''==''*.(round(Lexp.mean[i.g,i.site],1))))
    if (i.site==1) legend(-0.5,0.4,c("mean","elast","CV","cor"),fill=pal,bty="n")
  }
  dev.off()
}


# Germination rates decreasing with elevation
elev <- c(930,1480,1980,2090,2590,2930)
names(elev) <- v.site 
germ <- 1-elev/3000 # Germination rate per site


# Initialize matrices
P <- mat.P
G <- mat.G
F <- mat.F
K <- mat.F
r <- Lexp.mean <- Lexp.sd <- array(NA,length(v.site))
nsd <- array(NA,c(k,num.years))
cnet <- array(NA,c(6, # Number of life-cycle components: seed germination, (survival, retro.growth, stasis, prog.growth, fec)
                   4, # Number of type of SLTRE contributions (mean, elast, cv, corr)
                   length(v.site)))


P1 <- S1 <- G1 <- F1 <- K1 <- list()
i.site <- 1
for (i.site in 1 : num.sites) {
  g <- germ[i.site] # Germination rate
  P1[[i.site]] <- S1[[i.site]] <- G1[[i.site]] <- F1[[i.site]] <- K1[[i.site]] <- list()
  i.year <- 1
  for (i.year in 1 : num.years) {
    
    surv <- colSums(P[[i.site]][[i.year]])
    nsd[,i.year] <- F[[i.site]][[i.year]][,1] / sum(F[[i.site]][[i.year]][,1])
    
    # Modify G
    g_sb <- c((1-g),g*nsd[,i.year])
    G1[[i.site]][[i.year]] <- rbind(rep(0,k),
                                    G[[i.site]][[i.year]])
    G1[[i.site]][[i.year]] <- cbind(g_sb,G1[[i.site]][[i.year]])
    
    # Modify F
    fj <- colSums(F[[i.site]][[i.year]]) / (g*(2-g))
    
    F1[[i.site]][[i.year]] <- nsd[,i.year] %*% t(fj) * g
    F1[[i.site]][[i.year]] <- rbind(fj * (1-g),
                                    F1[[i.site]][[i.year]])
    F1[[i.site]][[i.year]] <- cbind(rep(0,(k+1)),
                                    F1[[i.site]][[i.year]])
    
    # Modify S and K
    S1[[i.site]][[i.year]] <- matrix(c(1,surv),(k+1),(k+1),byrow=T)
    P1[[i.site]][[i.year]] <- S1[[i.site]][[i.year]]*G1[[i.site]][[i.year]]
    K1[[i.site]][[i.year]] <- P1[[i.site]][[i.year]] + F1[[i.site]][[i.year]]
  }
  
  r[i.site] <- stochGrowthRateSampleList(nRunIn=1000, tMax=11000, listIPMmatrix=K1[[i.site]])
  sle <- stoch.life.expectancy(P1[[i.site]])
  nsd_mean <- rowMeans(nsd)
  nsd_mean_1 <- c((1-g),g*nsd_mean)
  Lexp.mean[i.site] <- weighted.mean(sle$mean, nsd_mean_1)
  Lexp.sd[i.site] <- weighted.mean(sle$sd, nsd_mean_1)
  
}

cnet <- calc.SLTRE.sb(K1,
                      G1,
                      P1,
                      F1)


# Examine contributions of SLTRE
dimnames(cnet) <- list(par=c("S","G-","G=","G+","F","germ"),
                      stat=c("mu","e","cc","rho"),
                      site = v.site)

pal <- brewer.pal(4, "Accent")
png(paste0(getwd(),"/Figures/Fig_SLTRE_BySite_sb_varElev.png"),width=13,height=9,units="cm",res=300)
par(mfrow=c(2,3),mar=c(3,3,2,1))
for (i.site in 1 : num.sites) {
  c1 <- c2 <- c3 <- cnet[,,i.site]
  c2[c2>0] <- 0 
  c3[c3<0] <- 0 
  #myRange <- c( min(rowSums(c2)), max(rowSums(c3)) )
  myRange <- c(-0.2,0.45)
  barplot(t(c2),col=pal,space=0,border=NA,axisnames=FALSE, ylim=myRange,main=v.site[i.site])
  barplot(t(c3),col=pal,space=0,border=NA,axisnames=FALSE, add=T)
  axis(1, at=c(0:5)+0.5, lab= row.names(cnet), las=2, cex.axis=1)
  text(4,0.4,bquote(log~lambda*''==''*.(round(r[i.site],2))))
  text(4,0.3,bquote(italic(L)[exp]*''==''*.(round(Lexp.mean[i.site],1))))
  if (i.site==1) legend(-0.5,0.4,c("mean","elast","CV","cor"),fill=pal,bty="n")
}
dev.off()
