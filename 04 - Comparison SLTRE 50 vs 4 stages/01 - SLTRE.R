# SLTRE Arabis Davison et al. 2013 Am Nat
# Stochastic Life-Table Response Experiment (SLTRE) (Davison et al. 2013 Am Nat)
# considering F0, F1 and F2 and not F

rm(list=ls())

library(fields)
library(popbio)
library(abind)
library(RColorBrewer)
library(corrgram)
library(ggplot2)
library(ggpubr)
source("00 - Functions_SLTRE.R")

red.mat <- F
v.site <- c("BRU","CHA","VIL","LAU","GAL","PIC") # Needed ?
num.sites <- length(v.site)
load("Results_IPM_OneRepl_perSLTRE_50vs4.RData") # Results obtained by running the IPM model for one replicate only and saving pfruit.all, log.nfruits.all and nsd

# Size of the (reduced) matrices, num.sites and num.years
if (red.mat) k <- 4 else k <- dim(matrices.G[[1]][[1]])[1]
num.sites <- length(matrices.G)
num.years <- length(matrices.G[[1]])


# Transform lists into arrays and reduce matrices to k stages if red.mat=T
G <- P <- array(NA,c(k, k, num.sites, num.years)) 
F0 <- F1 <- F2 <- array(NA,c(k, num.sites, num.years))
for (i.site in 1 : num.sites){
  for(i.year in 1 : num.years) {
    if(red.mat) {
      G[,,i.site,i.year] <- reduce.mat(matrices.G[[i.site]][[i.year]])
      P[,,i.site,i.year] <- reduce.mat(matrices.P[[i.site]][[i.year]])
      #P[,,i.site,i.year] <- reduce.mat(matrices.P[[i.site]][[i.year]])
      F0[,i.site,i.year] <- reduce.F(pfruit.all[[i.site]][,i.year])
      F1[,i.site,i.year] <- reduce.F(log.nfruits.all[[i.site]][,i.year])
      F2[,i.site,i.year] <- reduce.F2(nsd[[i.site]][,i.year])
    } else {
      G[,,i.site,i.year] <- matrices.G[[i.site]][[i.year]]
      P[,,i.site,i.year] <- matrices.P[[i.site]][[i.year]]
      F0[,i.site,i.year] <- pfruit.all[[i.site]][,i.year]
      F1[,i.site,i.year] <- log.nfruits.all[[i.site]][,i.year]
      F2[,i.site,i.year] <- nsd[[i.site]][,i.year]
    }
  }
}
rm(matrices.G, matrices.P, pfruit.all, log.nfruits.all, nsd) # anche log.nfruits etc?


# Calculate parameters
surv <- array(NA,c(k,num.sites,num.years)) # k, year, site
g <- array(NA,c((k*k),num.sites,num.years))
for (i.site in 1 : num.sites){
  for(i.year in 1 : num.years) {
    surv[,i.site,i.year] <- colSums(P[,,i.site,i.year])
    g[,i.site,i.year] <- as.vector(G[,,i.site,i.year]) # g11, g21, g31, g41, g12, ... g44
  }
}
num.par <- k + (k*k) + k + k + k # Stage-specific surv(j), gamma(i,j), F0(j), F1(j), F2(i)
### erase matrices (A, G, P, F) here?

# Arrange parameters in the same array
park <- abind(surv, g, F0, F1, F2, along=1)

# Calculate parameters of the reference site (mean site)
survR <- apply(surv, c(1,3), mean)
gR <- apply(g, c(1,3), mean)
F0R <- apply(F0, c(1,3), mean)
F1R <- apply(F1, c(1,3), mean)
F2R <- apply(F2, c(1,3), mean)
parkR <- array(abind(survR, gR, F0R, F1R, F2R, along=1), c(num.par, 1, num.years))

# Add it to the array
park <- abind(park, parkR, along=2)
dimnames(park) <- list(par=c(1:num.par),site=c(1:(num.sites+1)),year=c(1:num.years))
rm(surv,g,F0,F1,F2,survR,gR,F0R,F1R,F2R,parkR)

# Calculate means over years
parkmean <- apply(park,c(1,2),mean)
dimnames(parkmean) <- list(par=c(1:num.par),site=c(1:(num.sites+1)))

# Calculate intrinsic growth rates (r) and elasticities in the mean population
r <- array(NA,(num.sites+1))
e <- array(NA,c(num.par,(num.sites+1)))

for (i.site in 1 : (num.sites+1)) {
  res <- calc.elast.vital.rates(surv   = parkmean[1:k,i.site],
                                growth = parkmean[(k+1) : (k+(k*k)), i.site],
                                F0 = parkmean[(k+(k*k)+1) : (k+(k*k)+k), i.site],
                                F1 = parkmean[(k+(k*k)+k+1) : (k+(k*k)+k+k), i.site],
                                F2 = parkmean[(k+(k*k)+k+k+1) : (k+(k*k)+k+k+k), i.site])
  e[,i.site] <- res$elast
  r[i.site] <- res$r
}
dimnames(e) <- dimnames(parkmean)

# Calculate coefficient of variation over years per par, per site 
cv <- apply(park, c(1,2), function(x) { sd(x) / mean(x) } )

# Calculate correlations rho between par, per site 
rho <- array(NA,c(num.par,num.par,(num.sites+1)))
for (i.site in 1 : (num.sites+1)) {
  rho[,,i.site] <- cor(t(park[,i.site,]))
}

# Matrices of ee, per site 
ee <- array(NA,c(num.par,num.par,(num.sites+1)))
for (i.site in 1 : (num.sites+1)) {
  ek <- e[,i.site] * matrix(1,num.par,num.par)
  el <- t(ek)
  ee[,,i.site] <- ek*el
}

# Matrices of cc per site 
cc <- array(NA,c(num.par,num.par,(num.sites+1)))
for (i.site in 1 : (num.sites+1)) {
  ck <- cv[,i.site] * matrix(1,num.par,num.par)
  cl <- t(ck)
  cc[,,i.site] <- ck*cl
}

# Matrices of ccrho, per site 
ccrho <- array(NA,c(num.par,num.par,(num.sites+1)))
for (i.site in 1 : (num.sites+1)) {
  ccrho[,,i.site] <- cc[,,i.site] * rho[,,i.site]
}



# Differences
diffpark <- array(NA,c(num.par,(num.sites+1)))
diffee <-
  meanee <- 
  diffcc <-
  meancc <-
  diffrho <-
  meanrho <-
  diffccrho <-
  meanccrho <- array(NA,c(num.par,num.par,(num.sites+1)))

for (i.site in 1 : num.sites) {
  diffpark[,i.site]    <- parkmean[,i.site] - parkmean[,(num.sites+1)]
  diffee[,,i.site]     <-  ee[,,i.site] - ee[,,(num.sites+1)]
  meanee[,,i.site]     <- (ee[,,i.site] + ee[,,(num.sites+1)]) / 2
  diffcc[,,i.site]     <-  cc[,,i.site] - cc[,,(num.sites+1)]
  meancc[,,i.site]     <- (cc[,,i.site] + cc[,,(num.sites+1)]) / 2
  diffrho[,,i.site]    <-  rho[,,i.site] - rho[,,(num.sites+1)]
  meanrho[,,i.site]    <- (rho[,,i.site] + rho[,,(num.sites+1)]) / 2
  diffccrho[,,i.site]  <-  ccrho[,,i.site] - ccrho[,,(num.sites+1)]
  meanccrho[,,i.site]  <- (ccrho[,,i.site] + ccrho[,,(num.sites+1)]) / 2
}

# Calculating contributions (C_x) and net stochastic contributions (C_net_e)
C_mu <- array(NA,c(num.par,num.sites))
C_e     <-
  C_ccrho <-
  C_cc    <-
  C_rho   <- array(NA,c(num.par,num.par,num.sites))
C_net_e     <- 
  C_net_ccrho <-
  C_net_cc    <-
  C_net_rho   <- array(NA,c(num.par,num.sites))

for (i.site in 1 : num.sites) {
  # Contributions
  C_mu[,i.site] <- diffpark[,i.site] / parkmean[,(num.sites+1)] * e[,i.site]            # Davison et al. 2013, eq B4
  C_e[,,i.site]     <- -0.5 *  meanccrho[,,i.site] * diffee[,,i.site] 
  C_ccrho[,,i.site] <- -0.5 * meanee[,,i.site] * diffccrho[,,i.site]            
  C_cc[,,i.site]    <- -0.5 * meanee[,,i.site] * meanrho[,,i.site] * diffcc[,,i.site]
  C_rho[,,i.site]   <- -0.5 * meanee[,,i.site] * meancc[,,i.site] * diffrho[,,i.site]
  
  # Net stochastic contributions
  C_net_e[,i.site]     <- colSums(C_e[,,i.site])
  C_net_ccrho[,i.site] <- colSums(C_ccrho[,,i.site])
  C_net_cc[,i.site]    <-colSums(C_cc[,,i.site])
  C_net_rho[,i.site]   <-colSums(C_rho[,,i.site])
}

#Name dimensions
dimnames(C_mu)        <-
  dimnames(C_net_e)     <-
  dimnames(C_net_ccrho) <-
  dimnames(C_net_cc)    <-
  dimnames(C_net_rho)  <- list(par=c(1:num.par), site=v.site)

# Bind the four types of contributions into the same array
C_net <- abind(C_mu, C_net_e, C_net_cc, C_net_rho, along=3)
dimnames(C_net)  <- list(par=c(1:num.par), site=v.site, stat=c("mu","e","cc","rho"))
cnet <- aperm(C_net,c(1,3,2))

# Save
if (red.mat) {
  save(cnet, file="SLTRE_Reduce.RData")
} else {
  save(cnet, file="SLTRE_noReduce.RData")
}


## COMPARE REDUCE VS NOREDUCE
# Compare the contributions of the different statistics in the different sites

# Total effects of each statistic as percent of summed absolute contributions
tes.p <- array(NA,dim=c(4,6,2)) # 4 statistics, 6 sites, 2 versions
load("SLTRE_noReduce.RData")
tes <- apply(abs(cnet),c(2,3),mean) # Total effect of statistic in each site
tes.p[,,1] <- apply(tes, 2, function(x) {x/sum(x)})
load("SLTRE_Reduce.RData")
tes <- apply(abs(cnet),c(2,3),mean)
tes.p[,,2] <- apply(tes, 2, function(x) {x/sum(x)})
dimnames(tes.p) <- list(stat=c("mu","e","cc","rho"),
                       site = v.site,
                       matrix = c("50 classes","4 classes"))
data <- data.frame(m50=as.vector(tes.p[,,1]),
                   m4=as.vector(tes.p[,,2]),
                   stat=dimnames(tes.p)[[1]],
                   site=factor(rep(dimnames(tes.p)[[2]],each=4)))
plot(m4~m50, data=data, col=stat, pch=16)

# To plot stat with the same colours as in the main Figure
data$stat <- as.character(data$stat)
data$stat <- factor(data$stat,levels=c("mu","e","cc","rho"))

# Set themes for plotting
theme_set(theme_classic())
theme_update(axis.text.x = element_text(size=8, color="black"),
             axis.text.y = element_text(size=8, color="black"),
             axis.title.x = element_text(size=8, color="black"),
             axis.title.y = element_text(size=8, color="black"),
             plot.title = element_text(size = 8, color = "black", hjust=0.5),
             legend.position = "bottom")

png("Fig_comparison50vs4.png",width=7,height=8,units="cm",res=200)
ggplot(data,aes(x=m50, y=m4, col=stat)) +
  geom_point() +
  geom_abline(intercept=0,slope=1,linetype=2) +
  scale_colour_brewer(palette="Accent", name=NULL, labels= c(expression(italic("\u03BC")),
                                                             expression(italic(E)),
                                                             expression(italic(CV)),
                                                             expression(rho)) ) +
  scale_x_log10() +
  scale_y_log10() +
  stat_cor(aes(x=m50, y=m4),inherit.aes=F,size=3) + 
  labs(x="fifty-stages matrix",y="four-stages matrix")
dev.off()
  
cor.test(data$m50, data$m4, method="pe")

legend(0,-20,c(expression(italic("\u03BC")),expression(italic(E)),
               expression(italic(CV)),expression(rho)),fill=pal,bty="n",ncol=4)
