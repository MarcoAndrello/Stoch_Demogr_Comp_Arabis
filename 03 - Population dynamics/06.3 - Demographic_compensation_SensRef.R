# Demographic compensation
# Marco Andrello
# Last modification: 23/01/2020

rm(list=ls())

library(tidyverse)
library(abind)
library(corrgram)
library(fields)
library(gridExtra)
library(RColorBrewer)
source("00.1 - Functions_MPM.R")

v.scen <- c(1,5)
i.scen <- 2
v.site <- c("BRU", "CHA", "VIL", "LAU", "GAL", "PIC")
v.year <- c(2008:2013)

load(paste0(getwd(),"/Results SensRef/Results_MPM_SensRef_scen_",v.scen[i.scen],"_.RData"))
cnet_2000 <- cnet
log.lambda.R_2000 <- log.lambda.R
nrepl <- dim(cnet_2000)[4]
# First rearrange parameters so as to have them in columns
data <- abind(cnet_2000[,1,,],cnet_2000[,2,,],cnet_2000[,3,,],cnet_2000[,4,,],along=1)
data <- aperm(data,c(2,1,3))
data.names <- expand.grid(dimnames(cnet_2000)[[1]], dimnames(cnet_2000)[[2]])
dimnames(data) <- list(Site = v.site,
                       Var = paste0(data.names$Var1,"_",data.names$Var2),
                       Repl = c(1:dim(data)[3]))


# CALCULATING CORRELATIONS
cor.mat.obs <- array(NA,c(28,28,2000))
for (i.repl in 1 : nrepl) {
    data.repl <- data[,,i.repl]
    cor.mat.obs[,,i.repl] <- cor(data.repl,method="spearman")
}
dimnames(cor.mat.obs) <- list(var1 = colnames(data), var2 = colnames(data), repl = c(1:nrepl) )
med.obs  <- apply(cor.mat.obs,c(1,2),mean)
low.obs  <- apply(cor.mat.obs,c(1,2),quantile,0.025)
high.obs <- apply(cor.mat.obs,c(1,2),quantile,0.975)
#neg.corr.par.vec <- unique(rownames(which(high.obs<0, arr.ind=T))) # Still useful ?
# save(low.obs,med.obs,high.obs,neg.corr.par.vec, file=paste0("Observed_negative_correlations_SensRef_",v.scen[i.scen],"_.RData"))
save(low.obs,med.obs,high.obs, file=paste0("Observed_negative_correlations_SensRef_",v.scen[i.scen],"_.RData"))
# Plot the correlation matrix
data.plot <- data[,,2000]
colnames(data.plot) <- rep("",dim(data.plot)[2])
png(paste0("Figure_3_corrplot_SensRef_scen_",v.scen[i.scen],"_.png"),height=15,width=15,units="cm",res=600)
corrgram(data.plot,cor.method="spearman",lower.panel=NULL, upper.panel=panel.fct, diag.panel=NULL,
         cex.labels=0.6, gap=0.15, col.regions = colorRampPalette(c("red", "salmon","white", "lightblue", "royalblue")),
         data1=data.plot, med=med.obs, low=low.obs, high=high.obs)

dev.off()
# NUMBER OF CORRELATIONS
med.triang.obs <- med.obs # Define triangular matrices of correlations to count things only once
med.triang.obs[lower.tri(med.triang.obs,diag=T)] <- NA
low.triang.obs <- low.obs
low.triang.obs[lower.tri(low.triang.obs,diag=T)] <- NA
high.triang.obs <- high.obs
high.triang.obs[lower.tri(high.triang.obs,diag=T)] <- NA
num.pos.obs <- dim(which(low.triang.obs>0,arr.ind=T))[1] # Calculate number of significant negative and positive correlations
num.neg.obs <- dim(which(high.triang.obs<0,arr.ind=T))[1]
# Permutation test
nrand <- 1000
num.pos.rand <- num.neg.rand <- rep(NA,nrand)
for (i.rand in 1: nrand) {
    cat(i.rand,"\n"); flush.console()
    # Permute each parameter (all the 2000 bootstrap samples) across sites
    data.rand <- array(NA,dim(data))
    for (i.par in 1 : ncol(data)) {
        id.rand <- sample(6)
        data.rand[,i.par,] <- data[id.rand,i.par,]
    }
    # Set up correlation matrix
    cor.mat.rand <- array(NA,c(ncol(data),ncol(data),2000))
    # Calculate Spearman correlation for each replicate
    for (i.repl in 1 : nrepl) { 
        data.rand.repl <- data.rand[,,i.repl]
        cor.mat.rand[,,i.repl] <- cor(data.rand.repl,method="spearman")
    }
    cor.mat.rand[is.na(cor.mat.rand)] <- 0 # SEED BANK: replace the NA with 0
    # Calculate lower and upper limits of correlations
    low.rand <- apply(cor.mat.rand,c(1,2),quantile,0.025)
    high.rand <- apply(cor.mat.rand,c(1,2),quantile,0.975)
    # Define triangular matrices of correlations to count things only once
    low.triang.rand <- low.rand
    low.triang.rand[lower.tri(low.triang.rand, diag = T)] <- NA
    high.triang.rand <- high.rand
    high.triang.rand[lower.tri(high.triang.rand, diag = T)] <- NA
    # Calculate number of significant negative and positive correlations
    num.pos.rand[i.rand] <- dim(which(low.triang.rand > 0, arr.ind = T))[1]
    num.neg.rand[i.rand] <- dim(which(high.triang.rand < 0, arr.ind = T))[1]
}
save(num.pos.obs,num.neg.obs,num.pos.rand,num.neg.rand,file=paste0(getwd(),"/Results SensRef/Res.numberCorr_SensRef_",v.scen[i.scen],"_.RData"))

# Plot the histogram
png(paste0(getwd(),"/Results SensRef/Figure_3_Histograms_SensRef_scen_",v.scen[i.scen],"_.png"),height=15,width=7.5,units="cm",res=300)
par(mfrow=c(2,1),mar=c(4,4,3,1),cex.main=0.85)
hist(num.pos.rand,xlab="",ylab="",main="Number of\npositive correlations",col="gray",border="darkgray",
     xlim=c(0,30),ylim=c(0,250),las=1)
abline(v=num.pos.obs,lwd=1,lty=2)
hist(num.neg.rand,xlab="",ylab="",main="Number of\n negative correlations",col="gray",border="darkgray",
     xlim=c(0,30),ylim=c(0,250),las=1)
abline(v=num.neg.obs,lwd=1,lty=2)
dev.off()
#



# EFFECTIVENESS OF DEMOGRAPHIC COMPENSATION PER PARAMETER
load(paste0(getwd(),"/Results SensRef/Observed_negative_correlations_SensRef_",v.scen[i.scen],"_.RData"))
data.mean <- apply(data,c(1,2),mean)
log.lambda.R <- mean(log.lambda.R_2000)
log.lambda.obs <- log.lambda.R + rowSums(data.mean)
var_log.lambda.obs <- var(log.lambda.obs) # which is equal to sum(rowSums(data)^2)/5
# Permutation test
nrand <- 1000
var_log.lambda.rand <- array(NA,c(28,nrand))
for (i.par in 1 : 28 ){
    cat("Permutation parameter",i.par,"\n")
    data.rand <- data.mean
    for (i.rand in 1: nrand) {
        data.rand[,i.par] <- data.mean[sample(6),i.par] # Permute the parameter across sites
        log.lambda.rand <- log.lambda.R + rowSums(data.rand) #Variance in log.lambda between sites
        var_log.lambda.rand[i.par,i.rand] <- var(log.lambda.rand)
    }
}
data.var <- var_log.lambda.obs / var_log.lambda.rand
data.var <- var_log.lambda.rand
data.var <- data.frame(
    low = apply(data.var,1,quantile,0.025),
    med = apply(data.var,1,mean),
    high = apply(data.var,1,quantile,0.975)
)
data.var$lcc <- factor(c("S","G-","G=","G+","F0","F1","F2"),levels=c("S","G-","G=","G+","F0","F1","F2"))
data.var$descr <- factor(rep(c("Mean","Elasticity","CV","Temporal correlations"),each=7),levels=c("Mean","Elasticity","CV","Temporal correlations"))
summary(data.var)
# Theme for plotting
theme_set(theme_classic())
theme_update(axis.text.x = element_text(size=10, color="black"),
             axis.text.y = element_text(size=10, color="black"),
             axis.title.x = element_text(size=10, color="black"),
             axis.title.y = element_text(size=10, color="black"),
             panel.spacing.x = unit(0.45, "cm"),
             plot.margin = margin(0.5,0.2,0,0, unit="cm"),
             strip.text.x = element_text(size = 6.5),
             strip.background.x = element_rect(linetype = "blank"),
             legend.position = "bottom",
             legend.direction = "horizontal",
             legend.box.spacing = unit(0, "cm"),
             legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),
             legend.title=element_text(colour="white"))

png(filename=paste0("Figure_4_errorbars_SensRef_scen_",v.scen[i.scen],"_.png"),width=11,height=5.5,units="cm",res=300)
ggplot(data.var,aes(x=lcc, y=med, ymin=low, ymax=high, col=descr)) +
    geom_point(size=1.5,position=position_dodge(width=1))+
    geom_errorbar(width=0.75,size=0.5,position=position_dodge(width=1)) +
    #facet_wrap(~descr,scale="free_x") +
    geom_hline(yintercept=var_log.lambda.obs,linetype="dashed") +
    labs(y=bquote(sigma[log~lambda[s]]^2),x="") +
    scale_color_brewer(type="qual",palette=1,
                       labels=c(bquote(mu), bquote(e),bquote(CV),bquote(rho))) +
    scale_x_discrete(labels=c("S",expression(G^'-'),expression(G^'='),expression(G^'+'),expression(F[0]),expression(F[1]),expression(F[2]))) +
    geom_vline(xintercept=c(1:6)+0.5,linetype="dotted",size=0.25)

dev.off()


par.to.permute <- which(data.var$med>var_log.lambda.obs) # Which are the parameters whose permutation increases Var(log_lambda)
# Permutation of those
nrand <- 1000
var_log.lambda.rand.overall <- rep(NA,nrand)
for (i.rand in 1: nrand) {
    # Permute the (parameters to permute) across sites
    data.rand <- data.mean
    for (i.par in 1 : 28) {
        if (any(i.par == par.to.permute)) {
            id.rand <- sample(6)
            data.rand[,i.par] <- data.mean[id.rand,i.par]
        }
    }
    #Variance in log.lambda between sites
    log.lambda.rand <- log.lambda.R + rowSums(data.rand)
    var_log.lambda.rand.overall[i.rand] <- var(log.lambda.rand)
}

png(filename=paste0("Figure_4_histogram_SensRef_scen_",v.scen[i.scen],"_.png"),width=5.5,height=5.5,units="cm",res=300)
par(mar=c(4,3,1,1),cex=0.8)
hist(var_log.lambda.rand.overall,xlab=bquote(sigma[log~lambda[s]]^2),ylab="",main="",col="gray",border="darkgray",
     las=1)
abline(v=var_log.lambda.obs,lwd=1,lty=2)
dev.off()

var_log.lambda.obs
median(var_log.lambda.rand.overall)
var_log.lambda.obs / median(var_log.lambda.rand.overall)
length(which(var_log.lambda.rand.overall<var_log.lambda.obs)) / nrand


