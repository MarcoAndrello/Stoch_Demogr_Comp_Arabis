# Analyse results
rm(list=ls())

library(tidyverse)
library(gridExtra)
library(fields)
library(RColorBrewer)
library(abind)
library(corrgram)
library(shades)

load("Results_IPM.RData")

v.site <- c("BRU", "CHA", "VIL", "LAU", "GAL", "PIC")
v.year <- c(2008:2013)
v.stage <- c(1,5,10,15)

# Function to calculate coefficient of variation
cv <- function(x) { sd(x) / mean(x) }

# Set themes for plotting
theme_set(theme_classic())
theme_update(axis.text.x = element_text(size=15, color="black"),
             axis.text.y = element_text(size=15, color="black"),
             axis.title.x = element_text(size=15, color="black"),
             axis.title.y = element_text(size=15, color="black"),
             plot.title = element_text(size = 15, color = "black", hjust=0.5))


#####################################################
# Add dimnames
#####################################################

# Vital rates
dimnames(nfruits) <- 
  dimnames(pfruit) <- 
  dimnames(surv) <- 
  dimnames(progr) <- 
  dimnames(retro) <- 
  dimnames(stasis) <-  list(stage = v.stage,
                            site = v.site,
                            year = v.year,
                            repl=c(1:dim(mean.fec)[3]))

# Life-cycle components
dimnames(mean.nfruits.net) <- 
dimnames(mean.nfruits) <- 
dimnames(mean.pfruit) <- 
dimnames(mean.recr.size) <- 
dimnames(mean.surv) <- 
dimnames(mean.progr) <- 
dimnames(mean.retro) <- 
dimnames(mean.stasis) <- 
dimnames(mean.fec) <- list(site = v.site,
                             year = v.year,
                             repl=c(1:dim(mean.fec)[3]))
# Stage distribution
dimnames(w) <- list(stage= y,
                    site = v.site,
                    year = v.year,
                    repl = c(1:dim(mean.fec)[3]))
dimnames(ws) <- list(stage= y,
                     site = v.site,
                     repl = c(1:dim(mean.fec)[3]))
# Other variables
row.names(Lexp.mean) <- 
  row.names(Lexp.sd)   <- 
  row.names(r)         <- v.site
dimnames(lambda.det) <- list(site = v.site,
                             year = v.year,
                             repl=c(1:dim(mean.fec)[3]))

elev <- c(930,1480,1980,2090,2590,2930)
names(elev) <- v.site 
cnet_2000 <- cnet

#########################################################
#########################################################
#
# 1. Elevational patterns. 
#
#########################################################
#########################################################

# 1.1 Life-cycle components and elasticities

# 1.1.1 Combine observed values into a data.frame

names.var <- c("mean.surv","mean.retro","mean.stasis","mean.progr","mean.pfruit","mean.nfruits","mean.recr.size",
               "mean.surv","mean.retro","mean.stasis","mean.progr","mean.pfruit","mean.nfruits","mean.recr.size",
               dimnames(elast)[[2]][seq(1,15,2)][-5],
               dimnames(elast)[[2]][seq(2,16,2)][-5])
names.rate <- c("Mean S","Mean G-","Mean G=","Mean G+","Mean F0","Mean F1","Mean F2",
                "CV S","CV G-","CV G=","CV G+","CV F0","CV F1","CV F2",
                "Elast mu S","Elast mu G-","Elast mu G=","Elast mu G+","Elast mu F0","Elast mu F1","Elast mu F2",
                "Elast sigma S","Elast sigma G-","Elast sigma G=","Elast sigma G+","Elast sigma F0","Elast sigma F1","Elast sigma F2")

i <- 1
data <- 
  data.frame(
    Site  = row.names(get(names.var[i])),
    med   = get(names.var[i]) %>% apply(1,mean,na.rm=T),
    low   = get(names.var[i]) %>% apply(1,quantile,0.025,na.rm=T),
    high  = get(names.var[i]) %>% apply(1,quantile,0.975,na.rm=T),
    Rate  = names.rate[i]
  )

for (i in 2 : 28) {
  typ.var <- substr(names.rate[i],1,2)
  if (typ.var == "Me") {
    data <- rbind(data,
                  data.frame(
                    Site  = row.names(get(names.var[i])),
                    med   = get(names.var[i]) %>% apply(1,mean,na.rm=T),
                    low   = get(names.var[i]) %>% apply(1,quantile,0.025,na.rm=T),
                    high  = get(names.var[i]) %>% apply(1,quantile,0.975,na.rm=T),
                    Rate  = names.rate[i]
                  ))
    next
  }
  if (typ.var == "CV") {
    data <- rbind(data,
                  data.frame(
                    Site = row.names(get(names.var[i])),
                    med  = get(names.var[i]) %>% apply(c(1,3),cv) %>% apply(1,mean,na.rm=T),
                    low  = get(names.var[i]) %>% apply(c(1,3),cv) %>% apply(1,quantile,0.025,na.rm=T),
                    high = get(names.var[i]) %>% apply(c(1,3),cv) %>% apply(1,quantile,0.975,na.rm=T),
                    Rate = names.rate[i]
                  ))
    next
  }
  if (typ.var == "El") {
    data <- rbind(data,
                  data.frame(
                    Site  = row.names(elast),
                    med   = elast[,names.var[i],] %>% apply(1,mean,na.rm=T),
                    low   = elast[,names.var[i],] %>% apply(1,quantile,0.025,na.rm=T),
                    high  = elast[,names.var[i],] %>% apply(1,quantile,0.975,na.rm=T),
                    Rate  = names.rate[i]
                  ))
    next
  }
}

# Summary for descriptors and lcc averaged over sites
data %>% group_by(Rate) %>% summarise(mean=mean(med)) -> summary.lcc


# 1.1.2 Predict values from linear regressions

newdat <- data.frame(elev=seq(900,3000,100))
nrepl <- 2000

# Means
names.var <- c("mean.surv","mean.retro","mean.stasis","mean.progr","mean.pfruit","mean.nfruits","mean.recr.size")
names.rate <- c("Mean S","Mean G-","Mean G=","Mean G+","Mean F0","Mean F1","Mean F2")
ypred <- array(NA,c(dim(newdat)[1],nrepl*6,7))
coef.b <- array(NA,c(nrepl*6,7))
for (i.var in 1 : 7) {
  i.seq <- 1
  for (i.repl in 1 : nrepl) {
    if(i.repl %% 100 == 0) {
      cat(names.rate[i.var],i.repl,"\n")
      flush.console()
    }
    for (i.year in 1 : 6) {
      mod <- lm(get(names.var[i.var])[,i.year,i.repl]~elev)
      coef.b[i.seq,i.var] <- coef(mod)[2]
      ypred[,i.seq,i.var] <- predict(mod,newdat)
      i.seq <- i.seq + 1
    }
  }
}
data.pred.mean <- data.frame(
  elev = newdat$elev,
  med = as.vector(apply(ypred,c(1,3),mean)),
  low = as.vector(apply(ypred,c(1,3),quantile,0.025)),
  high = as.vector(apply(ypred,c(1,3),quantile,0.975)),
  Rate = rep(names.rate,each=dim(newdat)[1])
)
data.coef.mean <- data.frame(
  med = as.vector(apply(coef.b,2,mean)),
  low = as.vector(apply(coef.b,2,quantile,0.025)),
  high = as.vector(apply(coef.b,2,quantile,0.975)),
  Rate = names.rate
)


# CV
names.var <- c("mean.surv","mean.retro","mean.stasis","mean.progr","mean.pfruit","mean.nfruits","mean.recr.size")
names.rate <- c("CV S","CV G-","CV G=","CV G+","CV F0","CV F1","CV F2")
ypred <- array(NA,c(dim(newdat)[1],nrepl,7))
coef.b <- array(NA,c(nrepl,7))
for (i.var in 1 : 7) {
  for (i.repl in 1 : nrepl) {
    if(i.repl %% 100 == 0) {
      cat(names.rate[i.var],i.repl,"\n")
      flush.console()
    }
    yobs <- get(names.var[i.var])[,,i.repl] %>% apply(1,cv)
    mod <- lm(yobs~elev)
    coef.b[i.repl,i.var] <- coef(mod)[2]
    ypred[,i.repl,i.var] <- predict(mod,newdat)
    i.seq <- i.seq + 1
  }
}
data.pred.cv <- data.frame(
  elev = newdat$elev,
  med = as.vector(apply(ypred,c(1,3),mean)),
  low = as.vector(apply(ypred,c(1,3),quantile,0.025)),
  high = as.vector(apply(ypred,c(1,3),quantile,0.975)),
  Rate = rep(names.rate,each=dim(newdat)[1])
)
data.coef.cv <- data.frame(
  med = as.vector(apply(coef.b,2,mean)),
  low = as.vector(apply(coef.b,2,quantile,0.025)),
  high = as.vector(apply(coef.b,2,quantile,0.975)),
  Rate = names.rate
)

# Elasticities
names.var <- c(dimnames(elast)[[2]][seq(1,15,2)][-5],
               dimnames(elast)[[2]][seq(2,16,2)][-5])
names.rate <- c("Elast mu S","Elast mu G-","Elast mu G=","Elast mu G+","Elast mu F0","Elast mu F1","Elast mu F2",
                "Elast sigma S","Elast sigma G-","Elast sigma G=","Elast sigma G+","Elast sigma F0","Elast sigma F1","Elast sigma F2")
ypred <- array(NA,c(dim(newdat)[1],nrepl,14))
coef.b <- array(NA,c(nrepl,14))
for (i.var in 1 : 14) {
  for (i.repl in 1 : nrepl) {
    if(i.repl %% 100 == 0) {
      cat(names.rate[i.var],i.repl,"\n")
      flush.console()
    }
    yobs <- elast[,names.var[i.var],i.repl]
    mod <- lm(yobs~elev)
    coef.b[i.repl,i.var] <- coef(mod)[2]
    ypred[,i.repl,i.var] <- predict(mod,newdat)
    i.seq <- i.seq + 1
  }
}
data.pred.elast <- data.frame(
  elev = newdat$elev,
  med = as.vector(apply(ypred,c(1,3),mean)),
  low = as.vector(apply(ypred,c(1,3),quantile,0.025)),
  high = as.vector(apply(ypred,c(1,3),quantile,0.975)),
  Rate = rep(names.rate,each=dim(newdat)[1])
)
data.coef.elast <- data.frame(
  med = as.vector(apply(coef.b,2,mean)),
  low = as.vector(apply(coef.b,2,quantile,0.025)),
  high = as.vector(apply(coef.b,2,quantile,0.975)),
  Rate = names.rate
)


# 1.1.3 Combine all predicted values and coefficients

data.pred <- rbind(
  data.pred.mean,
  data.pred.cv,
  data.pred.elast)
data.coef <- rbind(
  data.coef.mean,
  data.coef.cv,
  data.coef.elast)
data.coef <- mutate(data.coef,significant=sign(low)*sign(high))
data.coef$linetype <- "solid"
data.coef$linetype[data.coef$significant==-1] <- "dotted"

# Reorder rates 
data.pred$Rate <- as.character(data.pred$Rate)
data.pred$Rate <- factor(data.pred$Rate,
                         levels=c("Mean S","Mean G-","Mean G=","Mean G+","Mean F0","Mean F1","Mean F2",
                                  "CV S","CV G-","CV G=","CV G+","CV F0","CV F1","CV F2",
                                  "Elast mu S","Elast mu G-","Elast mu G=","Elast mu G+","Elast mu F0","Elast mu F1","Elast mu F2",
                                  "Elast sigma S","Elast sigma G-","Elast sigma G=","Elast sigma G+","Elast sigma F0","Elast sigma F1","Elast sigma F2")
)


save(data,file="Data_for_Figure_1.RData")  

# 1.1.4 Plot everything
load("Data_for_Figure_1.RData")
data$elev <- elev[match(data$Site,names(elev))]

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


dose.labs <- c("D0.5", "D1", "D2")
names(dose.labs) <- c("0.5", "1", "2")

png("Figure 1.png",width=11,height=20,units="cm",res=600)
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



# 1.2 Stochastic growth rates and life-expectancies

# 1.2.1 Combine observed values into a data.frame

data <- rbind(
    data.frame(
        Site  = row.names(r),
        med   = r %>% apply(1,mean,na.rm=T),
        low   = r %>% apply(1,quantile,0.025,na.rm=T),
        high  = r %>% apply(1,quantile,0.975,na.rm=T),
        Rate  = "log lambda"
    ),
    data.frame(
        Site  = row.names(Lexp.mean),
        med   = Lexp.mean %>% apply(1,mean,na.rm=T),
        low   = Lexp.mean %>% apply(1,quantile,0.025,na.rm=T),
        high  = Lexp.mean %>% apply(1,quantile,0.975,na.rm=T),
        Rate  = "Life expectancy"
    )
)

# 1.2.2 Predict values from linear regressions

newdat <- data.frame(elev=seq(900,3000,100))
nrepl <- 2000

names.var <- c("r","Lexp.mean")
names.rate <- c("log lambda","Life expectancy")
ypred <- array(NA,c(dim(newdat)[1],nrepl,length(names.var)))
coef.b <- array(NA,c(nrepl,length(names.var)))
for (i.var in 1 : length(names.var)) {
    for (i.repl in 1 : nrepl) {
        if(i.repl %% 100 == 0) {
            cat(names.rate[i.var],i.repl,"\n")
            flush.console()
        }
        mod <- lm(get(names.var[i.var])[,i.repl]~elev)
        coef.b[i.repl,i.var] <- coef(mod)[2]
        ypred[,i.repl,i.var] <- predict(mod,newdat)
    }
}

data.pred <- data.frame(
    elev = newdat$elev,
    med = as.vector(apply(ypred,c(1,3),mean)),
    low = as.vector(apply(ypred,c(1,3),quantile,0.025)),
    high = as.vector(apply(ypred,c(1,3),quantile,0.975)),
    Rate = rep(names.rate,each=dim(newdat)[1])
)
data.coef <- data.frame(
    med = as.vector(apply(coef.b,2,mean)),
    low = as.vector(apply(coef.b,2,quantile,0.025)),
    high = as.vector(apply(coef.b,2,quantile,0.975)),
    Rate = names.rate
)

# 1.2.3 Combine all predicted values and coefficients

data.coef <- mutate(data.coef,significant=sign(low)*sign(high))
data.coef$linetype <- "solid"
data.coef$linetype[data.coef$significant==-1] <- "dotted"

data.pred$Rate <- as.character(data.pred$Rate)
data.pred$Rate <- factor(data.pred$Rate,
                         levels=c("log lambda","Life expectancy"))

# 1.2.4 Plot everything

data$elev <- elev[match(data$Site,names(elev))]

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

png("Figure S11.png",width=11,height=5.7,units="cm",res=600)
ggplot(data,aes(x=elev,y=med, ymin=low, ymax=high)) +
    geom_ribbon(data=data.pred,fill="gray") +
    geom_line(data=data.pred,colour="blue",aes(linetype=Rate),size=1) +
    scale_linetype_manual(values=data.coef$linetype) +
    
    geom_errorbar(width=0,size=0.5) +
    geom_point(size=1.5) + 
    labs(y="",x="") +
    scale_x_continuous(breaks=c(1000,2000,3000)) +
    labs(y="",x="Elevation (m)") +
    facet_wrap(~Rate,scales="free_y",ncol=2,dir="v")
dev.off()



#########################################################
#########################################################
#
# Temporal covariation in life-cycle components
# 
#########################################################
#########################################################

nrepl <- dim(mean.surv)[3]
ncor <- (7*8)/2 - 7
corr <- array(NA,c(ncor,6,nrepl))
for (i.site in 1 : length(v.site)) {
  for (i.repl in 1 : nrepl) {
    dat <- cbind(
      mean.surv[i.site,,i.repl],
    mean.retro[i.site,,i.repl],
    mean.stasis[i.site,,i.repl],
    mean.progr[i.site,,i.repl],
    mean.pfruit[i.site,,i.repl],
    mean.nfruits[i.site,,i.repl],
    mean.recr.size[i.site,,i.repl]
    )
    colnames(dat) <- c("S","G-","G=","G+","F0","F1","F2")
    temp <- cor(dat)
    temp[upper.tri(temp,diag=T)] <- NA
    temp <- as.data.frame(as.table(temp))
    temp <- temp[!is.na(temp$Freq),]
    corr[,i.site,i.repl] <- temp$Freq
  }
}
dimnames(corr) <- list(pair=paste0("(",temp$Var1,",",temp$Var2,")"),
                       site=v.site,
                       repl=c(1:nrepl))

# Calculate mean, lower and upper confidence limits and store them in data
data.frame(
  cor = paste(temp$Var1,temp$Var2,sep=" vs "),
  med = as.vector(apply(corr,c(1,2),mean)),
  low = as.vector(apply(corr,c(1,2),quantile,0.025)),
  high = as.vector(apply(corr,c(1,2),quantile,0.975)),
  site = rep(v.site,each=ncor)
) %>% as_tibble -> data


# Plot correlations per site
# Update theme to have vertical x text
theme_update(axis.text.x = element_text(size=11, color="black", angle=90, vjust=0.5),
             legend.position="none")
plots <- list()

for (i.site in 1 : 6) {
  
  # Select values of this site only and arrange them in increasing order
  data %>%
    filter(site==v.site[i.site]) %>%
    arrange(med) %>%
    mutate(sign.lo = high<0,
           sign.hi = low > 0,
           sign = "NS") ->
    data.site
  
  # Set values for negative, ns and pos
  data.site$sign[data.site$sign.lo] <- "Neg"
  data.site$sign[data.site$sign.hi] <- "Pos"
  data.site$sign <- factor(data.site$sign, levels=c("Neg","NS","Pos"))
  
  # Keep order for correlations
  data.site$cor <- factor(data.site$cor,levels=as.character(data.site$cor))
  
  plots[[i.site]] <- ggplot(data.site,aes(x=cor,y=med,col=sign)) + 
    geom_errorbar(data=data.site,mapping=aes(x=cor, ymin=low, ymax=high),
                  width=0.5) +
    geom_point() + 
    geom_hline(yintercept=0,linetype=2) + 
    labs(title=v.site[i.site],y="",x="") +
    scale_color_manual(values=c("red", "gray", "blue"), drop=FALSE)
}

png("Fig_TempCorr.png", width=20, height=30, units="cm", res=300)
marrangeGrob(plots,nrow=3,ncol=2,top="")
dev.off()

# Test correlations between (correlation coefficients between life-cycle components) and (elevation)
elev <- c(930,1480,1980,2090,2590,2930)
corel <- array(NA,c(ncor,nrepl))
for (i.pair in 1 : ncor) {
  corel[i.pair,] <- apply(corr[i.pair,,],2,cor,elev)
}
dimnames(corel) <- list(pair=paste0("(",temp$Var1,",",temp$Var2,")"),
                        repl=c(1:nrepl))
tibble(
  Pair = dimnames(corel)[[1]],
  med = apply(corel,1,mean,na.rm=T),
  low = apply(corel,1,quantile,0.025,na.rm=T),
  high = apply(corel,1,quantile,0.975,na.rm=T)
) -> data


# To have them in the order specified above (instead of alphabetic order)
data$Pair <-factor(data$Pair,levels=data$Pair)

#png("Fig_CorrElev_TempCor_.png", height=25, width = 15, units="cm", res=300)

  ggplot(data,aes(x=med,y=Pair)) +
    geom_errorbarh(data=data,mapping=aes(y=Pair,xmin=low,xmax=high),height=0.5) +
    geom_point(size=3) +
    geom_vline(xintercept=0,linetype="dashed") +
    labs(title="Correlation with elevation",x="Pearson's correlation",x="")

dev.off()



#########################################################
#########################################################
#
# stable stage distribution
# 
#########################################################
#########################################################

data <- list()
for (i.site in 1 : length(v.site)) {
  data[[i.site]] <- data.frame(
    stage = y,
    me = apply(ws[,i.site,],1,mean),
    lo = apply(ws[,i.site,],1,quantile,0.025),
    hi = apply(ws[,i.site,],1,quantile,0.975)
  )
}

# Find stage (size) corresponding to median and 90% distribution
size50 <- size90 <- vector()
for (i.site in 1 : length(v.site)){
  size50[i.site] <- which(cumsum(data[[i.site]]$me)>0.5)[1]
  size90[i.site] <- which(cumsum(data[[i.site]]$me)>0.90)[1]
}

theme_set(theme_classic(base_size=9, base_line_size = 0.25, base_rect_size = 0.25))
theme_update(axis.text.x = element_text(size=15, color="black"),
             axis.text.y = element_text(size=15, color="black"),
             axis.title.x = element_text(size=15, color="black"),
             plot.title = element_text(size = 15, color = "black", hjust=0.5),
             plot.margin = margin(.5, .5, 0, 0, "cm"))

plots <- list()
for (i.site in 1 : length(v.site)){
  plots[[i.site]] <- ggplot(data[[i.site]], aes(x=stage, y=me, ymin=lo, ymax=hi)) +
    geom_ribbon(fill="gray") +
    geom_line() + 
    labs(title=v.site[i.site],y="",x="") + 
    geom_vline(xintercept=size50[i.site], linetype=3) + 
    geom_vline(xintercept=size90[i.site], linetype=2) + 
    coord_cartesian(xlim = c(1,20), ylim=c(0,0.35))
}

# Plot
png("Fig_ws.png",width=20,height=30,units="cm",res=300)
marrangeGrob(plots,nrow=3,ncol=2)
dev.off()




#########################################################
#########################################################
#
# SLTRE (with F0, F1, F2)
# 
#########################################################
#########################################################



#######################################################################
# First we take the mean of cnet
#######################################################################

cnet <- apply(cnet_2000,c(1,2,3),mean)
# Reminder: dimensions of cnet are:
# par = c("S","G-","G=","G+","F0","F1","F2"),
# stat =c("mu","e","cc","rho"),
# site = v.site

# Calculate Relative magnitude of contributions for each site
# Total effects of each life-cycle component as percent of summed absolute contributions
tel <- apply(abs(cnet),c(1,3),mean) # Total (absolute) effect of each life-cycle component in each site
tel.p <- apply(tel, 2, function(x) {x/sum(x)*100}) # Total (absolute) effect of each life-cycle component in each site in percent
write.csv(t(tel.p),file="Total_effect_l.csv")
# Total effects of each statistic as percent of summed absolute contributions
tes <- apply(abs(cnet),c(2,3),mean) # Total (absolute) effect of each statistic in each site
tes.p <- apply(tes, 2, function(x) {x/sum(x)*100})
write.csv(t(tes.p),file="Total_effect_s.csv")

png("Fig_SLTRE_Effects_LifeCycleComp_1.png",width=8,height=8,units="cm",res=300)
par(xpd=NA,cex=0.7)
pal <- brewer.pal(7,"Set1")
barplot(tel.p,col=pal,border=NA,las=1)
legend(0,-15,c("S","G-","G=","G+","F0","F1","F2"),fill=pal,bty="n",ncol=4)
dev.off()

png("Fig_SLTRE_Effects_stat_1.png",width=8,height=8,units="cm",res=300)
par(xpd=NA,cex=0.7)
pal <- brewer.pal(4,"Accent")
barplot(tes.p,col=pal,border=NA,las=1)
legend(0,-20,c(expression(italic("\u03BC")),expression(italic(E)),
               expression(italic(CV)),expression(rho)),fill=pal,bty="n",ncol=4)
dev.off()

# Plot
r.plot <- r %>% apply(1,mean,na.rm=T)
Lexp.plot <- Lexp.mean %>% apply(1,mean,na.rm=T)
v.site.elev <- c("BRU (930m)", "CHA (1480m)", "VIL (1980m)", "LAU (2090m)", "GAL (2590m)", "PIC (2930m)")
png("Fig_SLTRE_BySite.png",width=15,height=9,units="cm",res=300)
pal <- brewer.pal(4, "Accent")
par(mfrow=c(2,3), mar=c(3,3,2,1), cex.axis=0.75)
for (i.site in 1 : length(v.site)) {
  c1 <- c2 <- c3 <- cnet[,,i.site]
  c2[c2>0] <- 0 
  c3[c3<0] <- 0 
  myRange <- c( min(rowSums(c2)), max(rowSums(c3)) )
  myRange <- c(-0.2,0.45)
  barplot(t(c2),col=pal,space=0,border=NA,axisnames=FALSE, ylim=myRange,main=v.site.elev[i.site],las=1)
  barplot(t(c3),col=pal,space=0,border=NA,axisnames=FALSE, add=T, axes=F)
  axis(1, at=c(0:6)+0.5, lab= row.names(cnet), las=1)
  text(5.5,0.4,bquote(log~lambda[s]*''==''*.(round(r.plot[i.site],2))))
  #text(5.5,0.3,bquote(italic(L)[exp]*''==''*.(round(Lexp.plot[i.site],1))))

}
dev.off()



#####################################################################
# Correlations of SLTRE contributions  with elevation
# Contributions of means of each life-cycle component
# This is the old way (correlations) while I should use the new way (regression on bootstrapped) values
#####################################################################
names.descr <- dimnames(cnet_2000)[[2]]
for (i in 1 : 4) {
  cnet_descr <- cnet_2000[,i,,]
  corel.descr <- apply(cnet_descr,c(1,3),cor,elev)
  tibble(
    Variable = dimnames(corel.descr)[[1]],
    med = apply(corel.descr,1,mean,na.rm=T),
    low = apply(corel.descr,1,quantile,0.025,na.rm=T),
    high = apply(corel.descr,1,quantile,0.975,na.rm=T)
  )-> data.corel
  # To have them in the order specified above (instead of alphabetic order)
  data.corel$Variable <-factor(data.corel$Variable,levels=data.corel$Variable)
  png(paste0("Fig_SLTRE_corrElev_",i,".png"), height=10, width = 10, units="cm", res=300)
  print(
    ggplot(data.corel,aes(x=med,y=Variable)) +
      geom_errorbarh(data=data.corel,mapping=aes(y=Variable,xmin=low,xmax=high),height=0.5) +
      geom_point(size=3) +
      geom_vline(xintercept=0,linetype="dashed") +
      labs(title=paste("Contributions",names.descr[i]),x="Pearson's correlation",x="")
  )
  dev.off()
}



############################################################################################################
############################################################################################################
#
# Test of demographic compensation
#
############################################################################################################
############################################################################################################

# In order to test the significance of the correlations, we repeat the correlation test over the 2000 resmapled datasets

cor.mat.sp <- cor.mat.pe  <- array(NA,c(28,28,2000))
for (i.repl in 1 : 2000) { 
  cnet <- cnet_2000[,,,i.repl]
  temp <- aperm(cnet,c(3,1,2))
  # Keeping all the four descriptors
  data <- abind(temp[,,1],temp[,,2],temp[,,3],temp[,,4]) 
  colnames(data) <- c(paste0(colnames(temp),".me"),
                      paste0(colnames(temp),".e"),
                      paste0(colnames(temp),".cv"),
                      paste0(colnames(temp),".rho"))
  cor.mat.sp[,,i.repl] <- cor(data,method="spearman")
  }
dimnames(cor.mat.sp) <- dimnames(cor.mat.pe) <- list(var1 = colnames(data), var2 = colnames(data), repl = c(1:2000) )

cor.mat <- cor.mat.sp

med <- apply(cor.mat,c(1,2),mean)
low <- apply(cor.mat,c(1,2),quantile,0.025)
high <- apply(cor.mat,c(1,2),quantile,0.975)

# Customized panel function:
# data contains one replicate of the bootstrap, but it is not used to calculate the correlation
# It is only used to find the [i,j] indices of the couple of variable being considered, by comparing 
# the values contained in "data" with the ones extracted by corrgram and passed to the panel function as "x" and "y" arguments
# the correlation is already stored in "med". its "significance" is assessed using "low" and "high"
panel.fct <- function(x, y, corr = NULL, col.regions, cor.method, digits = 2, cex.cor, data1, med, low, high, ...) {
  # Find i and j
  i <- which(apply(data,2,function(z){all(x==z)}))
  j <- which(apply(data,2,function(z){all(y==z)}))
  corr <- med[i,j] #a$estimate
  signif <- 0; if (low[i,j]>0 | high[i,j] <0) signif <- 1
  auto <- missing(cex.cor)
  usr <- par("usr")
  ncol <- 14
  pal <- col.regions(ncol)
  col.ind <- as.numeric(cut(corr, breaks = seq(from = -1, to = 1, length.out = ncol + 1), include.lowest = TRUE))
  abscorr <- formatC(abs(corr), digits = digits, format = "f")
  med <- formatC(med, digits = digits, format = "f")
  if (signif==1) col.border <- "black" else col.border <- "lightgray"
  if (auto) cex.cor <- 0.7/strwidth(abscorr)
  rect(usr[1], usr[3], usr[2], usr[4], col = pal[col.ind], border = NA)
  box(col = col.border)
  on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  if (signif==1) text(0.5, 0.5, "*", cex = 2.5, col = "black")
}
colnames(data) <- rep("",dim(data)[2])
png("Figure_3_corrplot.png",height=15,width=15,units="cm",res=300)
corrgram(data,cor.method="spearman",lower.panel=NULL, upper.panel=panel.fct, diag.panel=NULL,
         cex.labels=0.6, gap=0.15, col.regions = colorRampPalette(c("red", "salmon","white", "lightblue", "royalblue")),
         data1=data,med=med, low=low, high=high)
dev.off()

# Draw a colorbar
col.regions <- colorRampPalette(c("red", "salmon","white", "lightblue", "royalblue"))
ncol <- 14
pal <- col.regions(ncol)
breaks = seq(from = -1, to = 1, length.out = ncol + 1)
png("Figure_3_legend.png",height=15,width=15,units="cm",res=300)
image.plot(med, horizontal=T, col=pal,breaks=breaks,legend.only=T)
dev.off()

# Number of positive and negative correlations
num.pos <- length(which(low[upper.tri(low,diag=F)]>0))
num.neg <- length(which(high[upper.tri(high,diag=F)]<0))

# Is the number of positive correlations lower than expetced by chance?
# To derive the distribution, we need to permute the mean SLTRE contributions 1000 times across sites
# (using the mean contributions over the 2000 replicates)
cnet <- apply(cnet_2000,c(1,2,3),mean)
data <- abind(temp[,,1],temp[,,2],temp[,,3],temp[,,4]) 
colnames(data) <- c(paste0(colnames(temp),".me"),
                    paste0(colnames(temp),".e"),
                    paste0(colnames(temp),".cv"),
                    paste0(colnames(temp),".rho"))
nrand <- 1000
nvar <- dim(data)[2]
ncomp <- nvar*(nvar+1)/2-nvar # Number of comparisons: n * (n + 1) / 2 - n
cmr.rand <- cmp.rand <- array(NA,ncomp)
num.neg.rand <- num.pos.rand <- vector()
for (repl in 1: nrand) {
  data.rand <- apply(data,2,sample)
  k <- 1
  for (i in 1 : dim(data.rand)[2]) {
    for (j in i : dim(data.rand)[2]) {
      if (i == j) next
      b <- cor.test(data.rand[,i],data.rand[,j],method="sp")
      cmr.rand[k] <- b$estimate
      cmp.rand[k] <- b$p.value
      k <- k + 1
    }
  }
  id.neg.rand <- which(cmr.rand<0)
  num.neg.rand[repl] <- length(which(cmp.rand[id.neg.rand]<0.05))
  id.pos.rand <- which(cmr.rand>0)
  num.pos.rand[repl] <- length(which(cmp.rand[id.pos.rand]<0.05))
}
save(data, nvar, ncomp, num.pos, num.neg, num.neg.rand, num.pos.rand,
     file="Results_DemogrComp_all4_3F.RData")

par(mfrow=c(1,2))
hist(num.pos.rand,breaks=c(0:30),main="Num positive corr")
abline(v=num.pos)
quantile(num.pos.rand,0.975)
hist(num.neg.rand,breaks=c(0:30),main="Num negative corr")
abline(v=num.neg)
quantile(num.neg.rand,0.975)




#########################################################
#########################################################
#
# Compare correlations between life-cycle components and correlations between SLTRE contributions
#
#########################################################
#########################################################

# We do this on the pairs of life-cycle components showing significant negative SLTRE contributions


#########################
# Mean S and mean F1
#########################
i1 <- 1 # S
i2 <- 6 # F1
j <- 1 # means
# 1. SLTRE contributions
# 1.1 Calculate correlations
cor.sp <- array(NA,2000)
a <- cnet_2000[c(i1,i2),j,,]
for (i.repl in 1 : 2000) { 
  cnet <- t( cnet_2000[c(i1,i2),j,,i.repl] )
  cor.sp[i.repl] <- cor(cnet[,1],cnet[,2],method="spearman")
}
cor.stat <- quantile(cor.sp,c(0.025,0.5,0.975))
# 1.2 Arrange contributions in a dataframe
data.frame( t(apply(a,c(1,2),mean)),
            t(apply(a,c(1,2),quantile,0.025)),
            t(apply(a,c(1,2),quantile,0.975)),
            dimnames(a)[[2]]
) %>% as_tibble -> data
colnames(data) <- c("med.x","med.y","low.x","low.y","high.x","high.y","Site")
data$Site <- as.character(data$Site)
data$Site <- factor(data$Site, levels=v.site)
# 1.3 PLot     
png("Fig_Corr_meanS_meanF1_Contr.png", height=7, width = 7, units="cm", res=300)
ggplot(data,aes(x=med.x, xmin=low.x, xmax=high.x,
                y=med.y, ymin=low.y, ymax=high.y,
                col=Site)) +
  geom_errorbar(width=0,size=0.5) +
  geom_errorbarh(height=0,size=0.5) +
  geom_point(size=1.5) + 
  labs(y="Contribution mean F1",x="Contribution mean S") +
  annotate("text",0.2,0.03,
           label=paste0("Spearman r = ",round(cor.stat[2],2),
                        "\n(95% CI: ",round(cor.stat[1],2),
                        ", ",round(cor.stat[3],2),")"),
           size=2.5)
dev.off()

# 2. Life-cycle components
# 2.1 Calculate correlations
cor.sp <- array(NA,2000)
meanS <- apply(mean.surv,c(1,3),mean)
meanF1 <- apply(mean.nfruits,c(1,3),mean)
for (i.repl in 1 : 2000) { 
    cor.sp[i.repl] <- cor(meanS[,i.repl],meanF1[,i.repl],method="spearman")
}
cor.stat <- quantile(cor.sp,c(0.025,0.5,0.975))
# 2.2 Arrange life-cycle components in a dataframe
data.frame( med.x = apply(meanS,1,mean),
            med.y = apply(meanF1,1,mean),
            low.x = apply(meanS,1,quantile,0.025),
            low.y = apply(meanF1,1,quantile,0.025),
            high.x = apply(meanS,1,quantile,0.975),
            high.y = apply(meanF1,1,quantile,0.975),
            Site=dimnames(meanS)[[1]]
) %>% as_tibble -> data
data$Site <- as.character(data$Site)
data$Site <- factor(data$Site, levels=v.site)
# 2.3 Plot     
png("Fig_Corr_meanS_meanF1_LCC.png", height=7, width = 7, units="cm", res=300)
ggplot(data,aes(x=med.x, xmin=low.x, xmax=high.x,
                y=med.y, ymin=low.y, ymax=high.y,
                col=Site)) +
    geom_errorbar(width=0,size=0.5) +
    geom_errorbarh(height=0,size=0.5) +
    geom_point(size=1.5) + 
    labs(x="Mean S",y="Mean F1") +
    annotate("text",0.6,20,
             label=paste0("Spearman r = ",round(cor.stat[2],2),
                          "\n(95% CI: ",round(cor.stat[1],2),
                          ", ",round(cor.stat[3],2),")"),
             size=2.5)
dev.off()

#########################
# CV of S and mean F0
#########################
i1 <- 1 # S
j1 <- 3 # CV
i2 <- 5 # F0
j2 <- 1 # means
# 1. SLTRE contributions
# 1.1 Calculate correlations
cor.sp <- array(NA,2000)
cnet_cvS <- cnet_2000[i1,j1,,]
cnet_meanF0 <- cnet_2000[i2,j2,,]
for (i.repl in 1 : 2000) { 
    cor.sp[i.repl] <- cor(cnet_cvS[,i.repl],cnet_meanF0[,i.repl],method="spearman")
}
cor.stat <- quantile(cor.sp,c(0.025,0.5,0.975))
# 1.2 Arrange contributions in a dataframe
data.frame( med.x = apply(cnet_cvS,1,mean),
            med.y = apply(cnet_meanF0,1,mean),
            low.x = apply(cnet_cvS,1,quantile,0.025),
            low.y = apply(cnet_meanF0,1,quantile,0.025),
            high.x = apply(cnet_cvS,1,quantile,0.975),
            high.y = apply(cnet_meanF0,1,quantile,0.975),
            Site = dimnames(cnet_cvS)[[1]]
) %>% as_tibble -> data
data$Site <- as.character(data$Site)
data$Site <- factor(data$Site, levels=v.site)
# 1.3 PLot     
png("Fig_Corr_cvS_meanF0_Contr.png", height=7, width = 7, units="cm", res=300)
ggplot(data,aes(x=med.x, xmin=low.x, xmax=high.x,
                y=med.y, ymin=low.y, ymax=high.y,
                col=Site)) +
    geom_errorbar(width=0,size=0.5) +
    geom_errorbarh(height=0,size=0.5) +
    geom_point(size=1.5) + 
    labs(x="Contribution CV S",y="Contribution mean F0") +
    annotate("text",-0.04,0.03,
             label=paste0("Spearman r = ",round(cor.stat[2],2),
                          "\n(95% CI: ",round(cor.stat[1],2),
                          ", ",round(cor.stat[3],2),")"),
             size=2.5)
dev.off()

# 2. Life-cycle components
# 2.1 Calculate correlations
cor.sp <- array(NA,2000)
cvS <- mean.surv %>% apply(c(1,3),cv)
meanF0 <- apply(mean.pfruit,c(1,3),mean)
for (i.repl in 1 : 2000) { 
    cor.sp[i.repl] <- cor(cvS[,i.repl],meanF0[,i.repl],method="spearman")
}
cor.stat <- quantile(cor.sp,c(0.025,0.5,0.975))
# 2.2 Arrange life-cycle components in a dataframe
data.frame( med.x = apply(cvS,1,mean),
            med.y = apply(meanF0,1,mean),
            low.x = apply(cvS,1,quantile,0.025),
            low.y = apply(meanF0,1,quantile,0.025),
            high.x = apply(cvS,1,quantile,0.975),
            high.y = apply(meanF0,1,quantile,0.975),
            Site=dimnames(cvS)[[1]]
) %>% as_tibble -> data
data$Site <- as.character(data$Site)
data$Site <- factor(data$Site, levels=v.site)
# 2.3 Plot     
png("Fig_Corr_cvS_meanF0_LCC.png", height=7, width = 7, units="cm", res=300)
ggplot(data,aes(x=med.x, xmin=low.x, xmax=high.x,
                y=med.y, ymin=low.y, ymax=high.y,
                col=Site)) +
    geom_errorbar(width=0,size=0.5) +
    geom_errorbarh(height=0,size=0.5) +
    geom_point(size=1.5) + 
    labs(x="CV S",y="Mean F0") +
    annotate("text",0.4,0.45,
             label=paste0("Spearman r = ",round(cor.stat[2],2),
                          "\n(95% CI: ",round(cor.stat[1],2),
                          ", ",round(cor.stat[3],2),")"),
             size=2.5)
dev.off()


#########################################################
#########################################################
#
# Save site- and year-specific matrices for data availability
#
#########################################################
#########################################################

names(mat.G) <-
names(mat.P) <-
names(mat.F) <-
names(mat.K) <- v.site
for (i.site in 1 : 6) {
    names(mat.G[[i.site]]) <-
    names(mat.P[[i.site]]) <-
    names(mat.F[[i.site]]) <-
    names(mat.K[[i.site]]) <- v.year
}

save(mat.G, mat.P, mat.F, mat.K, file="Transition_matrices.RData")

