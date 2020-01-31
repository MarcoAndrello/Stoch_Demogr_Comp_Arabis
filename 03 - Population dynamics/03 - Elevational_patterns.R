# Elevational patterns
# Marco Andrello
# 13/01/2020

rm(list=ls())

library(tidyverse)
library(gridExtra)
library(abind)

load("Results_MPM.RData")

v.site <- c("BRU", "CHA", "VIL", "LAU", "GAL", "PIC")
v.year <- c(2008:2013)
v.stage <- c(1,5,10,15)
elev <- c(930,1480,1980,2090,2590,2930)
names(elev) <- v.site 

# Function to calculate coefficient of variation
cv <- function(x) { sd(x) / mean(x) }

# Theme for plotting
theme_set(theme_classic())
theme_update(axis.text.x = element_text(size=7, color="black"),
             axis.text.y = element_text(size=7, color="black"),
             axis.title.x = element_text(size=7, color="black", margin = margin(0.5,0,0,0, unit="cm")),
             axis.title.y = element_text(size=7, color="black"),
             panel.spacing.x = unit(0.45, "cm"),
             plot.margin = margin(0,0.5,0,0, unit="cm"),
             strip.text.x = element_text(size = 0),
             strip.background.x = element_rect(linetype = "blank"),
             legend.position = "none")


# 1 Life-cycle components and elasticities

# 1.1 Combine observed values into a data.frame

names.var <- c("mean.surv","mean.retro","mean.stasis","mean.progr","mean.pfruit","mean.nfruits","mean.recr.size",
               "mean.surv","mean.retro","mean.stasis","mean.progr","mean.pfruit","mean.nfruits","mean.recr.size",
               dimnames(elast)[[2]][seq(1,13,2)],
               dimnames(elast)[[2]][seq(2,14,2)])
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

# Summary for descriptors and lcc averaged over sites. 
data %>% group_by(Rate) %>% summarise(mean=mean(med)) -> summary.lcc


# 1.2 Predict values from linear regressions

newdat <- data.frame(elev=seq(900,3000,100))
nrepl <- 2000

# Means. 
# Before it fitted year-specific regression, resulting in non-significant relationships except F2
# Now relationships are significant for all variables except F0
names.var <- c("mean.surv","mean.retro","mean.stasis","mean.progr","mean.pfruit","mean.nfruits","mean.recr.size")
names.rate <- c("Mean S","Mean G-","Mean G=","Mean G+","Mean F0","Mean F1","Mean F2")
ypred <- array(NA,c(dim(newdat)[1],nrepl,7))
coef.b <- array(NA,c(nrepl,7))
# For each variable, fit a linear model for each replicate
for (i.var in 1 : 7) {
    for (i.repl in 1 : nrepl) {
        if (i.repl %% 100 == 0) cat("Var:",names.rate[i.var],"Repl:",i.repl,"\n"); flush.console()
        yobs <- rowMeans(get(names.var[i.var])[,,i.repl])
        mod <- lm(yobs~elev)
        coef.b[i.repl,i.var] <- coef(mod)[2]
        ypred[,i.repl,i.var] <- predict(mod,newdat)
    }
}


# Calculate mean and 95%CI of the coefficients
data.coef.mean <- data.frame(
  med = as.vector(apply(coef.b,2,mean)),
  low = as.vector(apply(coef.b,2,quantile,0.025)),
  high = as.vector(apply(coef.b,2,quantile,0.975)),
  Rate = names.rate
)
# Calculate mean and 95%CI of the predicted values
data.pred.mean <- data.frame(
    elev = newdat$elev,
    med = as.vector(apply(ypred,c(1,3),mean)),
    low = as.vector(apply(ypred,c(1,3),quantile,0.025)),
    high = as.vector(apply(ypred,c(1,3),quantile,0.975)),
    Rate = rep(names.rate,each=dim(newdat)[1])
)


# CV
names.var <- c("mean.surv","mean.retro","mean.stasis","mean.progr","mean.pfruit","mean.nfruits","mean.recr.size")
names.rate <- c("CV S","CV G-","CV G=","CV G+","CV F0","CV F1","CV F2")
ypred <- array(NA,c(dim(newdat)[1],nrepl,7))
coef.b <- array(NA,c(nrepl,7))
for (i.var in 1 : 7) {
  for (i.repl in 1 : nrepl) {
    if (i.repl %% 100 == 0) cat("Var:",names.rate[i.var],"Repl:",i.repl,"\n"); flush.console()
    yobs <- get(names.var[i.var])[,,i.repl] %>% apply(1,cv)
    mod <- lm(yobs~elev)
    coef.b[i.repl,i.var] <- coef(mod)[2]
    ypred[,i.repl,i.var] <- predict(mod,newdat)
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
names.var <- c(dimnames(elast)[[2]][seq(1,13,2)],
               dimnames(elast)[[2]][seq(2,14,2)])
names.rate <- c("Elast mu S","Elast mu G-","Elast mu G=","Elast mu G+","Elast mu F0","Elast mu F1","Elast mu F2",
                "Elast sigma S","Elast sigma G-","Elast sigma G=","Elast sigma G+","Elast sigma F0","Elast sigma F1","Elast sigma F2")
ypred <- array(NA,c(dim(newdat)[1],nrepl,14))
coef.b <- array(NA,c(nrepl,14))
for (i.var in 1 : 14) {
    for (i.repl in 1 : nrepl) {
        if(i.repl %% 100 == 0) cat("Var:",names.rate[i.var],"Repl:",i.repl,"\n"); flush.console()
        yobs <- elast[,names.var[i.var],i.repl]
        mod <- lm(yobs~elev)
        coef.b[i.repl,i.var] <- coef(mod)[2]
        ypred[,i.repl,i.var] <- predict(mod,newdat)
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


# 1.3 Combine all predicted values and coefficients

data.pred <- rbind(
  data.pred.mean,
  data.pred.cv,
  data.pred.elast)
data.coef <- rbind(
  data.coef.mean,
  data.coef.cv,
  data.coef.elast)
# Significant relationships are those with concordant signs among lower and upper bounds of CI
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


save(data,data.pred,data.coef,elev,file="Data_for_Figure_1.RData")  

# 1.4 Plot everything
load("Data_for_Figure_1.RData")
data$elev <- elev[match(data$Site,names(elev))]

png("Figure_1_unlabelled.png",width=11,height=20,units="cm",res=600)
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



# 2 Stochastic growth rates

rm(list=ls())

load("Results_MPM.RData")

v.site <- c("BRU", "CHA", "VIL", "LAU", "GAL", "PIC")
v.year <- c(2008:2013)
v.stage <- c(1,5,10,15)
elev <- c(930,1480,1980,2090,2590,2930)
names(elev) <- v.site 

# 2.1 Combine observed values into a data.frame

data <- 
    data.frame(
        Site  = row.names(r),
        med   = r %>% apply(1,mean,na.rm=T),
        low   = r %>% apply(1,quantile,0.025,na.rm=T),
        high  = r %>% apply(1,quantile,0.975,na.rm=T),
        Rate  = "log lambda"
    )

# Report log lambda in Table 1
data.round <- data
data.round$med <- round(data.round$med,1)
data.round$low <- round(data.round$low,1)
data.round$high <- round(data.round$high,1)
write.csv(data.round,file="Stoch_growth_rate.csv")
rm(data.round)

# 2.2 Predict values from linear regressions

newdat <- data.frame(elev=seq(900,3000,100))
nrepl <- 2000

names.var <- c("r")
names.rate <- c("log lambda")
ypred <- array(NA,c(dim(newdat)[1],nrepl,length(names.var)))
coef.b <- array(NA,c(nrepl,length(names.var)))
for (i.var in 1 : length(names.var)) {
    for (i.repl in 1 : nrepl) {
        if (i.repl %% 100 == 0) cat("Var:",names.rate[i.var],"Repl:",i.repl,"\n"); flush.console()
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

# 2.3 Combine all predicted values and coefficients

data.coef <- mutate(data.coef,significant=sign(low)*sign(high))
data.coef$linetype <- "solid"
data.coef$linetype[data.coef$significant==-1] <- "dotted"

data.pred$Rate <- as.character(data.pred$Rate)
data.pred$Rate <- factor(data.pred$Rate,
                         levels=c("log lambda"))

# 2.4 Plot everything

data$elev <- elev[match(data$Site,names(elev))]

png("Figure S5.png",width=5.7,height=5.7,units="cm",res=600)
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

