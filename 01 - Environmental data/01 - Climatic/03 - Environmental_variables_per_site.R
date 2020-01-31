# Environmental variables per site
# Marco Andrello
# 05/06/2019

rm(list=ls())

library(tidyverse)
library(ggpmisc)

v.site <- c("BRU","CHA","VIL","LAU","GAL","PIC")
elev <- c(930,1480,1980,2090,2590,2930)
names(elev) <- v.site 

# Load soilveg data
load(paste0(getwd(),"/../09 - PCA Soil and vegetation/PCA_soil_veg.RData")) 
soilveg$SiteQuad <- row.names(soilveg)
soilveg <- 
  as_tibble(soilveg) %>%
  add_column(Site=substr(.[["SiteQuad"]],1,3)) %>%
  add_column(Quad=paste0("Q",substr(.[["SiteQuad"]],5,5))) %>%
  # Reorder sites
  mutate(Site = factor(Site,levels=c("BRU","CHA","VIL","LAU","GAL","PIC"))) %>%
  mutate(SiteQuad=paste0(Site, "-", Quad)) %>%
  arrange(Site)

data1 <- soilveg %>% group_by(Site) %>% summarise(low = min(N_mean), med = mean(N_mean), high=max(N_mean), Variable="Nitrogen %")
data2 <- soilveg %>% group_by(Site) %>% summarise(low = min(C_mean), med = mean(C_mean), high=max(C_mean), Variable="Carbon %")
data3 <- soilveg %>% group_by(Site) %>% summarise(low = min(pH), med = mean(pH), high=max(pH), Variable="pH")
data4 <- soilveg %>% group_by(Site) %>% summarise(low = min(SLA), med = mean(SLA), high=max(SLA), Variable="SLA")
data5 <- soilveg %>% group_by(Site) %>% summarise(low = min(VegHeight), med = mean(VegHeight), high=max(VegHeight), Variable="Height")

data.soilveg <- rbind(data1,data2,data3,data4,data5)

# Load average climatic data 
load("reconstruct_monthly.RData")
monthly$Pop <- factor(monthly$Pop, levels=c("BRU","CHA","VIL","LAU","GAL","PIC"))
monthly <- monthly %>% select(-NbFreeze2) %>% arrange(Pop) %>% as.data.frame()

## Set outlier temperature to NA (correct aberrant prediction for VIL 2008
monthly$AvgMaxTemp[which(monthly$AvgMaxTemp>32)] <- NA
monthly$AvgAmp[which(monthly$AvgAmp>32)] <- NA


names.var <- c("T.mean","T.min","T.max","T.range","Freeze.n")
i.var <- 1
ii.var <- c(4,6,8,10,12)[i.var]
clim <- 
  data.frame(
    Site=levels(monthly$Pop),
    low = as.vector(by(monthly[,ii.var],monthly$Pop,min,na.rm=T)),
    med = as.vector(by(monthly[,ii.var],monthly$Pop,mean,na.rm=T)),
    high = as.vector(by(monthly[,ii.var],monthly$Pop,max,na.rm=T)),
    Variable = names.var[i.var]
  )
for (i.var in 2 : 5) {
  ii.var <- c(4,6,8,10,12)[i.var]
  clim <- rbind(
    clim,
    data.frame(
      Site=levels(monthly$Pop),
      low = as.vector(by(monthly[,ii.var],monthly$Pop,min,na.rm=T)),
      med = as.vector(by(monthly[,ii.var],monthly$Pop,mean,na.rm=T)),
      high = as.vector(by(monthly[,ii.var],monthly$Pop,max,na.rm=T)),
      Variable = names.var[i.var]
    )
  )
}


# Combine soilveg with clim
data <- rbind(data.soilveg,clim)


############################################################################
# Produce Figure S3
############################################################################

# Add elevation
data$elev <- elev[match(data$Site,names(elev))]

# Order Variable in the order for plotting
data$Variable <- factor(data$Variable,
                        levels=c("T.mean","T.min","T.max","T.range","Freeze.n",
                                 "Nitrogen %","Carbon %","pH","SLA","Height")
)

# Actual plotting
theme_set(theme_classic())
theme_update(axis.text.x = element_text(size=7, color="black"),
             axis.text.y = element_text(size=7, color="black"),
             axis.title.y = element_text(size=7, color="black"),
             plot.title = element_text(size = 7, color = "black", hjust=0.5),
             panel.spacing.x = unit(0.45, "cm"),
             plot.margin = margin(0,0.5,0,0, unit="cm"),
             strip.text.x = element_text(size = 6.5),
             strip.background.x = element_rect(linetype = "blank"),
             legend.position = "none")

ltype <- data.frame(Variable = levels(data$Variable), linetype="dotted")
ltype$linetype <- as.character(ltype$linetype)
ltype$linetype[c(1,2,6,7)] <- "solid"


png("Figure S2.png",width=13.75,height=10,units="cm",res=600)
ggplot(data,aes(x=elev,y=med, ymin=low, ymax=high)) +
  geom_smooth(method="lm",colour="blue",aes(linetype=Variable)) + 
  geom_errorbar(width=0,size=0.5) +
  geom_point(size=1.5) + 
  labs(y="",x="") +
  scale_linetype_manual(values=ltype$linetype) +
  stat_fit_glance(aes(label = sprintf('italic(P)~"="~%.3f',
                                      stat(p.value))),
                  size=2.5,
                  parse = TRUE) + 
  scale_x_continuous(breaks=c(1000,2000,3000)) +
  labs(y="",x="Elevation (m)") +
  facet_wrap(~Variable,scales="free_y",ncol=5,dir="h")
dev.off()


