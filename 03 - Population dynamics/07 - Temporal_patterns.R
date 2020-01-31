# Tenmporal patterns in deterministic lambda in each site
# Marco Andrello
# 08/01/2020

rm(list=ls())

library(tidyverse)

load("Results_MPM.RData")

v.site <- c("BRU", "CHA", "VIL", "LAU", "GAL", "PIC")
v.year <- c(2008:2013)

names.var <- c("mean.surv","mean.retro","mean.stasis","mean.progr","mean.pfruit","mean.nfruits","mean.recr.size","lambda.det")
names.rate <- c("S","G-","G=","G+","F0","F1","F2","lambda_d")

i.site <- 1
i.year <- 1
data <- 
    data.frame(
        Site  = v.site[i.site],
        Year  = v.year[i.year], 
        med   = mean(lambda.det[i.site,i.year,],na.rm=T),
        low   = quantile(lambda.det[i.site,i.year,],0.025,na.rm=T),
        high  = quantile(lambda.det[i.site,i.year,],0.975,na.rm=T),
        Rate = "lambda_d"
    )


for (i.var in 1 : 8) {
    for (i.site in 1 : 6) {
        for (i.year in 1 : 6){
            data <- rbind(data,
                          data.frame(
                              Site  = v.site[i.site],
                              Year  = v.year[i.year], 
                              med   = get(names.var[i.var])[i.site,i.year,] %>% mean(na.rm=T),
                              low   = get(names.var[i.var])[i.site,i.year,] %>%quantile(0.025,na.rm=T),
                              high  = get(names.var[i.var])[i.site,i.year,] %>%quantile(0.975,na.rm=T),
                              Rate  = names.rate[i.var]
                          )
            )
            
        }
    }
}

data <- data[-1,]
data$Year <- factor(data$Year)
rownames(data) <- NULL

# Label the "Rate" factor so that the lables can be parsed
data$Rate <- factor(data$Rate, levels=levels(data$Rate), labels=c("lambda","S",bquote(G^'-'),bquote(G^'='),bquote(G^'+'),bquote(F[0]),bquote(F[1]),bquote(F[2])))

# # To draw horizontal line at 1 only for lambda
# data$yintercept <- 0
# data$yintercept[data$Rate=="lambda_d"] <- 1
# data$colour <- "white"
# data$colour[data$Rate=="lambda_d"] <- "black"

# Theme for plotting
theme_set(theme_classic())
theme_update(axis.text.x = element_text(size=7, color="black"),
             axis.text.y = element_text(size=7, color="black"),
             axis.title.x = element_text(size=7, color="black"),
             axis.title.y = element_text(size=7, color="black"),
             panel.spacing.x = unit(0.45, "cm"),
             plot.margin = margin(0,0.5,0,0, unit="cm"),
             strip.text.x = element_text(size = 9),
             strip.background.x = element_rect(linetype = "blank"),
             legend.position = "right",
             legend.direction = "vertical")

png("Figure_S4.png",width=17.5,height=15,units="cm",res=600)
ggplot(data,aes(x=Year,y=med, ymin=low, ymax=high, col=Site)) +
    geom_point(size=1.5,position=position_dodge(width=0.75)) +
    geom_errorbar(width=0,size=0.5,position=position_dodge(width=0.75)) +
    geom_vline(xintercept=c(1:5)+0.5,linetype="dashed") +
    facet_wrap(~Rate,scales="free",nrow=4,dir="v",labeller=label_parsed) +
    scale_colour_brewer(type="div",palette=7)
    # geom_hline(aes(yintercept=data$yintercept, colour=data$colour),linetype="dotted")
dev.off()        



