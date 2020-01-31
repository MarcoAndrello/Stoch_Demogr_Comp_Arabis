# Temporal correlations between life-cycle components
# Marco Andrello
# 24/12/2019

rm(list=ls())

library(tidyverse)
library(gridExtra)

load("Results_MPM.RData")

v.site <- c("BRU", "CHA", "VIL", "LAU", "GAL", "PIC")
v.year <- c(2008:2013)
v.stage <- c(1,5,10,15)
elev <- c(930,1480,1980,2090,2590,2930)
names(elev) <- v.site 

# # Set themes for plotting
theme_set(theme_classic())
theme_update(axis.text.x = element_text(size=11, color="black", angle=90, vjust=0.5),
             axis.text.y = element_text(size=15, color="black"),
             axis.title.x = element_text(size=15, color="black"),
             axis.title.y = element_text(size=15, color="black"),
             plot.title = element_text(size = 15, color = "black", hjust=0.5),
             legend.position="none")



# For each site and repl, calculate temporal correlations between life-cycle components
nrepl <- dim(mean.surv)[3]
ncor <- (7*8)/2 - 7
corr <- array(NA,c(ncor,length(v.site),nrepl))
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
  cor = paste0("(",temp$Var1,", ",temp$Var2,")"),
  # paste(temp$Var1,temp$Var2,sep=" vs "),
  med = as.vector(apply(corr,c(1,2),mean)),
  low = as.vector(apply(corr,c(1,2),quantile,0.025)),
  high = as.vector(apply(corr,c(1,2),quantile,0.975)),
  site = rep(v.site,each=ncor)
) %>% as_tibble -> data


# Plot temporal correlations per site

plots <- list()
for (i.site in 1 : length(v.site)) {
  
  # Select values of this site only and arrange them in increasing order
  data %>%
    filter(site==v.site[i.site]) %>%
    arrange(med) %>%
    mutate(sign.lo = high < 0,
           sign.hi = low > 0,
           sign = "NS") ->
    data.site
  
  # Set values for negative, ns and pos
  data.site$sign[data.site$sign.lo] <- "Neg"
  data.site$sign[data.site$sign.hi] <- "Pos"
  data.site$sign <- factor(data.site$sign, levels=c("Neg","NS","Pos"))
  
  # Keep order for correlations
  data.site$cor <- factor(data.site$cor,levels=as.character(data.site$cor))
  
  i.plot <- c(1,3,5,2,4,6)[i.site] # in this order: BRU, LAU, CHA, GAL, VIL, PIC to have them in the right order in the figure
  cat(i.plot,"\n")
  plots[[i.plot]] <- ggplot(data.site,aes(x=cor,y=med,col=sign)) + 
    geom_errorbar(data=data.site,mapping=aes(x=cor, ymin=low, ymax=high),
                  width=0.5) +
    geom_point() + 
    geom_hline(yintercept=0,linetype=2) + 
    labs(title=v.site[i.site],y="",x="") +
    scale_color_manual(values=c("red", "gray", "blue"), drop=FALSE)
}

png("Figure S5.png", width=30, height=20, units="cm", res=600)
marrangeGrob(plots,nrow=2,ncol=3,top="")
dev.off()


# Total number of positive and negative temporal correlations:
# Count them in the graphs: 21 positive, 21 negative
21 / (ncor*length(v.site))

# Test correlations between (correlation coefficients between life-cycle components) and (elevation)
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

ggplot(data,aes(x=med,y=Pair)) +
    geom_errorbarh(data=data,mapping=aes(y=Pair,xmin=low,xmax=high),height=0.5) +
    geom_point(size=3) +
    geom_vline(xintercept=0,linetype="dashed") +
    labs(title="Correlation with elevation",x="Pearson's correlation",x="")
# None of the correlation coefficients is significantly different from 0




