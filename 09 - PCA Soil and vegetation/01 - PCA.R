
rm(list=ls())
library(ade4)
library(corrgram)
library(factoextra)
# Example plots:
# http://www.sthda.com/english/articles/31-principal-component-methods-in-r-practical-guide/119-pca-in-r-using-ade4-quick-scripts/


# Read environmental variables
load(paste0(getwd(),"/Soil and fecundity/Soil.statistics.RData"))
load(paste0(getwd(),"/Vegetation/Arabis.community.data.new.RData"))
veg.stat <- dat.noA
rm(androsace1,ComSpList,dat,dat.noA,DAtraits.mean,DAtraits.meanH,divgrass,spxplot,spxplotH,spxplotHR,spxplotT)

# Add site and SiteQuad columns to veg.stat
veg.stat$Site <- factor(substr(rownames(veg.stat),1,3))
veg.stat$SiteQuad <- row.names(veg.stat)

# Combine them in a single dataframe
all(soil.stat$SiteQuad==veg.stat$SiteQuad)
data <- cbind(soil.stat,veg.stat)

# Add massif variable
Massif <- rep("Verc",dim(data)[1])
Massif[data$Site=="LAU"] <- "Galib"; Massif[data$Site=="GAL"] <- "Galib"; Massif[data$Site=="PIC"] <- "Galib"
data$Massif <- Massif
data.or <- data[,c(16,1,2,3,5,14,12,15)]
names(data.or)[6:7] <- c("VegHeight","SLA")
rm(data)

# PCA without Seed mass
data <- data.or[,c(3:7)]
names(data) <- c("N%", "C%","pH","Height","SLA")
a <- dudi.pca(data, scannf = FALSE, nf = 3)
a$eig / sum(a$eig)
cumsum(a$eig/sum(a$eig))
# Reverse the sign to ease interpretation
a$li <- a$li * -1
a$c1 <- a$c1 * -1
a$co <- a$co * -1


m <- 3 # Multiplier for the figure size
png("PCA_SoilVeg_Biplot.png",width=5.2*m,height=2.9*m,units="cm",res=300)
fviz_pca_biplot(a, repel = TRUE, col.var = "black", col.ind = data.or$Site,
                label="var", mean.point=F, pch=16,
                legend.title="Site", title="")
dev.off()

# Results for Variables
res.var <- get_pca_var(a)
res.var$coord          # Coordinates
res.var$contrib        # Contributions to the PCs
res.var$cos2           # Quality of representation 

# Results for individuals
res.ind <- get_pca_ind(a)
res.ind$coord          # Coordinates
res.ind$contrib        # Contributions to the PCs
res.ind$cos2           # Quality of representation 

data <- cbind(data,a$li)

corrgram(data,upper.panel=panel.conf)
soilveg <- data
save(soilveg,file="PCA_soil_veg.RData")
data <- soilveg

# For Soilveg1 and SoilVeg2: calculate mean values per site by averaging on plots
load("PCA_soil_veg.RData")
soilveg$Site <- factor(substr(rownames(soilveg),1,3))
data.site <- data.frame(Site=levels(soilveg$Site))
data.site$Axis1 <- as.vector(by(soilveg$Axis1,soilveg$Site,mean))
data.site$Axis2 <- as.vector(by(soilveg$Axis2,soilveg$Site,mean))
data.site.soilveg <- data.site
save(data.site.soilveg,file="Data.site.soilveg.RData")
#








# #####################################
# 
# # PCA with Seed mass
# data <- data.or[,c(3:8)]
# a <- dudi.pca(data, scannf = FALSE, nf = 3)
# a$eig / sum(a$eig)
# cumsum(a$eig/sum(a$eig))
# # Reverse the sign to ease interpretation
# a$li <- a$li * -1
# a$c1 <- a$c1 * -1
# a$co <- a$co * -1
# 
# s.class(a$li, data.or$Site, cpoint = 1)
# s.class(a$li[,c(2,3)], data.or$Site, cpoint = 1)
# s.arrow(a$c1, lab = names(a$tab))
# s.corcircle(a$co, lab = names(a$tab), box = TRUE)
# 
# data <- cbind(data,a$li)
# 
# corrgram(data,upper.panel=panel.conf)
# soilveg <- data
# save(soilveg,file="PCA_soil_veg_withSeedMass.RData")
# data <- soilveg
# 
# ###














