# Look up vegetation data
rm(list=ls())

load("Arabis.Community.Data.RData")

dat.noA$SiteQuad <- factor(rownames(dat.noA))
dat.noA$Site <- factor(substr(dat.noA$SiteQuad,1,3))

veg.stat <- dat.noA
save(veg.stat,file="Vegetation.statistics.RData")




