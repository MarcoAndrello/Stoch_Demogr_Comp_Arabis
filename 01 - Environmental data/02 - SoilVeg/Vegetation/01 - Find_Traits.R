#########################################################
#	Find Traits for Arabis alpina quadrats
#########################################################
# Marta Carboni
# January 2019

##### Set directory
setwd("~/Dropbox/0_leca/projects/Marco Arabis") #SY Mac
#setwd("C:/Users/Marta Carboni/Dropbox/0_leca/projects/Marco Arabis") #hamlet


##### Packages
library(reshape2)		  # dcast
library(Taxonstand)		  # find synonyms
library(FD)				  # CWMs and FDs
library(R.utils)		  #capitalize

##### Load Data

#final data
#load("Arabis.Community.Data.RData")

# Communities with Arabis ------------------------------------

com<-read.table("Species_checklist_Arabis.txt", sep="\t", head=T) 
head(com) ; dim(com)

# Complete species list

ComSpList <- read.table(file="Arabis_Community_SpList_corrected.txt", sep="\t") 
head(ComSpList) ; dim(ComSpList)

# Trait databases ---------------------------------------------

# Traits from Androsace (Database di Wilfried/Julien)

# version 1 - (data of all species in French alps available in 2010, as used by Isa B.)
androsace1 = read.table("Traits/All_Traits_allSpecies_29Sept2010.txt", h=T, sep="\t", dec=".")
head(androsace1)
androsace1$Species_nosub = sapply(strsplit(as.character(androsace1$Genus_species), split="_"),function(x) paste(x[1:2], collapse=" "))
#androsace1$nom_sp= paste("X", androsace1$No_Taxon_CBNA, sep="") 

# # version 2 (only Ecrins?)
# androsace2<-read.table("Traits/AllTraits_quanti.txt", head=T,sep = "\t") #no sps names, only cbna code
# head(androsace2)
# # load the list of taxa in Alps with CBNA code to match with species names
# txList <- read.table('Traits/txList_EcrinsPeriph.txt',h=T, sep="\t")
# head(txList) ; dim(txList)
# setdiff(androsace2$androsace2,androsace1$No_Taxon_CBNA) #no extra species! 

# version 3 - (2018 data from Wilfried)
androsace3 = read.csv("Traits/traits.Marco.csv", h=T, sep="\t") #, dec="."
androsace3$value2 = as.numeric(as.character(androsace3$value))
head(androsace3)
#unique(androsace3$Species.original) 
#unique(androsace3$source) 

#by source and subspecies
androsace3W <- dcast(androsace3[,c("Species.original","libelle","source","CODE","value2")], Species.original + libelle + source ~  CODE, value.var="value2") #, mean, na.rm=T
androsace3W_GF <- dcast(androsace3[,c("Species.original","libelle","source","CODE","value")], Species.original + libelle + source ~  CODE, value.var="value") #, mean, na.rm=T
androsace3W$GROWTHF <- androsace3W_GF$GROWTHF
head(androsace3W) ; dim(androsace3W)

#averages aross sources and subspecies
androsace3W_means <- dcast(androsace3[,c("Species.original","libelle","source","CODE","value2")], Species.original ~  CODE, mean, na.rm=T, value.var="value2") #
androsace3W_means = merge(androsace3W_means[,-2:-3], androsace3W_GF[!is.na(androsace3W_GF$GROWTHF),c("Species.original", "GROWTHF")], by="Species.original", all.x=T)
head(androsace3W_means) ; dim(androsace3W_means)

#----------------------------------------------

# Traits from DIVGRASS (originally from TRY)

divgrass=read.csv("Traits/mean_TRAITS_TAXREF_SHORT_20140130b.csv",sep=";",h=T)
head(divgrass)
divgrass$SPECIES<-as.character(divgrass$SPECIES)
divgrass$Species_nosub <-sub(" X ", " ", divgrass$SPECIES)
divgrass$Species_nosub <-sapply(strsplit(tolower(divgrass$Species_nosub), " "),function(x) paste(capitalize(x[1]),x[2], collapse=" "))
divgrass$Species_nosub <- sub(" NA", " sp", divgrass$Species_nosub)
divgrass$Genus_species <- sapply(strsplit(as.character(divgrass$Species_nosub), split=" "),function(x) paste(x[1:2], collapse="_"))

#----------------------------------------------

# Traits from TRY

Traits_Try.mean<-read.table(file="Traits/Traits_Try.mean.txt", sep="\t")

summary(Traits_Try.mean[,c("VegHeight" , "LNC_T14" , "SM_T26" , "SLA_try")]) ; dim(Traits_Try.mean)
# 65 herbaceous species with trait info in NEW androsace + DIVGRASS out of 66 possible
# 65 with height
# 51 sps with SLA, 46 sps with seed mass, 42 sps with LNC


##################################################################
#	Prepare and recover data 
##################################################################

#####	Fix community data ---------------------------------------

# Add simple species names and genus
Genus_species<-sapply(strsplit(as.character(com$Species_id), split=" "),function(x) paste(x[1:2], collapse="_"))
com<-transform(com, Genus_species = Genus_species)
com<-transform(com, Genus = sapply(strsplit(as.character(Genus_species), split="_"),function(x) paste(x[1])))

# Rescale abundance in percentages (+/- median value of range)
com <-transform(com, abund=NA)
com$abund[com$Braun.Blanquet=="r"]<-0.5
com$abund[com$Braun.Blanquet=="+"]<-1
com$abund[com$Braun.Blanquet=="1"]<-3
com$abund[com$Braun.Blanquet=="2"]<-15
com$abund[com$Braun.Blanquet=="3"]<-37.5
com$abund[com$Braun.Blanquet=="4"]<-62.5
com$abund[com$Braun.Blanquet=="5"]<-87.5
head(com) ; dim(com)

# Determine the species list
#ComSpList <- unique(com[,c("Species", "Species_id", "Genus_species", "Genus")]) #unique species list
ComSpList <- unique(com[,c("Species_id","Genus_species", "Genus")]) #unique species list
ComSpList$Species_id.subsp <- as.character(ComSpList$Genus_species) 
ComSpList$Species_id.subsp[ComSpList$Species_id.subsp=="Cerastium_arvense"] <- "Cerastium_arvense subsp.strictum"
ComSpList$Species_nosub <- sapply(strsplit(as.character(ComSpList$Genus_species), split="_"),function(x) paste(x[1:2], collapse=" "))
head(ComSpList); dim(ComSpList)
# 81 unique species after fixing typos


#####	Find accepted names for all species -------------------------------

# In community data

# add name from the plant list (Taxonstand): TAKES LONG TIME!!! (1 min for 82 sps)
#c1 <- TPL(ComSpList$Species_id, corr=TRUE)
#c2 <- TPL(ComSpList$Species_nosub, corr=TRUE)
#c1$New.Species[c1$New.Species=="NA"]<-NA
#head(c1);str(c1)
#c1[,c(1,10,12:14,16,21:22)]
#save(c1, file="Traits/taxonstand_c1")

# # load(file="Traits/taxonstand_c1")
# ComSpList$Species.taxonstand<-paste(c1$New.Genus,c1$New.Species, sep="_")
# head(ComSpList); dim(ComSpList) 
# ComSpList[,c(1:2,6)]
com1 = merge(com, ComSpList[,c('Genus_species', 'Species.taxonstand')], by="Genus_species", all.x=T )
dim(com1); head(com1)

#save species list with corrected names
#write.table(ComSpList, file="Arabis_Community_SpList_corrected.txt", sep="\t")
ComSpList <- read.table(file="Arabis_Community_SpList_corrected.txt", sep="\t") 
names(ComSpList)[1]<-"Species.original"

# In Androsaceae database

# # add name from the plant list (Taxonstand): TAKES LONG TIME!!! (ca. 10 min)
# r2 <- TPL(androsace1$Species_nosub, corr=TRUE)
# r2$New.Species[r2$New.Species=="NA"]<-NA
# head(r2); str(r2)
# save(r2, file="Traits/taxonstand_r2")

load(file="Traits/taxonstand_r2")
androsace1$Species.taxonstand<-paste(r2$New.Genus,r2$New.Species, sep="_")
head(androsace1)


# In DIVGRASS database

# # add name from the plant list (Taxonstand): TAKES LONG TIME!!! (ca. 60 min)
# d2 <- TPL(divgrass$Species_nosub, corr=TRUE)   
# d2$New.Species[r2$New.Species=="NA"]<-NA
# head(d2); str(d2)
# save(d2, file="Traits/taxonstand_d2")

load(file="Traits/taxonstand_d2")
divgrass$Species.taxonstand<-paste(d2$New.Genus,d2$New.Species, sep="_")
head(divgrass)


#####	Create species by plot matrix -------------------------------

# Create species by plot matrix
spxplot<- data.frame(dcast(com1[,c("Species.taxonstand","Id","abund")], Id ~  Species.taxonstand, value.var="abund"), row.names=1)
head(spxplot) ; dim(spxplot)
# 82 sp by 18 sites. / 78 after fixing typos
# Note that A. alpina is missing in site PIC-2



#####	Check which traits are available -------------------------------

##### Androsace OLD

#traits = merge(ComSpList, androsace1[,c("Genus_species.subsp", 'Genus_species', "Plant.type", "SLAall" ,"HautVegall")], by.x="Genus_species", by.y="Genus_species")
#traits2 = merge(ComSpList, androsace1[,c("Genus_species.subsp", 'Genus_species', "Plant.type", "SLAall" ,"HautVegall")], by.x="Species.taxonstand", by.y="Genus_species")
traits3 = merge(ComSpList, androsace1[,c("Genus_species.subsp", 'Species.taxonstand', "Plant.type", "SLAall" ,"HautVegall")], by="Species.taxonstand", all.x=T)
traits3 = traits3[!traits3$Genus_species %in% "Cerastium_arvense" | traits3$Genus_species.subsp %in% "Cerastium_arvense subsp.strictum",]
#dim(traits) ; head(traits) #68 rows
#dim(traits2) ; head(traits2) #79 rows
dim(traits3) ; head(traits3) #88 rows with traits / 100 tot
#setdiff(traits$Species.taxonstand,traits2$Species.taxonstand)
#setdiff(traits2$Species.taxonstand,traits$Species.taxonstand)
traits3[,c(1,3,7:10)]

# How many species with just genus?
traits3$Species.taxonstand[grep("_sp", traits3$Species.taxonstand)]
length(grep("_sp", traits3$Species.taxonstand)) #8 (of which 2 mosses)

#Fix Plant.type, add mosses
traits3$Plant.type<-as.character(traits3$Plant.type)
traits3[is.na(traits3$Plant.type),]
traits3[traits3$Genus=="Alchemilla", "Plant.type"] <- "forb"
traits3[traits3$Genus=="Galium", "Plant.type"] <- "forb"
traits3[traits3$Genus=="Geranium", "Plant.type"] <- "forb"
traits3[traits3$Genus=="Hyeracium", "Plant.type"] <- "forb"
traits3[traits3$Genus=="Moehringia", "Plant.type"] <- "forb"
traits3[traits3$Genus=="Ranunculus", "Plant.type"] <- "forb"
traits3[traits3$Genus_species =="Saxifraga_exarata", "Plant.type"] <- "forb"
traits3[traits3$Genus=="Festuca", "Plant.type"] <- "grass"
traits3[traits3$Genus=="Phalaris", "Plant.type"] <- "grass"
traits3[traits3$Genus=="Dicronum", "Plant.type"] <- "moss"
traits3[traits3$Genus=="Moss", "Plant.type"] <- "moss"
traits3[traits3$Genus=="Mnium", "Plant.type"] <- "moss"

#take average for subspecies
tmp1<-unique(traits3[,c("Species.taxonstand","Genus", "Plant.type")]) #qualitative info
tmp2<-aggregate(traits3[,c("SLAall" , "HautVegall")], by=list(traits3$Species.taxonstand),FUN="mean", na.rm=T) #average traits
traits.mean0 = merge(tmp1[complete.cases(tmp1),], tmp2, by.x="Species.taxonstand", by.y="Group.1") #regroup
dim(traits.mean0) ; head(traits.mean0) 
#55 spp without taxonstand 
#59 with taxonstand just for com 
#69 with taxonstand for both com and androsace / corrected typos (SLA and height)
#78 (with plant type)
#save(traits.mean0, file="Androsaceae_Traits_old.RData")


#Keep only species that have trait values and are not trees
traits.meanH<-traits.mean0[traits.mean0$Plant.type != "tree" & traits.mean0$Plant.type != "moss" ,]
dim(traits.meanH) #72 herbaceous species and no mosses
traits.mean<-traits.meanH[!(is.na(traits.meanH$SLAall)&is.na(traits.meanH$HautVegall)) ,]
traits.mean<-data.frame(traits.mean[,c("Species.taxonstand","SLAall", "HautVegall")], row.names=1)
dim(traits.mean) ; head(traits.mean)
# 45 herbaceous species with trait info in androsace
# 51 with taxonstand for com
# 54 with taxonstand for both 
# out of 66 possible, because: 3 mosses, 3 trees, 6 species only at genus level (78 - 3 -3-6 = 66 species total for which could find trait values)


##### Androsace NEW

Atraits = merge(ComSpList, androsace3W_means, by="Species.original", all.x=T)
Atraits.mean0=Atraits[,c("Species.original","Species.taxonstand","LCCp", "LNC", "PL_VEG_H", "SEEDM",  "SLA", "GROWTHF")]
Atraits.mean0=Atraits.mean0[rowSums(apply(Atraits.mean0,2,function(x)is.na(x)))<6,]       #remove rows that have no trait data
dim(Atraits) ; head(Atraits) #81 rows tot
dim(Atraits.mean0) ; head(Atraits.mean0) #58 rows with some traits

#Add Plant.type, add mosses and traits from androsace
dim(traits.mean0[rowSums(apply(traits.mean0,2,function(x)is.na(x)))<2,])
Atraits.mean0 <- merge(Atraits.mean0, traits.mean0[,c("Species.taxonstand","Plant.type")], by="Species.taxonstand", all.x=T)
Atraits.mean0<-Atraits.mean0[,c("Species.original","Species.taxonstand", "LNC", "SLA", "PL_VEG_H", "SEEDM", "Plant.type")]
dim(Atraits.mean0) ; head(Atraits.mean0)           #58 sps with data from androsaceae (from wil)
summary(Atraits.mean0[,-1:-2])
#58 species with height and plant type info (but still includes trees!)
#51 sps with SLA, 47 sps with seed mass, 36 sps with LNC
# out of 69 possible in sps list, because: 3 mosses, 6 species only at genus level (78 -3-6 = 69 species total for which could find trait values)


#Keep only species that have trait values and are not trees
Atraits.meanH<-Atraits.mean0[Atraits.mean0$Plant.type != "tree" & Atraits.mean0$Plant.type != "moss" ,]
dim(Atraits.meanH) #55 herbaceous species and no mosses
Atraits.mean <-data.frame(Atraits.meanH[,c("Species.taxonstand","Plant.type","SLA","LNC", "PL_VEG_H", "SEEDM")], row.names=1)
dim(Atraits.mean) ; head(Atraits.mean)
summary(Atraits.mean)
# 55 herbaceous species with trait info in NEW androsace 
# out of 66 possible, because: 3 mosses, 3 trees, 6 species only at genus level (78 - 3 -3-6 = 66 species total for which could find trait values)
# 55 sps with height, 49 sps with SLA, 44 sps with seed mass, 33 sps with LNC

##### Divgrass 

# #Dtraits1 = merge(ComSpList, divgrass[,c("Genus_species", "SLA" ,"PH")], by.x="Species.taxonstand", by.y="Genus_species", all.x=T)
# #dim(Dtraits1) ; head(Dtraits1) #81 rows 

# Dtraits = merge(ComSpList, divgrass[,c("Species.taxonstand", "SLA" ,"PH", "SM", "LNC_m")], by="Species.taxonstand", all.x=T)
# dim(Dtraits) ; head(Dtraits) #85 rows 

# #Add Plant.type, add mosses and traits from androsace
# Dtraits <- merge(Dtraits, Atraits.mean0[,c("Species.taxonstand","Plant.type",   "SLA" , "PL_VEG_H", "SEEDM","LNC" )], by="Species.taxonstand", all.x=T)
# dim(Dtraits) ; head(Dtraits) #85 rows 
# Dtraits[,c(-2:-3,-5:-6)]
# #merge(Dtraits[,c(-2:-3,-5:-6)],Dtraits1[,c(1,7:8)], by="Species.taxonstand", all.x=T)

# #take average for subspecies
# Dtmp1<-unique(Dtraits[,c("Species.taxonstand","Genus")]) #qualitative info
# Dtmp2<-aggregate(Dtraits[,c("SLA.x" , "PH","SM","LNC_m")], by=list(Dtraits$Species.taxonstand),FUN="mean", na.rm=T) #average traits
# Dtraits.mean0 = merge(Dtmp1[complete.cases(tmp1),], Dtmp2, by.x="Species.taxonstand", by.y="Group.1") #regroup
# dim(Dtraits.mean0) ; head(Dtraits.mean0) 
# #78 tot

# #Combine with Androsaceae data
# DAtraits.mean0 = merge(Atraits.mean0[,c("Species.taxonstand","Plant.type",   "SLA" , "PL_VEG_H", "SEEDM","LNC" )], Dtraits.mean0[,-2], by="Species.taxonstand", all.x=T, all.y=T) 
# #average SLA values from androsace and divgrass (very similar and same scale!), 
# #DAtraits.mean0$SLAave<-rowMeans(DAtraits.mean0[,c("SLA","SLA.x")],na.rm=T)
# #SLA, keep only androsaceae when available
# DAtraits.mean0$SLA2<-DAtraits.mean0$SLA
# DAtraits.mean0$SLA2[is.na(DAtraits.mean0$SLA)]<-DAtraits.mean0$SLA.x[is.na(DAtraits.mean0$SLA)]
# #height is different scale (cm vs. m, and also not consistent), keep only androsaceae when available
# DAtraits.mean0$Height<-DAtraits.mean0$PL_VEG_H
# DAtraits.mean0$Height[is.na(DAtraits.mean0$PL_VEG_H)]<-DAtraits.mean0$PH[is.na(DAtraits.mean0$PL_VEG_H)]*100
# #Seed mass, keep only androsaceae when available
# DAtraits.mean0$SeedMass<-DAtraits.mean0$SEEDM
# DAtraits.mean0$SeedMass[is.na(DAtraits.mean0$SEEDM)]<-DAtraits.mean0$SM[is.na(DAtraits.mean0$SEEDM)]
# dim(DAtraits.mean0) ; head(DAtraits.mean0) 
# #LNC, keep only androsaceae when available
# DAtraits.mean0$LNC2<-DAtraits.mean0$LNC
# DAtraits.mean0$LNC2[is.na(DAtraits.mean0$LNC)]<-DAtraits.mean0$LNC_m[is.na(DAtraits.mean0$LNC)]
# dim(DAtraits.mean0) ; head(DAtraits.mean0) 


# #Keep only species that have trait values and are not trees
# DAtraits.meanH<-DAtraits.mean0[!DAtraits.mean0$Plant.type %in% c("tree","moss") ,]
# #dim(DAtraits.meanH) #75 herbaceous species and no mosses
# DAtraits.meanH <- DAtraits.meanH[!sapply(strsplit(as.character(DAtraits.meanH$Species.taxonstand), split="_"), function(x) x[2] ) %in% c("sp","spp") ,]
# #dim(DAtraits.meanH) #67 herbaceous species and no mosses/trees and no sps with just genus
# #DAtraits.mean<-DAtraits.meanH[!(is.na(DAtraits.meanH$SLA2)&is.na(DAtraits.meanH$Height)) ,]
# DAtraits.meanH=DAtraits.meanH[rowSums(apply(DAtraits.meanH,2,function(x)is.na(x)))<13,]
# #dim(DAtraits.meanH) #64 herbaceous species and no mosses/trees and no sps with just genus that have some data
# DAtraits.mean <- data.frame(DAtraits.meanH[,c("Species.taxonstand","SLA2","LNC2", "Height", "SeedMass")], row.names=1)
# dim(DAtraits.mean) ; head(DAtraits.mean)
# # 58 species with trait info in divgrass with taxonstand for com
# # 64 species with trait info in divgrass and androsaceae combined (without trees/mosses)
# # out of 66 possible, because: 3 mosses, 3 trees, 6 species only at genus level (78 - 3 -3-6 = 66 species total for which could find trait values)
# summary(DAtraits.mean)
# # 64 herbaceous species with trait info in NEW androsace + DIVGRASS out of 66 possible
# # 64 with height
# # 58 sps with SLA, 54 sps with seed mass, 51 sps with LNC

##### TRY ---------------------------------

head(Traits_Try.mean[,c("Species.taxonstand","VegHeight" , "LNC_T14" , "SM_T26" , "SLA_try")])

#Combine with Androsaceae data
TAtraits.mean0 = merge(Atraits.mean0[,c("Species.taxonstand","Plant.type",   "SLA" , "PL_VEG_H", "SEEDM","LNC" )], 
					   Traits_Try.mean[,c("Species.taxonstand","VegHeight" , "LNC_T14" , "SM_T26" , "SLA_try")], 
					   by="Species.taxonstand", all.x=T, all.y=T) 
#SLA, keep only androsaceae when available
TAtraits.mean0$SLA2<-TAtraits.mean0$SLA
TAtraits.mean0$SLA2[is.na(TAtraits.mean0$SLA)]<-TAtraits.mean0$SLA_try[is.na(TAtraits.mean0$SLA)]
#height is already on same scale (cm), keep only androsaceae when available
TAtraits.mean0$Height<-TAtraits.mean0$PL_VEG_H
TAtraits.mean0$Height[is.na(TAtraits.mean0$PL_VEG_H)]<-TAtraits.mean0$VegHeight[is.na(TAtraits.mean0$PL_VEG_H)]
#Seed mass, keep only androsaceae when available
TAtraits.mean0$SeedMass<-TAtraits.mean0$SEEDM
TAtraits.mean0$SeedMass[is.na(TAtraits.mean0$SEEDM)]<-TAtraits.mean0$SM_T26[is.na(TAtraits.mean0$SEEDM)]
dim(TAtraits.mean0) ; head(TAtraits.mean0) 
#LNC, keep only androsaceae when available
TAtraits.mean0$LNC2<-TAtraits.mean0$LNC
TAtraits.mean0$LNC2[is.na(TAtraits.mean0$LNC)]<-TAtraits.mean0$LNC_T14[is.na(TAtraits.mean0$LNC)]
dim(TAtraits.mean0) ; head(TAtraits.mean0) 


#Keep only species that have trait values and are not trees
TAtraits.meanH<-TAtraits.mean0[!TAtraits.mean0$Plant.type %in% c("tree","moss") ,]
#dim(TAtraits.meanH) #66 herbaceous species and no mosses
TAtraits.meanH <- TAtraits.meanH[!sapply(strsplit(as.character(TAtraits.meanH$Species.taxonstand), split="_"), function(x) x[2] ) %in% c("sp","spp") ,]
#dim(TAtraits.meanH) #66 herbaceous species and no mosses/trees and no sps with just genus
TAtraits.meanH=TAtraits.meanH[rowSums(apply(TAtraits.meanH,2,function(x)is.na(x)))<13,]
#dim(TAtraits.meanH) #65 herbaceous species and no mosses/trees and no sps with just genus that have some data
TAtraits.mean <- data.frame(TAtraits.meanH[,c("Species.taxonstand","SLA2","LNC2", "Height", "SeedMass")], row.names=1)
dim(TAtraits.mean) ; head(TAtraits.mean)
# 65 species with trait info in divgrass and androsaceae combined (without trees/mosses)
# out of 66 possible, because: 3 mosses, 3 trees, 6 species only at genus level (78 - 3 -3-6 = 66 species total for which could find trait values)
summary(TAtraits.mean)
# 65 herbaceous species with trait info in NEW androsace + TRY out of 66 possible
# 65 with height
# 60 sps with SLA, 53 sps with seed mass, 51 sps with LNC


##################################################################
#	Calculate CWMs and community metrics 
##################################################################

# traits.mean0		# all species with plant type infor (from androsaceae old) 
# Atraits.mean0 	# all sps that have Androsaceae new traits 
# Atraits.meanH 	# Androsaceae new traits, but no trees/mosses
# TAtraits.meanH	# Androsaceae new traits, plus TRY new, no trees/mosses

#####	Prepare data -------------------------------

setdiff(traits.mean0$Species.taxonstand,colnames(spxplot))	
setdiff(Atraits.mean0$Species.taxonstand,colnames(spxplot))							# fix "Chenopodium_bonus-henricus"
setdiff(TAtraits.meanH$Species.taxonstand,colnames(spxplot))						# fix "Chenopodium_bonus-henricus"
traits.mean0$Species.taxonstand <- sub("-",".",traits.mean0$Species.taxonstand)	# fix "Chenopodium_bonus-henricus"
Atraits.mean0$Species.taxonstand <- sub("-",".",Atraits.mean0$Species.taxonstand)	# fix "Chenopodium_bonus-henricus"
Atraits.meanH$Species.taxonstand <- sub("-",".",Atraits.meanH$Species.taxonstand)
traits.meanH$Species.taxonstand <- sub("-",".",traits.meanH$Species.taxonstand)
row.names(TAtraits.mean)<-sub("-",".",row.names(TAtraits.mean))						# fix "Chenopodium_bonus-henricus"
TAtraits.meanH$Species.taxonstand<-sub("-",".",TAtraits.meanH$Species.taxonstand)	# fix "Chenopodium_bonus-henricus"


# only herbaceous species in comunity
spxplotH <- spxplot[,traits.meanH$Species.taxonstand] #only herbaceous species (no mosses)
dim(spxplotH)  #72 herb sps

# Total cover in community
Cover <- rowSums(spxplot, na.rm=T) #all species
CoverH <- rowSums(spxplotH, na.rm=T) #only herbaceous species (no mosses)
Cover.noA <- rowSums(spxplot[,setdiff(colnames(spxplot),"Arabis_alpina")], na.rm=T) #Without Arabis
CoverH.noA <- rowSums(spxplot[,setdiff(traits.meanH$Species.taxonstand,"Arabis_alpina")], na.rm=T) #Without Arabis

# Crop matrix for species with trait data (arabis + TRY) / only herbs
spxplotT <- spxplot[,rownames(TAtraits.mean)]
head(spxplotT) ; dim(spxplotT) 	# 65 species x 18 plots

# Cover of species with traits (for height)
CoverT <- rowSums(spxplotT, na.rm=T)
CoverT/Cover  	# low representation in BRU site (but always above 50%)
CoverT/CoverH  	# but better if not counting trees and mosses (65 sps with traits out of 72 herb species)

# Cover of species with SLA
CoverTsla <- rowSums(spxplotT[,rownames(TAtraits.mean[!is.na(TAtraits.mean$SLA2),])], na.rm=T)
CoverTsla/Cover  	# very low representation in BRU site (20% cover), but ok in other sites
CoverTsla/CoverH  	# 37% in BRU site cover if not counting trees and mosses (60 sps with SLA out of 72 herb species)

# Cover of species with LNC
CoverTlnc <- rowSums(spxplotT[,rownames(TAtraits.mean[!is.na(TAtraits.mean$LNC2),])], na.rm=T)
CoverTlnc/Cover  	# very low representation in BRU site (12% cover), but ok in other sites
CoverTlnc/CoverH  	# 21% in BRU site cover if not counting trees and mosses (51 sps with LNC out of 72 herb species)

#Relative abundance matrix
spxplotR <- spxplot/Cover  #Full matrix
rowSums(spxplotR, na.rm=T) 
spxplotHR <- spxplotH/CoverH  #Herbaceous sps matrix
rowSums(spxplotHR, na.rm=T)
spxplotRT <- spxplotR[,row.names(TAtraits.mean)]  #Trait sub-matrix
rowSums(spxplotRT, na.rm=T) 
spxplotHRT <- spxplotHR[,row.names(TAtraits.mean)]  #Trait sub-matrix
rowSums(spxplotHRT, na.rm=T) 


#####	Calculate metrics -----------------------------

# Calculate CWMs

#with Arabis
CWM <- functcomp(as.matrix(TAtraits.mean), as.matrix(spxplotHRT))
head(CWM)
#without Arabis
CWM.noA <- functcomp(TAtraits.mean[setdiff(rownames(TAtraits.mean),"Arabis_alpina"),], 
					as.matrix(spxplotHRT)[,setdiff(rownames(TAtraits.mean),"Arabis_alpina")])
head(CWM.noA)

# Calculate Sps richness

Richness <- rowSums(spxplot>0, na.rm=T)
RichnessH <- rowSums(spxplotH>0, na.rm=T) #(no trees,no mosses)
Richness.noA <- rowSums(spxplot[,setdiff(colnames(spxplot),"Arabis_alpina")]>0, na.rm=T)
RichnessH.noA <- rowSums(spxplotH[,setdiff(colnames(spxplotH),"Arabis_alpina")]>0, na.rm=T)

# Functional diversity metrics interesting? Species average dissimilarities interesting?
# FD<-dbFD(DAtraits.mean, as.matrix(spxplotRT),calc.CWM = F)            

#####	Put metrics together  -----------------------------

#with Arabis
dat<-cbind(Arabis_alpina = spxplot$Arabis_alpina, TotCover=Cover, HerbCover=CoverH, TotRichness= Richness, HerbRichness= RichnessH,CWM)
head(dat)

#without Arabis
dat.noA<-cbind(Arabis_alpina = spxplot$Arabis_alpina, TotCover=Cover.noA, HerbCover=CoverH.noA, TotRichness= Richness.noA, HerbRichness= RichnessH.noA, CWM.noA)
head(dat.noA)

##################################################################
#	SAVE DATA FOR ANALYSIS 
##################################################################


save(list=c("androsace3", "Traits_Try.mean", "TAtraits.mean",						#trait data
			"ComSpList", "spxplot", "spxplotH", "spxplotHR", "spxplotT",		#all species and plot data
			"dat.noA", "dat"),                     								#community metrics
             file="Arabis.Community.Data.new.RData")

# "androsace3" 		: Original database of traits from alpine species from Wilfried
# "Traits_Try.mean"
# "TAtraits.mean" 	: Species x traits matrix (species with trait data available in communities)
# "ComSpList"  		: Species list
# "spxplot"  		: Sp x plot abundance matrix (braun blanquet transformed into %)
# "spxplotR" 		: Sp x plot abundance matrix, relative abundances
# "spxplotT" 		: Sp x plot abundance matrix, matching traits
# "spxplotH" 		: Sp x plot abundance matrix, no trees and mosses (only vascular herbaceous sps)
# "dat.noA", "dat" 	: Community and CWM metrics x plot matrix, with and without Arabis in calculation 



load("Arabis.Community.Data.RData")









