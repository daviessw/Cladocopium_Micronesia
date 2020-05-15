#Cladocopium Discriminant function on binary data from Micronesia
#Author: Sarah W. Davies March 2020
##Note C3 wqas later desingated as C21 so please use C21 files
library(vegan)
library(plotrix)
library(adegenet)
library(ggplot2)
library(tidyr)
library(tidyverse)

setwd("~/Dropbox/Texas/CladeC/2019_MEC_Take2/SUBMISSION_March2019/Resubmission_2019/New-All-Figs")

dat <- read.table("Cladocopium_Binary_AllData.txt", sep="\t", header=TRUE, row.names=1)
dat$island=as.factor(dat$island)
dat$site=as.factor(dat$site)
dat$spp=as.factor(dat$spp)
head(dat)
#reorder levels of islands
dat$site <- ordered(dat$site, levels = c("West Channel", "Lighthouse Reef", "Ngulu Atoll", "Goofnuw Channel", "South Tip", "Pago Bay", "Tanguisson", "West Polle", "South East Pass", "Ant Atoll", "Roj", "Coral Garden", "Hiroshi Point"))
dat$island <- ordered(dat$island, levels = c("Palau", "Ngulu", "Yap", "Guam", "Chuuk", "Pohnpei", "Kosrae"))
dat$spp <- ordered(dat$spp, levels = c("A. digitifera", "A. hyacinthus"))

gendat=dat[,-c(1:3)]
gendat=gendat[,colSums(gendat) > 0]
head(gendat)
ncol(gendat) #56 alleles total
obj <- genind(gendat, ploidy = 1, type = "PA")
obj

grp=find.clusters(obj, max.n.clust=40)
#18 PCs and 2 clusters saw the strongest drop in BIC ay 2 clusters
names(grp)
head(grp$grp, 2)
grp$size #174 in group 1 and 394 in group 2
dapc1=dapc(obj, grp$grp) #18PC 1 DF

pdf('All-Data-Ghost.pdf', width = 4, height = 4)
scatter(dapc1, posi.da="topleft", bg="white", cex=1.5, legend=TRUE, posi.leg="topright", col=c("chartreuse2", "darkolivegreen"))
dev.off()

compoplot(dapc1, posi="bottomright", txt.leg=paste("Cluster", 1:2), lab="", xlab="individuals", col=c("chartreuse2", "darkolivegreen"))
pred <- predict(dapc1)
names(pred)
new=cbind.data.frame(pred$posterior, dat$spp, dat$island, dat$site)
head(new)
colnames(new)=c("assign1", "assign2", "spp", "island", "site")
head(new)

new2 <- new %>%
  rownames_to_column(var = 'indiv') %>%
  gather(assignment, proportion, assign1:assign2)

head(new2)
ggplot(new2, aes(x = indiv, y = proportion, fill = assignment)) +
  geom_bar(stat="identity", width = 1) +
  facet_wrap( ~ island, scales = 'free_x', nrow = 1)+
  scale_x_discrete(NULL) +
  scale_y_continuous("Proportion",
                     expand = c(0,0)) +
  scale_fill_manual(labels = c("C3", "C40"),
                    values=c("chartreuse2", "darkolivegreen")) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )


########   Create df of only C3 or C40, (>.9 prob) from original df
ass1 <- new %>%
  rownames_to_column(var = "sample") %>%
  filter(assign1 > .9) %>%
  pull(sample) %>%
  as.vector()
ass1
#172 samples

ass2 <- new %>%
  rownames_to_column(var = "sample") %>%
  filter(assign2 > .9) %>%
  pull(sample) %>%
  as.vector()
ass2
#388 samples

dat %>%
  rownames_to_column(var = "sample") %>%
  select(sample, everything()) %>%
  filter(sample %in% ass1) %>%
  write_csv("C3_alldata.csv")

dat %>%
  rownames_to_column(var = "sample") %>%
  select(sample, everything()) %>%
  filter(sample %in% ass2) %>%
  write_csv("C40_alldata.csv")

##############CLONAL ANALYSIS
#Written by Kelsey Moreland
# ADD IT IN HERE

#####now looking only within C40 individuals with only unique samples from each site, species, and group combination retained 
dat_c40 <- read.table("C40_noclones_withinsppsite.csv", sep=",", header=TRUE)
head(dat_c40)
nrow(dat_c40)
#127 samples left / 172 (45 samples removed)
summary(dat_c40)
dat_c40$island=as.factor(dat_c40$island)
dat_c40$site=as.factor(dat_c40$site)
dat_c40$spp=as.factor(dat_c40$spp)
dat_c40$site
dat_c40$combo=paste(dat_c40$spp, dat_c40$site)
dat_c40$combo=as.factor(dat_c40$combo)
length(which(dat_c40$spp == "A. digitifera"))
# number after clones removed/number assigned to C40
# A. digitifera West Channel: 16/25 
# A. hyacinthus West Channel: 13/23
# A. digitifera Lighthouse Reef: 19/24
# A. hyacinthus Lighthouse Reef: 18/25
#A. hyacinthus Ngulu Atoll: 28/42
#A. digitifera South Tip: 1/1
#A. hyacinthus South Tip: 0/0
#A. digitifera Goofnuw Channel: 17/17
#A. hyacinthus Goofnuw Channel: 0/0
#A. digitifera Pago Bay: 0/0
#A. digitifera Tanguisson: 0/0
#A. digitifera West Polle: 2/2
#A. hyacinthus West Polle: 1/1
# A. digitifera South East Pass: 1/1
# A. hyacinthus South East Pass: 2/2
# A. digitifera Ant Atoll: 0/0
# A. hyacinthus Ant Atoll: 3/3
# A. digitifera Roj: 0/0
# A. hyacinthus Roj: 1/1
# A. digitifera Coral Garden: 1/1
# A. hyacinthus Coral Garden: 2/2
# A. digitifera Hiroshi Point: 2/2
# A. hyacinthus Hiroshi Point: 0/0
#127/172 total, 68/99 A hya and 59/73 A digi

#looking at all data by spp
head(dat_c40)
colnames(dat_c40)
genall=dat_c40[,-c(1:4, 61)]
ncol(genall) #56 alleles total to start
genall=genall[,colSums(genall) > 0]
head(genall)
ncol(genall)
#44 alleles total in C40
nrow(genall)
#127 indiv assigned to C40
#allelic diversity in A. hya
hya =subset(dat_c40, spp=="A. hyacinthus")
genhya=hya[,-c(1:4, 61)]
genhya=genhya[,colSums(genhya) > 0]
ncol(genhya)
nrow(genhya)
#44 alleles total in C40 hosted by 68 A. hyacinthus indivs

digi =subset(dat_c40, spp=="A. digitifera")
gendigi=digi[,-c(1:4, 61)]
gendigi=gendigi[,colSums(gendigi) > 0]
ncol(gendigi)
nrow(gendigi)
#29 alleles total in C40 hosted by 59 A. digitifera indivs

obj <- genind(genall, ploidy = 1, type = "PA")
obj
truenames(obj)
pop(obj)=dat_c40$spp
obj2 <- genind2genpop(obj)
library(ade4)
pca1 <- dudi.pca(obj2, scannf=FALSE, scale=FALSE)
dapcall <- dapc(obj, n.da=50, n.pca=50)
temp <- optim.a.score(dapcall)
#tells us that for all dataset we keep only 16 PCs
dapcspp <- dapc(obj, dat_c40$spp, n.pca=16, n.da=1)
dapcspp
summary(dapcspp)
# $var (proportion of conserved variance): 0.914
# $assign.prop
# [1] 0.8110236
# $assign.per.pop
# A. digitifera A. hyacinthus 
    # 0.8135593     0.8088235  
mycol=c("purple4", "orchid2")
pdf('c40_hostspecies.pdf', width = 4, height = 4)
scatter(dapcspp, posi.da="topright", bg="white", col=mycol, cex=1.5, legend=TRUE, posi.leg="topright")
dev.off()

#looking at C40 by island
head(genall)
obj <- genind(genall, ploidy = 1, type = "PA")
obj
truenames(obj)
pop(obj)=dat_c40$island
obj2 <- genind2genpop(obj)
pca1 <- dudi.pca(obj2, scannf=FALSE, scale=FALSE)
dapcall <- dapc(obj, n.da=50, n.pca=50)
temp <- optim.a.score(dapcall)
#tells us that for all dataset we keep only 17 PCs
dapcis <- dapc(obj, dat_c40$island, n.pca=17, n.da=5)
dapcis
summary(dapcis)
# $var (proportion of conserved variance): 0.924
# $assign.prop
# [1] 0.8503937
# $assign.per.pop
    # Chuuk    Kosrae     Ngulu     Palau   Pohnpei       Yap 
# 0.8333333 0.8000000 0.8928571 0.9090909 0.7500000 0.6111111  

# Island colors: symbols
# Palau: red2, 8
# Ngulu: tomato1, 15
# Guam: aquamarine3, 16
# Yap: goldenrod1, 17
# Chuuk: steelblue2, 18
# Pohnpei: royalblue4, 3
# Kosrae: brown4, 4

mypch=c(18, 4, 15, 8, 3, 17)
myCol <- c("steelblue2","brown4","tomato1","red2","royalblue4","goldenrod1") 

pdf('c40_island.pdf', width = 4, height = 4)
scatter(dapcis, posi.da="topright", bg="white", col=myCol, pch=mypch, cex=1.5, legend=TRUE, posi.leg="bottomright")
dev.off()


#diversity of C40 across islands
pal =subset(dat_c40, island=="Palau")
genpal=pal[,-c(1:4, 61)]
genpal=genpal[,colSums(genpal) > 0]
ncol(genpal)
nrow(genpal)
#30 alleles total in C40 on Palau in 66 indivs

ngu =subset(dat_c40, island=="Ngulu")
genngu=ngu[,-c(1:4, 61)]
genngu=genngu[,colSums(genngu) > 0]
ncol(genngu)
nrow(genngu)
#23 alleles total in C40 on Ngulu in 28 indivs

yap =subset(dat_c40, island=="Yap")
genyap=yap[,-c(1:4, 61)]
genyap=genyap[,colSums(genyap) > 0]
ncol(genyap)
nrow(genyap)
#16 alleles total in C40 on Yap in 18 indivs

gua =subset(dat_c40, island=="Guam")
gengua=gua[,-c(1:4, 61)]
gengua=gengua[,colSums(gengua) > 0]
ncol(gengua)
nrow(gengua)
#0 alleles total in C40 on Guam in 0 indivs

chu =subset(dat_c40, island=="Chuuk")
genchu=chu[,-c(1:4, 61)]
genchu=genchu[,colSums(genchu) > 0]
ncol(genchu)
nrow(genchu)
#19 alleles total in C40 on Chuuk in 6 indivs

poh =subset(dat_c40, island=="Pohnpei")
genpoh=poh[,-c(1:4, 61)]
genpoh=genpoh[,colSums(genpoh) > 0]
ncol(genpoh)
nrow(genpoh)
#15 alleles total in C40 on Pohnpei in 4 indivs

kos =subset(dat_c40, island=="Kosrae")
genkos=kos[,-c(1:4, 61)]
genkos=genkos[,colSums(genkos) > 0]
ncol(genkos)
nrow(genkos)
#20 alleles total in C40 on Kosrae in 5 indivs

#looking C40 data within Palau with clonal individuals removed from  
head(dat_c40)
pal <- subset(dat_c40, island=="Palau")
ncol(pal)
head(pal)
genall=pal[,-c(1:4, 61)]
ncol(genall) #56 total in all data
genall=genall[,colSums(genall) > 0]
head(genall)
ncol(genall)
#30 alleles total in C40 on Palau
nrow(genall)
#66/97 indiv assigned to C40 on Palau with clones removed

obj <- genind(genall, ploidy = 1, type = "PA")
obj
truenames(obj)
pal_sites <- as.vector(pal$combo)
pop(obj)=pal_sites
obj
obj2 <- genind2genpop(obj)
obj2
pca1 <- dudi.pca(obj2, scannf=FALSE, scale=FALSE)
dapcall <- dapc(obj, n.da=50, n.pca=50)
temp <- optim.a.score(dapcall)
#tells us that for all dataset we keep only 14 PCs
dapcpal <- dapc(obj, as.vector(pal$combo), n.pca=14, n.da=3)
dapcpal
summary(dapcpal)
# $var (proportion of conserved variance): 0.923
# $assign.prop
# [1] 0.6969697
# $assign.per.pop
# A. digitifera Lighthouse Reef    A. digitifera West Channel A. hyacinthus Lighthouse Reef 
                    # 0.8421053                     0.7500000                     0.6111111 
   # A. hyacinthus West Channel 
                    # 0.5384615 
pdf('c40_Palau.pdf', width = 5, height = 5)
scatter(dapcpal, posi.da="bottomleft", bg="white", , pch=c(17, 15, 17, 15), cex=1.5, legend=TRUE, posi.leg="topleft", col=c("purple4",   "purple4", "orchid2", "orchid2"))
dev.off()


##diversity within Palau
pal =subset(dat_c40, island=="Palau")
hyapal=subset(pal, spp=="A. hyacinthus")
hyapalLH=subset(hyapal, site=="Lighthouse Reef")
genhyapalLH=hyapalLH[,-c(1:4, 61)]
genhyapalLH=genhyapalLH[,colSums(genhyapalLH) > 0]
ncol(genhyapalLH)
nrow(genhyapalLH)
#15/30 alleles in 18 indivs

hyapalWC=subset(hyapal, site=="West Channel")
genhyapalWC=hyapalWC[,-c(1:4, 61)]
genhyapalWC=genhyapalWC[,colSums(genhyapalWC) > 0]
ncol(genhyapalWC)
nrow(genhyapalWC)
#21/30 alleles in 13 indivs

digipal=subset(pal, spp=="A. digitifera")
digipalLH=subset(digipal, site=="Lighthouse Reef")
gendigipalLH=digipalLH[,-c(1:4, 61)]
gendigipalLH=gendigipalLH[,colSums(gendigipalLH) > 0]
ncol(gendigipalLH)
nrow(gendigipalLH)
#19/30 alleles in 19 indivs

digipalWC=subset(digipal, site=="West Channel")
gendigipalWC=digipalWC[,-c(1:4, 61)]
gendigipalWC=gendigipalWC[,colSums(gendigipalWC) > 0]
ncol(gendigipalWC)
nrow(gendigipalWC)
#14/30 alleles in 16 indivs

####################################################################
#################now looking only within C3 individuals with only unique samples from each site, species, and group combination retained 
dat_c3 <- read.table("C3_noclones_withinsppsite.csv", sep=",", header=TRUE)
head(dat_c3)
nrow(dat_c3)
#328 samples/ original 388 (lost 60 due to clones)
summary(dat_c3)
dat_c3$island=as.factor(dat_c3$island)
dat_c3$site=as.factor(dat_c3$site)
dat_c3$spp=as.factor(dat_c3$spp)
dat_c3$site
dat_c3$combo=paste(dat_c3$spp, dat_c3$site)
dat_c3$combo=as.factor(dat_c3$combo)
length(which(dat_c3$combo == "A. hyacinthus Hiroshi Point"))
# number after clones removed/number assigned to C40
# A. digitifera West Channel: 0/0
# A. hyacinthus West Channel: 1/1
# A. digitifera Lighthouse Reef: 0/0
# A. hyacinthus Lighthouse Reef: 0/0
#A. hyacinthus Ngulu Atoll: 0/0
#A. digitifera South Tip: 23/23
#A. hyacinthus South Tip: 20/25
#A. digitifera Goofnuw Channel: 5/7
#A. hyacinthus Goofnuw Channel: 24/25
#A. digitifera Pago Bay: 20/26
#A. digitifera Tanguisson: 17/20
#A. digitifera West Polle: 11/13
#A. hyacinthus West Polle: 22/22
# A. digitifera South East Pass: 19/20
# A. hyacinthus South East Pass: 13/21
# A. digitifera Ant Atoll: 17/24
# A. hyacinthus Ant Atoll: 19/20
# A. digitifera Roj: 21/24
# A. hyacinthus Roj: 21/23
# A. digitifera Coral Garden: 19/23
# A. hyacinthus Coral Garden: 20/23
# A. digitifera Hiroshi Point: 14/23
# A. hyacinthus Hiroshi Point: 22/25
#328/388 total, 162/185 A hya and 166/203 A digi

#looking at all data by spp
head(dat_c3)
colnames(dat_c3)
genall=dat_c3[,-c(1:4, 61)]
ncol(genall) #56 total
genall=genall[,colSums(genall) > 0]
head(genall)
ncol(genall)
#49 alleles total in C3
nrow(genall)
#328 indiv assigned to C3

obj <- genind(genall, ploidy = 1, type = "PA")
obj
truenames(obj)
pop(obj)=dat_c3$spp
obj2 <- genind2genpop(obj)
pca1 <- dudi.pca(obj2, scannf=FALSE, scale=FALSE)
dapcall <- dapc(obj, n.da=50, n.pca=50)
temp <- optim.a.score(dapcall)
#tells us that for all dataset we keep only 17 PCs
dapcspp <- dapc(obj, dat_c3$spp, n.pca=17, n.da=1)
dapcspp
summary(dapcspp)
# $var (proportion of conserved variance): 0.915
# $assign.prop
# [1] 0.8384146
# $assign.per.pop
# A. digitifera A. hyacinthus 
    # 0.9277108     0.7469136 
mycol=c("purple4", "orchid2")

pdf('c3_hostspecies.pdf', width = 4, height = 4)
scatter(dapcspp, posi.da="bottomleft", bg="white", col=mycol, cex=1.5, legend=TRUE, posi.leg="topright")
dev.off()

#looking at C3 by island
head(genall)
obj <- genind(genall, ploidy = 1, type = "PA")
obj
truenames(obj)
pop(obj)=dat_c3$island
obj2 <- genind2genpop(obj)
pca1 <- dudi.pca(obj2, scannf=FALSE, scale=FALSE)
dapcall <- dapc(obj, n.da=50, n.pca=50)
temp <- optim.a.score(dapcall)
#tells us that for all dataset we keep only 11 PCs
dapcis <- dapc(obj, dat_c3$island, n.pca=11, n.da=5)
dapcis
summary(dapcis)
# $var (proportion of conserved variance): 0.804
# $assign.prop
# [1] 0.7042683
# $assign.per.pop
    # Chuuk      Guam    Kosrae     Palau   Pohnpei       Yap 
# 0.5384615 0.7567568 0.8800000 1.0000000 0.7179487 0.6250000 

# Island colors: symbols
# Palau: red2, 8
# Ngulu: tomato1, 15
# Guam: aquamarine3, 16
# Yap: goldenrod1, 17
# Chuuk: steelblue2, 18
# Pohnpei: royalblue4, 3
# Kosrae: brown4, 4

mypch=c(18, 16, 4, 8, 3, 17)
myCol <- c("steelblue2","aquamarine3", "brown4","red2","royalblue4","goldenrod1") 

pdf('c3_island.pdf', width = 4, height = 4)
scatter(dapcis, posi.da="bottomright", bg="white", col=myCol, pch=mypch, cex=1.5, legend=TRUE, posi.leg="topright")
dev.off()

#looking C3 data within Yap
head(dat_c3)
yap <- subset(dat_c3, island=="Yap")
genall=yap[,-c(1:4, 61)]
genall=genall[,colSums(genall) > 0]
head(genall)
ncol(genall)
#27 alleles total in C3 on Yap
nrow(genall)
#72/80 indiv assigned to C3 on Yap

obj <- genind(genall, ploidy = 1, type = "PA")
obj
truenames(obj)
yap_sites <- as.vector(yap$combo)
pop(obj)=yap_sites
obj
obj2 <- genind2genpop(obj)
obj2
pca1 <- dudi.pca(obj2, scannf=FALSE, scale=FALSE)
dapcall <- dapc(obj, n.da=50, n.pca=50)
temp <- optim.a.score(dapcall)
#tells us that for all dataset we keep only 16 PCs
dapcyap <- dapc(obj, as.vector(yap$combo), n.pca=16, n.da=3)
dapcyap
summary(dapcyap)
# $var (proportion of conserved variance): 0.961
# $assign.prop
# [1] 0.8333333

# $assign.per.pop
# A. digitifera Goofnuw Channel       A. digitifera South Tip A. hyacinthus Goofnuw Channel 
                    # 0.6000000                     0.9130435                     0.9166667 
      # A. hyacinthus South Tip 
                    # 0.7000000 
pdf('c3_Yap.pdf', width = 5, height = 5)
scatter(dapcyap, posi.da="bottomleft", bg="white", , pch=c(17, 15, 17, 15), cex=1.5, legend=TRUE, posi.leg="topleft", col=c("purple4",   "purple4", "orchid2", "orchid2"))
dev.off()

#looking C3 data within Guam
head(dat_c3)
guam <- subset(dat_c3, island=="Guam")
genall=guam[,-c(1:4, 61)]
genall=genall[,colSums(genall) > 0]
head(genall)
ncol(genall)
#25 alleles total in C3 on Guam
nrow(genall)
#37/46 indiv assigned to C3 on Guam

obj <- genind(genall, ploidy = 1, type = "PA")
obj
truenames(obj)
guam_sites <- as.vector(guam$combo)
pop(obj)=guam_sites
obj
obj2 <- genind2genpop(obj)
obj2
pca1 <- dudi.pca(obj2, scannf=FALSE, scale=FALSE)
dapcall <- dapc(obj, n.da=50, n.pca=50)
temp <- optim.a.score(dapcall)
#tells us that for all dataset we keep only 3 PCs
dapcguam <- dapc(obj, as.vector(guam$combo), n.pca=3, n.da=3)
dapcguam
summary(dapcguam)
# $var (proportion of conserved variance): 0.522
# $assign.prop
# [1] 0.8108108
# $assign.per.pop
  # A. digitifera Pago Bay A. digitifera Tanguisson 
               # 0.9000000                0.7058824 
pdf('c3_Guam.pdf', width = 5, height = 5)
scatter(dapcguam, posi.da="bottomleft", bg="white", cex=1.5, legend=TRUE, posi.leg="topright", col=c("purple4",   "mediumpurple1"))
dev.off()

#looking C3 data within Chuuk
head(dat_c3)
chuuk <- subset(dat_c3, island=="Chuuk")
genall=chuuk[,-c(1:4, 61)]
genall=genall[,colSums(genall) > 0]
head(genall)
ncol(genall)
#35 alleles total in C3 on Chuuk
nrow(genall)
#65/76 indiv assigned to C3 on Chuuk

obj <- genind(genall, ploidy = 1, type = "PA")
obj
truenames(obj)
chuuk_sites <- as.vector(chuuk$combo)
pop(obj)=chuuk_sites
obj
obj2 <- genind2genpop(obj)
obj2
pca1 <- dudi.pca(obj2, scannf=FALSE, scale=FALSE)
dapcall <- dapc(obj, n.da=50, n.pca=50)
temp <- optim.a.score(dapcall)
#tells us that for all dataset we keep only 10 PCs
dapcchuuk <- dapc(obj, as.vector(chuuk$combo), n.pca=10, n.da=3)
dapcchuuk
summary(dapcchuuk)
# $var (proportion of conserved variance): 0.823
# $assign.prop
# [1] 0.8153846
# # $assign.per.pop
# A. digitifera South East Pass      A. digitifera West Polle A. hyacinthus South East Pass 
                    # 0.7894737                     0.7272727                     1.0000000 
     # A. hyacinthus West Polle 
                    # 0.7727273 
pdf('c3_Chuuk.pdf', width = 5, height = 5)
scatter(dapcchuuk, posi.da="bottomleft", bg="white", , pch=c(17, 15, 17, 15), cex=1.5, legend=TRUE, posi.leg="topleft", col=c("purple4",   "purple4", "orchid2", "orchid2"))
dev.off()

#looking C3 data within Pohnpei
head(dat_c3)
poh <- subset(dat_c3, island=="Pohnpei")
genall=poh[,-c(1:4, 61)]
genall=genall[,colSums(genall) > 0]
head(genall)
ncol(genall)
#31 alleles total in C3 on Pohnpei
nrow(genall)
#78/91 indiv assigned to C3 on Pohnpei

obj <- genind(genall, ploidy = 1, type = "PA")
obj
truenames(obj)
poh_sites <- as.vector(poh$combo)
pop(obj)=poh_sites
obj
obj2 <- genind2genpop(obj)
obj2
pca1 <- dudi.pca(obj2, scannf=FALSE, scale=FALSE)
dapcall <- dapc(obj, n.da=50, n.pca=50)
temp <- optim.a.score(dapcall)
#tells us that for all dataset we keep only 9 PCs
dapcpoh <- dapc(obj, as.vector(poh$combo), n.pca=9, n.da=3)
dapcpoh
summary(dapcpoh)
# $var (proportion of conserved variance): 0.814
# $assign.prop
# [1] 0.8461538
# $assign.per.pop
# A. digitifera Ant Atoll       A. digitifera Roj A. hyacinthus Ant Atoll       A. hyacinthus Roj 
              # 0.9411765               0.9047619               0.7894737               0.7619048 
pdf('c3_Pohnpei.pdf', width = 5, height = 5)
scatter(dapcpoh, posi.da="topright", bg="white", , pch=c(17, 15, 17, 15), cex=1.5, legend=TRUE, posi.leg="topleft", col=c("purple4",   "purple4", "orchid2", "orchid2"))
dev.off()

#looking C3 data within Kosrae
head(dat_c3)
kos <- subset(dat_c3, island=="Kosrae")
genall=kos[,-c(1:4, 61)]
genall=genall[,colSums(genall) > 0]
head(genall)
ncol(genall)
#22 alleles total in C3 on Kosrae
nrow(genall)
#75/94 indiv assigned to C3 on Kosrae

obj <- genind(genall, ploidy = 1, type = "PA")
obj
truenames(obj)
kos_sites <- as.vector(kos$combo)
pop(obj)=kos_sites
obj
obj2 <- genind2genpop(obj)
obj2
pca1 <- dudi.pca(obj2, scannf=FALSE, scale=FALSE)
dapcall <- dapc(obj, n.da=50, n.pca=50)
temp <- optim.a.score(dapcall)
#tells us that for all dataset we keep only 7 PCs
dapckos <- dapc(obj, as.vector(kos$combo), n.pca=7, n.da=3)
dapckos
summary(dapckos)
# $var (proportion of conserved variance): 0.792
# $assign.prop
# [1] 0.84
# $assign.per.pop
 # A. digitifera Coral Garden A. digitifera Hiroshi Point  A. hyacinthus Coral Garden 
                  # 0.8947368                   0.5714286                   1.0000000 
# A. hyacinthus Hiroshi Point 
                  # 0.8181818 

pdf('c3_Kosrae.pdf', width = 5, height = 5)
scatter(dapckos, posi.da="bottomright", bg="white", , pch=c(17, 15, 17, 15), cex=1.5, legend=TRUE, posi.leg="topleft", col=c("purple4",   "purple4", "orchid2", "orchid2"))
dev.off()

#################################################################
#################################################################
#####PCA analysis
library(vegan)
library(ggplot2)
library(ggrepel)
library(tidyverse)
library(ade4)
library(ape)

####First on all C40 with clones removed
dat <- read.table("C40_noclones_withinsppsite.csv", sep=",", header=TRUE, row.names=1)
head(dat)
ncol(dat)
###all data visualizaed

# PCoA on a percentage difference (Bray-Curtis) dissimilarity matrix of species
dist <- vegdist(dat[,4:59], method="bray")
pc=pcoa(dist)

pc1=as.data.frame(pc$vectors[,1:2])
head(pc1)
pc1$island=dat$island
pc1$site=dat$site
pc1$spp=dat$spp
head(pc1)

# Island colors: symbols
# Palau: red2, 8
# Ngulu: tomato1, 15
# Guam: aquamarine3, 16
# Yap: goldenrod1, 17
# Chuuk: steelblue2, 18
# Pohnpei: royalblue4, 3
# Kosrae: brown4, 4
pc1$island

pdf('c40_Islands_PCA.pdf', width = 5.5, height = 4)
cbPalette <- c("steelblue2",  "brown4", "tomato1", "red2", "royalblue4", "goldenrod1")
ggplot(pc1, aes(Axis.1, Axis.2, color = island, pch = island)) +
  geom_point(size=3) +
  # geom_text_repel(aes(label=Samples)) +
  scale_colour_manual(values=cbPalette)+
  theme_bw() +
  stat_ellipse()+
  xlab("PC1") +
  ylab("PC2") +
  xlim(-0.4, 0.52)+
  ylim(-0.55, 0.46)
dev.off()

adonis(dist ~ island, data = dat, method='bray')
           # Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# island      5    5.9520 1.19040  16.779 0.40945  0.001 ***
# Residuals 121    8.5845 0.07095         0.59055           
# Total     126   14.5365                 1.00000           

mycol=c("purple4", "orchid2")
pdf('c40_spp_PCA.pdf', width = 5.5, height = 4)
ggplot(pc1, aes(Axis.1, Axis.2, color = spp, pch = spp)) +
  geom_point(size=3) +
  # geom_text_repel(aes(label=Samples)) +
  scale_colour_manual(values=mycol)+
  theme_bw() +
  stat_ellipse()+
  xlab("PC1") +
  ylab("PC2") 
dev.off()

adonis(dist ~ spp, data = dat, method='bray')
           # Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# spp         1    1.0803 1.08031  10.035 0.07432  0.001 ***
# Residuals 125   13.4562 0.10765         0.92568           
# Total     126   14.5365                 1.00000 

###now just Palau
head(dat)
ncol(dat)
pal=subset(dat, island=="Palau")
pal$combo=paste(pal$site,pal$spp)
head(pal)
dist <- vegdist(pal[,4:59], method="bray")
pc=pcoa(dist)

pc1=as.data.frame(pc$vectors[,1:2])
head(pc1)
pc1$site=pal$site
pc1$spp=pal$spp
pc1$combo=pal$combo
head(pc1)

mycol=c("purple4","orchid2",  "purple4","orchid2")
mypch=c(17,17,15,15)
pdf('c40_Palau_PCA.pdf', width = 6, height = 4)
ggplot(pc1, aes(Axis.1, Axis.2, color = combo, pch = combo)) +
  geom_point(size=3) +
  # geom_text_repel(aes(label=Samples)) +
  scale_colour_manual(values=mycol)+
  scale_shape_manual(values=mypch)+
theme_bw() +
  stat_ellipse()+
  xlab("PC1") +
  ylab("PC2") 
dev.off()
     
adonis(dist ~ site+spp, data = pal, method='bray')
          # Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)    
# site       1    0.2860 0.286043  5.6456 0.07611  0.001 ***
# spp        1    0.2800 0.280037  5.5271 0.07452  0.002 ** 
# Residuals 63    3.1920 0.050667         0.84937           
# Total     65    3.7581                  1.00000 

####First on all C3 with clones removed
dat <- read.table("C3_noclones_withinsppsite.csv", sep=",", header=TRUE, row.names=1)
head(dat)
ncol(dat)

# PCoA on a percentage difference (Bray-Curtis) dissimilarity matrix of species
dist <- vegdist(dat[,4:59], method="bray")
pc=pcoa(dist)

pc1=as.data.frame(pc$vectors[,1:2])
head(pc1)
pc1$island=dat$island
pc1$site=dat$site
pc1$spp=dat$spp
head(pc1)

# Island colors: symbols
# Palau: red2, 8
# Ngulu: tomato1, 15
# Guam: aquamarine3, 16
# Yap: goldenrod1, 17
# Chuuk: steelblue2, 18
# Pohnpei: royalblue4, 3
# Kosrae: brown4, 4
pc1$island

pdf('c3_Islands_PCA.pdf', width = 5.5, height = 4)
cbPalette <- c("steelblue2", "aquamarine3", "brown4", "red2", "royalblue4", "goldenrod1")
ggplot(pc1, aes(Axis.1, Axis.2, color = island, pch = island)) +
  geom_point(size=3) +
  # geom_text_repel(aes(label=Samples)) +
  scale_colour_manual(values=cbPalette)+
  theme_bw() +
  stat_ellipse()+
  xlab("PC1") +
  ylab("PC2") 
dev.off()

# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# island      5    6.1455 1.22910  18.324 0.22151  0.001 ***
# Residuals 322   21.5985 0.06708         0.77849           
# Total     327   27.7441                 1.00000

mycol=c("purple4", "orchid2")
pdf('c3_spp_PCA.pdf', width = 5.5, height = 4)
ggplot(pc1, aes(Axis.1, Axis.2, color = spp, pch = spp)) +
  geom_point(size=3) +
  # geom_text_repel(aes(label=Samples)) +
  scale_colour_manual(values=mycol)+
  theme_bw() +
  stat_ellipse()+
  xlab("PC1") +
  ylab("PC2") 
dev.off()

adonis(dist ~ spp, data = dat, method='bray')
           # Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# spp         1    1.4826 1.48263  18.405 0.05344  0.001 ***
# Residuals 326   26.2614 0.08056         0.94656           
# Total     327   27.7441                 1.00000           

###now just Yap
head(dat)
ncol(dat)
yap=subset(dat, island=="Yap")
yap$combo=paste(yap$site,yap$spp)
head(yap)
dist <- vegdist(yap[,4:59], method="bray")
pc=pcoa(dist)

pc1=as.data.frame(pc$vectors[,1:2])
head(pc1)
pc1$site=yap$site
pc1$spp=yap$spp
pc1$combo=yap$combo
head(pc1)

mycol=c("purple4","orchid2",  "purple4","orchid2")
mypch=c(17,17,15,15)
pdf('c3_Yap_PCA.pdf', width = 6, height = 4)
ggplot(pc1, aes(Axis.1, Axis.2, color = combo, pch = combo)) +
  geom_point(size=3) +
  # geom_text_repel(aes(label=Samples)) +
  scale_colour_manual(values=mycol)+
  scale_shape_manual(values=mypch)+
theme_bw() +
  stat_ellipse()+
  xlab("PC1") +
  ylab("PC2") 
dev.off()

adonis(dist ~ spp+site, data = yap, method='bray')
          # Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# spp        1    0.4019 0.40188  8.3888 0.10452  0.001 ***
# site       1    0.1375 0.13745  2.8691 0.03575  0.016 *  
# Residuals 69    3.3056 0.04791         0.85973           
# Total     71    3.8449                 1.00000 

###now just Guam
head(dat)
ncol(dat)
gua=subset(dat, island=="Guam")
gua$combo=paste(gua$site,gua$spp)
head(gua)
dist <- vegdist(gua[,4:59], method="bray")
pc=pcoa(dist)

pc1=as.data.frame(pc$vectors[,1:2])
head(pc1)
pc1$site=gua$site
head(pc1)

mycol=c("purple4","purple4")
mypch=c(17,15)
pdf('c3_Guam_PCA.pdf', width = 4.75, height = 4)
ggplot(pc1, aes(Axis.1, Axis.2, color = site, pch = site)) +
  geom_point(size=3) +
  # geom_text_repel(aes(label=Samples)) +
  scale_colour_manual(values=mycol)+
  scale_shape_manual(values=mypch)+
theme_bw() +
  stat_ellipse()+
  xlab("PC1") +
  ylab("PC2") 
dev.off()

adonis(dist ~ site, data = gua, method='bray')
          # Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)  
# site       1   0.17631 0.176308  3.2511 0.08499  0.018 *
# Residuals 35   1.89808 0.054231         0.91501         
# Total     36   2.07439                  1.00000  

###now just Chuuk
head(dat)
ncol(dat)
chu=subset(dat, island=="Chuuk")
chu$combo=paste(chu$site,chu$spp)
head(chu)
dist <- vegdist(chu[,4:59], method="bray")
pc=pcoa(dist)

pc1=as.data.frame(pc$vectors[,1:2])
head(pc1)
pc1$site=chu$site
pc1$spp=chu$spp
pc1$combo=chu$combo
head(pc1)

mycol=c("purple4","orchid2",  "purple4","orchid2")
mypch=c(17,17,15,15)
pdf('c3_Chuuk_PCA.pdf', width = 6, height = 4)
ggplot(pc1, aes(Axis.1, Axis.2, color = combo, pch = combo)) +
  geom_point(size=3) +
  # geom_text_repel(aes(label=Samples)) +
  scale_colour_manual(values=mycol)+
  scale_shape_manual(values=mypch)+
theme_bw() +
  stat_ellipse()+
  xlab("PC1") +
  ylab("PC2") 
dev.off()

adonis(dist ~ site+spp, data = chu, method='bray')
          # Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# site       1    0.6217 0.62170  8.3386 0.10437  0.001 ***
# spp        1    0.7126 0.71262  9.5581 0.11963  0.001 ***
# Residuals 62    4.6225 0.07456         0.77600           
# Total     64    5.9569                 1.00000  

###now just Pohnpei
head(dat)
ncol(dat)
poh=subset(dat, island=="Pohnpei")
poh$combo=paste(poh$site,poh$spp)
head(poh)
dist <- vegdist(poh[,4:59], method="bray")
pc=pcoa(dist)

pc1=as.data.frame(pc$vectors[,1:2])
head(pc1)
pc1$site=poh$site
pc1$spp=poh$spp
pc1$combo=poh$combo
head(pc1)

mycol=c("purple4","orchid2",  "purple4","orchid2")
mypch=c(17,17,15,15)
pdf('c3_Pohnpei_PCA.pdf', width = 5.5, height = 4)
ggplot(pc1, aes(Axis.1, Axis.2, color = combo, pch = combo)) +
  geom_point(size=3) +
  # geom_text_repel(aes(label=Samples)) +
  scale_colour_manual(values=mycol)+
  scale_shape_manual(values=mypch)+
theme_bw() +
  stat_ellipse()+
  xlab("PC1") +
  ylab("PC2") 
dev.off()

adonis(dist ~ spp+site, data = poh, method='bray')
          # Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# spp        1    0.4952 0.49519  11.334 0.10692  0.001 ***
# site       1    0.8593 0.85934  19.668 0.18555  0.001 ***
# Residuals 75    3.2769 0.04369         0.70753           
# Total     77    4.6314                 1.00000      

###now just Kosrae
head(dat)
ncol(dat)
kos=subset(dat, island=="Kosrae")
kos$combo=paste(kos$site,kos$spp)
head(kos)
dist <- vegdist(kos[,4:59], method="bray")
pc=pcoa(dist)

pc1=as.data.frame(pc$vectors[,1:2])
head(pc1)
pc1$site=kos$site
pc1$spp=kos$spp
pc1$combo=kos$combo
head(pc1)

mycol=c("purple4","orchid2",  "purple4","orchid2")
mypch=c(17,17,15,15)
pdf('c3_Kosrae_PCA.pdf', width = 5.75, height = 4)
ggplot(pc1, aes(Axis.1, Axis.2, color = combo, pch = combo)) +
  geom_point(size=3) +
  # geom_text_repel(aes(label=Samples)) +
  scale_colour_manual(values=mycol)+
  scale_shape_manual(values=mypch)+
theme_bw() +
  stat_ellipse()+
  xlab("PC1") +
  ylab("PC2") 
dev.off()

adonis(dist ~ site+spp, data = kos, method='bray')
          Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
site       1    1.9671 1.96710  59.019 0.38639  0.001 ***
spp        1    0.7241 0.72411  21.725 0.14223  0.001 ***
Residuals 72    2.3998 0.03333         0.47138           
Total     74    5.0910                 1.00000 



END
####################################################################################################################################################################################################################################


##not using
dat <- read.table("C3_AllData.csv", sep=",", header=TRUE, row.names=1)
summary(dat)
head(dat)
nrow(dat)
gendat=dat[,-c(1:3)]
head(gendat)
plot(hclust(vegdist(gendat,binary=T, method="manhattan")), cex=0.2)
hc=hclust(vegdist(gendat,binary=T, method="manhattan"))

library(plotrix) # for color.scale()
library(ape)
pops=paste(dat$spp, dat$island)
pops=paste(dat$spp)
cols=WGCNA::labels2colors(pops)
#plot(as.phylo(hc),cex=0.5,tip.color=color.scale(log(meta$H),c(1,0,1),c(1,1,0),alpha=1))
plot(as.phylo(hc),type="unrooted", cex=0.3,tip.color=cols)
plot(as.phylo(hc), cex=0.3,tip.color=cols)

kos=subset(dat, island=="Kosrae")
gendat=kos[,-c(1:3)]
nrow(gendat)
head(gendat)
pops=paste(kos$spp, kos$site)
cols=WGCNA::labels2colors(pops)
#plot(as.phylo(hc),cex=0.5,tip.color=color.scale(log(meta$H),c(1,0,1),c(1,1,0),alpha=1))
hc=hclust(vegdist(gendat,binary=T, method="manhattan"))
plot(as.phylo(hc),type="unrooted", cex=0.3,tip.color=cols)
plot(as.phylo(hc), cex=0.3,tip.color=cols)
