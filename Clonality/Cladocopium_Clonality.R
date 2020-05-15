##Clonality analysis on Cladocopium C40 and C3 presence/absense allele data
##Davies et a., 2020 Molecular Ecology

library(vegan)
library(plotrix)
library(adegenet)
library(ggplot2)
library(ape)
library(WGCNA)
library(GO.db)
library(dplyr)
library(tidyr)
library(plyr)

########################################
######ASSIGNING GENOTYPE GROUPS#########
########################################
###Note C3 was later designated as C21
c3dat <- read.csv("C3_alldata.csv")
c3gendat = c3dat[,-c(2:4)]
plot(hclust(vegdist(c3gendat, binary=T, method="manhattan")), cex=0.2)

plot(hclust(vegdist(c3gendat,binary=T, method="manhattan")), cex=0.2)
c3hc=hclust(vegdist(c3gendat,binary=T, method="manhattan"))

#basic
c3pops=paste(c3dat$spp, c3dat$island)
cols=WGCNA::labels2colors(c3pops)

#to look specifically at other tip labels
#c3spp=paste(c3dat$spp)
#c3site=paste(c3dat$site)
#c3isl=paste(c3dat$island)
#colspp=WGCNA::labels2colors(c3spp)
#colsite=WGCNA::labels2colors(c3site)
#colisl=WGCNA::labels2colors(c3isl)

# tree
plot(as.phylo(c3hc),type="fan", cex=0.3,tip.color=cols)

c3cut <- cutree(c3hc, h=.2) #group assignment

#make dataframes
C3_alldata_frame <- data.frame(C3_alldata)
C3_alldata_frame$x1 <- c3cut #add group assignment to dataframe
C3_alldata_frame <- subset(C3_alldata_frame, select=c(sample:site,x1,X92:X279))
#head(C3_alldata_frame)

#write.table(C3_alldata_frame, file = "C3_workingtable.csv") 

#repeat for C40#
c40dat <- read.csv("C40_alldata.csv")
c40gendat = c40dat[,-c(2:4)]

plot(hclust(vegdist(c40gendat,binary=T, method="manhattan")), cex=0.2)
c40hc=hclust(vegdist(c40gendat,binary=T, method="manhattan"))
#c40pops=paste(c40dat$spp, c40dat$island)
cols=WGCNA::labels2colors(c40pops)

#good plot
plot(as.phylo(c40hc),type="fan", cex=0.3,tip.color=cols)
tree<-as.phylo(c40hc)

c40cut <- cutree(c40hc, h=.2) #group assignment

C40_alldata_frame <- data.frame(C40_alldata)
C40_alldata_frame$x1 <- c40cut #add group assignment to dataframe
C40_alldata_frame <- subset(C40_alldata_frame, select=c(sample:site,x1,X92:X279))
#head(C40_alldata_frame)

#write.table(C40_alldata_frame, file = "C40_workingtable.csv")


###########################################
############SIMULATION#####################
###########################################

###remove alleles that only happen once####
#C40#
c40dat <- read.csv("C40_alldata.csv")
C40alleles <- subset(c40dat, select=c(X92:X279))

rpc40alleles<-C40alleles

for (i in 1:ncol(rpc40alleles)){
  if(sum(rpc40alleles[,i])==1){
    rpc40alleles[,-i]->rpc40alleles
  }
}

#simulation#
simC40<-as.matrix(rpc40alleles)
sim_c40<-matrix(NA, ncol = ncol(simC40), nrow=100000)

for (i in 1:ncol(simC40)){
  x <- vector()
  for (j in 1:100000){
    x[j] <- sample(simC40[,i], 1)
  }
  sim_c40[,i] <- x
}

#C3#
c3dat <- read.csv("C3_alldata.csv")
C3alleles <- subset(c3dat, select=c(X92:X279)

#remove non-repeating alleles#
rpc3alleles<-C3alleles

for (i in 1:ncol(rpc3alleles)){
  if(sum(rpc3alleles[,i])==1){
    rpc3alleles[,-i]->rpc3alleles
  }
}

#simulation#
simC3<-as.matrix(rpc3alleles)
sim_c3<-matrix(NA, ncol = ncol(simC3), nrow=100000)

for (i in 1:ncol(simC3)){
  x <- vector()
  for (j in 1:100000){
    x[j] <- sample(simC3[,i], 1)
  }
  sim_c3[,i] <- x
}

sim_C3_artificial <- as.data.frame(sim_c3)
colnames(sim_C3_artificial) <- colnames(simC3)

colnames(simC3) -> colnamesC3
colnames(simC40) -> colnamesC40

unique_C3alleles <- subset(C3_alldata_frame, select=c(x1,X92:X279))
unique_C3alleles <- unique(unique_C3alleles)

unique_C40alleles <- subset(C40_alldata_frame, select=c(x1, X92:X279))
unique_C40alleles <- unique(unique_C40alleles)

unique_C40alleles %>% select(x1,colnamesC40) -> workc40
unique_C3alleles %>% select(x1,colnamesC3) -> workc3

#count matches 
library(plyr)
workc3 %>% filter(x1==293) -> of_interest
nrow(match_df(sim_C3_artificial, of_interest))
#record number of matches in excel


###############################################
##########make stacked barplots################
#y is counts, x is site, stratify for species##
###############################################

#one of these at a time
plot<-C40_alldata_frame
plot<-C3_alldata_frame

#then go through the rest
plot$group = plot$x1
plot %>%
  add_count(x1) %>%
  mutate(group = ifelse(n==1, "not clones", "clones")) %>% #either yes or no
  dplyr::select(-n) %>%
  add_count(site, group,spp) %>% select(-X92:-X279)-> plot 
head(plot)
plot$abrv<-plot$site

plot %>% 
  mutate(abrv = ifelse(site == "South East Pass","SEP", abrv)) %>%
  mutate(abrv = ifelse(site == "West Polle","WP", abrv)) %>%
  mutate(abrv = ifelse(site == "Pago Bay","PB", abrv)) %>%
  mutate(abrv = ifelse(site == "Tanguisson","T", abrv)) %>%
  mutate(abrv = ifelse(site == "Coral Garden","CG", abrv)) %>%
  mutate(abrv = ifelse(site == "Hiroshi Point","HP", abrv)) %>%
  mutate(abrv = ifelse(site == "Lighthouse Reef","LR", abrv)) %>%
  mutate(abrv = ifelse(site == "West Channel","WC", abrv)) %>%
  mutate(abrv = ifelse(site == "Ant Atoll","AA", abrv)) %>%
  mutate(abrv = ifelse(site == "Roj","R", abrv)) %>%
  mutate(abrv = ifelse(site == "Goofnuw Channel","GC", abrv)) %>%
  mutate(abrv = ifelse(site == "South Tip","ST", abrv)) -> plot


plot<-transform(plot,
                abrv=factor(abrv,levels=c("SEP","WP",
                                    "PB", "T",
                                    "CG", "HP",
                                    "LR","WC",
                                    "AA", "R",
                                    "GC", "ST")))

plot$species<-plot$spp
plot%>% mutate(species=ifelse(spp=="A. digitifera", "d", species))%>%
  mutate(species=ifelse(spp=="A. hyacinthus", "h", species))->plot

ggplot(plot, aes(fill=group,x=species)) +
  geom_bar(position="fill") +
  geom_text(aes(label = ifelse(group == "clones", n, "")), y = 0.95, size=2) +
  geom_text(aes(label = ifelse(group == "not clones", n, "")), y = 0.05, size =2) +
  facet_wrap(~ abrv, scales = "free_x", ncol = 11) +
  scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 5)) +
  theme_bw() +
  labs(x="", y="Proportion of Clones", title = "C40") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values = c("gray57", "lightgrey"))
  #scale_color_manual(values = c25, na.value = "gray57")+
  #scale_fill_manual(values = c25, na.value = "gray57")

###############################################################
##############genotype frequency of MLG for each site##########
###############################################################

#only one of these at a time
gfdata <- C3_alldata_frame
gfdata <- C40_alldata_frame

#then continue with the rest
#this counts the MLGs at each site and gives clonal fraction
#input renamed dataframe
gfdata %>% 
  group_by(site) %>% 
  dplyr::summarise(N = n(), MLG = length(unique(x1))) %>%
  mutate(clonal_fraction = 1-(MLG/N))

##this calculates genotype diversity by 1-(sum of squared MLG f)##
#############
##by c3/c40##
#############

#only one of these at a time
gfdata <- C3_alldata_frame
gfdata <- C40_alldata_frame

#then continue with the rest

gfdata %>% 
  group_by(site) %>%
  dplyr::mutate(total_N = n()) %>%
  ungroup(site) %>%
  group_by(site,x1,total_N) %>%
  dplyr::summarise(N=n()) %>%
  filter(total_N >= 3) %>%
  #filter(site == "Roj") #%>% #check total_N in specific site
  mutate(p_square = (N/total_N)^2) %>%
  group_by(site) %>%
  #dplyr::summarise(total_N = median(total_N)) %>% #check total individuals in each leftover group
  dplyr::summarise(sum = sum(p_square)) %>% 
  mutate(gd = 1-sum) %>%
  dplyr::select(-sum) -> c40gddata #change the name for c3/c40!

#make groups
c3gddata %>%
  mutate(group = "C3") -> c3gddata 
c40gddata %>%
  mutate(group = "C40") -> c40gddata
gddata<- merge(c3gddata, c40gddata, all = TRUE)

#make violin##
library(ggplot2)
colors<-c('#b2df6eff','#b9c3adff')
ggplot(gddata, aes(x = group, y = gd, y, fill=group)) +
  geom_violin(trim = FALSE,
              #fill = "darkgrey", 
              #width=1, 
              position=position_dodge(),
              scale = "width",
              aes(alpha=.9)) +
  geom_jitter(shape=16, position=position_jitter(0.2)) +
  geom_boxplot(width=0.045, fill="white"
               #aes(alpha=0.9)
               ) +
  theme_bw() +
  scale_fill_manual(values = colors)+
  ylim(0.4, 1) +
  labs(x="", y="Genetic Diversity") +
  #scale_x_discrete(limits = c("A. digitifera",
  # "A. hyacinthus"), expand=c(0.8, 0.1)) +
  theme(legend.position="none") #+
#ylim(0,1)

##########
##by spp##
##########
library(ggplot2)
library(dplyr)

#only one of these at a time#
gfdata <- C40_alldata_frame
gfdata <- C3_alldata_frame

#continue with the rest
gfdata %>% dplyr::filter(spp=="A. digitifera") -> gfdigi
gfdata %>% dplyr::filter(spp=="A. hyacinthus") -> gfhya

gfhya %>% 
  group_by(site) %>%
  dplyr::mutate(total_N = n()) %>% #how many samples at each site
  ungroup(site) %>%
  group_by(site,x1,total_N) %>%
  dplyr::summarise(N=n()) %>% #how many of site and group at each site
  filter(total_N >= 3) %>% 
  #filter(site == "Roj") #%>% #check total_N in specific site
  mutate(p_square = (N/total_N)^2) %>% 
  group_by(site) %>%
  #dplyr::summarise(total_N = median(total_N)) %>% #check total individuals in each leftover group
  dplyr::summarise(sum = sum(p_square)) %>%
  mutate(gd = 1-sum) %>%
  dplyr::select(-sum) -> gdhya

gfdigi %>% 
  group_by(site) %>%
  dplyr::mutate(total_N = n()) %>% #how many samples at each site
  ungroup(site) %>%
  group_by(site,x1,total_N) %>%
  dplyr::summarise(N=n()) %>% #how many of site and group at each site
  filter(total_N >= 3) %>% 
  #filter(site == "Roj") #%>% #check total_N in specific site
  mutate(p_square = (N/total_N)^2) %>% 
  group_by(site) %>%
  #dplyr::summarise(total_N = median(total_N)) %>% #check total individuals in each leftover group
  dplyr::summarise(sum = sum(p_square)) %>%
  mutate(gd = 1-sum) %>%
  dplyr::select(-sum) -> gddigi

#add group labels
gddigi %>%
  mutate(group = "A. digitifera") -> gddigi
gdhya %>%
  mutate(group = "A. hyacinthus") -> gdhya
gddata<- merge(gddigi, gdhya, all = TRUE)

#make violin##
colors<-c("orchid2","purple4")
ggplot(gddata, aes(x = group, y = gd, y, fill=group)) +
  geom_violin(trim = FALSE,
              #fill = "darkgrey", 
              #width=1, 
              position=position_dodge(width=.5),
              scale="width",
              aes(alpha=0.5)) +
  geom_jitter(shape=16, position=position_jitter(0.2)) +
  geom_boxplot(width=0.045, fill="white", aes(alpha=0.5)) +
  theme_bw() +
  scale_fill_manual(values = colors)+
  ylim(0.4, 1) +
  labs(x="Species", y="Genetic Diversity") +
  #scale_x_discrete(limits = c("A. digitifera",
                             # "A. hyacinthus"), expand=c(0.8, 0.1)) +
  theme(legend.position="none") #+
  #ylim(0,1)

#####################
#####histogram#######
#####################
library(ggplot2)
colours=c("#b2df6eff", "#b9c3adff")
##barplots2 was generated in excel##
ggplot(barplots2, aes(x=count, fill=group)) +
  geom_histogram(
    #position="identity",
    position="dodge",
    colour="grey40", 
    #alpha=.4, 
    bins = 13) +
  scale_color_manual(values = colours)+
  scale_fill_manual(values = colours)+
  theme_bw()+
  labs(x="Group Size", y="Counts")

#f barplot
library(ggplot2)
library(dplyr)
colours=c("#b2df6eff", "#b9c3adff")
##barplots was generated in excel##
barplots$total <- barplots$yes+barplots$no
barplots$f <- barplots$yes/barplots$total
barplots %>% select(-yes, -no) -> bar
ggplot(barplots, aes(x=clones, y=f, fill=group)) +
  geom_bar(color="grey40", alpha=1,
           #position="dodge",
           #position="identity",
           stat="identity") +
  scale_color_manual(values = colours)+
  scale_fill_manual(values = colours)+
  geom_text(aes(label = total), 
            position=position_dodge(width=0.9), 
            vjust=-0.25,size=3.0)+
 # geom_text(aes(label = ifelse(group == "C40", yes, "")),
           # size=2.5, vjust=-.25) +
  #geom_text(aes(label = ifelse(group == "C3", yes, "")),
           # y = 0.05, size =2.5) 
  facet_wrap(~group)+
  theme_bw()+
  labs(x="Group Size", y="Frequency")
  
###############################
#####distance scatterplots#####
###############################
library(dplyr)

##creating distance matrix
C3_alldata_frame %>% select(-X92:-X279) -> c3info
c3island %>% select(-num_groups:-gd) -> c3coord
merge(c3info, c3coord) -> c3dist
c3dist<- subset(c3dist, select=c(site, lat, long))
c3distuniq<-unique(c3dist)
row.names(c3distuniq) <- c3distuniq$site
c3distuniq<-subset(c3distuniq, select=c(lat, long))

#converting dms to decimal degrees
dms2dec <- function(dms, separators = c("º", "°", "\'", "\"")) {
  # version 1.0 (25 Sep 3013)
  # dms: a vector (or column) of latitude or longitude in degrees-minutes-seconds-hemisfere, e.g. 41° 34' 10.956" N (with or without spaces)
  # separators: the characters that are separating degrees, minutes and seconds in dms
  
  dms <- as.character(dms)
  dms <- gsub(pattern = " ", replacement = "", x = dms)
  for (s in separators) dms <- gsub(pattern = s, replacement = "_splitHere_", x = dms)
  
  splits <- strsplit(dms, split = "_splitHere_")
  n <- length(dms)
  deg <- min <- sec <- hem <- vector("character", n)
  
  for (i in 1:n) {
    deg[i] <- splits[[i]][1]
    min[i] <- splits[[i]][2]
    sec[i] <- splits[[i]][3]
    hem[i] <- splits[[i]][4]
  }
  
  dec <- as.numeric(deg) + (as.numeric(min) / 60) + (as.numeric(sec) / 3600)
  sign <- ifelse (hem %in% c("N", "E"), 1, -1)
  dec <- sign * dec
  return(dec)
}  # end dms2dec function,  Zanolla et al. (2018)



c3distuniq$lat.dec <- dms2dec(c3distuniq$lat)
c3distuniq$long.dec <- dms2dec(c3distuniq$long)

c3tomat<-subset(c3distuniq, select = c(lat.dec, long.dec))
as.matrix(dist(c3tomat)) -> c3mat


##finding mismatched sites in a group
c3info %>% 
  ungroup(x1) %>% group_by(x1, site) %>% summarise(N=n()) %>% 
  group_by(x1) %>% filter(n()>1) %>% #now only clones with mismatched sites
  group_by(x1) ->c3sum

c3info %>% group_by(x1) %>% 
  mutate(n=n()) %>% ungroup(x1) %>% 
  select(x1, n) %>% filter(n>1)-> c3counts

merge(c3sum, c3counts) -> c3sum
c3sum <- unique(c3sum)

c3sum %>% select(x1, n) %>% unique() -> c3counts

#how many members in the group
c3info %>% filter(x1=="243") %>% nrow() #put in group number of interest
#find the sites the group spans
c3sum %>% filter(x1=="243") #put in group number of interest, generates sites of interest
c3sum %>% filter(site=="Coral Garden"|site=="Hiroshi Point") %>% 
  group_by(x1) %>% filter(n()>1) #check if any other groups also span these sites

##recording largest distance in those groups to excel
c3mat["Coral Garden", "Hiroshi Point"] #replace with sites of interest


####C40#####
##creating distance matrix
C40_alldata_frame %>% select(-X92:-X279) -> c40info
c40island %>% select(island:total) -> c40coord
merge(c40info, c40coord) -> c40dist
#c3dist$site <- paste(c3dist$island,"-",c3dist$site)
c40dist<- subset(c40dist, select=c(site, lat, long))
c40distuniq<-unique(c40dist)
row.names(c40distuniq) <- c40distuniq$site
c40distuniq<-subset(c40distuniq, select=c(lat, long))

dms2dec <- function(dms, separators = c("º", "°", "\'", "\"")) {
  # version 1.0 (25 Sep 3013)
  # dms: a vector (or column) of latitude or longitude in degrees-minutes-seconds-hemisfere, e.g. 41° 34' 10.956" N (with or without spaces)
  # separators: the characters that are separating degrees, minutes and seconds in dms
  
  dms <- as.character(dms)
  dms <- gsub(pattern = " ", replacement = "", x = dms)
  for (s in separators) dms <- gsub(pattern = s, replacement = "_splitHere_", x = dms)
  
  splits <- strsplit(dms, split = "_splitHere_")
  n <- length(dms)
  deg <- min <- sec <- hem <- vector("character", n)
  
  for (i in 1:n) {
    deg[i] <- splits[[i]][1]
    min[i] <- splits[[i]][2]
    sec[i] <- splits[[i]][3]
    hem[i] <- splits[[i]][4]
  }
  
  dec <- as.numeric(deg) + (as.numeric(min) / 60) + (as.numeric(sec) / 3600)
  sign <- ifelse (hem %in% c("N", "E"), 1, -1)
  dec <- sign * dec
  return(dec)
}  # end dms2dec function,  Zanolla et al. (2018)



c40distuniq$lat.dec <- dms2dec(c40distuniq$lat)
c40distuniq$long.dec <- dms2dec(c40distuniq$long)

c40tomat<-subset(c40distuniq, select = c(lat.dec, long.dec))
as.matrix(dist(c40tomat)) -> c40mat


##finding mismatched sites in a group
c40info %>% group_by(x1) %>% mutate(n=n()) %>% filter(n>1) %>% #now only clone groups
  ungroup(x1) %>% group_by(x1, site) %>% mutate(N=n()) %>% group_by(x1) %>% filter(n()>1) #now only mismatched sites
  ->c40sum

c40info %>% 
  ungroup(x1) %>% group_by(x1, site) %>% summarise(N=n()) %>% 
  group_by(x1) %>% filter(n()>1) %>% #now only clones with mismatched sites
  group_by(x1) ->c40sum

c40info %>% group_by(x1) %>% 
  mutate(n=n()) %>% ungroup(x1) %>% 
  select(x1, n) %>% filter(n>1)-> c40counts
  
merge(c40sum, c40counts) -> c40sum
c40sum <- unique(c40sum)

c40sum %>% select(x1, n) %>% unique() -> c40counts

c40sum %>% filter(x1=="2")
c40sum %>% filter(site=="Lighthouse Reef"|site=="Goofnuw Channel") %>% 
  group_by(x1) %>% filter(n()>1)
##recording largest distance in those groups from the following
c40mat["West Channel", "Goofnuw Channel"]

#append c3/c40counts to groups file
c40ddata %>% filter(dist!="NA") %>% select(group:dist) -> c40groups
c40groups$n<-c40groups$count
c40groups$x1 <- c40groups$group
c40groups %>% select(-count, -group) -> c40groups
merge(c40counts, c40groups, all = TRUE) -> c40scat

c3ddata %>% filter(dist!="NA") %>% select(repeating_group:dist) -> c3groups
c3groups$n<-c3groups$count
c3groups$x1 <- c3groups$repeating_group
c3groups %>% select(-count, -repeating_group) -> c3groups
merge(c3counts, c3groups, all = TRUE) -> c3scat

c40scat$group<-"C40"
c3scat$group<-"C3"
merge(c3scat, c40scat, all=TRUE) -> mergetest


#####plot scatterplots#####
ggplot(c3scat, aes(x=dist, y=n)) +
  geom_point(color="#b2df6eff")+
  theme_bw()+
    ylim(0,15)

ggplot(c40scat, aes(x=dist, y=n)) +
  geom_point()+
  theme_bw()+
  ylim(0,15)

colours=c("#b2df6eff", "#b9c3adff")
ggplot(mergetest, aes(x=dist, y=n, color=group)) +
  geom_point()+
  theme_bw()+
  scale_color_manual(values = colours)+
  scale_fill_manual(values = colours)+
  ylim(0,15)+
  facet_wrap(~group)+
  labs(x="Distance", y="Group Size")

##for each row in the x1 col of test, take the value and 

##########only have one of each genotype group###############
##########create dataframe for re-analysis###################
C3_lessdata <- C3_alldata_frame
C3_lessdata %>% 
  group_by(x1) %>% #one rep with independent features
  filter(row_number() == 1) %>%
  ungroup(x1) %>%
  dplyr::select(-x1)-> C3_lessdata
write.table(C3_lessdata, file = "C3_lessdata.csv")
#head(C3_lessdata) #used in other analysis

C40_lessdata <- C40_alldata_frame
C40_lessdata %>%
  group_by(x1) %>%
  filter(row_number() == 1) %>%
  ungroup(x1) %>%
  dplyr::select(-x1) -> C40_lessdata
#write.table(C40_lessdata, file = "C40_lessdata.csv") #used in other analysis

#########################
#########TREES###########
#########################

#make unlabelled trees
c3dat <- read.csv("C3_alldata.csv")
c3gendat = c3dat[,-c(2:4)]
c3hc=hclust(vegdist(c3gendat,binary=T, method="manhattan"))
plot(as.phylo(c3hc),type="fan", cex=0.3, show.tip.label=FALSE)

c40dat <- read.csv("C40_alldata.csv")
c40gendat = c40dat[,-c(2:4)]
c40hc=hclust(vegdist(c40gendat,binary=T, method="manhattan"))
plot(as.phylo(c40hc),type="fan", cex=0.3, show.tip.label=FALSE)


##tree set-up##
#C3#
c3pops=paste(c3dat$spp, c3dat$island)
c3pops=paste(c3dat$site)
cols=WGCNA::labels2colors((c3pops))
c3hc$order -> c3order
as.phylo(c3hc) -> phylo
plot(phylo,type="fan", cex=0.3,tip.color=cols) #check labels

labels(phylo)<-c3pops[c3order]

#C40#
c40pops=paste(c40dat$spp, c40dat$island)
c40pops=paste(c40dat$site)
cols=WGCNA::labels2colors((c40pops))
c40hc$order -> c40order
as.phylo(c40hc) -> phylo
plot(phylo,type="fan", cex=0.3,tip.color=cols) #check labels

labels(phylo)<-c40pops[c40order]

######################################
####tree island color wheel labels####
######################################
#C3#
c3island=paste(c3dat$island)
c3islorder<-c3island[c3order]
colisl=WGCNA::labels2colors(c3island)
library(ape)
plot(as.phylo(c3hc),type="fan", cex=0.3,tip.color=colisl) -> isl
isl$tip.color -> isltips #extracted character string of colors
isltips[c3order] -> c3islcolorder

gsub("brown","brown4",c3islcolorder)->c3islcolorder #kosrae
gsub("blue","aquamarine3",c3islcolorder)->c3islcolorder #guam
gsub("turquoise","steelblue2",c3islcolorder)->c3islcolorder #chuuk
gsub("red","goldenrod1",c3islcolorder)->c3islcolorder #yap
gsub("yellow","red2",c3islcolorder)->c3islcolorder #palau
gsub("green","royalblue4",c3islcolorder)->c3islcolorder #pohnpei

#do.call(rbind, Map(data.frame, A=c3islorder, B=c3islcolorder)) -> c3islmerge
#unique(c3islmerge) ->c3islunique ##used this as a visual reference

#creating species color labels#
c3spp=paste(c3dat$spp)
colspp=WGCNA::labels2colors(c3spp)
par(mar = c(1, 1, 1, 1))
library(ape)
plot(as.phylo(c3hc),type="fan", cex=0.3,
     tip.color=colspp, show.tip.label=FALSE) -> spp
spp$tip.color -> spptips
spptips[c3order] -> spporder #extracted character string of colors
gsub("turquoise", "purple4", spporder) -> c3colspporder #get right colors
gsub("blue", "orchid2", c3colspporder) -> c3colspporder #get right colors

#C40#
c40island=paste(c40dat$island)
c40islorder<-c40island[c40order]
colisl=WGCNA::labels2colors(c40island)
library(ape)
plot(as.phylo(c40hc),type="fan", cex=0.3,tip.color=colisl) -> isl
isl$tip.color -> isltips #extracted character string of colors
isltips[c40order] -> c40islcolorder

gsub("brown","tomato1",c40islcolorder)->c40islcolorder #ngulu
gsub("blue","brown4",c40islcolorder)->c40islcolorder #kosrae
gsub("turquoise","steelblue2",c40islcolorder)->c40islcolorder #chuuk
gsub("red","goldenrod1",c40islcolorder)->c40islcolorder #yap
gsub("yellow","red2",c40islcolorder)->c40islcolorder #palau
gsub("green","royalblue4",c40islcolorder)->c40islcolorder #pohnpei

#do.call(rbind, Map(data.frame, A=c40islorder, B=c40islcolorder)) -> c40islmerge
#unique(c40islmerge) ->c40islunique #visual check

###############################
#creating species color labels#
###############################
c40spp=paste(c40dat$spp)
colspp=WGCNA::labels2colors(c40spp)
par(mar = c(1, 1, 1, 1))
as.phylo(c40hc) -> spp
plot(spp,type="fan", cex=0.3,
     tip.color=colspp#, show.tip.label=FALSE
     ) -> spp
plot(spp,type="fan", cex=0.3,tip.color=colspp) 
labels(spp)<-c40spp[c40order]
spp$tip.color -> spptips
spptips[c40order] -> spporder #extracted character string of colors
gsub("turquoise", "purple4", spporder) -> colspporder #digi
gsub("blue", "orchid2", colspporder) -> colspporder #hya


###make a color wheel with amount of nodes
#c40 has 172 tips 
#c3 has 388 tips

# Simple Pie Chart

###C3 island wheel##
slices <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
            1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
            1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
            1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
            1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
            1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
            1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
            1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)

colors<-c3islcolorder #C3 island wheel
pie(slices, radius =1, labels=NA, density=NULL, col=colors)

###C3 spp wheel##
slices <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
            1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
            1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
            1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
            1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
            1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
            1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
            1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)

colors<-c3colspporder #C3 species wheel
pie(slices, radius =1, labels=NA, density=NULL, col=colors)

##C40 island wheel##
slices <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
            1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
            1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
            1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1
            )

colors<-c40islcolorder #C40 island wheel
pie(slices, radius =1, labels=NA, density=NULL, col=colors)

##C40 species wheel##
slices <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
            1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
            1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
            1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1
            )

colors<-colspporder #C40 species wheel
pie(slices, radius =1, labels=NA, density=NULL, col=colors)

#colors[1]="red"


######all images were exported as an .svg file and edited in Inkscape to fit pub.
#for tree, the wheel and unlabeled trees were aligned
#for barplots, site color labels were added to the facet labels and species colors were also added