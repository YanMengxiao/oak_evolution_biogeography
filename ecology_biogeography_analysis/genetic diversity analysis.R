# 多样性统计分析


library(geosphere)
library(ape)
library(phytools)
library(vegan)

setwd("E:/research/analysis/cpdna biogeography quercus/climate/mydata")

############################
###        ILEX          ###
############################

clim.ilex<-read.csv("ilex_res2.5.csv",as.is=T)
clim.ilex<-clim.ilex[-grep("DM6572|DM6938|DM3688|DM8227|DM8229",clim.ilex$sequence),]
# colnames(pcoa.ilex)<-c("sequence","axis1")
# clim.ilex$axis1<-pcoa.ilex$axis1[match(clim.ilex$tip.label,pcoa.ilex$sequence)]
# write.csv(clim.ilex,"clim.ilex.csv")
##多重线性分析，检测相关气候因子
# library(psych)
# corr.test.ilex<-corr.test(clim.ilex[,4:22], use="complete")
# write.csv(corr.test.ilex[[1]],"corr.ilex.csv")
# write.csv(corr.test.ilex[[4]],"corr.test.p-value.ilex.csv")
# cor(clim.ilex[,4:22],method = "pearson")
# 矩阵计算与读入
# dist.clim.ilex<-as.vector(vegdist(clim.ilex[,4:22],method="manhattan",upper = T,diag=T))
# write.csv(dist.clim.ilex,"ilex.clim2.csv")
# 遗传距离计算 Pairwise Distances From A Phylogenetic Tree
# library(ape)
# tr.ilex<-read.tree("ilex.newick")
# tr.ilex<-drop.tip(tr.ilex,grep("Litho|Q_aquifolioides_KP340971|Q_aquifolioides_NC026913|DM6572|DM6938|DM3688|DM8227|DM8229",tr.ilex$tip.label))
# # tr.ilex<-drop.tip(tr.ilex,setdiff(tr.ilex$tip.label,clim.ilex$tip.label))
# # dist.ilex.gentic<-cophenetic(tr.ilex)

library(vegan)
dist.clim.ilex<-vegdist(clim.ilex[,6:24],method="gower",upper = T,diag=T)
select.ilex<-paste("bio",c(1:4,7:9,14:15,19),sep="")
dist.clim.selc.ilex<-vegdist(clim.ilex[,select.ilex],method="gower",upper = T,diag=T)
# dist.clim.ilex.m<-vegdist(clim.ilex[,4:22],method="manhattan",upper = T,diag=T)
# select.ilex<-paste("bio",c(1:4,7:9,14:15,19),sep="")
# dist.clim.selc.ilex.m<-vegdist(clim.ilex[,select.ilex],method="manhattan",upper = T,diag=T)
# 之前遗传距离计算自genalex导出
dist.ilex.gentic<-read.csv("ilex.genetic.distance.csv",head=F)
dimnames(dist.ilex.gentic)<-list(1:119,1:119)
ilex.geo.dist<-read.csv("ilex.geo.dist.csv", head=F)
dimnames(ilex.geo.dist)<-list(1:119,1:119)

##glm
ilex.dist.clim<-data.frame(clim=as.vector(dist.clim.ilex),
                           selc=as.vector(dist.clim.selc.ilex))
write.csv(ilex.dist.clim,"ilex.dist.clim.csv")
ilex.all<-read.csv("ilex.dist.clim.csv")
plot(density(ilex.all$genetic))
plot(density(ilex.all$clim,na.rm=T))
glm.ilex.GE<-glm(genetic~clim,quasi(link = "identity", variance = "constant"),data=ilex.all)
summary(glm.ilex.GE)
glm.ilex.GG<-glm(genetic~geo,quasi(link = "identity", variance = "constant"),data=ilex.all)
summary(glm.ilex.GG)
glm.ilex.GGE<-glm(genetic~clim + geo,quasi(link = "identity", variance = "constant"),data=ilex.all)
summary(glm.ilex.GGE)
glm.ilex.GGES<-glm(genetic~selc + geo,quasi(link = "identity", variance = "constant"),data=ilex.all)
summary(glm.ilex.GGES)
glm.ilex.GGEM<-glm(genetic~clim.m + geo,quasi(link = "identity", variance = "constant"),data=ilex.all)
summary(glm.ilex.GGEM)
glm.ilex.GGESM<-glm(genetic~selc.m + geo,quasi(link = "identity", variance = "constant"),data=ilex.all)
summary(glm.ilex.GGESM)



# library(vegan)
# mantel.ilex.GE<-mantel(dist.ilex.gentic, dist.ilex, method="pearson", permutations=999)
# mantel.ilex.GG<-mantel(dist.ilex.gentic, ilex.geo.dist, method="pearson", permutations=999)
library(phytools)
mantel.ilex.GE<-multi.mantel(dist.ilex.gentic, dist.clim.ilex, nperm=999)
mantel.ilex.GE[1:6]
mantel.ilex.GES<-multi.mantel(dist.ilex.gentic, dist.clim.selc.ilex, nperm=999)
mantel.ilex.GES[1:6]
mantel.ilex.GG<-multi.mantel(dist.ilex.gentic, as.matrix(ilex.geo.dist), nperm=999)
mantel.ilex.GG[1:6]
mantel.ilex.GGE<-multi.mantel(dist.ilex.gentic, list(dist.clim.ilex,as.matrix(ilex.geo.dist)), nperm=999)
mantel.ilex.GGE[1:6]
mantel.ilex.GGES<-multi.mantel(dist.ilex.gentic, list(dist.clim.selc.ilex,as.matrix(ilex.geo.dist)), nperm=999)
mantel.ilex.GGES[1:6]
# library(vegan)
mantel.ilex.partial<-mantel.partial(dist.ilex.gentic, dist.clim.selc.ilex,ilex.geo.dist, permutations = 999)
mantel.ilex.partial


############################
###        CYC          ###
############################

clim.cyc<-read.csv("cyc_res2.5.csv",as.is=T)
# pcoa.cyc<-read.csv("cyc.pcoa.csv")
# colnames(pcoa.cyc)<-c("sequence","axis1")
# clim.cyc$axis1<-pcoa.cyc$axis1[match(clim.cyc$tip.label,pcoa.cyc$sequence)]
# write.csv(clim.cyc,"clim.cyc.csv")
##多重线性分析，检测相关气候因子
library(psych)
corr.test.cyc<-corr.test(clim.cyc[,5:23], use="complete")
write.csv(corr.test.cyc[[1]],"corr.cyc.csv")
# write.csv(corr.test.cyc[[4]],"corr.test.p-value.cyc.csv")
# cor(clim.cyc[,4:22],method = "pearson")
dist.clim.cyc<-vegdist(clim.cyc[,4:22],method="gower",upper = T,diag=T)
select.cyc<-paste("bio",c(1:5,7,12,14),sep="")
dist.clim.selc.cyc<-vegdist(clim.cyc[,select.cyc],method="gower",upper = T,diag=T)
# write.csv(as.vector(dist.clim.cyc),"cyc.clim.dist.csv")
dist.cyc.gentic<-read.csv("cyc.genetic.dist.csv",head=F)
dimnames(dist.cyc.gentic)<-list(1:147,1:147)
cyc.geo.dist<-read.csv("cyc.geo.dist.csv", head=F)
dimnames(cyc.geo.dist)<-list(1:147,1:147)


##glm
cyc.dist.clim<-data.frame(clim=as.vector(dist.clim.cyc),
                          selc=as.vector(dist.clim.selc.cyc))
write.csv(cyc.dist.clim,"cyc.dist.clim.csv")
cyc.all<-read.csv("cyc.all.dist.csv")
plot(density(cyc.all$genetic))
plot(density(cyc.all$clim,na.rm=T))
glm.cyc.GE<-glm(genetic~clim,gaussian(link = "identity"),data=cyc.all)
summary(glm.cyc.GE)
glm.cyc.GG<-glm(genetic~geo,gaussian(link = "identity"),data=cyc.all)
summary(glm.cyc.GG)
glm.cyc.GGE<-glm(genetic~clim + geo, gaussian(link = "identity"),data=cyc.all)
summary(glm.cyc.GGE)
glm.cyc.GGES<-glm(genetic~selc + geo, gaussian(link = "identity"),data=cyc.all)
summary(glm.cyc.GGES)
glm.cyc.GGEM<-glm(genetic~clim.m + geo, gaussian(link = "identity"),data=cyc.all)
summary(glm.cyc.GGEM)
glm.cyc.GGESM<-glm(genetic~selc.m + geo, gaussian(link = "identity"),data=cyc.all)
summary(glm.cyc.GGESM)

# library(vegan)
# mantel.cyc.GE<-mantel(dist.cyc.gentic, dist.cyc, method="pearson", permutations=999)
# mantel.cyc.GG<-mantel(dist.cyc.gentic, cyc.geo.dist, method="pearson", permutations=999)
# library(phytools)
mantel.cyc.GE<-multi.mantel(dist.cyc.gentic, dist.clim.cyc, nperm=999)
mantel.cyc.GE[1:6]
mantel.cyc.GES<-multi.mantel(dist.cyc.gentic, dist.clim.selc.cyc, nperm=999)
mantel.cyc.GES[1:6]
mantel.cyc.GG<-multi.mantel(dist.cyc.gentic, as.matrix(cyc.geo.dist), nperm=999)
mantel.cyc.GG[1:6]
mantel.cyc.GGE<-multi.mantel(dist.cyc.gentic, list(dist.clim.cyc,as.matrix(cyc.geo.dist)), nperm=999)
mantel.cyc.GGE[1:6]
mantel.cyc.GGES<-multi.mantel(dist.cyc.gentic, list(dist.clim.selc.cyc,as.matrix(cyc.geo.dist)), nperm=999)
mantel.cyc.GGES[1:6]
# library(vegan)
mantel.cyc.partial<-mantel.partial(dist.cyc.gentic, dist.clim.selc.cyc, cyc.geo.dist, permutations = 999)
mantel.cyc.partial

############################
###        PQ          ###
############################

clim.pq<-read.csv("pq_res2.5.csv",as.is=T)
clim.pq<-clim.pq[-grep("chrysolepis|tomentella|vacciniifolia|greggii",clim.pq$tip.label),]
# pcoa.pq<-read.csv("pq.pcoa.csv")
# colnames(pcoa.pq)<-c("sequence","axis1")
# clim.pq$axis1<-pcoa.pq$axis1[match(clim.pq$tip.label,pcoa.pq$sequence)]
# write.csv(clim.pq,"clim.pq.csv")
##多重线性分析，检测相关气候因子
# library(psych)
# corr.test.pq<-corr.test(clim.pq[,4:22], use="complete")
# write.csv(corr.test.pq[[1]],"corr.pq.csv")
# write.csv(corr.test.pq[[4]],"corr.test.p-value.pq.csv")
# # cor(clim.pq[,4:22],method = "pearson")
dist.clim.pq<-vegdist(clim.pq[,6:24],method="gower")
select.pq<-paste("bio",c(1:2,5,7:8,12:15,18:19),sep="")
dist.clim.selc.pq<-vegdist(clim.pq[,select.pq],method="gower")
# dist.clim.pq.m<-vegdist(clim.pq[,4:22],method="manhattan",upper = T,diag=T)
# select.pq<-paste("bio",c(1:2,4:5,8,12:15,18:19),sep="")
# dist.clim.selc.pq.m<-vegdist(clim.pq[,select.pq],method="manhattan",upper = T,diag=T)
dist.pq.gentic<-read.csv("pq.genetic.dist2.csv",head=F)
dimnames(dist.pq.gentic)<-list(1:89,1:89)
pq.geo.dist<-read.csv("pq.geo.dist2.csv", head=F)
dimnames(pq.geo.dist)<-list(1:89,1:89)

library(ape)
tr.pq<-read.tree("pq.newick")
tr.pq<-drop.tip(tr.pq,grep("Litho",tr.pq$tip.label))
tr.pq<-drop.tip(tr.pq,setdiff(tr.pq$tip.label,clim.pq$tip.label))
dist.pq.gentic<-cophenetic(tr.pq)


##glm
pq.dist.clim<-data.frame(clim=as.vector(dist.clim.pq),
                         selc=as.vector(dist.clim.selc.pq))
write.csv(pq.dist.clim,"pq.dist.clim2.csv")
pq.all<-read.csv("pq.dist.clim2.csv")
plot(density(pq.all$genetic,na.rm=T))
plot(density(pq.all$clim,na.rm=T))
glm.pq.GE<-glm(genetic~clim,quasi(link = "identity", variance = "constant"),data=pq.all)
summary(glm.pq.GE)
glm.pq.GG<-glm(genetic~geo,quasi(link = "identity", variance = "constant"),data=pq.all)
summary(glm.pq.GG)
glm.pq.GGE<-glm(genetic~clim + geo,quasi(link = "identity", variance = "constant"),data=pq.all)
summary(glm.pq.GGE)
glm.pq.GGES<-glm(genetic~selc + geo,quasi(link = "identity", variance = "constant"),data=pq.all)
summary(glm.pq.GGES)
glm.pq.GGEM<-glm(genetic~clim.m + geo,quasi(link = "identity", variance = "constant"),data=pq.all)
summary(glm.pq.GGEM)
glm.pq.GGESM<-glm(genetic~selc.m + geo,quasi(link = "identity", variance = "constant"),data=pq.all)
summary(glm.pq.GGESM)

# library(vegan)
# mantel.pq.GE<-mantel(dist.pq.gentic, dist.pq, method="pearson", permutations=999)
# mantel.pq.GG<-mantel(dist.pq.gentic, pq.geo.dist, method="pearson", permutations=999)
# library(phytools)
mantel.pq.GE<-multi.mantel(dist.pq.gentic, dist.clim.pq, nperm=999)
mantel.pq.GE[1:6]
mantel.pq.GES<-multi.mantel(dist.pq.gentic, dist.clim.selc.pq, nperm=999)
mantel.pq.GES[1:6]
mantel.pq.GG<-multi.mantel(dist.pq.gentic, as.matrix(pq.geo.dist), nperm=999)
mantel.pq.GG[1:6]
mantel.pq.GGE<-multi.mantel(dist.pq.gentic, list(dist.clim.pq,as.matrix(pq.geo.dist)), nperm=999)
mantel.pq.GGE[1:6]
mantel.pq.GGES<-multi.mantel(dist.pq.gentic, list(dist.clim.selc.pq,as.matrix(pq.geo.dist)), nperm=999)
mantel.pq.GGES[1:6]
# library(vegan)
mantel.pq.partial<-mantel.partial(dist.pq.gentic, dist.clim.selc.pq, pq.geo.dist, permutations = 999)
mantel.pq.partial

#多重线性回归
fit.ilex<-lm(axis1 ~ bio1 + bio2 + bio3 + bio4 + bio7 + bio8 + bio9 + bio14 + bio15 + bio19,data=clim.ilex)
summary(fit.ilex)
# library(cat)
# cat(print(summary(fit.ilex)),"fit.ilex.txt")
# confint(fit.ilex)
# 结果检测
pdf("fit.ilex.SELECTBIO.pdf")
par(mfrow=c(2 ,2))
plot(fit.ilex)
dev.off()
clim.ilex[,2]<-as.numeric(clim.ilex[,2])



row.names(gps.ilex)<-gps.ilex$tip.label

geodis.ilex<-distHaversine(c(110,33), c(100,26), r=6378137)

geodis.ilex<- apply(gps.ilex[,2:3],1,distHaversine)
  
distHaversine(gps.ilex$longitude, gps.ilex$latitude, r=6378137)

geodis.ilex <-data.frame()
for (i in 1:length(gps.ilex)) {geodis.ilex[i]<-distHaversine(gps.ilex[i,2:3], r=6378137)}


haversine.matrix <- function(x, lat.long.labels = c('Lat', 'Long'), name.column = NA, ...) {
  if(is.na(name.column[1])) nameVector <- row.names(x)
  else nameVector <- x[, name.column]
  x <- apply(x[, lat.long.labels], 1:2, as.numeric)
  out <- matrix(NA, nrow = dim(x)[1], ncol = dim(x)[1])
  for(i in 1:dim(x)[1]) {
    for(j in 1:i) {
      out[i, j] <- haversine.default(lat1 = x[i, 1], long1 = x[i, 2], lat2 = x[j, 1], long2 = x[j, 2])
    }}
  dimnames(out) <- list(nameVector, nameVector)
  out <- as.dist(out, ...)
  return(out)
}
geodis.ilex<-haversine.matrix(gps.ilex[,2:3], lat.long.labels = c('Lat', 'Long')) 


haversine <- function(x, lat.long.labels = c('Lat', 'Long'), name.column = NA, ...) UseMethod('haversine')

haversine.default <-
  ## Arguments:
  ##  lat1 and long1: latitude and longitude for site one, in decimal format (e.g., N40deg 30', W 90deg 45' = 40.500, -90.750)
  ##  lat2 and long2: latitude and longitude for site two
  ##  r = radius of the earth in the units of interest. The default value is in Kilometers. This particular version of the formula assumes a spherical earth
  ## Andrew Hipp (ahipp@mortonarb.org), January 2008
  function(lat1, long1, lat2, long2, r = 6372.795) {
    lat1 = lat1 * pi / 180
    long1 = long1 * pi / 180
    lat2 = lat2 * pi / 180
    long2 = long2 * pi / 180
    deltaLong = long2 - long1
    deltaLat = lat2 - lat1
    a = sin(deltaLat/2)^2 + cos(lat1)*cos(lat2)*sin(deltaLong/2)^2
    if(length(dim(a) == 2)) a <- apply(a, c(1,2), min, 1)
    else a <- min(a,1)
    c = 2 * asin(sqrt(a))
    d = r * c
    return(d) }

haversine.data.frame <- function(x, lat.long.labels = c('Lat', 'Long'), name.column = NA, ...) haversine(as.matrix(x), lat.long.labels, name.column, ...)

haversine.matrix <- function(x, lat.long.labels = c('Lat', 'Long'), name.column = NA, ...) {
  if(is.na(name.column[1])) nameVector <- row.names(x)
  else nameVector <- x[, name.column]
  x <- apply(x[, lat.long.labels], 1:2, as.numeric)
  out <- matrix(NA, nrow = dim(x)[1], ncol = dim(x)[1])
  for(i in 1:dim(x)[1]) {
    for(j in 1:i) {
      out[i, j] <- haversine.default(lat1 = x[i, 1], long1 = x[i, 2], lat2 = x[j, 1], long2 = x[j, 2])
    }}
  dimnames(out) <- list(nameVector, nameVector)
  out <- as.dist(out, ...)
  return(out)
}

pq<-read.csv("111.csv",header=F)
fastaFile <- readDNAStringSet("E:/research/analysis/cpdna biogeography quercus/quercus/4ampy.quercus.protobalanus.outgroup.edit.fas")
seq<- data.frame(name = names(fastaFile), sequence = paste(fastaFile))

simple.seq<-data.frame(name=seq$name[match(pq$V1,seq$name)],
                       sequence=seq$sequence[match(pq$V1,seq$name)])

hap.seq<-paste(">",simple.seq$name,sep="")
hap.seq<-paste(hap.seq,simple.seq$sequence,sep="###")
write.csv(hap.seq,"E:/research/analysis/cpdna biogeography quercus/quercus/pq.csv")
