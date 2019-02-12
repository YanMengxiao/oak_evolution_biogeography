library(ape)
library(phytools)
library(vegan)

setwd("E:/research/analysis/cpdna biogeography quercus/climate/mydata")

## cyc

clim.cyc<-read.csv("cyc_res2.5.csv",as.is=T)
rownames(clim.cyc)<-clim.cyc$tip.label
clim.cyc<-clim.cyc[green.cyc,]
# clim.cyc<-clim.cyc[-which(is.na(clim.cyc$bio1)),]

dist.clim.cyc<-vegdist(clim.cyc[,5:23],method="gower")
select.cyc<-paste("bio",c(1:2,5,7:8,12:15,18:19),sep="")
dist.clim.selc.cyc<-vegdist(clim.cyc[,select.cyc],method="gower",upper = T,diag=T)

dist.cyc.gentic<-read.csv("cyc.EA.genetic.csv",head=F)
dimnames(dist.cyc.gentic)<-list(1:132,1:132)
cyc.geo.dist<-read.csv("cyc.EA.geo.csv", head=F)
dimnames(cyc.geo.dist)<-list(1:132,1:132)

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
mantel.cyc.partial<-mantel.partial(dist.cyc.gentic, cyc.geo.dist,dist.clim.selc.cyc, permutations = 999)
mantel.cyc.partial

cyc.dist.clim<-data.frame(clim=as.vector(dist.clim.cyc),
                          selc=as.vector(dist.clim.selc.cyc))
write.csv(cyc.dist.clim,"cyc.EA.dist.csv")
cyc.all<-read.csv("cyc.EA.dist.csv")
plot(density(cyc.all$genetic))
plot(density(cyc.all$clim,na.rm=T))
glm.cyc.GE<-glm(genetic~clim,quasi(link = "identity", variance = "constant"),data=cyc.all)
summary(glm.cyc.GE)
glm.cyc.GG<-glm(genetic~geo,quasi(link = "identity", variance = "constant"),data=cyc.all)
summary(glm.cyc.GG)
glm.cyc.GGE<-glm(genetic~clim + geo,gaussian(link = "identity"),data=cyc.all)
summary(glm.cyc.GGE)
glm.cyc.GGES<-glm(genetic~selc + geo,gaussian(link = "identity"),data=cyc.all)
summary(glm.cyc.GGES)

## quercus
clim.que<-read.csv("pq_res2.5.csv",as.is=T)
clim.que<-clim.que[-grep("chrysolepis|tomentella|vacciniifolia|greggii",clim.que$tip.label),]
rownames(clim.que)<-clim.que$tip.label
clim.que<-clim.que[usa.pq,]
# clim.que<-clim.que[-which(is.na(clim.que$bio1)),]

dist.clim.que<-vegdist(clim.que[,5:23],method="gower")
select.que<-paste("bio",c(1:2,5,7:8,12:15,18:19),sep="")
dist.clim.selc.que<-vegdist(clim.que[,select.que],method="gower",upper = T,diag=T)

dist.que.gentic<-read.csv("que.NA.genetic.csv",head=F)
dimnames(dist.que.gentic)<-list(1:37,1:37)
que.geo.dist<-read.csv("que.NA.geo.csv", head=F)
dimnames(que.geo.dist)<-list(1:37,1:37)

mantel.que.GE<-multi.mantel(dist.que.gentic, dist.clim.que, nperm=999)
mantel.que.GE[1:6]
mantel.que.GES<-multi.mantel(dist.que.gentic, dist.clim.selc.que, nperm=999)
mantel.que.GES[1:6]
mantel.que.GG<-multi.mantel(dist.que.gentic, as.matrix(que.geo.dist), nperm=999)
mantel.que.GG[1:6]
mantel.que.GGE<-multi.mantel(dist.que.gentic, list(dist.clim.que,as.matrix(que.geo.dist)), nperm=999)
mantel.que.GGE[1:6]
mantel.que.GGES<-multi.mantel(dist.que.gentic, list(dist.clim.selc.que,as.matrix(que.geo.dist)), nperm=999)
mantel.que.GGES[1:6]
mantel.que.partial<-mantel.partial(dist.que.gentic, que.geo.dist,dist.clim.selc.que, permutations = 999)
mantel.que.partial

que.dist.clim<-data.frame(clim=as.vector(dist.clim.que),
                          selc=as.vector(dist.clim.selc.que))
write.csv(que.dist.clim,"que.NA.dist.csv")
que.all<-read.csv("que.NA.dist.csv")
plot(density(que.all$genetic))
plot(density(que.all$clim,na.rm=T))
glm.que.GE<-glm(genetic~clim,quasi(link = "identity", variance = "constant"),data=que.all)
summary(glm.que.GE)
glm.que.GG<-glm(genetic~geo,quasi(link = "identity", variance = "constant"),data=que.all)
summary(glm.que.GG)
glm.que.GGE<-glm(genetic~clim + geo,gaussian(link = "identity"),data=que.all)
summary(glm.que.GGE)
glm.que.GGES<-glm(genetic~selc + geo,gaussian(link = "identity"),data=que.all)
summary(glm.que.GGES)







clim.ilex<-read.csv("ilex_res2.5.csv",as.is=T)
clim.ilex<-clim.ilex[-grep("DM6572|DM6938|DM3688|DM8227|DM8229",clim.ilex$sequence),]
rownames(clim.ilex)<-clim.ilex$sequence
clim.ilex<-clim.ilex[blue.ilex,]
# clim.ilex<-clim.ilex[-which(is.na(clim.ilex$bio1)),]

dist.clim.ilex<-vegdist(clim.ilex[,6:24],method="gower")
select.ilex<-paste("bio",c(1:4,7:9,14:15,19),sep="")
dist.clim.selc.ilex<-vegdist(clim.ilex[,select.ilex],method="gower",upper = T,diag=T)

dist.ilex.gentic<-read.csv("ilex.SWC.genetic.csv",head=F)
dimnames(dist.ilex.gentic)<-list(1:69,1:69)
ilex.geo.dist<-read.csv("ilex.SWC.geo.csv", head=F)
dimnames(ilex.geo.dist)<-list(1:69,1:69)

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
mantel.ilex.partial<-mantel.partial(dist.ilex.gentic, ilex.geo.dist,dist.clim.selc.ilex, permutations = 999)
mantel.ilex.partial

ilex.dist.clim<-data.frame(clim=as.vector(dist.clim.ilex),
                           selc=as.vector(dist.clim.selc.ilex))
write.csv(ilex.dist.clim,"ilex.SWC.dist.csv")
ilex.all<-read.csv("ilex.SWC.dist.csv")
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

