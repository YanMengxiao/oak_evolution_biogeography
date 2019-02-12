source("https://bioconductor.org/biocLite.R")
biocLite(c("ggtree"))

rm(list=ls())

library(ape)
library(phytools)
library(phangorn)

setwd("E:/research/analysis/cpdna biogeography quercus/phylo_to_map")


######################################################################################################
########################################## CYCLOBALANOPSIS PLOT ##############################################
######################################################################################################

##tree import & manipulate
tr.cyc <- read.nexus("4ampy_cyc_litho_edit.nex.con.tre")
tr.cyc <- ladderize(root(tr.cyc, grep('Litho', tr.cyc$tip.label, value = T)))
tr.cyc <-drop.tip(tr.cyc,grep('Litho', tr.cyc$tip.label, value = T))
# tr.cyc <-drop.tip(tr.cyc,which(is.na(info$longitude[match(tr.cyc$tip.label,info$sequence)])))

##geographic data
info<-read.csv("E:/research/analysis/cpdna biogeography quercus/combination/add.outgroup/info_edit.csv",as.is=T)
gps.cyc <- info[match(tr.cyc$tip.label,info$sequence),16:17]
row.names(gps.cyc)<- info$sequence[match(tr.cyc$tip.label,info$sequence)]
# gps.cyc$tip.label<- info$sequence[match(tr.cyc$tip.label,info$sequence)]
# write.csv(gps.cyc,"gps.cyc.sequence.csv")
colnames(gps.cyc) <- c("latitude","longitude")
gps.cyc$longitude <- as.numeric(gps.cyc$longitude)
#tr.cyc$tip.label <- info$tip.label[match(tr.cyc$tip.label,info$sequence)]
#gps.cyc <- info[grep("Cyclobalanopsis",info$section),16:17]
#gps.cyc$name <- info[grep("Cyclobalanopsis",info$section),13]
#row.names(gps.cyc) <- info[grep("Cyclobalanopsis",info$section),13]
#info$sequence<- gsub("(Cyc_.*[0-9])_[0-9]","\\1",info$sequence,perl=TRUE)

## plot phyto to map
# 获取所有点的经纬度范围
# lat.range.cyc<-min(gps.cyc$latitude,na.rm=TRUE)
# lat.range.cyc<-c(16,53)
# log.range.cyc<-min(gps.cyc$longitude,na.rm=TRUE)
# log.range.cyc<-c(log.range.cyc,max(gps.cyc$longitude,na.rm=TRUE))

tr.cyc2 <- multi2di(tr.cyc)
#获得所有投点的点及线的，组成2列矩阵
green.cyc<-tr.cyc2$tip.label[Descendants(tr.cyc2,153,"tips")[[1]]] #Descendants,library(phangorn)
blue.cyc<-setdiff(tr.cyc2$tip.label,green.cyc)
colors.cyc<-matrix(NA,nrow(gps.cyc),2,dimnames=list(rownames(gps.cyc)))
for(i in 1:length(green.cyc))
  colors.cyc[green.cyc[i],1:2]<-c("palegreen3","palegreen3")
for(i in 1:length(blue.cyc))
  colors.cyc[blue.cyc[i],1:2]<-c("steelblue2","steelblue2")

pdf("cyc.phylotomap.node.pdf")
map.cyc <- phylo.to.map(tr.cyc2, gps.cyc, plot=F)
plot(map.cyc,ftype="off",lwd=1,colors=colors.cyc,ylim=c(16,53),adj=c(0.5,0.5), xlim=c(50,150),lty="dashed",ftype="off",fsize=0.3,psize=1)
nodelabels(round(as.numeric(tr.cyc$node.label),2),adj=c(x=0,y=1),frame='n',cex=0.8)
dev.off()

######################################################################################################
########################################## ILEX PLOT ##############################################
######################################################################################################
tr.ilex <- read.nexus("4ampy_ilex_edit.nex.con.tre")
tr.ilex <- root(tr.ilex, grep('Litho', tr.ilex$tip.label, value = T))
tr.ilex <-drop.tip(tr.ilex,grep('Litho', tr.ilex$tip.label, value = T))
tr.ilex <-drop.tip(tr.ilex,grep('DM6572|DM6938|DM3688', tr.ilex$tip.label, value = T))
# ilex DM6572| coccocifera DM6938 经BLAST 可能序列有误；baronii DM3688 测序可能有问题
# tr.ilex <-drop.tip(tr.ilex,which(!info$gps[match(tr.ilex$tip.label,info$sequence)])) 去除原始gps错误者，现已更正
# tr.ilex <-drop.tip(tr.ilex,which(is.na(info$longitude[match(tr.ilex$tip.label,info$sequence)])))
##ladderize/rotate tree都没用，plot出来都是反的
# tr2<-tr.ilex
# for(i in length(tr2$tip)+1:tr2$Nnode) tr2<-rotate(tr2,i)
# plot(tr2,cex=0.6)
# tr.ilex<-tr2
##geographic data
info<-read.csv("E:/research/analysis/cpdna biogeography quercus/combination/add.outgroup/info_edit.csv",as.is=T)
gps.ilex <- info[match(tr.ilex$tip.label,info$sequence),c(17,16)]
row.names(gps.ilex)<- info$sequence[match(tr.ilex$tip.label,info$sequence)]
# gps.ilex$tip.label<- info$sequence[match(tr.ilex$tip.label,info$sequence)]
# write.csv(gps.ilex,"gps.ilex.sequence.csv")
# colnames(gps.ilex) <- c("latitude","longitude")
#tr.ilex$tip.label <- info$tip.label[match(tr.ilex$tip.label,info$sequence)]

## plot phyto to map
# 获取所有点的经纬度范围
# lat.range.ilex<-min(gps.ilex$latitude,na.rm=TRUE)
# lat.range.ilex<-c(lat.range.ilex,max(gps.ilex$latitude,na.rm=TRUE))
# log.range.ilex<-min(gps.ilex$longitude,na.rm=TRUE)
# log.range.ilex<-c(log.range.ilex,max(gps.ilex$longitude,na.rm=TRUE))

# 多歧分枝转化为二歧分枝
tr.ilex2 <- multi2di(tr.ilex)
plot(tr.ilex2,cex=0.6)
nodelabels(cex=0.8)
# #获得所有投点的点及线的，组成2列矩阵
# phylo.to.map(tr.ilex2,gps.ilex,ftype="off")
# plot(map.ilex,ftype="off",lwd=1,colors=colors,ylim=c(16,53), xlim=c(0,130),lty="dashed",ftype="off",fsize=0.3,psize=2)
# plot(tr.ilex2,cex=0.6)

green.ilex<-tr.ilex2$tip.label[Descendants(tr.ilex2,212,"tips")[[1]]] #Descendants,library(phangorn)
blue.ilex<-tr.ilex2$tip.label[Descendants(tr.ilex2,144,"tips")[[1]]]
purple.ilex<-tr.ilex2$tip.label[Descendants(tr.ilex2,125,"tips")[[1]]]
colors.ilex<-matrix(NA,nrow(gps.ilex),2,dimnames=list(rownames(gps.ilex)))
for(i in 1:length(green.ilex))
  colors.ilex[green.ilex[i],1:2]<-c("palegreen3","palegreen3")
for(i in 1:length(blue.ilex))
  colors.ilex[blue.ilex[i],1:2]<-c("steelblue2","steelblue2")
for(i in 1:length(purple.ilex))
  colors.ilex[purple.ilex[i],1:2]<-c("mediumorchid2","mediumorchid2")


# green.ilex<-tr.ilex2$tip.label[Descendants(tr.ilex2,129,"tips")[[1]]] #Descendants,library(phangorn)
# blue.ilex<-tr.ilex2$tip.label[Descendants(tr.ilex2,164,"tips")[[1]]]
# purple.ilex<-tr.ilex2$tip.label[Descendants(tr.ilex2,233,"tips")[[1]]]
# colors.ilex<-matrix(NA,nrow(gps.ilex),2,dimnames=list(rownames(gps.ilex)))
# for(i in 1:length(green.ilex))
#   colors.ilex[green.ilex[i],1:2]<-c("palegreen3","palegreen3")
# for(i in 1:length(blue.ilex))
#   colors.ilex[blue.ilex[i],1:2]<-c("steelblue2","steelblue2")
# for(i in 1:length(purple.ilex))
#   colors.ilex[purple.ilex[i],1:2]<-c("mediumorchid2","mediumorchid2")

pdf("ilex.phylotomap.pdf")
map.ilex <- phylo.to.map(tr.ilex2, gps.ilex, plot=F,rotate=F)
plot(map.ilex,ftype="off",lwd=1,colors=colors.ilex,ylim=c(16,53), xlim=c(-10,140),lty="dashed",ftype="off",fsize=0.3,psize=1)
# nodelabels(round(as.numeric(tr.ilex2$node.label,2)),adj=c(x=0,y=1),frame="n",bg="white")
dev.off()


######################################################################################################
########################################## QUERCUS PLOT ##############################################
######################################################################################################
tr.pq <- read.nexus("4ampy.quercus.protobalanus.outgroup.edit.nex.con.tre")
tr.pq <- ladderize(root(tr.pq, grep('Litho', tr.pq$tip.label, value = T)))
tr.pq <-drop.tip(tr.pq,grep('Litho', tr.pq$tip.label, value = T))
tr.pq <-drop.tip(tr.pq,grep('greggii', tr.pq$tip.label, value = T))
tr.pq <-drop.tip(tr.pq,grep('chrysolepis|tomentella|vacciniifolia', tr.pq$tip.label, value = T))
# tr.pq <-drop.tip(tr.pq,which(is.na(info$longitude[match(tr.pq$tip.label,info$sequence)])))
## plot phyto to map
tr.pq2 <- multi2di(tr.pq)
plot(tr.pq2,cex=0.8)
nodelabels(cex=0.8)
## drop protobalanus
usa.pq<-tr.pq2$tip.label[Descendants(tr.pq2,152,"tips")[[1]]] #Descendants,library(phangorn)
purple.pq<-tr.pq2$tip.label[Descendants(tr.pq2,97,"tips")[[1]]]
colors.pq<-matrix(NA,nrow(gps.pq),2,dimnames=list(rownames(gps.pq)))
for(i in 1:length(usa.pq))
  colors.pq[usa.pq[i],1:2]<-c("orange","orange")
# for(i in 1:length(grey.pq))
#   colors.pq[grey.pq[i],1:2]<-c("black","black")
for(i in 1:length(purple.pq))
  colors.pq[purple.pq[i],1:2]<-c("mediumorchid2","mediumorchid2")

# ## 保留 protobalanus
# usa.pq<-tr.pq2$tip.label[Descendants(tr.pq2,157,"tips")[[1]]] #Descendants,library(phangorn)
# purple.pq<-tr.pq2$tip.label[Descendants(tr.pq2,102,"tips")[[1]]]
# proto.pq<-info$sequence[grep("Protobalanus",info$section)]
# colors.pq<-matrix(NA,nrow(gps.pq),2,dimnames=list(rownames(gps.pq)))
# for(i in 1:length(usa.pq))
#   colors.pq[usa.pq[i],1:2]<-c("orange","orange")
# for(i in 1:length(purple.pq))
#   colors.pq[purple.pq[i],1:2]<-c("mediumorchid2","mediumorchid2")
# colors.pq[intersect(row.names(colors.pq),proto.pq),1:2]<-"firebrick2"
# # row.names(gps.pq)<- info$sequence[match(tr.pq$tip.label,info$sequence)]

##geographic data
#info<-read.csv("E:/research/analysis/cpdna biogeography quercus/combination/add.outgroup/info_edit.csv",as.is=T)
# gps.pq <- info[match(tr.pq2$tip.label,info$sequence),16:17]
# row.names(gps.pq)<- info$sequence[match(tr.pq2$tip.label,info$sequence)]
# # gps.pq$tip.label<- info$sequence[match(tr.pq$tip.label,info$sequence)]
# # write.csv(gps.pq,"gps.pq.sequence.csv")
# colnames(gps.pq) <- c("latitude","longitude")
# gps.pq$latitude <- as.numeric(gps.pq$latitude)
# #tr.pq$tip.label <- info$tip.label[match(tr.pq$tip.label,info$sequence)]
# # 获取所有点的经纬度范围
# # lat.range.pq<-min(gps.pq$latitude,na.rm=TRUE)
# # lat.range.pq<-c(lat.range.pq,max(gps.pq$latitude,na.rm=TRUE))
# # log.range.pq<-min(gps.pq$longitude,na.rm=TRUE)
# # log.range.pq<-c(log.range.pq,max(gps.pq$longitude,na.rm=TRUE))
# #获得所有投点的点及线的，组成2列矩阵

##geographic data
gps.pq<-data.frame(name=tr.pq2$tip.label)
gps.pq$latitude<-as.numeric(info$latitude[match(gps.pq$name,info$sequence)])
gps.pq$longitude<-info$longitude[match(gps.pq$name,info$sequence)]
rownames(gps.pq)<-gps.pq$name
##投影
pdf("pq.drop.proto.phylotomap.node.pdf")
map.pq <- phylo.to.map(tr.pq2, gps.pq[,2:3], plot=F)
plot(map.pq,ftype="off",lwd=1,colors=colors.pq,ylim=c(10,65), xlim=c(-130,140),lty="dashed",ftype="off",fsize=0.3,psize=1)
nodelabels(round(as.numeric(tr.pq$node.label,2)),adj=c(x=0,y=1),frame="n",bg="white")
dev.off()
##note
# plot(direction="rightwards"#树变垂直方向)

# 调整地图边界，去除小型岛屿
# obj<-phylo.to.map(tree,X,regions="Chile",direction="rightwards", plot=FALSE)
# ii<-which(map.cyc$map$x<(-77))
# obj$map$x<-obj$map$x[-ii]
# obj$map$y<-obj$map$y[-ii]
# obj$map$range[1]<--81
#tip名字查重
#tr.cyc$tip.label[which(duplicated(tr.cyc$tip.label))]
#taxa.cyc <- setdiff(row.names(gps.cyc), tr.cyc$tip.label)#gps.cyc$name[match(tr.cyc$tip.label,gps.cyc$name)]



#taxa <- intersect(row.names(gps.cyc), tr.cyc$tip.label) 
#dat <- dat[taxa, ]
#tr.cyc <- drop.tip(tr, which(!tr.cyc$tip.label %in% gps.cyc$name))


row.names(dat) <- gsub("[|-]", "_", row.names(dat)) # fixes mismatch between row names and tip labels
taxa <- intersect(row.names(dat), tr$tip.label) # only 92 names match... you may still need to clean the label names
dat <- dat[taxa, ]
tr <- drop.tip(tr, which(!tr$tip.label %in% taxa))
tr2 <- multi2di(tr)
phylo.to.map(tr2, dat[,2:3],ftype = 'off')



match(info$sequence,tr.cyc$tip.label)

tr.cyc$tip.label[99:103]


tr.que <- root(tr.que, grep('Litho', tr.que$tip.label, value = T))
tr.que <-ladderize(drop.tip(tr.que,grep('Litho', tr.que$tip.label, value = T)))




row.names(dat) <- gsub("[|-]", "_", row.names(dat)) # fixes mismatch between row names and tip labels
taxa <- intersect(row.names(dat), tr$tip.label) # only 92 names match... you may still need to clean the label names
dat <- dat[taxa, ]
tr <- drop.tip(tr, which(!tr$tip.label %in% taxa))
tr2 <- multi2di(tr)
phylo.to.map(tr2, dat[,2:3],ftype = 'off')

tr.que <- read.nexus("E:/research/analysis/cpdna biogeography quercus/quercus/4ampy.quercus.protobalanus.outgroup.nex.con.tre")
tr.que <- root(tr.que, grep('Litho', tr.que$tip.label, value = T))
#tr.que <-drop.tip(tr.que,grep('Litho', tr.que$tip.label, value = T))
dat <- read.csv('E:/research/analysis/cpdna biogeography quercus/quercus/haplotype/gps.que.csv', as.is = T, row.names=1)
row.names(dat) <- gsub("[|-]", "_", row.names(dat)) # fixes mismatch between row names and tip labels
taxa <- intersect(row.names(dat), tr$tip.label) # only 92 names match... you may still need to clean the label names
dat <- dat[taxa, ]
tr <- drop.tip(tr, which(!tr$tip.label %in% taxa))
tr2 <- multi2di(tr)
phylo.to.map(tr2, dat[,2:3],ftype = 'off')

tr <- read.nexus("E:/research/analysis/cpdna biogeography quercus/quercus/4ampy.quercus.protobalanus.outgroup.nex.con.tre")
tr <- root(tr, grep('Litho', tr$tip.label, value = T))
#tr <-drop.tip(tr,grep('Litho', tr$tip.label, value = T))
dat <- read.csv('E:/research/analysis/cpdna biogeography quercus/quercus/haplotype/gps.que.csv', as.is = T, row.names=1)
row.names(dat) <- gsub("[|-]", "_", row.names(dat)) # fixes mismatch between row names and tip labels
taxa <- intersect(row.names(dat), tr$tip.label) 
dat <- dat[taxa, ]
tr <- drop.tip(tr, which(!tr$tip.label %in% taxa))
tr2 <- multi2di(tr)
phylo.to.map(tr2, dat[1:5,2:3],ftype = 'off')

info<-read.csv("E:/research/analysis/cpdna biogeography quercus/combination/add.outgroup/info_edit.csv",as.is=T)





plot(ladderize(tree.que.root),cex=0.8)
tree.que.root <- multi2di(tree.que.root)
tree.que.root$tip.label<-1:99
gps.que<-read.csv("E:/research/analysis/cpdna biogeography quercus/quercus/gps.que2.csv")
rownames(gps.que)<-1:99

phylo.to.map(tree.que.root,gps.que[,2:3],plot=T)
plot(xx,type="phylogram")

ilex<-read.tree("E:/research/analysis/cpdna biogeography quercus/quercus/4ampy.quercus.protobalanus.outgroup.nwk")
ilex$tip.label<-1:99
phylo.to.map(ilex,gps.que[1:99,2:3])

l<-subtrees(tree.que)

### plot all the subtrees
for (i in 1:11) plot(l[[i]], sub=paste("Node", l[[i]]$node.label[1]))
title(main="",adj=0)
par(mfrow=c(1,1))
# }
sub5<-l[[7]]
rm(l)
tree<-rtree(n=99)
tree$tip.label<-1:99
tree.map<-phylo.to.map(tree,gps.que[1:99,2:3])
plot(tree.map,fsize=0.3)

plot(tree)

tree$tip.label<-1:10
tree.map<-phylo.to.map(tree,gps.que[1:10,2:3])
plot(tree.map,fsize=1)


data(bird.families)
is.binary(bird.families)
is.binary(multi2di(bird.families))
all.equal(di2multi(multi2di(bird.families)), bird.families)
### To see the results of randomly resolving a trichotomy:
tr <- read.tree(text = "(a:1,b:1,c:1);")
layout(matrix(1:4, 2, 2))
for (i in 1:4)
  plot(multi2di(tr), use.edge.length = FALSE, cex = 1.5)
layout(1)



##GPS数据格式修正
gps.que<-read.csv("E:/research/analysis/cpdna biogeography quercus/quercus/Pham KK et al.2017-gps.csv")
colnames(gps.que)<-c("sample","latitude","longitude")
#latitude
gps.que$latitude<-gsub("unknown", " N W",gps.que$latitude)
gps.que$latitude<-gsub(" N.*W", "",gps.que$latitude)
gps.que$latitude<-gsub(" N.*E", "",gps.que$latitude)
#longitude
gps.que$longitude<-gsub("unknown", " N W",gps.que$longitude)
gps.que$longitude<-gsub(".*N", "",gps.que$longitude)
gps.que$longitude<-gsub(",", "",gps.que$longitude)
write.csv(gps.que,"E:/research/analysis/cpdna biogeography quercus/quercus/gps.que.csv")##excel中去除西经的负号
gps.que<-read.csv("E:/research/analysis/cpdna biogeography quercus/quercus/gps.que.csv")
gps.que$longitude<-gsub("([0-9].*) W", "-\\1",gps.que$longitude)
##利用替换将西经改为负数，gsub匹配多次，sub只匹配一次
##()内模式编组，\\1代表()内模式
gps.que$longitude<-gsub(" W", "",gps.que$longitude)
gps.que$longitude<-gsub(" E", "",gps.que$longitude)
write.csv(gps.que,"E:/research/analysis/cpdna biogeography quercus/quercus/gps.que.csv")
gps.que<-read.csv("E:/research/analysis/cpdna biogeography quercus/quercus/gps.que.csv")
gps.que.gps<-gps.que[,3:4]
row.names(gps.que.gps)<-gps.que$sample
tree.que$tip.label<-1:99
row.names(gps.que)<-1:99
old<-phylo.to.map(tree.que,gps.que[,2:3],plot=T,rotate=F,fsize=0.3)
plot(old,fsize=0.3)

gps.que.gps<-as.matrix(gps.que.gps)
as.numeric(gps.que.gps[,1:2])

tree<-pbtree(n=99,scale=100)
tree$tip.label <- 1:99
lat<-fastBM(tree,sig2=10,bounds=c(-90,90))
long<-fastBM(tree,sig2=80,bounds=c(-180,180))
phylo.to.map(tree,cbind(lat,long),plot=T)
layout(matrix(c(1,1),1,1))


tree<-pbtree(n=26,scale=100)
tree$tip.label<-LETTERS[26:1]
lat<-fastBM(tree.que,sig2=10,bounds=c(-90,90))
long<-fastBM(tree.que,sig2=80,bounds=c(-180,180))
# now plot projection
xx<-phylo.to.map(tree.que,cbind(lat,long),plot=T,rotate=F)
layout(matrix(c(1,2),2,1),heights=c(0.61,0.39))
plot(xx,type="phylogram",asp=1.3,mar=c(0.1,0.5,3.1,0.1))
title(main="(a)",adj=0)
plot(xx,type="direct",asp=1.3,mar=c(0.1,0.5,3.1,0.1))
title(main="(b)",adj=0)






# now plot projection
phylo.to.map(tree.que,cbind(lat,long),plot=T,rotate=F)
plot(xx,type="phylogram")
cw$tip.label
lat<-fastBM(tree.que3,sig2=10,bounds=c(-90,90))
long<-fastBM(tree.que3,sig2=80,bounds=c(-180,180))
# now plot projection
xx<-phylo.to.map(tree.que3,cbind(lat,long),tip.labels=FALSE, plot=TRUE)
#tiplabels(tree.que3$tip.label,fsize=0.0)
plot(xx,type="phylogram",show.tip.label = F)
ggtree(xx)
tree.que3 <- drop.tip(tree.que3,"99")
PLOT
