library(phytools)
library(ape)
library(vegan)
library(phangorn)
setwd("E:/research/analysis/cpdna biogeography quercus/climate/mydata")
info<-read.csv("E:/research/analysis/cpdna biogeography quercus/combination/add.outgroup/info_edit.csv",as.is=T)


############################
###        ILEX          ###
############################

ilex.clim<-read.csv("ilex_res2.5.csv")
ilex.clim<-ilex.clim[-grep("DM6572|DM6938|DM3688",ilex.clim$sequence),]
rownames(ilex.clim)<- ilex.clim$sequence
# select.ilex<-paste("bio",c(1:4,7:9,14:15,19),sep="_")
# ilex.clim.selc<-ilex.clim[,select.ilex]
#NMDS k=2 计算
mds.ilex<-metaMDS(ilex.clim[,grep("bio",colnames(ilex.clim))], 'gower', k=2)
mds.ilex.all<-metaMDS(ilex.clim[,4:24], 'gower', k=2)

# mds.ilex.selc<-metaMDS(ilex.clim.selc, 'gower', k=2)
#检查拟合度
stressplot(mds.ilex.geo,main="ilex stressplot")
#NMDS k=1:10 计算
mds.ilex.multi <- lapply(1:10, function(x) {metaMDS(ilex.clim[,grep("bio",colnames(ilex.clim))], 'gower', k=x)})
# mds.ilex.multi.selc <- lapply(1:10, function(x) {metaMDS(ilex.clim.selc, 'gower', k=x)})
# k=1:10 stress变化趋势
stress.ilex<-c()
for (i in 1:10) {stress.ilex<-c(stress.ilex,mds.ilex.multi[[i]]$stress)}
plot(stress.ilex, type = 'l')
# 计算各因子贡献率，选取投影的因子
factor.ilex<-envfit(mds.ilex.all, ilex.clim[,grep("bio",colnames(ilex.clim))])
factor.ilex.all<-envfit(mds.ilex.all, ilex.clim[,4:24])

# 系统树预处理
tr.ilex <- read.nexus("E:/research/analysis/cpdna biogeography quercus/phylo_to_map/4ampy_ilex_edit.nex.con.tre")
tr.ilex <- root(tr.ilex, grep('Litho', tr.ilex$tip.label, value = T))
tr.ilex <-drop.tip(tr.ilex,grep('Litho', tr.ilex$tip.label, value = T))
tr.ilex <-drop.tip(tr.ilex,grep('DM3688', tr.ilex$tip.label, value = T))
tr.ilex <-drop.tip(tr.ilex,setdiff(tr.ilex$tip.label,tr.ilex$tip.label[match(ilex.clim$sequence,tr.ilex$tip.label)]))
tr.ilex <-drop.tip(tr.ilex,grep('DM6572|DM6938|DM3688', tr.ilex$tip.label, value = T))


# tr.ilex <-drop.tip(tr.ilex,which(is.na(info$longitude[match(tr.ilex$tip.label,info$sequence)])))
#tr.ilex <-drop.tip(tr.ilex,which(!info$gps[match(tr.ilex$tip.label,info$sequence)]))
# 规定各分支的tip，并指定颜色，注意tip顺序需同NMDS结果
plot(tr.ilex,cex=0.5)
nodelabels(cex=0.8)
green.ilex<-tr.ilex$tip.label[Descendants(tr.ilex,170,"tips")[[1]]] #Descendants,library(phangorn)
blue.ilex<-tr.ilex$tip.label[Descendants(tr.ilex,138,"tips")[[1]]]
purple.ilex<-tr.ilex$tip.label[Descendants(tr.ilex,123,"tips")[[1]]]
color.bg.ilex<-c()
color.bg.ilex[green.ilex]<-rep("palegreen3",length(green.ilex))
color.bg.ilex[blue.ilex]<-rep("steelblue2",length(blue.ilex))
color.bg.ilex[purple.ilex]<-rep("mediumorchid2",length(purple.ilex))
color.bg.ilex<-color.bg.ilex[match(row.names(mds.ilex$point),names(color.bg.ilex))]

# ilex.pop<-read.csv("E:/research/analysis/cpdna biogeography quercus/climate/mydata/genalex/ilex.pop.csv")
# ilex.pop$Pop<-color.bg.ilex[match(ilex.pop$Sample,names(color.bg.ilex))]
# write.csv(ilex.pop,"ilex.pop2.csv")

## bio 1:19
pdf("ilex.NMDS.gower.bio15.pdf")
# 环境变量投影
ordisurf(mds.ilex, ilex.clim[,'bio15'],
         xlim = range(mds.ilex$points[, 1]),
         lwd.cl = 0.5, col = 'gray80', main = '', cex = 0, add = FALSE, axes = F)
# 系统树投影在二维坐标
phylomorphospace(tr.ilex, mds.ilex$points[,1:2],
                 lwd = 0.75, xlab = '', ylab = '', axes = F, node.size=0,
                 xlim = range(mds.ilex$points[, 1]),
                 ylim = range(c(-10, mds.ilex$points[, 2])), add = TRUE,label="off")
box()
axis(side = 1, cex.axis = 0.5)
axis(side = 2, cex.axis = 0.5)
points(mds.ilex$points[,1:2], pch = 21, col = 'black', bg = color.bg.ilex, cex = 2)
text(-10.5,10, substitute(paste("section ", italic('Ilex'))), cex = 0.6, pos = 4)
legend("bottomleft", c("East Asia", "SW China","Himilaya-Mediterranean"),
       pch = 21, col = 'black',
       pt.bg = c("palegreen3","steelblue2","mediumorchid2"),
       bty = 'n', cex = 1)
title(ylab="MDS axis 2", outer=T, line=2)
title(xlab="MDS axis 1", outer=T, line=1)
dev.off()

## geo+bio
pdf("ilex.NMDS.gower.geobio.bio15.pdf")
# 环境变量投影
ordisurf(mds.ilex.all, ilex.clim[,'bio15'],
         xlim = range(mds.ilex.all$points[, 1]),
         lwd.cl = 0.5, col = 'gray80', main = '', cex = 0, add = FALSE, axes = F)
# 系统树投影在二维坐标
phylomorphospace(tr.ilex, mds.ilex.all$points[,1:2],
                 lwd = 0.75, xlab = '', ylab = '', axes = F, node.size=0,
                 xlim = range(mds.ilex.all$points[, 1]),
                 ylim = range(c(-10, mds.ilex.all$points[, 2])), add = TRUE,label="off")
box()
axis(side = 1, cex.axis = 0.5)
axis(side = 2, cex.axis = 0.5)
points(mds.ilex.all$points[,1:2], pch = 21, col = 'black', bg = color.bg.ilex, cex = 2)
text(-10.5,10, substitute(paste("section ", italic('Ilex'))), cex = 0.6, pos = 4)
legend("bottomleft", c("East Asia", "SW China","Himilaya-Mediterranean"),
       pch = 21, col = 'black',
       pt.bg = c("palegreen3","steelblue2","mediumorchid2"),
       bty = 'n', cex = 1)
title(ylab="MDS axis 2", outer=T, line=2)
title(xlab="MDS axis 1", outer=T, line=1)
dev.off()

############################
###        CYC          ###
############################

cyc.clim<-read.csv("cyc_res2.5.csv")
rownames(cyc.clim)<- cyc.clim$tip.label
mds.cyc<-metaMDS(cyc.clim[,grep("bio",colnames(cyc.clim))], 'gower', k=2)
mds.cyc.all<-metaMDS(cyc.clim[,3:23], 'gower', k=2)

mds.cyc.multi <- lapply(1:10, function(x) {metaMDS(cyc.clim[,grep("bio",colnames(cyc.clim))], 'gower', k=x)})
stressplot(mds.cyc,main="Cyc stressplot")
stress.cyc<-c()
for (i in 1:10) {stress.cyc<-c(stress.cyc,mds.cyc.multi[[i]]$stress)}
plot(stress.cyc, type = 'l')
factor.cyc.bio<-envfit(mds.cyc.all, cyc.clim[,grep("bio",colnames(cyc.clim))])
factor.cyc.gb<-envfit(mds.cyc.all, cyc.clim[,3:23])

tr.cyc <- read.nexus("E:/research/analysis/cpdna biogeography quercus/phylo_to_map/4ampy_cyc_litho_edit.nex.con.tre")
tr.cyc <- root(tr.cyc, grep('Litho', tr.cyc$tip.label, value = T))
tr.cyc <-drop.tip(tr.cyc,grep('Litho', tr.cyc$tip.label, value = T))
tr.cyc <-drop.tip(tr.cyc,which(is.na(info$longitude[match(tr.cyc$tip.label,info$sequence)])))
green.cyc<-tr.cyc$tip.label[Descendants(tr.cyc,160,"tips")[[1]]] #Descendants,library(phangorn)
blue.cyc<-setdiff(tr.cyc$tip.label,green.cyc)
color.bg.cyc<-c()
color.bg.cyc[green.cyc]<-rep("palegreen3",length(green.cyc))
color.bg.cyc[blue.cyc]<-rep("steelblue2",length(blue.cyc))
color.bg.cyc<-color.bg.cyc[match(row.names(mds.cyc$point),names(color.bg.cyc))]

## geo+bio
pdf("cyc.NMDS.gower.geo.bio.bio12.pdf")
ordisurf(mds.cyc.all, cyc.clim[, 'bio12'],
         xlim = range(mds.cyc.all$points[, 1]),
         lwd.cl = 0.5, col = 'gray80', main = '', cex = 0, add = FALSE, axes = F)
phylomorphospace(tr.cyc, mds.cyc.all$points[,1:2],
                 lwd = 0.75, xlab = '', ylab = '', axes = F, node.size=0,
                 xlim = range(mds.cyc.all$points[, 1]),
                 ylim = range(c(-10, mds.cyc.all$points[, 2])), add = TRUE,label="off")
box()
axis(side = 1, cex.axis = 0.5)
axis(side = 2, cex.axis = 0.5)
points(mds.cyc.all$points[,1:2], pch = 21, col = 'black', bg = color.bg.cyc, cex = 2)
text(-10.5,10, substitute(paste("section ", italic('Cyclobalanopsis'))), cex = 0.6, pos = 4)
legend("bottomleft", c("East Asia", "SW China and Vietnam"),
       pch = 21, col = 'black',
       pt.bg = c("palegreen3","steelblue2"),
       bty = 'n', cex = 1)
title(ylab="MDS axis 2", outer=T, line=2)
title(xlab="MDS axis 1", outer=T, line=1)
dev.off()

## bio 1:19  
pdf("cyc.NMDS.gower.bio12.pdf")
ordisurf(mds.cyc, cyc.clim[, 'bio12'],
         xlim = range(mds.cyc$points[, 1]),
         lwd.cl = 0.5, col = 'gray80', main = '', cex = 0, add = FALSE, axes = F)
phylomorphospace(tr.cyc, mds.cyc$points[,1:2],
                 lwd = 0.75, xlab = '', ylab = '', axes = F, node.size=0,
                 xlim = range(mds.cyc$points[, 1]),
                 ylim = range(c(-10, mds.cyc$points[, 2])), add = TRUE,label="off")
box()
axis(side = 1, cex.axis = 0.5)
axis(side = 2, cex.axis = 0.5)
points(mds.cyc$points[,1:2], pch = 21, col = 'black', bg = color.bg.cyc, cex = 2)
text(-10.5,10, substitute(paste("section ", italic('Cyclobalanopsis'))), cex = 0.6, pos = 4)
legend("bottomleft", c("East Asia", "SW China and Vietnam"),
       pch = 21, col = 'black',
       pt.bg = c("palegreen3","steelblue2"),
       bty = 'n', cex = 1)
title(ylab="MDS axis 2", outer=T, line=2)
title(xlab="MDS axis 1", outer=T, line=1)
dev.off()

############################
 ###        PQ          ###
############################

pq.clim<-read.csv("pq_res2.5.csv")
pq.clim<-pq.clim[-grep("chrysolepis|tomentella|vacciniifolia|greggii",pq.clim$tip.label),]
rownames(pq.clim)<- pq.clim$tip.label
mds.pq<-metaMDS(pq.clim[,grep("bio",colnames(pq.clim))], 'gower', k=2)
mds.pq.all<-metaMDS(pq.clim[,3:23], 'gower', k=2)

mds.pq.multi <- lapply(1:10, function(x) {metaMDS(pq.clim[,grep("bio",colnames(pq.clim))], 'gower', k=x)})
stressplot(mds.pq,main="P&Q stressplot")
stress.pq<-c()
for (i in 1:10) {stress.pq<-c(stress.pq,mds.pq.multi[[i]]$stress)}
plot(stress.pq, type = 'l')
factor.pq<-envfit(mds.pq.all, pq.clim[,grep("bio",colnames(pq.clim))])
factor.pq.all<-envfit(mds.pq.all, pq.clim[,3:23])

tr.pq <- read.nexus("E:/research/analysis/cpdna biogeography quercus/phylo_to_map/4ampy.quercus.protobalanus.outgroup.edit.nex.con.tre")
tr.pq <- root(tr.pq, grep('Litho', tr.pq$tip.label, value = T))
tr.pq <-drop.tip(tr.pq,grep('Litho', tr.pq$tip.label, value = T))
tr.pq <-drop.tip(tr.pq,grep('chrysolepis|tomentella|vacciniifolia|greggii', tr.pq$tip.label, value = T))
tr.pq <-drop.tip(tr.pq,which(is.na(info$longitude[match(tr.pq$tip.label,info$sequence)])))
# plot(tr.pq,cex=0.6)
# nodelabels(cex=0.8)
usa.pq<-tr.pq$tip.label[Descendants(tr.pq,91,"tips")[[1]]] #Descendants,library(phangorn)
purple.pq<-tr.pq$tip.label[Descendants(tr.pq,107,"tips")[[1]]]
# proto.pq<-info$sequence[grep("Protobalanus",info$section)]

color.bg.pq<-c()
color.bg.pq[purple.pq]<-rep("mediumorchid2",length(purple.pq))
color.bg.pq[usa.pq]<-rep("orange",length(usa.pq))
# color.bg.pq[proto.pq]<-rep("firebrick2",length(proto.pq))
color.bg.pq<-color.bg.pq[match(row.names(mds.pq.all$point),names(color.bg.pq))]

# names(color.bg.pq)<-gsub("___","_|_",names(color.bg.pq))
# pq.pop<-read.csv("E:/research/analysis/cpdna biogeography quercus/climate/mydata/genalex/pq.pop.csv")
# pq.pop$Pop<-color.bg.pq[match(pq.pop$Sample,names(color.bg.pq))]
# write.csv(pq.pop,"pq.pop2.csv")

## geo+bio
pdf("pq.NMDS.gower.geo.bio.bio11.pdf")
ordisurf(mds.pq.all, pq.clim[, 'bio11'],
         xlim = range(mds.pq.all$points[, 1]),
         lwd.cl = 0.5, col = 'gray80', main = '', cex = 0, add = FALSE, axes = F)
phylomorphospace(tr.pq, mds.pq.all$points[,1:2],
                 lwd = 0.75, xlab = '', ylab = '', axes = F,node.size=0,
                 xlim = range(mds.pq.all$points[, 1]),
                 ylim = range(c(-10, mds.pq.all$points[, 2])), add = TRUE,label="off")
box()
axis(side = 1, cex.axis = 0.5)
axis(side = 2, cex.axis = 0.5)
points(mds.pq.all$points[,1:2], pch = 21, col = 'black', bg = color.bg.pq, cex = 2)
text(-10.5,10, substitute(paste("section ", italic('Quercus'))), cex = 0.6, pos = 4)
legend("bottomleft", c("Eurasia", "America"),
       pch = 21, col = 'black',
       pt.bg = c("mediumorchid2","orange"),
       bty = 'n', cex = 1)
title(ylab="MDS axis 2", outer=T, line=2)
title(xlab="MDS axis 1", outer=T, line=1)
dev.off()

## bio 1;19
pdf("pq.NMDS.gower.new.bio11.pdf")
ordisurf(mds.pq, pq.clim[, 'bio11'],
         xlim = range(mds.pq$points[, 1]),
         lwd.cl = 0.5, col = 'gray80', main = '', cex = 0, add = FALSE, axes = F)
phylomorphospace(tr.pq, mds.pq$points[,1:2],
                 lwd = 0.75, xlab = '', ylab = '', axes = F,node.size=0,
                 xlim = range(mds.pq$points[, 1]),
                 ylim = range(c(-10, mds.pq$points[, 2])), add = TRUE,label="off")
box()
axis(side = 1, cex.axis = 0.5)
axis(side = 2, cex.axis = 0.5)
points(mds.pq$points[,1:2], pch = 21, col = 'black', bg = color.bg.pq, cex = 2)
text(-10.5,10, substitute(paste("section ", italic('Quercus'))), cex = 0.6, pos = 4)
legend("bottomleft", c("Eurasia", "America"),
       pch = 21, col = 'black',
       pt.bg = c("mediumorchid2","orange"),
       bty = 'n', cex = 1)
title(ylab="MDS axis 2", outer=T, line=2)
title(xlab="MDS axis 1", outer=T, line=1)
dev.off()


############################
###        ALL          ###
############################

clim<-read.csv("all_res2.5.csv")
clim<-clim[-grep("DM6572|DM6938|DM3688|chrysolepis|tomentella|vacciniifolia|greggii",clim$sequence),]
rownames(clim)<-clim$sequence
mds.all.bio<-metaMDS(clim[,grep("bio",colnames(clim))], 'gower', k=2)
mds.all.gb<-metaMDS(clim[,2:22], 'gower', k=2)

factor.all.bio<-envfit(mds.all.gb, clim[,grep("bio",colnames(clim))])
factor.all.gb<-envfit(mds.all.gb, clim[,2:22])

info<-info[match(clim$sequence,info$sequence),]
taxa.tip<-list(quercus.tip=info$sequence[grep("Quercus",info$section)],
               cyc.tip=info$sequence[grep("Cyclobalanopsis",info$section)],
               ilex.tip=info$sequence[grep("Ilex",info$section)])

color.bg<-c()
color.bg[taxa.tip$quercus.tip]<-rep("mediumorchid2",length(taxa.tip$quercus.tip))
color.bg[taxa.tip$cyc.tip]<-rep("palegreen3",length(taxa.tip$cyc.tip))
color.bg[taxa.tip$ilex.tip]<-rep("steelblue2",length(taxa.tip$ilex.tip))
# color.bg[proto.pq]<-rep("firebrick2",length(proto.pq))
color.bg<-color.bg[match(row.names(mds.all.bio$point),names(color.bg))]

## bio 1;19
pdf("all.NMDS.gower.bio.bio11.pdf")
ordisurf(mds.all.bio, clim[, 'bio11'],
         xlim = range(mds.all.bio$points[, 1]),
         lwd.cl = 0.5, col = 'gray80', main = '', cex = 0, add = FALSE, axes = F)
# phylomorphospace(tr.pq, mds.all.bio$points[,1:2],
#                  lwd = 0.75, xlab = '', ylab = '', axes = F,node.size=0,
#                  xlim = range(mds.all.bio$points[, 1]),
#                  ylim = range(c(-10, mds.all.bio$points[, 2])), add = TRUE,label="off")
box()
axis(side = 1, cex.axis = 0.5)
axis(side = 2, cex.axis = 0.5)
points(mds.all.bio$points[,1:2], pch = 21, col = 'black', bg = color.bg, cex = 2)
# text(-10.5,10, substitute(paste("section ", italic('Quercus'))), cex = 0.6, pos = 4)
legend("bottomleft", c("Section cyclobalanopsis","Section Ilex", "Section Quercus"),
       pch = 21, col = 'black',
       pt.bg = c("palegreen3","steelblue2","mediumorchid2"),
       bty = 'n', cex = 1)
title(ylab="MDS axis 2", outer=T, line=2)
title(xlab="MDS axis 1", outer=T, line=1)
dev.off()


## bio+geo
pdf("all.NMDS.gower.GB.bio11.pdf")
ordisurf(mds.all.gb, clim[, 'bio11'],
         xlim = range(mds.all.gb$points[, 1]),
         lwd.cl = 0.5, col = 'gray80', main = '', cex = 0, add = FALSE, axes = F)
# phylomorphospace(tr.pq, mds.all.gb$points[,1:2],
#                  lwd = 0.75, xlab = '', ylab = '', axes = F,node.size=0,
#                  xlim = range(mds.all.gb$points[, 1]),
#                  ylim = range(c(-10, mds.all.gb$points[, 2])), add = TRUE,label="off")
box()
axis(side = 1, cex.axis = 0.5)
axis(side = 2, cex.axis = 0.5)
points(mds.all.gb$points[,1:2], pch = 21, col = 'black', bg = color.bg, cex = 2)
legend("bottomleft", c("Section cyclobalanopsis","Section Ilex", "Section Quercus"),
       pch = 21, col = 'black',
       pt.bg = c("palegreen3","steelblue2","mediumorchid2"),
       bty = 'n', cex = 1)
title(ylab="MDS axis 2", outer=T, line=2)
title(xlab="MDS axis 1", outer=T, line=1)
dev.off()

##投点与环境值
# pdf('out/SUPPLEMENT.eco.mds.ordisurf.v2.pdf', 12, 18)
# todo <- c('bio10', 'bio11', 'bio4', 'bio12', 'latitude', 'longitude')
# layout(matrix(1:6, 3, 2, byrow = TRUE))
# eco.mds.surface <- vector('list', length(todo))
# for(i in todo) {
#   plot(eco.mds.multi[[2]]$all, main = i)
#   eco.mds.surface[[i]] <- ordisurf(eco.mds.multi[[2]]$all, eco.means$all[, i], add = TRUE)
#   text(-10, 10, paste('r2 = ', round(eco.mds.fitted$k2$vectors$r[i], 3)))
# }c
# dev.off()






# ordisurf(mds.ilex.multi[[2]], ilex.clim[, 'bio11']/10,
#          xlim = range(mds.ilex$points[, 1]),
#          lwd.cl = 0.5, col = 'gray80', main = '', cex = 0, add = FALSE, axes = F)
# # phylomorphospace(tr.ilex, mds.ilex$points[,1:2],control=list(col.node=color.bg))
# # phylomorphospace(tree,XX,control=list(col.node=cols),
# #                  xlab="X1",ylab="X2")                 
# phylomorphospace(tr.ilex, mds.ilex$points[,1:2],
#                  lwd = 0.75, xlab = '', ylab = '', axes = F,
#                  control=list(col.node=color.bg),
#                  xlim = range(mds.ilex$points[, 1]),
#                  ylim = range(c(-10, mds.ilex$points[, 2])), add = TRUE,label="off")
# box()
# axis(side = 2, cex.axis = 0.5)
# # points(mds.ilex$points[,1:2], pch = 21, col = 'black',cex = 1)
# points(mds.ilex$points[,1:2], pch = 21, col = 'black', bg = color.bg, cex = 2)
# dev.off()
# 
# 
# box()
# phylomorphospace(drop.tip(tr.eco, which(!tr.eco$tip.label %in% taxa$ro)), eco.mds[["all"]]$points[ro.points,mdsDims],
#                  label = 'off', node.size = 0, lwd = 0.75, xlab = '', ylab = '', axes = F,
#                  xlim = range(eco.mds[["all"]]$points[c(wo.points, ro.points), mdsDims[1]]),
#                  ylim = range(c(-10, eco.mds[["all"]]$points[c(wo.points, ro.points), mdsDims[2]])),
#                  add = TRUE)
# axis(side = 1, cex.axis = 0.5)
# axis(side = 2, cex.axis = 0.5)
# points(eco.mds[["all"]]$points[ro.points,mdsDims], pch = 21, col = 'black',
#        bg = regionColors.figs[sect.species.translate[row.names(eco.mds$all$points)[ro.points], 'subclade']], cex = 1)
# text(-10.5,10, substitute(paste("section ", italic('Lobatae'))), cex = 0.6, pos = 4)

