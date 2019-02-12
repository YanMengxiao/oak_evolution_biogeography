##lineages through time
library(ape)
library(phytools)

setwd("E:/research/analysis/cpdna biogeography quercus/LTT")

trees<-read.nexus("beast2.BD.1node.4-8.trees")
mcc<-read.nexus("beast2.BD.1node.4-8.tre")
mcc<-drop.tip(mcc,grep("Litho",mcc$tip.label))
mcc_copy<-mcc
#mcc<-ladderize(mcc)
info<-read.csv("E:/research/analysis/cpdna biogeography quercus/combination/add.outgroup/info_edit.csv",as.is=T)

##抽树集
trees1000<-sample(trees,size=1000)

num<-sample(5000:6000,100)
trees1000<-trees[num]

trees1000<-lapply(trees1000,drop.tip,tip=trees1000[[1]]$tip.label[grep("Litho",trees1000[[1]]$tip.label)])
class(trees1000)<-"multiPhylo"
trees1000_copy <- trees1000

##展示
info$tip.label<-gsub("[|]","_",info$tip.label)
info$tip.label<-gsub("-","_",info$tip.label)
info$tip.label<-gsub(" ","_",info$tip.label)
taxa.tip<-list(quercus.tip=info$tip.label[grep("Quercus",info$section)],
               cyc.tip=info$tip.label[grep("Cyclobalanopsis",info$section)],
               ilex.tip=info$tip.label[grep("Ilex",info$section)],
               proto.tip=info$tip.label[grep("Protobalanus",info$section)],
               cer.tip=info$tip.label[grep("Cerris",info$section)],
               lob.tip=info$tip.label[grep("Lobatae",info$section)],
               outgroup.tip=info$tip.label[grep("outgroup",info$section)])
##获得各section的edge
taxa.edge<-list(quercus.edge=which.edge(mcc,taxa.tip$quercus.tip),
                cyc.edge=which.edge(mcc,taxa.tip$cyc.tip),
                ilex.edge=which.edge(mcc,taxa.tip$ilex.tip),
                proto.edge=which.edge(mcc,taxa.tip$proto.tip),
                cer.edge=which.edge(mcc,taxa.tip$cer.tip),
                lob.edge=which.edge(mcc,taxa.tip$lob.tip),
                outgroup.edge=which.edge(mcc,taxa.tip$outgroup.tip))

##根据各section的edge颜色
color.edge<-rep("black",length(mcc$edge.length))
color.edge[taxa.edge$quercus.edge]<-rep("purple",length(taxa.edge$quercus.edge))
color.edge[taxa.edge$ilex.edge]<-rep("green4",length(taxa.edge$ilex.edge))
color.edge[taxa.edge$cyc.edge]<-rep("dodgerblue",length(taxa.edge$cyc.edge))
color.edge[taxa.edge$proto.edge]<-rep("red3",length(taxa.edge$proto.edge))
color.edge[taxa.edge$cer.edge]<-rep("gold",length(taxa.edge$cer.edge))
color.edge[taxa.edge$lob.edge]<-rep("orangered",length(taxa.edge$lob.edge))

plot(mcc,"p",cex=0.4,edge.color=color.edge)

###############ltt plot###########################
##获得各section的tip
info$tip.label<-gsub("[|]","_",info$tip.label)
info$tip.label<-gsub("-","_",info$tip.label)
info$tip.label<-gsub(" ","_",info$tip.label)
taxa.tip<-list(quercus.tip=info$tip.label[grep("Quercus",info$section)],
               cyc.tip=info$tip.label[grep("Cyclobalanopsis",info$section)],
               ilex.tip=info$tip.label[grep("Ilex",info$section)],
               proto.tip=info$tip.label[grep("Protobalanus",info$section)],
               cer.tip=info$tip.label[grep("Cerris",info$section)],
               lob.tip=info$tip.label[grep("Lobatae",info$section)],
               outgroup.tip=info$tip.label[grep("outgroup",info$section)])
trees1000<-lapply(trees1000,drop.tip,tip=taxa.tip$lob.tip)
trees1000<-lapply(trees1000,drop.tip,tip=taxa.tip$cer.tip)
trees1000<-lapply(trees1000,drop.tip,tip=taxa.tip$outgroup.tip)
plot(trees1000[[1]],cex=0.3)
mcc<-drop.tip(mcc,taxa.tip$lob.tip)
mcc<-drop.tip(mcc,taxa.tip$cer.tip)
mcc<-drop.tip(mcc,taxa.tip$outgroup.tip)
plot(mcc,cex=0.3)

##cyc
trees.cyc<-trees1000
trees.cyc<-lapply(trees.cyc,drop.tip,tip=taxa.tip$ilex.tip)
trees.cyc<-lapply(trees.cyc,drop.tip,tip=taxa.tip$quercus.tip)
trees.cyc<-lapply(trees.cyc,drop.tip,tip=taxa.tip$proto.tip)
class(trees.cyc)<-"multiPhylo"
#plot(trees.cyc[[1]],cex=0.3)
write.tree(trees.cyc,"cyc.trees.trees")
mcc.cyc<-mcc
mcc.cyc<-drop.tip(mcc.cyc,taxa.tip$ilex.tip)
mcc.cyc<-drop.tip(mcc.cyc,taxa.tip$quercus.tip)
mcc.cyc<-ladderize(drop.tip(mcc.cyc,taxa.tip$proto.tip))
#plot(mcc.cyc,cex=0.3)
write.tree(mcc.cyc,"cyc.mcc.tree")

##ilex
trees.ilex<-trees1000
trees.ilex<-lapply(trees.ilex,drop.tip,tip=taxa.tip$cyc.tip)
trees.ilex<-lapply(trees.ilex,drop.tip,tip=taxa.tip$quercus.tip)
trees.ilex<-lapply(trees.ilex,drop.tip,tip=taxa.tip$proto.tip)
trees.ilex<-lapply(trees.ilex,drop.tip,tip=grep("DM12740|DM11859|DM9783|DM12432|DM6113|DM12995|DM11273|DM8349|DM3152|DM6576|DM13626",trees.ilex[[1]]$tip.label))
class(trees.ilex)<-"multiPhylo"
#plot(trees.ilex[[1]],cex=0.3)
write.tree(trees.ilex,"ilex.trees.trees")
mcc.ilex<-mcc
mcc.ilex<-drop.tip(mcc.ilex,taxa.tip$cyc.tip)
mcc.ilex<-drop.tip(mcc.ilex,taxa.tip$quercus.tip)
mcc.ilex<-ladderize(drop.tip(mcc.ilex,taxa.tip$proto.tip))
mcc.ilex<-ladderize(drop.tip(mcc.ilex,grep("DM12740|DM11859|DM9783|DM12432|DM6113|DM12995|DM11273|DM8349|DM3152|DM6576|DM13626",mcc.ilex$tip.label)))
#plot(mcc.ilex,cex=0.3)
write.tree(mcc.ilex,"ilex.mcc.tree")

##quercus&proto
trees.pq<-trees1000
trees.pq<-lapply(trees.pq,drop.tip,tip=taxa.tip$cyc.tip)
trees.pq<-lapply(trees.pq,drop.tip,tip=taxa.tip$ilex.tip)
trees.pq<-lapply(trees.pq,drop.tip,tip=trees.pq[[1]]$tip.label[grep("chrysolepis|tomentella|vacciniifolia|greggii|chrysolepnis",trees.pq[[1]]$tip.label)])
class(trees.pq)<-"multiPhylo"
#plot(trees.pq[[1]],cex=0.3)
write.tree(trees.pq,"que.trees.trees")

mcc.pq<-mcc
mcc.pq<-drop.tip(mcc.pq,taxa.tip$cyc.tip)
mcc.pq<-drop.tip(mcc.pq,taxa.tip$ilex.tip)
mcc.pq<-ladderize(drop.tip(mcc.pq,grep("chrysolepis|tomentella|vacciniifolia|greggii|chrysolepnis",mcc.pq$tip.label)))
#plot(mcc.pq,cex=0.3)
write.tree(mcc.pq,"que.mcc.tree")

#反转时间轴
plot(obj,xaxis="flipped")
plot(obj,xaxis="negative")

par(mfrow=c(4,3))
#opar <- par(no.readonly=T)
ltt_cyc <- ltt(trees.cyc,plot=T,col="gray")
#plot(ltt_cyc,col="gray",main="Cylobalanopsis")
ltt95_cyc<-ltt95(trees.cyc,xaxis=c("negative"))
mtext("Cylobalanopsis",side=3,line=1) ##设置line移动文本
ltt1_cyc <- ltt.plot(mcc.cyc,xlab="time",ylab="log(lineages)")
ltt_ilex <- ltt(trees.ilex,col="gray")
ltt95_ilex <-ltt95(trees.ilex,xaxis=c("negative"))
mtext("Ilex",side=3,line=1)
ltt1_ilex <- ltt.plot(mcc.ilex,xlab="time",ylab="log(lineages)")
ltt_pq <- ltt(trees.pq,col="gray")
ltt95_pq<-ltt95(trees.pq,xaxis=c("negative"))
mtext("Quercus and Protobalanus",side=3,line=1)
ltt1_pq <- ltt.plot(mcc.pq,xlab="time",ylab="log(lineages)")
ltt_all <- ltt(trees1000_copy,col="gray")
ltt95_all <-ltt95(trees1000_copy,xaxis=c("negative"))
mtext("All samples",side=3,line=1)
ltt1_allc <- ltt.plot(mcc_copy,xlab="time",ylab="log(lineages)")


#par(opar)

##
#max1 =c()
#for (i in 1:10) {max1<-c(max,max(ltt_cyc[[i]]$times))}


##ilex
trees100<-lapply(trees10,drop.tip,tip=taxa.tip$cyc.tip)
trees100<-lapply(trees10,drop.tip,tip=taxa.tip$quercus.tip)
trees100<-lapply(trees10,drop.tip,tip=taxa.tip$proto.tip)

##quercus&proto
trees100<-lapply(trees10,drop.tip,tip=taxa.tip$cyc.tip)
trees100<-lapply(trees10,drop.tip,tip=taxa.tip$ilex.tip)

plot(drop.tip(trees10[[1]],taxa.tip$cyc.tip))
class(tree)="multiphylo"
plot(trees10[[2]])
trees10<-lapply(trees10,drop.tip,tip=taxa.tip$cyc.tip)
tree2<-trees10[1]





tree1<-drop.tip(tree1,tip=taxa.tip$cyc.tip)
tree1<-drop.tip(tree1,tip=taxa.tip$cer.tip)
tree1<-drop.tip(tree1,tip=taxa.tip$quercus.tip)
tree1<-drop.tip(tree1,tip=taxa.tip$proto.tip)
tree1<-drop.tip(tree1,tip=taxa.tip$lob.tip)



plot(tree1,"p",cex=0.4,edge.color=color.edge)

trees100<-lapply(trees100,drop.tip,tip=trees100[[1]]$tip.label[grep("ilex",trees100[[3]]$tip.label)])


trees100<-lapply(trees100,drop.tip,tip=taxa.tip$lob.tip)
plot(trees100[1],cex=0.3)





