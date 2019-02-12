library(ape)
library(phytools)
setwd("E:/research/analysis/cpdna biogeography quercus/ancestral distribution")


##抽取100棵树
trees<-read.nexus("beast2.BD.1node.4-8.trees")
# info<-read.csv("E:/research/analysis/cpdna biogeography quercus/combination/add.outgroup/info_edit.csv",as.is=T)
# trees1000<-sample(trees,size=100)
num<-sample(5000:6000,100)
trees100<-trees[num]
trees100<-lapply(trees1000,drop.tip,tip=trees1000[[1]]$tip.label[grep("Litho",trees1000[[1]]$tip.label)])
class(trees100)<-"multiPhylo"
write.tree(trees100,"beast.trees100.trees")
# for (i in 1:100) {
#   trees100[[i]]<-ladderize(root(trees100[[i]],grep('Litho',trees100[[i]]$tip.label)))
# }
##不可循环批量改名#
##批量去掉tip，for循环里失败，提示droptip后tip数量不一致；完成后文件为list，修改为multiPhylo
# trees100<-lapply(trees100,drop.tip,tip=trees100[[1]]$tip.label[grep("Litho|Fagus|Castanea|Castanopsis|Notholitho",trees100[[3]]$tip.label)])
# class(trees100)<-"multiPhylo"
#plot(trees100[[55]],cex=0.3)

#trees200<-sapply(trees100, function(x) {drop.tip(x,grep("Fagus|Castanea|Castanopsis",x$tip.label))})

##consesus tree
mcc<-read.nexus("beast2.BD.1node.4-8.tre")
mcc<-drop.tip(mcc,grep("Litho",mcc$tip.label))
write.tree(mcc,"mcc.tree")

tree.con<-read.nexus("E:/research/analysis/cpdna biogeography quercus/combination/add.outgroup/all.section/haplotype/rasp/rasp.nex.con.tre")
tree.con<-ladderize(root(tree.con,grep('Litho',tree.con$tip.label)))
tree.con<-multi2di(tree.con)
tree.con<-drop.tip(tree.con,grep("Litho|Fagus|Castanea|Castanopsis|Notholitho",tree.con$tip.label))
write.tree(tree.con,"consensu.tre")

###distribution table
distribution<- data.frame(sequence=gsub("[|]","_",info$tip.label[which(info$rasp)]),
                          distribution=info$region[which(info$rasp)])
write.csv(distribution,"distribution.csv")


