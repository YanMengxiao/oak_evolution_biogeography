##dated tree with time scale

##安装所需的包
install.packages("strap",dependencies=TRUE)
#TRUE means to use c("Depends", "Imports", "LinkingTo", "Suggests") for pkgs and c("Depends", "Imports", "LinkingTo") for added dependencies: this installs all the packages needed to run pkgs, their examples, tests and vignettes (if the package author specified them correctly).
install.packages("phytools")
install.packages("devtools")
library("devtools")
install_url("http://www.christophheibl.de/phyloch_1.5-5.zip") #windows版本
install.packages("coda")

source("https://bioconductor.org/biocLite.R")
biocLite("ggtree")

library("phytools")
library("phyloch")
library("strap")
library("coda")
setwd("E:/research/analysis/cpdna biogeography quercus/LTT")
mcc<-read.beast("beast2.BD.1node.4-8.tre")
mcc<-ladderize(mcc,right=F)
mcc$root.time <- mcc$height[1]
info<-read.csv("E:/research/analysis/cpdna biogeography quercus/combination/add.outgroup/info_edit.csv",as.is=T)
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
taxa.edge<-list(quercus.edge=which.edge(mcc,taxa.tip$quercus.tip),
                cyc.edge=which.edge(mcc,taxa.tip$cyc.tip),
                ilex.edge=which.edge(mcc,taxa.tip$ilex.tip),
                proto.edge=which.edge(mcc,taxa.tip$proto.tip),
                cer.edge=which.edge(mcc,taxa.tip$cer.tip),
                lob.edge=which.edge(mcc,taxa.tip$lob.tip),
                outgroup.edge=which.edge(mcc,taxa.tip$outgroup.tip))
color.edge<-rep("black",length(mcc$edge.length))
color.edge[taxa.edge$quercus.edge]<-rep("purple",length(taxa.edge$quercus.edge))
color.edge[taxa.edge$ilex.edge]<-rep("green4",length(taxa.edge$ilex.edge))
color.edge[taxa.edge$cyc.edge]<-rep("dodgerblue",length(taxa.edge$cyc.edge))
color.edge[taxa.edge$proto.edge]<-rep("red3",length(taxa.edge$proto.edge))
color.edge[taxa.edge$cer.edge]<-rep("gold",length(taxa.edge$cer.edge))
color.edge[taxa.edge$lob.edge]<-rep("orangered",length(taxa.edge$lob.edge))

pdf("time.scale.tree.95.pdf",width=12,height=10)
geoscalePhylo(tree=mcc, boxes="Age", edge.color=color.edge,cex.tip=0.3,cex.age=1,cex.ts=1,label.offset=0,x.lim=c(-8,70),lwd=1.5)
nodelabels(paste(round(mcc$`height_95%_HPD_MIN`,2),round(mcc$`height_95%_HPD_MAX`,2),sep="-"),cex=1,frame="n",bg="white",adj=c(-0.2,0.2))
dev.off()

tr = rtree(5)
d=data.frame(node=1:9,color=sample(c("red","blue","green"),9,replace=T))
#colnames(d)=c("tip","color")
ggtree(tr) %<+% d +aes(color=I(color))

##ggtree
library(ggtree)
#file <- system.file("extdata/BEAST", "beast_mcc.tree", package="ggtree")
beast <- read.beast("beast2.BD.1node.ATGC.5-8.tre")
tree1 <- read.nexus("beast2.BD.1node.ATGC.5-8.tre")
# beast_data <- fortify(beast)
# head(beast_data)
info<-read.csv("E:/research/analysis/cpdna biogeography quercus/combination/add.outgroup/info_edit.csv",as.is=T)
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
taxa.edge<-list(quercus.edge=which.edge(beast@phylo,taxa.tip$quercus.tip),
                cyc.edge=which.edge(beast@phylo,taxa.tip$cyc.tip),
                ilex.edge=which.edge(beast@phylo,taxa.tip$ilex.tip),
                proto.edge=which.edge(beast@phylo,taxa.tip$proto.tip),
                cer.edge=which.edge(beast@phylo,taxa.tip$cer.tip),
                lob.edge=which.edge(beast@phylo,taxa.tip$lob.tip),
                outgroup.tip=which.edge(beast@phylo,taxa.tip$outgroup.tip))
color.edge<-rep("black",300)
color.edge[taxa.edge$quercus.edge]<-rep("purple",length(taxa.edge$quercus.edge))
color.edge[taxa.edge$ilex.edge]<-rep("green4",length(taxa.edge$ilex.edge))
color.edge[taxa.edge$cyc.edge]<-rep("dodgerblue",length(taxa.edge$cyc.edge))
color.edge[taxa.edge$proto.edge]<-rep("red3",length(taxa.edge$proto.edge))
color.edge[taxa.edge$cer.edge]<-rep("gold",length(taxa.edge$cer.edge))
color.edge[taxa.edge$lob.edge]<-rep("orangered",length(taxa.edge$lob.edge))
d=data.frame(node=1:300,
             color=color.edge)
#beast <- groupOTU(beast, taxa.tip)
plot.tree<-ggtree(beast,ladderize=FALSE) %<+% d +aes(color=I(color)) + 
  theme_tree2()  + 
  #geom_range(range='height_0.95_HPD', color='blue', alpha=.6, size=2) + 
  geom_tiplab(color="black",size=2)+
  xlim(0, 100)
revts(plot.tree) ##时间轴反向


##示范脚本
library("phytools")
library("phyloch")
library("strap")
library("coda")

setwd("E:/research/analysis/cpdna biogeography quercus/LTT")

t <- read.beast("beast2.BD.1node.4-8.trees")
t$root.time <- t$height[1]

log_data1 <- read.table("bearsDivtime_FBD.1.log",header=T)
log_data2 <- read.table("bearsDivtime_FBD.2.log",header=T)
nMCMC <- length(log_data1$originFBD)-1 ##MCMC代数
burnin <- 0.2 ##burnin比例
id1 <- as.integer(nMCMC*burnin)+1
id2 <- nMCMC+1
#comb_origin_data <- c(log_data1$originFBD[id1:id2],log_data2$originFBD[id1:id2]) 
#手动burnin，截取10001-50001代的log
#非FBD无此项参数 log_data1$originFBD
# origin_mcmc <-as.mcmc(comb_origin_data)
# origin_time <- mean(origin_mcmc)
# stem_length <- origin_time - t$root.time
# origin_HPD <- HPDinterval(origin_mcmc)

num_taxa <- length(t$tip.label)

##筛选tip里面有0的，即非化石的tip
# display_all_node_bars <- FALSE
# 
# names_list <-vector()
# for (name in t$tip){
#   v <- strsplit(name,"_")[[1]]
#   if(display_all_node_bars){
#     names_list = c(names_list,name)
#   }
#   else if(v[length(v)]=="0"){
#     names_list = c(names_list,name)
#   }
# }

##找祖先节点
# nids <- vector()
# pos <- 1
# len_nl <- length(names_list)
# for(n in names_list){
#   for(nn in names_list[pos:len_nl]){
#     if(n != nn){
#       m <- getMRCA(t,c(n,nn))
#       if(m %in% nids == FALSE){
#         nids <- c(nids,m)
#       }
#     }
#   }
#   pos<-pos+1
# }
##tip.label现生种去除0，化石种第三位标记X
# for(tp in 1:length(t$tip.label)){
#   v <- strsplit(t$tip.label,"_")[[tp]]
#   new_l <- v[1]
#   if(length(v) > 2){
#     new_l <- paste(v[1],v[2],sep="_")
#   }
#   if(tail(v,1) != "0"){
#     new_l <- paste(new_l, "X",sep=" ")
#   }
#   t$tip.label[tp] <- new_l
# }

root_max <- t$"height_95%_HPD_MAX"[1]
# x_max <- origin_HPD[2] * 0.1 + origin_HPD[2] 
#origin_HPD基于FBD参数计算出originFBD，origin_HPD <- HPDinterval(as.mcmc(c(log_data1$originFBD[id1:id2],log_data2$originFBD[id1:id2])))

########plot######
pdf("geoscaled_bears2.pdf", width=10, height=7)
geoscalePhylo(tree=ladderize(t,right=FALSE), boxes="Age", cex.tip=1.2,cex.age=1,
              cex.ts=1,label.offset=0,x.lim=c(-15,x_max),lwd=1.5)

lastPP<-get("last_plot.phylo",envir=.PlotPhyloEnv)

##origin bar&point
origin_xx <- c(lastPP$xx[num_taxa+1],-stem_length)
lines(origin_xx,c(lastPP$yy[num_taxa+1],lastPP$yy[num_taxa+1]),lwd=1.5,col="black")
bar_xx_o <- c(t$root.time-origin_HPD[1], t$root.time-origin_HPD[2])
lines(bar_xx_o,c(lastPP$yy[num_taxa+1],lastPP$yy[num_taxa+1]),col=rgb(1,0,0,alpha=0.3),lwd=8)
points(-stem_length,lastPP$yy[num_taxa+1],pch=15,col="red",cex=1.5)

##节点bar
for(nv in nids){
  bar_xx_a <- c(lastPP$xx[nv]+t$height[nv-num_taxa]-t$"height_95%_HPD_MIN"[nv-num_taxa], lastPP$xx[nv]-(t$"height_95%_HPD_MAX"[nv-num_taxa]-t$height[nv-num_taxa]))
  lines(bar_xx_a,c(lastPP$yy[nv],lastPP$yy[nv]),col=rgb(0,0,1,alpha=0.3),lwd=8)
}

t$node.label<-t$posterior
p <- character(length(t$node.label))
p[t$node.label >= 0.95] <- "black"
p[t$node.label < 0.95 & t$node.label >= 0.75] <- "gray"
p[t$node.label < 0.75] <- "white"
nodelabels(pch=21, cex=1.5, bg=p)
dev.off()

