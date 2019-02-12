## climate

setwd("E:/research/analysis/cpdna biogeography quercus/climate/section")
cyc.all<-read.csv("cyc GPS with climate data.csv")
ilex.all<-read.csv("ilex_clim.csv")
quercus.all<-read.csv("Quercus GPS with climate data.csv")

# cyc.all<-read.csv("cyc.clim.30.csv")
# ilex.all<-read.csv("ilex_clim.csv")
# quercus.all<-read.csv("que.clim.30.csv")

# 去除NA
cyc.all<-cyc.all[which(!is.na(cyc.all$bio1)),]
ilex.all<-ilex.all[which(!is.na(ilex.all$bio1)),]
quercus.all<-quercus.all[which(!is.na(quercus.all$bio1)),]

# 合并成总表
cyc.all$section<-"Cyclobalanopsis"
ilex.all$section<-"Ilex"
quercus.all$section<-"Quercus"
clim.all<-rbind(cyc.all,ilex.all,quercus.all)

地图投点
library(maps)
pdf("map.ilex.new.pdf",width = 96, height = 48)
map("world")
points(ilex.all$longitude,ilex.all$latitude,col="blue",cex=3)
dev.off()
# pdf("map.cyc.pdf",width = 96, height = 48)
# map("world")
# points(cyc.all$longitude,cyc.all$latitude,col="blue",cex=3)
# dev.off()
# pdf("map.quercus.pdf",width = 96, height = 48)
# map("world")
# points(quercus.all$longitude,quercus.all$latitude,col="blue",cex=3)
# dev.off()
# clim.all$section<-c(rep("Cyclobalanopsis",2241),rep("Ilex",4304),rep("Quercus",5788))

# 3者去除自相关重叠部分1,2,3,7,12,14,15
# maxant计算出的3者共同限制性因子1,3,7,12
# selc<-c(paste("bio",c(1,2,3,5,6,7,12,14,15,16,17),sep=""),"section")
selc<-c(paste("bio",c(1,3,7,12),sep=""),"section")
clim.selc<-clim.all[,selc]


ilex.all<-as.numeric(ilex.all$bio1)
ilex.all<-as.data.frame(apply(ilex.all,2,as.numeric))

library(vegan)
mds.ilex<-metaMDS(ilex.all[,selc], 'gower', k=2)
mds.cyc<-metaMDS(cyc.all[,selc], 'gower', k=2)
mds.quercus<-metaMDS(quercus.all[,selc], 'gower', k=2)

library(ggplot2)
clim.selc$section<-factor(clim.selc$section,levels=c("Cyclobalanopsis","Ilex","Quercus"))

mytheme<-theme(panel.background = element_rect(fill = "white", colour = "grey50"),
               panel.grid.major = element_line(colour = NA),
               panel.grid.minor = element_line(colour = NA))

bio1_p1<-ggplot(clim.selc,aes(bio1,colour=section,fill=section)) +   
  geom_density(alpha=.3) + 
  mytheme +
  labs(x="Annual Mean Temperature")
bio12_p2<-ggplot(clim.selc,aes(bio12,colour=section,fill=section)) +   
  geom_density(alpha=.3)+ 
  mytheme +
  labs(x="Annual Precipitation")
library(gridExtra)
all<-grid.arrange(bio1_p1,bio12_p2,nrow=1,ncol=2)
ggsave("climate.differentiation2.pdf",all,width=16, height=4)


## bio1
bio1_p1<-ggplot(clim.selc,aes(bio1,colour=section,fill=section)) +   
  geom_density(alpha=.3) 
bio1_p2<-ggplot(clim.selc,aes(section,bio1,fill=section)) +   
  geom_violin() + 
  geom_boxplot(width=.1,fill="white",alpha=.8)
bio1_p3<-ggplot(clim.selc,aes(section,bio1,fill=section)) +   
  geom_boxplot()

library(gridExtra)
bio1<-grid.arrange(bio1_p1,bio1_p2,bio1_p3,nrow=2,ncol=2)
ggsave("bio1.2.5.pdf",bio1,width=12, height=6)

## bio3
bio3_p1<-ggplot(clim.selc,aes(bio3,colour=section,fill=section)) +   
  geom_density(alpha=.3) 
bio3_p2<-ggplot(clim.selc,aes(section,bio3,fill=section)) +   
  geom_violin() + 
  geom_boxplot(width=.1,fill="white",alpha=.8)
bio3_p3<-ggplot(clim.selc,aes(section,bio3,fill=section)) +   
  geom_boxplot()

# library(gridExtra)
bio3<-grid.arrange(bio3_p1,bio3_p2,bio3_p3,nrow=2,ncol=2)
ggsave("bio3.2.5.pdf",bio3,width=12, height=6)

## bio7
bio7_p1<-ggplot(clim.selc,aes(bio7,colour=section,fill=section)) +   
  geom_density(alpha=.3) 
bio7_p2<-ggplot(clim.selc,aes(section,bio7,fill=section)) +   
  geom_violin() + 
  geom_boxplot(width=.1,fill="white",alpha=.8)
bio7_p3<-ggplot(clim.selc,aes(section,bio7,fill=section)) +   
  geom_boxplot()

# library(gridExtra)
bio7<-grid.arrange(bio7_p1,bio7_p2,bio7_p3,nrow=2,ncol=2)
ggsave("bio7.2.5.pdf",bio7,width=12, height=6)

## bio12
bio12_p1<-ggplot(clim.selc,aes(bio12,colour=section,fill=section)) +   
  geom_density(alpha=.3) 
bio12_p2<-ggplot(clim.selc,aes(section,bio12,fill=section)) +   
  geom_violin() + 
  geom_boxplot(width=.1,fill="white",alpha=.8)
bio12_p3<-ggplot(clim.selc,aes(section,bio12,fill=section)) +   
  geom_boxplot()

# library(gridExtra)
bio12<-grid.arrange(bio12_p1,bio12_p2,bio12_p3,nrow=2,ncol=2)
ggsave("bio12.2.5.pdf",bio12,width=12, height=6)

# ggplot(clim.selc,aes(bio17,colour=section,fill=section)) +   
#   geom_density(alpha=.3) 
# ggplot(clim.selc,aes(section,bio1,fill=section)) +   
#   geom_boxplot()
# ggplot(clim.selc,aes(section,bio1,fill=section)) +   
#   geom_violin() + 
#   geom_boxplot(width=.1,fill="white",alpha=.8)
# ggplot(clim.selc,aes(section,bio2,fill=section)) +   
#   geom_violin() + 
#   geom_boxplot(width=.1,fill="white",alpha=.8)

ggsave("p_dis_sections.pdf")
ggplot(clim.selc,aes(bio2,colour=section,fill=section)) +   
  geom_density(alpha=.3)  


ggplot(pcoa2,aes(axis1,colour=Locality)) + #用ggplot作图  
  geom_density()   

# plot(density(clim.selc[,1]))

##共线性分析
library(psych)
corr.test.ilex<-corr.test(cyc.all[,4:22], use="complete")
write.csv(corr.test.ilex[[1]],"corr.cyc.csv")
write.csv(corr.test.ilex[[4]],"corr.test.p-value.ilex.csv")


setwd("E:/research/analysis/cpdna biogeography quercus/genetic.diversity")
genalex<-read.csv("genalex.genetic.dist.csv")
# genalex$section[1:7140]<-"Ilex"
# genalex$section[7141:11145]<-"Quercus"
# genalex$section[11146:19923]<-"Cyclobalanopsis"
# write.csv(genalex,"genalex.genetic.dist.csv")

mytheme<-theme(panel.background = element_rect(fill = "white", colour = "grey50"),
               panel.grid.major = element_line(colour = NA),
               panel.grid.minor = element_line(colour = NA))

genalex$section<-factor(genalex$section,levels=c("Cyclobalanopsis","Ilex","Quercus"))
genalex.plot<-ggplot(genalex,aes(section,genetic.dist,fill=section)) +   
  geom_violin(alpha=.8) + 
  geom_boxplot(width=.1,fill="white",alpha=.8) +
  mytheme +
  labs(x="section",y="Genetic distance")
ggsave("genalex.genetic.dist.pdf",genalex.plot)

p.dist<-read.csv("p.dist.csv",as.is=T)

p.dist$section[3917:10937]<-"Ilex"
p.dist$section[1:3916]<-"Quercus"
p.dist$section[10938:21668]<-"Cyclobalanopsis"
write.csv(genalex,"p.dist2.csv")
p.dist$section<-factor(p.dist$section,levels=c("Cyclobalanopsis","Ilex","Quercus"))
p.dist.plot<-ggplot(p.dist,aes(section,genetic.dist,fill=section)) +   
  geom_violin(alpha=.8) + 
  geom_boxplot(width=.1,fill="white",alpha=.8) +
  mytheme +
  labs(x="section",y="Genetic distance")
ggsave("p.dist.pdf",p.dist.plot)
