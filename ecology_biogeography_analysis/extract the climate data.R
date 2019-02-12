library(raster)
library(dismo)
setwd("E:/research/analysis/cpdna biogeography quercus/climate")#设置工作路径

info<-read.csv("E:/research/analysis/cpdna biogeography quercus/combination/add.outgroup/info_edit.csv",as.is=T)
cyc.gps<-info$sequence[]

####cyc
dir <- list.files(path = "E:/research/analysis/cpdna biogeography quercus/climate/bio_2-5m_bil2",pattern = "bil",full.names = T)#导入气候变量一个或者多个均可
dir
clim <- stack(dir)
clim
occ_cyc_total <- read.csv("cyc.gps.final.csv",header = T,sep = ",")#导入分布点GPS
dim(occ_cyc_total)
occ_cyc_total <- occ_cyc_total[,-1]#删去第一列
climate_cyc_total<-extract(clim,occ_cyc_total) #提取气候因子
dim(climate_cyc_total)
climate_cyc_total<-cbind(occ_cyc_total,climate_cyc_total)
write.csv(climate_cyc_total,"cyc_total_climate.csv")

###ilex
# dir <- list.files(path = "E:/research/analysis/cpdna biogeography quercus/climate/bio_2-5m_bil",pattern = "bil",full.names = T)#导入气候变量一个或者多个均可
# dir
# clim <- stack(dir)
# clim
gps_ilex_total <- read.csv("ilex.distribution.csv",header = T,sep = ",")#导入分布点GPS
dim(gps_ilex_total)
gps_ilex_total <- gps_ilex_total[,-1]#删去第一列
climate_ilex_total<-extract(clim,gps_ilex_total)
dim(climate_ilex_total)
climate_ilex_total<-cbind(gps_ilex_total,climate_ilex_total)
write.csv(climate_ilex_total,"ilex_total_climate.csv")

###quercus
###quercus
dir <- list.files(path = "E:/research/analysis/cpdna biogeography quercus/climate/bio_2-5m_bil2",pattern = "bil",full.names = T)#导入气候变量一个或者多个均可
dir
clim <- stack(dir)
clim
gps_quercus_total <- read.csv("quercus.gps.final.csv",header = T,sep = ",")#导入分布点GPS
dim(gps_quercus_total)
# gps_quercus_total <- gps_quercus_total[,-1]#删去第一列
climate_quercus_total<-extract(clim,gps_quercus_total)
dim(climate_quercus_total)
climate_quercus_total<-cbind(gps_quercus_total,climate_quercus_total)
write.csv(climate_quercus_total,"quercus_total_climate.csv")

quercus.gps<-read.csv("E:/research/analysis/cpdna biogeography quercus/climate/cyc.distribution/cyc.combination.gps.csv",as.is=T)
quercus.gps$longtitude<-as.numeric(quercus.gps$longtitude)
quercus.gps$source[73197:151898]<-"gbif robur"
class(quercus.gps$source)
write.csv(quercus.gps,"E:/research/analysis/cpdna biogeography quercus/climate/quercus/quercus.qps.all2.csv")

quercus<-read.csv("andrew.csv")
gbif<-read.csv("quercus.gbif2.csv")
robur<-read.csv("gbif.robur.plot.csv")
yan<-read.csv("pq.gps.csv")
paper<-read.csv("paper.csv")
cvh<-read.csv("cvh.gps.csv")
mentisky<-read.csv("mentisky.csv")

library(maps)
pdf("map.cyc.pdf",width = 96, height = 48)
map("world")
# #points(quercus$longitude,quercus$latitude,col="red",cex=3,pch=16)
points(quercus.gps$longitude,quercus.gps$latitude,col="blue",cex=3,pch=16)
# #points(robur$decimallongitude,robur$decimallatitude,col="green",cex=3,pch=16)
#points(yan$longitude,yan$latitude,col="chocolate4",cex=3,pch=16)
# points(paper$longitude,paper$latitude,col="purple",cex=3,pch=16)
# points(cvh$longitude,cvh$latitude,col="deeppink1",cex=3,pch=16)
# points(mentisky$longitude,mentisky$latitude,col="deeppink1",cex=3,pch=16)
dev.off()

