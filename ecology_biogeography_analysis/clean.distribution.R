setwd("E:/research/analysis/cpdna biogeography quercus/cyc.distribution")

origin<-read.csv("cyc1.csv")
gps<-read.csv("cyc.orin.csv")

gps$species.name <- origin$species.name[match(gps$number,origin$number)]
gps$final.name <- origin$final.name[match(gps$number,origin$number)]

write.csv(gps,"cyc.collection.gps.csv")

setwd("E:/research/analysis/cpdna biogeography quercus/climate/quercus")
gbif<-read.csv("andrew.csv")
colnames(gbif)
gbif2<-data.frame(number=1:length(which(!is.na(gbif$decimallatitude))),
                  species=gbif$species[which(!is.na(gbif$decimallatitude))],
                  scientificname=gbif$scientificname[which(!is.na(gbif$decimallatitude))],
                  longitude=gbif$decimallongitude[which(!is.na(gbif$decimallongitude))],
                  latitude=gbif$decimallatitude[which(!is.na(gbif$decimallatitude))])
write.csv(gbif,"andrew.csv")

library(stringr)
content<-read.csv("NHM.csv") 
species<- unlist(str_match_all(content$content,"Quercus .*, Collector"))
species<- gsub(", [C|D].*", "", species)
lat<- unlist(str_match_all(content$content,"Lat: .*Lon"))
lat<-gsub("Lat: ([0-9].*[0-9]), Lon.*","\\1",lat,perl=T)
lat<-gsub("Lat: (-[0-9].*[0-9]), Lon.*","\\1",lat,perl=T)
lon<- unlist(str_match_all(content$content,"Lon:.*[0-9], "))
lon<-gsub("Lon: (.*[0-9]), .*","\\1",lon)
lon<-gsub(", .*","",lon)
locality<-unlist(str_match_all(content$content,"Lon.*"))
locality<- gsub("Lon: [0-9].*[0-9], ","", locality)
country<-
NHM <- data.frame(species=species, lon=lon, lat=lat, locality=locality)
write.csv(NHM,"NHM2.csv")
NHM_SEA<-NHM[grep("Laos|Malaysia|Thailand|Indonesia",NHM$locality),]
write.csv(NHM_SEA,"NHM_SEA.csv")
