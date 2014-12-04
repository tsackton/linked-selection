#this script will estimate area of each species range from use

#method: use alpah-hulls (which at the limit of alpha converge to convex hulls) to estimate range based on occurence data from GBIF, unless otherwise indicated.

source("ahull_to_polygon.R")

library(dismo)
library(maptools)
library(alphahull)
library(rgeos)
library(rgdal)
library(sp)
data(wrld_simpl)

cty.borders<-getData('countries')
oceans<-readOGR("ne_50m_ocean", "ne_50m_ocean") #Natural Earth oceans shapefile, this needs to be downloaded from http//www.naturalearthdata.com/download/50m/physical/ne_50m_ocean.zip and unzipped into this directory

#Polymorphism Species
#Anopheles gambiae (agam)
agam.use<-subset(agam.use, country != "Brazil") #gambiae was a brief invader to Brazil, but not part of its native range

agam.clean<-agam.use[,c("lon", "lat")]
agam.clean<-agam.clean[!duplicated(agam.clean),]
agam.ahull<-ahull(agam.clean$lon, agam.clean$lat, alpha=20); plot(agam.ahull)
agam.range<-ah2sp(agam.ahull); plot(agam.range)
proj4string(agam.range)<-CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
agam.range.clean<-gDifference(agam.range, oceans)
xlimits<-agam.range.clean@bbox[1,]
ylimits<-agam.range.clean@bbox[2,]
plot(wrld_simpl, axes=T, col="light yellow", xlim=xlimits, ylim=ylimits)
plot(agam.range.clean, lwd=2, axes=T, add=T, density=15, angle=45)
#project to Lambert Equal Area projection to calculate area in meters
agam.proj<-spTransform(agam.range.clean, CRS("+proj=cea +lon_0=0"))

#Apis mellifera scutellata
amel.use<-subset(amel.use, basisOfRecord!="unknown") #sufficient specimen data points to not require use of observation or unknown points
amel.use<-subset(amel.use, lon > -20) #native to old world
amel.use<-subset(amel.use, country != "Australia") #introduced to Australia
amel.use<-subset(amel.use, lon < 100) #introduced to east Asia
#remove wrong subspecies
amel.use<-subset(amel.use, !grepl("syriaca", species))
amel.use<-subset(amel.use, !grepl("ligustica", species))
amel.use<-subset(amel.use, !grepl("iberiensis", species))
amel.use<-subset(amel.use, !grepl("jemenitica", species))
#limit range to Africa since we used the African subspecies
amel.use<-subset(amel.use, lat < 20)

amel.clean<-amel.use[,c("lon", "lat")]
amel.clean<-amel.clean[!duplicated(amel.clean),]
amel.ahull<-ahull(amel.clean$lon, amel.clean$lat, alpha=1e6); plot(amel.ahull)
amel.range<-ah2sp(amel.ahull); plot(amel.range)
proj4string(amel.range)<-CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
amel.range.clean<-gDifference(amel.range, oceans)
xlimits<-amel.range.clean@bbox[1,]
ylimits<-amel.range.clean@bbox[2,]
plot(wrld_simpl, axes=T, col="light yellow", xlim=xlimits, ylim=ylimits)
plot(amel.range.clean, lwd=2, axes=T, add=T, density=15, angle=45)
#project to Lambert Equal Area projection to calculate area in meters
amel.proj<-spTransform(amel.range.clean, CRS("+proj=cea +lon_0=0"))

#Arabidopsis thaliana
atha.use=subset(atha.use, basisOfRecord == "observation" | basisOfRecord=="specimen")
atha.use=subset(atha.use, lon < 90 & lon > -30) #native to Eurasia
atha.use=subset(atha.use, lat < 100 & lat > 25) #native to Eurasia
atha.use=subset(atha.use, country != "Libya" & country != "Morocco") #native to Eurasia
atha.clean<-atha.use[,c("lon", "lat")]
atha.clean<-atha.clean[!duplicated(atha.clean),]
atha.ahull<-ahull(atha.clean$lon, atha.clean$lat, alpha=15); plot(atha.ahull)
atha.range<-ah2sp(atha.ahull, rnd=2); plot(atha.range)
proj4string(atha.range)<-CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
atha.range.clean<-gDifference(atha.range, oceans)
xlimits<-atha.range.clean@bbox[1,]
ylimits<-atha.range.clean@bbox[2,]
plot(wrld_simpl, axes=T, col="light yellow", xlim=xlimits, ylim=ylimits)
plot(atha.range.clean, lwd=2, axes=T, add=T, density=15, angle=45)
#project to Lambert Equal Area projection to calculate area in meters
atha.proj<-spTransform(atha.range.clean, CRS("+proj=cea +lon_0=0"))

#Brachypodium distachyon
#native range is circum-Mediterranean
bdis.use=subset(bdis.use, country != "Australia" & country != "Argentina" & country != "Colombia" & country != "Ecuador" & country != "Mexico" & country != "Namibia" & country != "Nepal" & country != "Peru" & country != "South Africa" & country != "United Kingdom" & country != "Ukraine" & country != "United States")
bdis.use=subset(bdis.use, basisOfRecord!="unknown")
#remove miscoded points in Africa with localities indicating country is in Europe
bdis.use=subset(bdis.use, lat > 25)
bdis.clean<-bdis.use[,c("lon", "lat")]
bdis.clean<-bdis.clean[!duplicated(bdis.clean),]
bdis.ahull<-ahull(bdis.clean$lon, bdis.clean$lat, alpha=20); plot(bdis.ahull)
bdis.range<-ah2sp(bdis.ahull); plot(bdis.range)
proj4string(bdis.range)<-CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
bdis.range.clean<-gDifference(bdis.range, oceans)
xlimits<-bdis.range.clean@bbox[1,]
ylimits<-bdis.range.clean@bbox[2,]
plot(wrld_simpl, axes=T, col="light yellow", xlim=xlimits, ylim=ylimits)
plot(bdis.range.clean, add=T, lwd=2, axes=T, density=15, angle=45)
#project to Lambert Equal Area projection to calculate area in meters
bdis.proj<-spTransform(bdis.range, CRS("+proj=cea +lon_0=0"))

#Bombyx mandarina
#load data from 40 genome resequencing paper which has geotagged samples sequenced
#reference: http://www.sciencemag.org/content/326/5951/433.short
bman.reseq<-data.frame(lon=c(104.08, 104.08, 104.08, 104.08, 112.57, 114.29, 118.77, 104.08, 112.98, 118.77, 106.55), lat=c(30.66, 30.66, 30.66, 30.66, 37.87, 30.57, 32.05, 30.66, 28.20, 32.05, 29.55))
bman.clean<-rbind(bman.reseq, bman.use[,c("lon", "lat")])
bman.clean<-bman.clean[!duplicated(bman.clean),]
bman.ahull<-ahull(bman.clean$lon, bman.clean$lat, alpha=25); plot(bman.ahull)
bman.range<-ah2sp(bman.ahull); plot(bman.range)
proj4string(bman.range)<-CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
bman.range.clean<-gDifference(bman.range, oceans)
xlimits<-bman.range.clean@bbox[1,]
ylimits<-bman.range.clean@bbox[2,]
plot(wrld_simpl, axes=T, col="light yellow", xlim=xlimits, ylim=ylimits)
plot(bman.range.clean, add=T, lwd=2, axes=T, density=15, angle=45)
#project to Lambert Equal Area projection to calculate area in meters
bman.proj<-spTransform(bman.range.clean, CRS("+proj=cea +lon_0=0"))

#Caenorhabditis briggsae - no occurrence data in GBIF
#get range data from:
#http://labs.eeb.utoronto.ca/cutter/pdf/Cutter-etal_2010_MolEcol.pdf
#http://www.genetics.org/content/173/4/2021.full.pdf (lookup longitude of city listed in table)
cbrig.gps<-read.table("cbrig.txt", header=F)
cbri.clean<-cbrig.gps[,c("V2", "V1")]
names(cbri.clean) = c("lon", "lat")
cbri.clean = cbri.clean[!duplicated(cbri.clean),]
cbri.ahull<-ahull(cbri.clean$lon, cbri.clean$lat, alpha=170); plot(cbri.ahull)
cbri.range<-ah2sp(cbri.ahull); plot(cbri.range)
proj4string(cbri.range)<-CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
cbri.range.clean<-gDifference(cbri.range, oceans)
xlimits<-cbri.range.clean@bbox[1,]
ylimits<-cbri.range.clean@bbox[2,]
plot(wrld_simpl, axes=T, col="light yellow", xlim=xlimits, ylim=ylimits)
plot(cbri.range.clean, add=T, lwd=2, axes=T, density=15, angle=45)
#project to Lambert Equal Area projection to calculate area in meters
cbri.proj<-spTransform(cbri.range.clean, CRS("+proj=cea +lon_0=0"))

#Caenorhabditis elegans - no occurrence data in GBIF
#use data from: http://www.nature.com.ezp-prod1.hul.harvard.edu/ng/journal/v44/n3/abs/ng.1050.html#supplementary-information
cele.gps<-read.table("cele-1.txt", header=F)
cele.clean<-cele.gps[,c("V2", "V1")]
names(cele.clean)=c("lon", "lat")
cele.clean=cele.clean[complete.cases(cele.clean),]
cele.clean=cele.clean[!duplicated(cele.clean),]
cele.ahull<-ahull(cele.clean$lon, cele.clean$lat, alpha=170); plot(cele.ahull)
cele.range<-ah2sp(cele.ahull); plot(cele.range)
proj4string(cele.range)<-CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
cele.range.clean<-gDifference(cele.range, oceans)
xlimits<-cele.range.clean@bbox[1,]
ylimits<-cele.range.clean@bbox[2,]
plot(wrld_simpl, axes=T, col="light yellow", xlim=xlimits, ylim=ylimits)
plot(cele.range.clean, add=T, lwd=2, axes=T, density=15, angle=45)
#project to Lambert Equal Area projection to calculate area in meters
cele.proj<-spTransform(cele.range.clean, CRS("+proj=cea +lon_0=0"))

#Citrullus lanatus lanatus
clan.use=subset(clan.use, basisOfRecord != "unknown")
clan.use=subset(clan.use, !grepl("vulgaris", species))
clan.use=subset(clan.use, !grepl("Colocynthis", species))
clan.use=subset(clan.use, !grepl("citroides", species))
clan.use=subset(clan.use, !grepl("colocynthoides", species))
clan.use=subset(clan.use, !grepl("Citroides", species))
clan.use=subset(clan.use, !grepl("colocynthis", species))
#native range is southern Africa
clan.use=subset(clan.use, lat < -10)
clan.use=subset(clan.use, lon > 0)
clan.use=subset(clan.use, lon < 40)
clan.clean=clan.use[,c("lon", "lat")]
clan.clean<-clan.clean[!duplicated(clan.clean),]
clan.ahull<-ahull(clan.clean$lon, clan.clean$lat, alpha=10); plot(clan.ahull)
clan.range<-ah2sp(clan.ahull); plot(clan.range)
proj4string(clan.range)<-CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
clan.range.clean<-gDifference(clan.range, oceans)
xlimits<-clan.range.clean@bbox[1,]
ylimits<-clan.range.clean@bbox[2,]
plot(wrld_simpl, axes=T, col="light yellow", xlim=xlimits, ylim=ylimits)
plot(clan.range.clean, add=T, lwd=2, axes=T, density=15, angle=45)
#project to Lambert Equal Area projection to calculate area in meters
clan.proj<-spTransform(clan.range.clean, CRS("+proj=cea +lon_0=0"))

#Canis lupus
clup.use=subset(clup.use, basisOfRecord=="observation" | basisOfRecord=="specimen")
clup.use=clup.use[!grepl("familiaris", clup.use$species, ignore.case=T),]
clup.use=clup.use[!grepl("dingo", clup.use$species, ignore.case=T),]
clup.use=clup.use[!grepl("rufus", clup.use$species, ignore.case=T),]
clup.use=subset(clup.use, country!="Australia") #remove dingos coded as Canis lupus
clup.use=subset(clup.use, lat > 10) # remove a few other random miscoded locations
clup.clean<-clup.use[,c("lon", "lat")]
clup.clean<-clup.clean[!duplicated(clup.clean),]
clup.ahull<-ahull(clup.clean$lon, clup.clean$lat, alpha=100); plot(clup.ahull)
clup.range<-ah2sp(clup.ahull); plot(clup.range)
proj4string(clup.range)<-CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
clup.range.clean<-gDifference(clup.range, oceans)
xlimits<-clup.range.clean@bbox[1,]
ylimits<-clup.range.clean@bbox[2,]
plot(wrld_simpl, axes=T, col="light yellow", xlim=xlimits, ylim=ylimits)
plot(clup.range.clean, add=T, lwd=2, axes=T, density=15, angle=45)
#project to Lambert Equal Area projection to calculate area in meters
clup.proj<-spTransform(clup.range.clean, CRS("+proj=cea +lon_0=0"))

#Capsella rubella
crub.use=subset(crub.use, grepl("Capsella\\s+rubella", species))
crub.use=subset(crub.use, lat > 10 & lon > -30 & lon < 100) #native range is Mediterranean     
crub.use=subset(crub.use, basisOfRecord != "unknown")
#add data from http://www.pnas.org/content/106/13/5246.full
crub.reseq<-data.frame(lon=c(-5.58,23.43,24.42,19.48,20.51,3.00,-16.34,-16.19,-0.06,24.42,11.02,14.37,26.83), lat=c(36.15,37.58,35.29,39.40,39.40,39.30,28.19,28.19,42.53,35.29,43.28,40.62,37.78))
crub.clean<-rbind(crub.reseq, crub.use[,c("lon", "lat")])
crub.clean<-crub.clean[!duplicated(crub.clean),]
crub.ahull<-ahull(crub.clean$lon, crub.clean$lat, alpha=20); plot(crub.ahull)
crub.range<-ah2sp(crub.ahull); plot(crub.range)
proj4string(crub.range)<-CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
crub.range.clean<-gDifference(crub.range, oceans)
xlimits<-crub.range.clean@bbox[1,]
ylimits<-crub.range.clean@bbox[2,]
plot(wrld_simpl, axes=T, col="light yellow", xlim=xlimits, ylim=ylimits)
plot(crub.range.clean, add=T, lwd=2, axes=T, density=15, angle=45)
#project to Lambert Equal Area projection to calculate area in meters
crub.proj<-spTransform(crub.range.clean, CRS("+proj=cea +lon_0=0"))

#Cucumis sativus var. hardwickii
csat.use=subset(csat.use, basisOfRecord != "unknown")
#too few records of just Hardwickii, so remove all wrong subspecies
csat.use=subset(csat.use, !grepl("antasiaticus", species))
csat.use=subset(csat.use, !grepl("Cucumber", species))
csat.use=subset(csat.use, !grepl("Gherkin", species))
csat.use=subset(csat.use, !grepl("subsp. sativus", species))
csat.use=subset(csat.use, !grepl("europeaus", species))
csat.use=subset(csat.use, !grepl("Melo", species))
csat.use=subset(csat.use, !grepl("gracilior", species))
csat.use=subset(csat.use, !grepl("izmir", species))
csat.use=subset(csat.use, !grepl("sikkimensis", species))
#native range
csat.use=subset(csat.use, country=="India" | country=="Nepal")
csat.use=subset(csat.use, lon > 76)
csat.clean<-csat.use[,c("lon", "lat")]
csat.clean<-csat.clean[!duplicated(csat.clean),]
csat.ahull<-ahull(csat.clean$lon, csat.clean$lat, alpha=20); plot(csat.ahull)
csat.range<-ah2sp(csat.ahull); plot(csat.range)
proj4string(csat.range)<-CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
csat.range.clean<-gDifference(csat.range, oceans)
xlimits<-csat.range.clean@bbox[1,]
ylimits<-csat.range.clean@bbox[2,]
plot(wrld_simpl, axes=T, col="light yellow", xlim=xlimits, ylim=ylimits)
plot(csat.range.clean, add=T, lwd=2, axes=T, density=15, angle=45)
#project to Lambert Equal Area projection to calculate area in meters
csat.proj<-spTransform(csat.range.clean, CRS("+proj=cea +lon_0=0"))

#Drosophila melanogaster - GBIF not useful, nothing from Africa. Need to get GPS locations of DPGP collections
#using data from bioRxiv publication and John Pool, personal communication
dmel.gps<-read.table("dmel-merge.txt", header=F)
dmel.clean<-dmel.gps[,c("V2", "V1")]
dmel.clean<-dmel.use[!duplicated(dmel.clean),]
names(dmel.clean)=c("lon", "lat")
dmel.ahull<-ahull(dmel.clean$lon, dmel.clean$lat, alpha=50); plot(dmel.ahull)
dmel.range<-ah2sp(dmel.ahull); plot(dmel.range)
proj4string(dmel.range)<-CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
dmel.range.clean<-gDifference(dmel.range, oceans)
xlimits<-dmel.range.clean@bbox[1,]
ylimits<-dmel.range.clean@bbox[2,]
plot(wrld_simpl, axes=T, col="light yellow", xlim=xlimits, ylim=ylimits)
plot(dmel.range.clean, add=T, lwd=2, axes=T, density=15, angle=45)
#project to Lambert Equal Area projection to calculate area in meters
dmel.proj<-spTransform(dmel.range.clean, CRS("+proj=cea +lon_0=0"))

#Danio rerio
drer.use=subset(drer.use, basisOfRecord!="unknown")
drer.use=subset(drer.use, country != "United States")
drer.use=subset(drer.use, lat > 20)
drer.clean<-drer.use[,c("lon", "lat")]
drer.clean<-drer.clean[!duplicated(drer.clean),]
drer.ahull<-ahull(drer.clean$lon, drer.clean$lat, alpha=10); plot(drer.ahull)
drer.range<-ah2sp(drer.ahull); plot(drer.range)
proj4string(drer.range)<-CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
drer.range.clean<-gDifference(drer.range, oceans)
xlimits<-drer.range.clean@bbox[1,]
ylimits<-drer.range.clean@bbox[2,]
plot(wrld_simpl, axes=T, col="light yellow", xlim=xlimits, ylim=ylimits)
plot(drer.range.clean, add=T, lwd=2, axes=T, density=15, angle=45)
#project to Lambert Equal Area projection to calculate area in meters
drer.proj<-spTransform(drer.range.clean, CRS("+proj=cea +lon_0=0"))

#Equus ferus przewalskii
efer.use<-subset(efer.use, basisOfRecord != "unknown")
efer.use<-subset(efer.use, lon > -30)
efer.clean<-efer.use[,c("lon", "lat")]
efer.clean<-efer.clean[!duplicated(efer.clean),]
efer.ahull<-ahull(efer.clean$lon, efer.clean$lat, alpha=500); plot(efer.ahull, ylim=c(30,70), xlim=c(-10,120))
efer.range<-ah2sp(efer.ahull); plot(efer.range, ylim=c(30,70), xlim=c(-10,120), axes=T)
proj4string(efer.range)<-CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
efer.range.clean<-gDifference(efer.range, oceans)
xlimits<-efer.range.clean@bbox[1,]
ylimits<-efer.range.clean@bbox[2,]
plot(wrld_simpl, axes=T, col="light yellow", xlim=xlimits, ylim=ylimits)
plot(efer.range.clean, add=T, lwd=2, axes=T, density=15, angle=45)
#project to Lambert Equal Area projection to calculate area in meters
efer.proj<-spTransform(efer.range.clean, CRS("+proj=cea +lon_0=0"))

#Gasterosteus aculeatus
gacu.use<-gacu.gbif
gacu.use = subset(gacu.use, basisOfRecord != "fossil" & basisOfRecord != "unknown")
coordinates(gacu.use)=c("lon", "lat")
proj4string(gacu.use)<-CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
gacu.clean=gIntersection(gacu.use, oceans)
gacu.ahull<-ahull(gacu.clean@coords, alpha=25); plot(gacu.ahull)
gacu.range<-ah2sp(gacu.ahull); plot(gacu.range)
proj4string(gacu.range)<-CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
gacu.range.clean<-gIntersection(gacu.range, oceans)
xlimits<-gacu.range.clean@bbox[1,]
ylimits<-gacu.range.clean@bbox[2,]
plot(wrld_simpl, axes=T, col="light yellow", xlim=xlimits, ylim=ylimits)
plot(gacu.range.clean, add=T, lwd=2, axes=T, density=15, angle=45)
#project to Lambert Equal Area projection to calculate area in meters
gacu.proj<-spTransform(gacu.range.clean, CRS("+proj=cea +lon_0=0"))

#Glycine soja
gmax.use<-subset(gmax.use, grepl("soja", gmax.use$species))
gmax.use<-subset(gmax.use, basisOfRecord != "unknown")
gmax.use$lat2=gmax.use$lon
gmax.use$lon2=gmax.use$lat
gmax.use$lat[gmax.use$lat > 100] = gmax.use$lat2[gmax.use$lat > 100]
gmax.use$lon[gmax.use$lat > 100] = gmax.use$lon2[gmax.use$lat > 100]
gmax.use<-subset(gmax.use, lon > 100)
gmax.clean<-gmax.use[,c("lon", "lat")]
gmax.clean<-gmax.clean[!duplicated(gmax.clean),]
gmax.ahull<-ahull(gmax.clean$lon, gmax.clean$lat, alpha=50); plot(gmax.ahull)
gmax.range<-ah2sp(gmax.ahull); plot(gmax.range)
proj4string(gmax.range)<-CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
gmax.range.clean<-gDifference(gmax.range, oceans)
xlimits<-gmax.range.clean@bbox[1,]
ylimits<-gmax.range.clean@bbox[2,]
plot(wrld_simpl, axes=T, col="light yellow", xlim=xlimits, ylim=ylimits)
plot(gmax.range.clean, add=T, lwd=2, axes=T, density=15, angle=45)
#project to Lambert Equal Area projection to calculate area in meters
gmax.proj<-spTransform(gmax.range.clean, CRS("+proj=cea +lon_0=0"))

#Gossypium raimondii
#not enough geotagged locations, but can geoencode based on location descriptors
grai.geo=subset(grai.use, is.na(lon) & !is.na(locality))
grai.geo.locs<-as.data.frame(unlist(rbind(c(-7.168912, -78.927403), c(-6.842775, -78.029056), c(-7.221396, -78.839573), c(-7.817045, -79.175495), c(-7.232938, -78.837665), c(-6.840499, -77.991026))))
names(grai.geo.locs)=c("lat", "lon")
grai.clean<-grai.use[,c("lon", "lat")]
grai.clean<-grai.clean[!duplicated(grai.clean),]
grai.clean<-grai.clean[complete.cases(grai.clean),]
grai.clean<-rbind(grai.clean, grai.geo.locs)
grai.ahull<-ahull(grai.clean$lon, grai.clean$lat, alpha=3; plot(grai.ahull)
grai.range<-ah2sp(grai.ahull); plot(grai.range)
proj4string(grai.range)<-CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
grai.range.clean<-gDifference(grai.range, oceans)
xlimits<-grai.range.clean@bbox[1,]
ylimits<-grai.range.clean@bbox[2,]
plot(wrld_simpl, axes=T, col="light yellow", xlim=xlimits, ylim=ylimits)
plot(grai.range.clean, add=T, lwd=2, axes=T, density=15, angle=45)
#project to Lambert Equal Area projection to calculate area in meters
grai.proj<-spTransform(grai.range.clean, CRS("+proj=cea +lon_0=0"))
                                    
#Heliconius melpomene
hmel.use<-subset(hmel.use, basisOfRecord!="unknown")
hmel.use<-subset(hmel.use, species=="Heliconius melpomene melpomene Linnaeus, 1758" | species=="Heliconius melpomene  (Linnaeus 1758)" |  species=="Heliconius melpomene Linnaeus, 1758")
hmel.clean<-hmel.use[,c("lon", "lat")]
hmel.clean<-hmel.clean[!duplicated(hmel.clean),]
hmel.ahull<-ahull(hmel.clean$lon, hmel.clean$lat, alpha=20); plot(hmel.ahull)
hmel.range<-ah2sp(hmel.ahull); plot(hmel.range)
proj4string(hmel.range)<-CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
hmel.range.clean<-gDifference(hmel.range, oceans)
xlimits<-hmel.range.clean@bbox[1,]
ylimits<-hmel.range.clean@bbox[2,]
plot(wrld_simpl, axes=T, col="light yellow", xlim=xlimits, ylim=ylimits)
plot(hmel.range.clean, add=T, lwd=2, axes=T, density=15, angle=45)
#project to Lambert Equal Area projection to calculate area in meters
hmel.proj<-spTransform(hmel.range.clean, CRS("+proj=cea +lon_0=0"))

#Lepisosteus oculatus
locu.use<-subset(locu.use, basisOfRecord != "unknown" & species != "Lepisosteus productus")
locu.use<-subset(locu.use, lon < 0) #native range USA
locu.clean<-locu.use[,c("lon", "lat")]
locu.clean<-locu.clean[!duplicated(locu.clean),]
locu.ahull<-ahull(locu.clean$lon, locu.clean$lat, alpha=10); plot(locu.ahull)
locu.range<-ah2sp(locu.ahull); plot(locu.range)
proj4string(locu.range)<-CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
locu.range.clean<-gDifference(locu.range, oceans)
xlimits<-locu.range.clean@bbox[1,]
ylimits<-locu.range.clean@bbox[2,]
plot(wrld_simpl, axes=T, col="light yellow", xlim=xlimits, ylim=ylimits)
plot(locu.range.clean, add=T, lwd=2, axes=T, density=15, angle=45)
#project to Lambert Equal Area projection to calculate area in meters
locu.proj<-spTransform(locu.range.clean, CRS("+proj=cea +lon_0=0"))

#Macaca mulatta
mmul.use=subset(mmul.use, basisOfRecord == "observation" | basisOfRecord == "specimen")
mmul.use=subset(mmul.use, lon > 0) #native range SE Asia
mmul.clean<-mmul.use[,c("lon", "lat")]
mmul.clean<-mmul.clean[!duplicated(mmul.clean),]
mmul.ahull<-ahull(mmul.clean$lon, mmul.clean$lat, alpha=50); plot(mmul.ahull)
mmul.range<-ah2sp(mmul.ahull); plot(mmul.range)
proj4string(mmul.range)<-CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
mmul.range.clean<-gDifference(mmul.range, oceans)
xlimits<-mmul.range.clean@bbox[1,]
ylimits<-mmul.range.clean@bbox[2,]
plot(wrld_simpl, axes=T, col="light yellow", xlim=xlimits, ylim=ylimits)
plot(mmul.range.clean, add=T, lwd=2, axes=T, density=15, angle=45)
#project to Lambert Equal Area projection to calculate area in meters
mmul.proj<-spTransform(mmul.range.clean, CRS("+proj=cea +lon_0=0"))

#Mus musculus castaneus
mcas.use<-subset(mcas.use, species=="Mus musculus castaneus" & basisOfRecord != "unknown")
mcas.use<-subset(mcas.use, lon > 55) #native range SE Asia
mcas.clean<-mcas.use[,c("lon", "lat")]
mcas.clean<-mcas.clean[!duplicated(mcas.clean),]
mcas.ahull<-ahull(mcas.clean$lon, mcas.clean$lat, alpha=60); plot(mcas.ahull)
mcas.range<-ah2sp(mcas.ahull); plot(mcas.range)
proj4string(mcas.range)<-CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
mcas.range.clean<-gDifference(mcas.range, oceans)
xlimits<-mcas.range.clean@bbox[1,]
ylimits<-mcas.range.clean@bbox[2,]
plot(wrld_simpl, axes=T, col="light yellow", xlim=xlimits, ylim=ylimits)
plot(mcas.range.clean, add=T, lwd=2, axes=T, density=15, angle=45)
#project to Lambert Equal Area projection to calculate area in meters
mcas.proj<-spTransform(mcas.range.clean, CRS("+proj=cea +lon_0=0"))

#Medicago truncatula
#native range med basin; keep "unknown" records given failure of "known" records to capture full range
mtru.use<-subset(mtru.use, lon < 60 & lat > 20)
mtru.use<-subset(mtru.use, lat < 47 & lon < 45)
mtru.clean<-mtru.use[,c("lon", "lat")]
mtru.clean<-mtru.clean[!duplicated(mtru.clean),]
mtru.ahull<-ahull(mtru.clean$lon, mtru.clean$lat, alpha=50); plot(mtru.ahull)
mtru.range<-ah2sp(mtru.ahull); plot(mtru.range)
proj4string(mtru.range)<-CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
mtru.range.clean<-gDifference(mtru.range, oceans)
xlimits<-mtru.range.clean@bbox[1,]
ylimits<-mtru.range.clean@bbox[2,]
plot(wrld_simpl, axes=T, col="light yellow", xlim=xlimits, ylim=ylimits)
plot(mtru.range.clean, add=T, lwd=2, axes=T, density=15, angle=45)
#project to Lambert Equal Area projection to calculate area in meters
mtru.proj<-spTransform(mtru.range.clean, CRS("+proj=cea +lon_0=0"))

#Ovis canadensis
ocan.use=subset(ocan.use, basisOfRecord == "observation" | basisOfRecord == "specimen")
ocan.use=subset(ocan.use, lon < -90) #native to Western NA
ocan.clean<-ocan.use[,c("lon", "lat")]
ocan.clean<-ocan.clean[!duplicated(ocan.clean),]
ocan.ahull<-ahull(ocan.clean$lon, ocan.clean$lat, alpha=20); plot(ocan.ahull)
ocan.range<-ah2sp(ocan.ahull); plot(ocan.range)
proj4string(ocan.range)<-CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
ocan.range.clean<-gDifference(ocan.range, oceans)
xlimits<-ocan.range.clean@bbox[1,]
ylimits<-ocan.range.clean@bbox[2,]
plot(wrld_simpl, axes=T, col="light yellow", xlim=xlimits, ylim=ylimits)
plot(ocan.range.clean, add=T, lwd=2, axes=T, density=15, angle=45)
#project to Lambert Equal Area projection to calculate area in meters
ocan.proj<-spTransform(ocan.range.clean, CRS("+proj=cea +lon_0=0"))

#Oryzias latipes
olat.use<-subset(olat.use, basisOfRecord != "unknown")
olat.use<-subset(olat.use, lon > 75) #native to SE Asia
olat.use<-subset(olat.use, country=="Japan") #sample is from a Japanese subspecies
olat.clean<-olat.use[,c("lon", "lat")]
olat.clean<-olat.clean[!duplicated(olat.clean),]
olat.ahull<-ahull(olat.clean$lon, olat.clean$lat, alpha=8); plot(olat.ahull)
olat.range<-ah2sp(olat.ahull); plot(olat.range)
proj4string(olat.range)<-CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
olat.range.clean<-gDifference(olat.range, oceans)
xlimits<-olat.range.clean@bbox[1,]
ylimits<-olat.range.clean@bbox[2,]
plot(wrld_simpl, axes=T, col="light yellow", xlim=xlimits, ylim=ylimits)
plot(olat.range.clean, add=T, lwd=2, axes=T, density=15, angle=45)
#project to Lambert Equal Area projection to calculate area in meters
olat.proj<-spTransform(olat.range.clean, CRS("+proj=cea +lon_0=0"))

#Oryza rufipogon
oruf.use<-subset(oruf.use, grepl("ufipogon", species))
oruf.use<-subset(oruf.use, species != "oryza nivara & rufipogon")
#keep unknown basis records do to limited sampling
#native range SE Asia & N Australia
oruf.use<-subset(oruf.use, lon > 0)
oruf.clean<-oruf.use[,c("lon", "lat")]
oruf.clean<-oruf.clean[!duplicated(oruf.clean),]
oruf.ahull<-ahull(oruf.clean$lon, oruf.clean$lat, alpha=40); plot(oruf.ahull)
oruf.range<-ah2sp(oruf.ahull); plot(oruf.range)
proj4string(oruf.range)<-CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
oruf.range.clean<-gDifference(oruf.range, oceans)
xlimits<-oruf.range.clean@bbox[1,]
ylimits<-oruf.range.clean@bbox[2,]
plot(wrld_simpl, axes=T, col="light yellow", xlim=xlimits, ylim=ylimits)
plot(oruf.range.clean, add=T, lwd=2, axes=T, density=15, angle=45)
#project to Lambert Equal Area projection to calculate area in meters
oruf.proj<-spTransform(oruf.range.clean, CRS("+proj=cea +lon_0=0"))

#Papio anubis
panu.use = subset(panu.use, basisOfRecord != "fossil")
panu.clean<-panu.use[,c("lon", "lat")]
panu.clean<-panu.clean[!duplicated(panu.clean),]
panu.ahull<-ahull(panu.clean$lon, panu.clean$lat, alpha=40); plot(panu.ahull)
panu.range<-ah2sp(panu.ahull); plot(panu.range)
proj4string(panu.range)<-CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
panu.range.clean<-gDifference(panu.range, oceans)
xlimits<-panu.range.clean@bbox[1,]
ylimits<-panu.range.clean@bbox[2,]
plot(wrld_simpl, axes=T, col="light yellow", xlim=xlimits, ylim=ylimits)
plot(panu.range.clean, add=T, lwd=2, axes=T, density=15, angle=45)
#project to Lambert Equal Area projection to calculate area in meters
panu.proj<-spTransform(panu.range.clean, CRS("+proj=cea +lon_0=0"))

#Prunus davidiana
pdav.use<-subset(pdav.use, basisOfRecord != "fossil" & pdav.use$country != "Canada")
#not enough geotagged locations, but can geoencode based on location descriptors
pdav.geo<-subset(pdav.use, !is.na(locality) & is.na(lon))
pdav.geo.locs<-as.data.frame(unlist(rbind(c(37.977587, 112.410736), c(31.013188, 103.609965), c(25.931874, 111.009431), c(39.855821, 114.984394), c(32.600766, 103.617217))))
names(pdav.geo.locs)=c("lat", "lon")
pdav.clean<-pdav.use[,c("lon", "lat")]
pdav.clean<-pdav.clean[!duplicated(pdav.clean),]
pdav.clean<-pdav.clean[complete.cases(pdav.clean),]
pdav.clean<-rbind(pdav.clean, pdav.geo.locs)
pdav.ahull<-ahull(pdav.clean$lon, pdav.clean$lat, alpha=8); plot(pdav.ahull)
#ahull doesn't work for Pdav, use gConvexHull (points are too sparsely sampled over too wide a region)
coordinates(pdav.clean)=c("lon", "lat")
pdav.range<-gConvexHull(pdav.clean); plot(pdav.range)
proj4string(pdav.range)<-CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
pdav.range.clean<-gDifference(pdav.range, oceans)
xlimits<-pdav.range.clean@bbox[1,]
ylimits<-pdav.range.clean@bbox[2,]
plot(wrld_simpl, axes=T, col="light yellow", xlim=xlimits, ylim=ylimits)
plot(pdav.range.clean, add=T, lwd=2, axes=T, density=15, angle=45)
#project to Lambert Equal Area projection to calculate area in meters
pdav.proj<-spTransform(pdav.range.clean, CRS("+proj=cea +lon_0=0"))
                  
#Populus trichocarpa
ptri.use<-subset(ptri.use, lon < (-90)) #native range is western US
ptri.clean<-ptri.use[,c("lon", "lat")]
ptri.clean<-ptri.clean[!duplicated(ptri.clean),]
ptri.ahull<-ahull(ptri.clean$lon, ptri.clean$lat, alpha=40); plot(ptri.ahull)
ptri.range<-ah2sp(ptri.ahull); plot(ptri.range)
proj4string(ptri.range)<-CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
ptri.range.clean<-gDifference(ptri.range, oceans)
xlimits<-ptri.range.clean@bbox[1,]
ylimits<-ptri.range.clean@bbox[2,]
plot(wrld_simpl, axes=T, col="light yellow", xlim=xlimits, ylim=ylimits)
plot(ptri.range.clean, add=T, lwd=2, axes=T, density=15, angle=45)
#project to Lambert Equal Area projection to calculate area in meters
ptri.proj<-spTransform(ptri.range.clean, CRS("+proj=cea +lon_0=0"))

#Sorghum bicolor subsp. verticilliflorum
sbic.use<-subset(sbic.use, grepl("verticilliflorum", species, ignore.case=T) | grepl("aethiopicum", species, ignore.case=T) | grepl("arundinaceum", species, ignore.case=T))
sbic.use<-subset(sbic.use, lon > -20 & lon < 70)
sbic.clean<-sbic.use[,c("lon", "lat")]
sbic.clean<-sbic.clean[!duplicated(sbic.clean),]
sbic.clean<-sbic.clean[order(sbic.clean$lon),]
sbic.ahull<-ahull(sbic.clean$lon, sbic.clean$lat, alpha=20); plot(sbic.ahull, xlim=c(0,40))
sbic.range<-ah2sp(sbic.ahull); plot(sbic.range)
proj4string(sbic.range)<-CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
sbic.range.clean<-gDifference(sbic.range, oceans)
xlimits<-sbic.range.clean@bbox[1,]
ylimits<-sbic.range.clean@bbox[2,]
plot(wrld_simpl, axes=T, col="light yellow", xlim=xlimits, ylim=ylimits)
plot(sbic.range.clean, add=T, lwd=2, axes=T, density=15, angle=45)
#project to Lambert Equal Area projection to calculate area in meters
sbic.proj<-spTransform(sbic.range.clean, CRS("+proj=cea +lon_0=0"))

#Sus scrofa
sscr.use<-subset(sscr.use, species == "Sus scrofa" | species == "SUS SCROFA" | species == "Sus scrofa (Linnaeus, 1758)" | species == "Sus scrofa LINNAEUS" | species == "Sus scrofa Linnaeus, 1758" | species == "Sus scrofa scrofa" | species == "Sus scrofa scrofa (LINNAEUS, 1758)" | species == "Sus scrofa scrofa Linnaeus, 1758")
sscr.use<-subset(sscr.use, basisOfRecord != "fossil")
sscr.use<-subset(sscr.use, country != "Australia")
sscr.use<-subset(sscr.use, lon > -25) #remove introduced pops in new world
sscr.use<-subset(sscr.use, lat > -20) #remove introduced pop in NZ plus mislabelled Australia pops
sscr.use<-subset(sscr.use, country != "Sao Tome and Principe")
sscr.clean<-sscr.use[,c("lon", "lat")]
sscr.clean<-sscr.clean[!duplicated(sscr.clean),]
sscr.ahull<-ahull(sscr.clean$lon, sscr.clean$lat, alpha=60); plot(sscr.ahull)
sscr.range<-ah2sp(sscr.ahull); plot(sscr.range)
proj4string(sscr.range)<-CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
sscr.range.clean<-gDifference(sscr.range, oceans)
xlimits<-sscr.range.clean@bbox[1,]
ylimits<-sscr.range.clean@bbox[2,]
plot(wrld_simpl, axes=T, col="light yellow", xlim=xlimits, ylim=ylimits)
plot(sscr.range.clean, add=T, lwd=2, axes=T, density=15, angle=45)
#project to Lambert Equal Area projection to calculate area in meters
sscr.proj<-spTransform(sscr.range.clean, CRS("+proj=cea +lon_0=0"))

#Zea mays ssp parviglumis
zmay.use=subset(zmay.use, grepl("parviglumis", species, ignore.case=T))
zmay.use=zmay.use[zmay.use$lon<(-96),]
zmay.clean<-zmay.use[,c("lon", "lat")]
zmay.clean<-zmay.clean[!duplicated(zmay.clean),]
zmay.ahull<-ahull(zmay.clean$lon, zmay.clean$lat, alpha=6); plot(zmay.ahull)
zmay.range<-ah2sp(zmay.ahull); plot(zmay.range)
proj4string(zmay.range)<-CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
zmay.range.clean<-gDifference(zmay.range, oceans)
xlimits<-zmay.range.clean@bbox[1,]
ylimits<-zmay.range.clean@bbox[2,]
plot(wrld_simpl, axes=T, col="light yellow", xlim=xlimits, ylim=ylimits)
plot(zmay.range.clean, add=T, lwd=2, axes=T, density=15, angle=45)
#project to Lambert Equal Area projection to calculate area in meters
zmay.proj<-spTransform(zmay.range.clean, CRS("+proj=cea +lon_0=0"))

#ficAlb
falb.use<-subset(falb.use, species=="Ficedula albicollis (Temminck 1815)" | species=="Ficedula albicollis (Temminck, 1815)" | species=="Ficedula albicollis albicollis" | species=="Ficedula albicollis albicollis (Temminck, 1815)" | species==" Ficedula albicollis subsp. albicollis")
#just breeding population in Europe
falb.use<-subset(falb.use, basisOfRecord != "unknown" & lat > 35)
falb.clean<-falb.use[,c("lon", "lat")]
falb.clean<-falb.clean[!duplicated(falb.clean),]
falb.ahull<-ahull(falb.clean$lon, falb.clean$lat, alpha=10); plot(falb.ahull)
falb.range<-ah2sp(falb.ahull); plot(falb.range)
proj4string(falb.range)<-CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
falb.range.clean<-gDifference(falb.range, oceans)
xlimits<-falb.range.clean@bbox[1,]
ylimits<-falb.range.clean@bbox[2,]
plot(wrld_simpl, axes=T, col="light yellow", xlim=xlimits, ylim=ylimits)
plot(falb.range.clean, add=T, lwd=2, axes=T, density=15, angle=45)
#project to Lambert Equal Area projection to calculate area in meters
falb.proj<-spTransform(falb.range.clean, CRS("+proj=cea +lon_0=0"))

#make range dataset
range<-data.frame(species=character(0), range=numeric(0))
range<-rbind(range, data.frame(species="agam", area=gArea(agam.proj)/1e6))
range<-rbind(range, data.frame(species="amel", area=gArea(amel.proj)/1e6))
range<-rbind(range, data.frame(species="atha", area=gArea(atha.proj)/1e6))
range<-rbind(range, data.frame(species="bdis", area=gArea(bdis.proj)/1e6))
range<-rbind(range, data.frame(species="bmor", area=gArea(bman.proj)/1e6))
range<-rbind(range, data.frame(species="cbri", area=gArea(cbri.proj)/1e6))
range<-rbind(range, data.frame(species="cele", area=gArea(cele.proj)/1e6))
range<-rbind(range, data.frame(species="clan", area=gArea(clan.proj)/1e6))
range<-rbind(range, data.frame(species="clup", area=gArea(clup.proj)/1e6))
range<-rbind(range, data.frame(species="crub", area=gArea(crub.proj)/1e6))
range<-rbind(range, data.frame(species="csat", area=gArea(csat.proj)/1e6))
range<-rbind(range, data.frame(species="dmel", area=gArea(dmel.proj)/1e6))
range<-rbind(range, data.frame(species="drer", area=gArea(drer.proj)/1e6))
range<-rbind(range, data.frame(species="ecab", area=gArea(efer.proj)/1e6))
range<-rbind(range, data.frame(species="gacu", area=gArea(gacu.proj)/1e6))
range<-rbind(range, data.frame(species="gmax", area=gArea(gmax.proj)/1e6))
range<-rbind(range, data.frame(species="grai", area=gArea(grai.proj)/1e6))
range<-rbind(range, data.frame(species="hmel", area=gArea(hmel.proj)/1e6))
range<-rbind(range, data.frame(species="locu", area=gArea(locu.proj)/1e6))
range<-rbind(range, data.frame(species="mmus", area=gArea(mcas.proj)/1e6))
range<-rbind(range, data.frame(species="mmul", area=gArea(mmul.proj)/1e6))
range<-rbind(range, data.frame(species="mtru", area=gArea(mtru.proj)/1e6))
range<-rbind(range, data.frame(species="oari", area=gArea(ocan.proj)/1e6))
range<-rbind(range, data.frame(species="olat", area=gArea(olat.proj)/1e6))
range<-rbind(range, data.frame(species="osat", area=gArea(oruf.proj)/1e6))
range<-rbind(range, data.frame(species="panu", area=gArea(panu.proj)/1e6))
range<-rbind(range, data.frame(species="pper", area=gArea(pdav.proj)/1e6))
range<-rbind(range, data.frame(species="ptri", area=gArea(ptri.proj)/1e6))
range<-rbind(range, data.frame(species="sbic", area=gArea(sbic.proj)/1e6))
range<-rbind(range, data.frame(species="sscr", area=gArea(sscr.proj)/1e6))
range<-rbind(range, data.frame(species="zmay", area=gArea(zmay.proj)/1e6))
range<-rbind(range, data.frame(species="falb", area=gArea(falb.proj)/1e6))


#Domesticated species
#Bos taurus - no way to accurately get range from occurrence data, unfortunately
#however, the reads we used are all from heritage Japanese breeds
#so we assume the range is roughly equivalent to the area of Japan
#Bos taurus == area of Japan from Google
range<-rbind(range, data.frame(species="btau", area=3.778e11/1e6))

#Citrus reticulata - domesticated species, use information on native range from http://www.ars-grin.gov/cgi-bin/npgs/html/taxon.pl?10778
#we assume range is equal to area of Vietnam, from Google
range<-rbind(range, data.frame(species="ccle", area=331210000000/1e6))

#Gallus gallus - domesticated, need to think about this
#data is from two Chinese heritage breeds, the Silkie and the Tawianese L2
#use data for area of eastern/central parts of China, excluding Inner Mongolia, Tibet, and other desert/mountain areas
range<-rbind(range, data.frame(species="ggal", area=(832028+793300+1556061+1014354+2365900-1183000-1228400)))

#Meleagris gallopavo - domesticated, need to think
#turkey orginated from South Mexico population and that source wild stock is included in our sample; use the range of current wild turkeys in Mexico as a proxy
mgal.use<-subset(mgal.use, country=="Mexico")
mgal.use<-subset(mgal.use, lon < -20)
mgal.clean<-mgal.use[,c("lon", "lat")]
mgal.clean<-mgal.clean[!duplicated(mgal.clean),]
mgal.ahull<-ahull(mgal.clean$lon, mgal.clean$lat, alpha=10); plot(mgal.ahull)
mgal.range<-ah2sp(mgal.ahull); plot(mgal.range)
proj4string(mgal.range)<-CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
mgal.range.clean<-gDifference(mgal.range, oceans)
xlimits<-mgal.range.clean@bbox[1,]
ylimits<-mgal.range.clean@bbox[2,]
plot(wrld_simpl, axes=T, col="light yellow", xlim=xlimits, ylim=ylimits)
plot(mgal.range.clean, add=T, lwd=2, axes=T, density=15, angle=45)
#project to Lambert Equal Area projection to calculate area in meters
mgal.proj<-spTransform(mgal.range.clean, CRS("+proj=cea +lon_0=0"))
range<-rbind(range, data.frame(species="mgal", area=gArea(mgal.proj)/1e6))

#Setaria italica -> domesticated
#use Chinese landrace locations from sequencing paper
#Table S1: http://www.nature.com/ng/journal/v45/n8/full/ng.2673.html
sita.gps<-read.table("sita-gps.txt", header=F)
names(sita.gps)=c("lon", "lat")
sita.clean<-sita.gps[,c("lon", "lat")]
sita.clean<-sita.clean[!duplicated(sita.clean),]
sita.clean<-sita.clean[complete.cases(sita.clean),]
sita.ahull<-ahull(sita.clean$lon, sita.clean$lat, alpha=10); plot(sita.ahull)
sita.range<-ah2sp(sita.ahull); plot(sita.range)
proj4string(sita.range)<-CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
sita.range.clean<-gDifference(sita.range, oceans)
xlimits<-sita.range.clean@bbox[1,]
ylimits<-sita.range.clean@bbox[2,]
plot(wrld_simpl, axes=T, col="light yellow", xlim=xlimits, ylim=ylimits)
plot(sita.range.clean, add=T, lwd=2, axes=T, density=15, angle=45)
#project to Lambert Equal Area projection to calculate area in meters
sita.proj<-spTransform(sita.range.clean, CRS("+proj=cea +lon_0=0"))
range<-rbind(range, data.frame(species="sita", area=gArea(sita.proj)/1e6))

#Homo sapiens -> obviously troublesome to get range...?
#there is really no option but to approximate as area of sub-Saharan Africa excluding Madagascar
#calculate based on area of Africa - area of Sahara - area of Madagascar
range<-rbind(range, data.frame(species="hsap", area=(30.22e6-587040-9400000)))

#Cynoglossus semilaevis - nothing in GBIF, not even occurence data that can be geoencoded
#FAO data and FishBase indicate range as East China Sea and Yellow Sea
#use areas of East China Sea and Yellow Sea from Encylopedia Britannica
range<-rbind(range, data.frame(species="csem", area=(750000+380000)))

#Drosophila pseudoobscura - nothing in GBIF, spreadsheet from Steve not sufficient
#use geocoded locations from Dobzhansky, Theodosius, and Carl Epling. "Taxonomy, geographic distribution, and ecology of Drosophila pseudoobscura and its relatives."Carnegie Inst. Washington Publ 554 (1944): 1-46.
#also Steve Schaffer, personal communication
dpse.gps<-read.table("dpse-final.txt", header=F)
names(dpse.gps)=c("lat", "lon")
dpse.clean<-dpse.gps[,c("lon", "lat")]
dpse.clean<-dpse.clean[!duplicated(dpse.clean),]
dpse.clean<-dpse.clean[complete.cases(dpse.clean),]
dpse.ahull<-ahull(dpse.clean$lon, dpse.clean$lat, alpha=20); plot(dpse.ahull)
dpse.range<-ah2sp(dpse.ahull); plot(dpse.range)
proj4string(dpse.range)<-CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
dpse.range.clean<-gDifference(dpse.range, oceans)
xlimits<-dpse.range.clean@bbox[1,]
ylimits<-dpse.range.clean@bbox[2,]
plot(wrld_simpl, axes=T, col="light yellow", xlim=xlimits, ylim=ylimits)
plot(dpse.range.clean, add=T, lwd=2, axes=T, density=15, angle=45)
#project to Lambert Equal Area projection to calculate area in meters
dpse.proj<-spTransform(dpse.range.clean, CRS("+proj=cea +lon_0=0"))
range<-rbind(range, data.frame(species="dpse", area=gArea(dpse.proj)/1e6))

#write final range data
write.table(file="range.final", range, sep="\t", col.names=T, row.names=F, quote=F)

#load IUCN polygons; these can be downloaded from IUCN if desired, otherwise ignore this section
#IUCN shapefiles are not used in the manuscript
iucn.area<-data.frame(species=character(0), iucn.area=numeric(0))
sp3746<-readOGR(dsn="IUCN/species_3746/", layer="species_3746")
proj4string(sp3746)<-CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
iucn.area<-rbind(iucn.area, data.frame(species="cfam", iucn.area=gArea(spTransform(sp3746, CRS("+proj=cea +lon_0=0")))/1e6))

sp12554<-readOGR(dsn="IUCN/species_12554/", layer="species_12554")
proj4string(sp12554)<-CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
iucn.area<-rbind(iucn.area, data.frame(species="mmul", iucn.area=gArea(spTransform(sp12554, CRS("+proj=cea +lon_0=0")))/1e6))

sp15735<-readOGR(dsn="IUCN/species_15735/", layer="species_15735")
proj4string(sp15735)<-CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
iucn.area<-rbind(iucn.area, data.frame(species="oair", iucn.area=gArea(spTransform(sp15735, CRS("+proj=cea +lon_0=0")))/1e6))

sp40647<-readOGR(dsn="IUCN/species_40647/", layer="species_40647")
proj4string(sp40647)<-CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
iucn.area<-rbind(iucn.area, data.frame(species="panu", iucn.area=gArea(spTransform(sp40647, CRS("+proj=cea +lon_0=0")))/1e6))

sp41775<-readOGR(dsn="IUCN/species_41775/", layer="species_41775")
proj4string(sp41775)<-CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
iucn.area<-rbind(iucn.area, data.frame(species="sscr", iucn.area=gArea(spTransform(sp41775, CRS("+proj=cea +lon_0=0")))/1e6))


