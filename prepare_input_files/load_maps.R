#load map data into R and process to clean individual species
#final version, updated 11/6/2014

#function to parse notes field in GFF files, needed for some maps
parse_gff_notes <- function(x) {
  x=as.character(x)
  com<-strsplit(x, split=c(";"), fixed=T)
  cm<-sub(pattern="Note\\s+", "", com[[1]][2], perl=T)
  cm<-as.numeric(sub(pattern="\\s+cM.*", "", cm, perl=T))
  mk<-sub(pattern="GMap\\s+", "", com[[1]][1], perl=T)
  return(c(cm, mk))
}

#make map info file
map.info<-data.frame(sp=character(0), num.chr=numeric(0), map.length=numeric(0), markers=numeric(0));

#agambiae
agam.epcr<-read.table("/Volumes/LaCie/Projects/Current/ne/final_data/maps/agam/agam.primer.map.run2")
agam.epcr$mb = ceiling((agam.epcr$V4 + agam.epcr$V5)/2)/1000000
agam<-read.table("/Volumes/LaCie/Projects/Current/ne/final_data/maps/agam/agam_map.txt", header=T)
map.info=rbind(map.info, data.frame(sp="agam", num.chr=length(unique(agam$chr)), map.length=sum(aggregate(agam$cm, by=list(chr=agam$chr), max)$x), markers=length(unique(agam$markerID))))
agam = merge(agam, agam.epcr, by.x="markerID", by.y="V1")
agam=agam[order(agam$chr, agam$cm),]
agam$mb[agam$V2=="3L"] = agam$mb[agam$V2=="3L"]+(53200684/1000000)
agam$mb[agam$V2=="2L"] = agam$mb[agam$V2=="2L"]+(61545105/1000000)
agam$V2 = sub("L", "", agam$V2)
agam$V2 = sub("R", "", agam$V2)
agam=subset(agam, agam$chr==agam$V2)
agam=agam[,c(2,1,3,11)]
names(agam) = c("chr", "marker", "cm", "mb")
agam$sp = "agam"

#amellifera
amel.epcr<-read.table("/Volumes/LaCie/Projects/Current/ne/final_data/maps/amel/amel.primer.map")
amel.epcr$mb = ceiling((amel.epcr$V4+amel.epcr$V5)/2)/1000000
amel.epcr$acc = sub("gi\\|\\d+\\|gb\\|", "", amel.epcr$V2, perl=T)
amel.epcr$acc = sub("|", "", amel.epcr$acc, fixed=T)
amel.epcr=amel.epcr[,c(1,9,10)]
names(amel.epcr)[1] = "marker"
amel.chrkey = read.table("/Volumes/LaCie/Projects/Current/ne/final_data/maps/amel/chr2acc", header=F)
amel.epcr = merge(amel.chrkey, amel.epcr, by.x="V2", by.y="acc")
names(amel.epcr) = c("acc", "chr.seq", "marker", "mb")
amel<-read.table("/Volumes/LaCie/Projects/Current/ne/final_data/maps/amel/amel_map.txt", header=T)
amel$cm = amel$cm/2 #correct for lack of recombination in males
map.info=rbind(map.info, data.frame(sp="amel", num.chr=length(unique(amel$chr)), map.length=sum(aggregate(amel$cm, by=list(chr=amel$chr), max)$x), markers=length(unique(amel$marker))))
amel = merge(amel, amel.epcr, by="marker")
amel=subset(amel, amel$chr == amel$chr.seq)
amel=amel[,c(2,1,3,6)]
amel$sp = "amel"
names(amel) = c("chr", "marker", "cm", "mb", "sp")

#athaliana
atha.map<-read.table("/Volumes/LaCie/Projects/Current/ne/final_data/maps/atha/singer_map.txt", header=T)
map.info=rbind(map.info, data.frame(sp="atha", num.chr=length(unique(atha.map$chr)), map.length=sum(aggregate(atha.map$xm, by=list(chr=atha.map$chr), max)$x), markers=length(unique(atha.map$marker))))
atha.gene<-read.table("/Volumes/LaCie/Projects/Current/ne/final_data/maps/atha/singer_markers.txt", header=T)
atha.pos<-read.table("/Volumes/LaCie/Projects/Current/ne/final_data/maps/atha/TAIR9_AGI_gene.data.data", header=F)
names(atha.pos) = c("id", "gene", "start", "end", "st", "seq.chr")
atha.pos$mb = ceiling((atha.pos$start+atha.pos$end)/2)/1000000
atha.pos = atha.pos[,c("gene", "mb", "seq.chr")]
atha.pos$gene = sub("\\.\\d*", "", atha.pos$gene, perl=T)
atha.pos = atha.pos[order(atha.pos$gene, atha.pos$mb),]
atha.pos = atha.pos[!duplicated(atha.pos$gene),]
atha.gene$gene=toupper(atha.gene$gene)
atha.pos = merge(atha.pos, atha.gene, all.x=F, all.y=T)
atha = merge(atha.map, atha.pos, by.x="marker", by.y="marker")
atha = subset(atha, atha$chr.x == atha$seq.chr)
atha = atha[,c("chr.x", "marker", "xm", "mb")]
atha$sp = "atha"
names(atha) = c("chr", "marker", "cm", "mb", "sp")

#bmori
bmori.epcr<-read.table("/Volumes/LaCie/Projects/Current/ne/final_data/maps/bmor/bmori.primer.map.run2")
bmori.epcr$mb = ceiling((bmori.epcr$V4+bmori.epcr$V5)/2)/1000000
bmor<-read.table("/Volumes/LaCie/Projects/Current/ne/final_data/maps/bmor/bmori_map.txt", header=T)
bmor = subset(bmor, bmor$chr!=1)
bmor$cm = bmor$cm / 2 #correct for lack of recombination in females
map.info=rbind(map.info, data.frame(sp="bmor", num.chr=length(unique(bmor$chr)), map.length=sum(aggregate(bmor$cm, by=list(chr=bmor$chr), max)$x), markers=length(unique(bmor$marker))))
bmor$marker = sub("^0", "", bmor$marker, perl=T)
bmor<-merge(bmor, bmori.epcr, by.x="marker", by.y="V1")
bmor$V2 = sub("chr", "", bmor$V2)
bmor=subset(bmor, bmor$chr==bmor$V2)
bmor=bmor[,c(2,1,3,11)]
names(bmor) = c("chr", "marker", "cm", "mb")
bmor$sp = "bmor"

#btaurus
btau<-read.table("/Volumes/LaCie/Projects/Current/ne/final_data/maps/btau/btau_map.txt", header=T, sep="\t") #skip lines that don't have 2 elements
map.info=rbind(map.info, data.frame(sp="btau", num.chr=length(unique(btau$chr)), map.length=sum(aggregate(btau$cm, by=list(chr=btau$chr), max)$x), markers=length(unique(btau$marker))))
btau.pos<-read.table("/Volumes/LaCie/Projects/Current/ne/final_data/maps/btau/hglft_genome_4491_23df20.bed.bed")
btau = merge(btau, btau.pos, by.x="marker", by.y="V4")
btau$mb = ceiling((btau$V2+btau$V3)/2)/1000000
btau$V1 = sub("chr", "", btau$V1)
btau = subset(btau, btau$chr == btau$V1)
btau = btau[,c(2,1,4,10)]
btau$sp = "btau"

#bdis
bdis<-read.table("/Volumes/LaCie/Projects/Current/ne/final_data/maps/bdis/bdis_map.txt", header=T, colClasses = c("character", "character", "numeric", "numeric", "numeric", "numeric"))
map.info=rbind(map.info, data.frame(sp="bdis", num.chr=length(unique(bdis$chromosome)), map.length=sum(aggregate(bdis$cm, by=list(chr=bdis$chromosome), max)$x), markers=length(unique(bdis$Locus))))
bdis$chromosome = sub("chr_", "", bdis$chromosome)
bdis = bdis[,c(2,1,6,5)]
bdis$sp = "bdis"
names(bdis) = c("chr", "marker", "cm", "mb", "sp")
bdis$mb = bdis$mb/1000000

#cbriggsae
cbri<-read.table("/Volumes/LaCie/Projects/Current/ne/final_data/maps/cbri/cbrig_map.txt", header=T)
cbri$cm = cbri$M * 100
cbri$mb = cbri$bp / 1000000
cbri$sp = "cbri"
cbri = cbri[,c(3,1,5,6,7)]
cbri = subset(cbri, cbri$chr != 6)

map.info=rbind(map.info, data.frame(sp="cbri", num.chr=length(unique(cbri$chr)), map.length=sum(aggregate(cbri$cm, by=list(chr=cbri$chr), max)$x), markers=length(unique(cbri$marker))))

#celegans
cele<-read.table("/Volumes/LaCie/Projects/Current/ne/final_data/maps/cele/c_elegans.current.genetic_limits.gff2", header=F, sep="\t", colClasses=c("character"))
cele = subset(cele, cele$V2=="absolute_pmap_position")
cele = cele[,c(1,4,5,9)]
cele$V4 = as.numeric(cele$V4)
cele$V5 = as.numeric(cele$V5)
cele$V9 = as.character(cele$V9)
cele = droplevels(subset(cele, cele$V1 != "X"))
cele$V1 = as.roman(cele$V1)
cele$chr = as.numeric(cele$V1)

cele$cm = apply(cele[,4, drop=F], 1, parse_gff_notes)[1,]
cele$marker = apply(cele[,4,drop=F], 1, parse_gff_notes)[2,]
cele$mb = ceiling((cele$V4+cele$V5)/2)/1000000
cele=cele[,c("chr", "marker", "cm", "mb")]
cele$sp = "cele"
cele$cm = as.numeric(as.character(cele$cm))

cele$cm[cele$chr==1] = cele$cm[cele$chr==1] + abs(min(cele$cm[cele$chr==1]))
cele$cm[cele$chr==2] = cele$cm[cele$chr==2] + abs(min(cele$cm[cele$chr==2]))
cele$cm[cele$chr==3] = cele$cm[cele$chr==3] + abs(min(cele$cm[cele$chr==3]))
cele$cm[cele$chr==4] = cele$cm[cele$chr==4] + abs(min(cele$cm[cele$chr==4]))
cele$cm[cele$chr==5] = cele$cm[cele$chr==5] + abs(min(cele$cm[cele$chr==5]))

map.info=rbind(map.info, data.frame(sp="cele", num.chr=length(unique(cele$chr)), map.length=sum(aggregate(cele$cm, by=list(chr=cele$chr), max)$x), markers=length(unique(cele$marker))))

#clupus
cf.gm<-read.table("/Volumes/LaCie/Projects/Current/ne/final_data/maps/clup/cm_positions_sexavg.txt", header=F)
cf.sm<-read.table("/Volumes/LaCie/Projects/Current/ne/final_data/maps/clup/marker_positions_canfam3.1.txt", header=F, fill=T)
names(cf.gm) = c("chr", "marker", "cm")
map.info=rbind(map.info, data.frame(sp="clup", num.chr=length(unique(cf.gm$chr)), map.length=sum(aggregate(cf.gm$cm, by=list(chr=cf.gm$chr), max)$x), markers=length(unique(cf.gm$marker))))
names(cf.sm) = c("marker", "seq.chr", "mb", "end")
clup<-merge(cf.sm, cf.gm, all.y=T)
clup$mb = clup$mb/1000000
clup$chr.agree = ifelse(clup$seq.chr == clup$chr, T, F)
clup=subset(clup, clup$chr.agree==T)
clup=clup[,c(5,1,6,3)]
clup$sp = "clup"

#capsella
crub<-read.table("/Volumes/LaCie/Projects/Current/ne/final_data/maps/crub/capsella_final.txt", header=T)
crub$marker = paste(crub$chr, crub$bp, sep="_")
crub$mb = crub$bp / 1000000
crub = crub[,c(1,4,3,5)]
crub$sp = "crub"
crub$chr = sub("scaffold_", "", crub$chr)
crub$chr = as.integer(crub$chr)
map.info=rbind(map.info, data.frame(sp="crub", num.chr=length(unique(crub$chr)), map.length=sum(aggregate(crub$cm, by=list(chr=crub$chr), max)$x), markers=length(unique(crub$marker))))

#clan
clan.pos<-read.table("/Volumes/LaCie/Projects/Current/ne/final_data/maps/clan/watermelon_scafs.txt", header=T)
clan.map<-read.table("/Volumes/LaCie/Projects/Current/ne/final_data/maps/clan/watermelon_map.txt", header=T)
clan.map$pos = ceiling((clan.map$start + clan.map$end)/2)

clan.pos$offset = ifelse(clan.pos$oriented=="R", clan.pos$end+1, clan.pos$start-1)
clan.pos$sign = ifelse(clan.pos$oriented=="R", -1, 1)
clan = merge(clan.map, clan.pos, by.x="scaf", by.y="scaffold", all.x=T, all.y=F)
clan$mb = clan$pos*clan$sign + clan$offset
clan = subset(clan, select=c("chromosome", "marker", "cm", "mb"))
clan$chromosome = sub("Chr", "", clan$chromosome)
clan$sp = "clan"
clan$mb = clan$mb / 1000000
names(clan) = c("chr", "marker", "cm", "mb", "sp")
map.info=rbind(map.info, data.frame(sp="clan", num.chr=length(unique(clan$chr)), map.length=sum(aggregate(clan$cm, by=list(chr=clan$chr), max)$x), markers=length(unique(clan$marker))))

#citrus - clementine
ccle.map<-read.table("/Volumes/LaCie/Projects/Current/ne/final_data/maps/ccle/clementine_avg_map.txt", header=T, na.strings="ABS")
map.info=rbind(map.info, data.frame(sp="ccle", num.chr=length(unique(ccle.map$chr)), map.length=sum(aggregate(ccle.map$cm, by=list(chr=ccle.map$chr), max)$x), markers=length(unique(ccle.map$marker))))
ccle.blast<-read.table("/Volumes/LaCie/Projects/Current/ne/final_data/maps/ccle/marker_ccle_forR.out", header=F)
ccle.blast = ccle.blast[,c(5,6,13,14,15)]
ccle.blast=ccle.blast[order(ccle.blast$V5, ccle.blast$V15),]
ccle.blast=ccle.blast[!duplicated(ccle.blast$V5),] #keep smallest e-value hit only
ccle<-merge(ccle.map, ccle.blast, by.x="marker", by.y="V5")
ccle$V6 = sub("scaffold_", "", ccle$V6)
ccle = subset(ccle, ccle$chr == ccle$V6, select=c("chr", "marker", "cm", "V13"))
names(ccle)[4] = "mb"
ccle$mb = ccle$mb/1000000
ccle$sp = "ccle"

#csat
csat.map<-read.table("/Volumes/LaCie/Projects/Current/ne/final_data/maps/csat/cuke_map.txt", header=F)
csat.map$V2 = sub("Chr", "", csat.map$V2, fixed=T)
names(csat.map) = c("marker", "chr", "cm")
map.info=rbind(map.info, data.frame(sp="csat", num.chr=length(unique(csat.map$chr)), map.length=sum(aggregate(csat.map$cm, by=list(chr=csat.map$chr), max)$x), markers=length(unique(csat.map$marker))))

csat.epcr<-read.table("/Volumes/LaCie/Projects/Current/ne/final_data/maps/csat/cuke.primer.map", header=F, comment.char="#")
csat.epcr = csat.epcr[,c(1,2,4,5)]
csat.epcr$mb = ceiling((csat.epcr$V4 + csat.epcr$V5)/2)/1000000
csat.epcr = csat.epcr[,c(1,2,5)]
names(csat.epcr) = c("marker", "seq.chr", "mb")
csat = merge(csat.map, csat.epcr)
csat$seq.chr = sub("Chr", "", csat$seq.chr, fixed=T)
csat = subset(csat, csat$seq.chr == csat$chr)
csat = csat[,c("chr", "marker", "cm", "mb")]
csat$sp = "csat"

#danio rerio
drer<-read.csv("/Volumes/LaCie/Projects/Current/ne/final_data/maps/drer/TableS1.csv", header=T)
map.info=rbind(map.info, data.frame(sp="drer", num.chr=length(unique(drer$Chr)), map.length=sum(aggregate(drer$cM.position, by=list(chr=drer$Chr), max)$x), markers=length(unique(drer$marker))))
drer = subset(drer, !is.na(drer$Zv9.chr))
drer = subset(drer, drer$Chr == drer$Zv9.chr)
drer = drer[,c(3,1,4,6)]
names(drer) = c("chr", "marker", "cm", "mb")
drer$mb = drer$mb / 1000000
drer$sp = 'drer'

#dmel
dmel<-read.table("/Volumes/LaCie/Projects/Current/ne/final_data/maps/dmel/comeron_map.txt", header=T)
dmel$marker = row.names(dmel)
dmel$mb = dmel$pos / 1000000
dmel=dmel[,c(1,5,4,6)]
dmel$sp = "dmel"
dmel$chr = sub("2L", 3, dmel$chr) #muller elements, coded as A=4, B=3, C=7, D=5, E=8, F=6
dmel$chr = sub("2R", 7, dmel$chr) 
dmel$chr = sub("3L", 5, dmel$chr)
dmel$chr = sub("3R", 8, dmel$chr)
dmel$cm = dmel$cm / 2 # correct of lack of recombination in males
map.info=rbind(map.info, data.frame(sp="dmel", num.chr=length(unique(dmel$chr)), map.length=sum(aggregate(dmel$cm, by=list(chr=dmel$chr), max)$x), markers=length(unique(dmel$marker))))

#pseudo
dpse<-read.table("/Volumes/LaCie/Projects/Current/ne/final_data/maps/dpse/dpse_map.txt", header=T)
dpse$chr = 2;
dpse$mb = dpse$bp / 1000000
dpse$sp = "dpse"
dpse = dpse[,c("chr", "marker", "cm", "mb", "sp")]
dpse$cm = dpse$cm / 2 #correct for lack of recombination in males
map.info=rbind(map.info, data.frame(sp="dpse", num.chr=length(unique(dpse$chr)), map.length=sum(aggregate(dpse$cm, by=list(chr=dpse$chr), max)$x), markers=length(unique(dpse$marker))))

#ecab
ecab<-read.table("/Volumes/LaCie/Projects/Current/ne/final_data/maps/ecab/fixed_horse_maps.txt", header=T)
ecab$mb = ecab$mb/1000000
map.info=rbind(map.info, data.frame(sp="ecab", num.chr=length(unique(ecab$chr)), map.length=sum(aggregate(ecab$cm, by=list(chr=ecab$chr), max)$x), markers=length(unique(ecab$marker))))
ecab.epcr = read.table("/Volumes/LaCie/Projects/Current/ne/final_data/maps/ecab/ecab.primer.map.run2", header=F)
ecab.epcr = ecab.epcr[,c(1,2,4,5)]
ecab.epcr$mb = ceiling((ecab.epcr$V4 + ecab.epcr$V5)/2)/1000000
ecab.epcr = ecab.epcr[,c(1,2,5)]
names(ecab.epcr) = c("marker", "chr.seq", "mb")
ecab = merge(ecab, ecab.epcr, by="marker", all=T)
ecab$mb = ecab$mb.x
ecab$mb[is.na(ecab$mb)] = ecab$mb.y[is.na(ecab$mb)]
ecab$seq.chr = as.character(ecab$chr.seq.x)
ecab$seq.chr[is.na(ecab$chr.seq.x)] = as.character(ecab$chr.seq.y[is.na(ecab$chr.seq.x)])
ecab = ecab[,c("chr", "marker", "cm", "mb", "seq.chr")]
ecab = subset(ecab, !is.na(ecab$mb) & ecab$seq.chr == ecab$chr)
ecab$cm[ecab$chr==25] = abs(ecab$cm[ecab$chr==25] - 48.7)
ecab$sp = "ecab"
ecab = ecab[,-5]

#ggal
ggal<-read.table("/Volumes/LaCie/Projects/Current/ne/final_data/maps/ggal/ggal_map.txt", header=T, sep="\t")
ggal<-ggal[,c(1,2,5,6)]
names(ggal) = c("chr", "marker", "cm", "mb")
ggal=droplevels(subset(ggal, chr != "Z"))
map.info=rbind(map.info, data.frame(sp="ggal", num.chr=length(unique(ggal$chr)), map.length=sum(aggregate(ggal$cm, by=list(chr=ggal$chr), max)$x), markers=length(unique(ggal$marker))))

ggal$sp = "ggal"
ggal = ggal[,c("marker", "cm", "sp")]
ggal.pos = read.table("/Volumes/LaCie/Projects/Current/ne/final_data/maps/ggal/hglft_genome_6f8e_23e4c0.bed.bed", header=F)
ggal.pos = ggal.pos[,c(4,1,3)]
ggal = merge(ggal, ggal.pos, by.x="marker", by.y="V4")
ggal = ggal[,c(4,1,2,5,3)]
names(ggal) = c("chr", "marker", "cm", "mb", "sp")
ggal$mb = ggal$mb/1000000
ggal$chr = sub("chr", "", ggal$chr)
ggal = subset(ggal, !grepl("random", ggal$chr))
ggal = subset(ggal, !grepl("Un", ggal$chr))

#gacu
gacu <- read.table("/Volumes/LaCie/Projects/Current/ne/final_data/maps/gacu/mec12322-sup-0004-AppendixS4.txt", header=T)
gacu = gacu[,c(5,2,7,6)]
gacu$sp = "gacu"
names(gacu) = c("chr", "marker", "cm", "mb", "sp")
gacu$mb = gacu$mb / 1000000
map.info=rbind(map.info, data.frame(sp="gacu", num.chr=length(unique(gacu$chr)), map.length=sum(aggregate(gacu$cm, by=list(chr=gacu$chr), max)$x), markers=length(unique(gacu$marker))))

#glycine max
gmax.map<-read.table("/Volumes/LaCie/Projects/Current/ne/final_data/maps/gmax/Soybean-GmConsensus40.txt", header=T, sep="\t")
gmax.map$cm = (gmax.map$feature_start + gmax.map$feature_stop)/2
map.info=rbind(map.info, data.frame(sp="gmax", num.chr=length(unique(gmax.map$map_name)), map.length=sum(aggregate(gmax.map$cm, by=list(chr=gmax.map$map_name), max)$x), markers=length(unique(gmax.map$feature_acc))))
gmax.mark<-read.table("/Volumes/LaCie/Projects/Current/ne/final_data/maps/gmax/marker_all.txt", header=T, sep="\t")
gmax.mark$Chromosome = sub("Gm0", "", gmax.mark$Chromosome, fixed=T)
gmax.mark$Chromosome = sub("Gm", "", gmax.mark$Chromosome, fixed=T)
gmax<-merge(gmax.map, gmax.mark, by.x="feature_name", by.y="Name")
gmax$mb = ceiling((gmax$Start+gmax$Stop)/2)/1000000
gmax=gmax[,c("Chromosome", "feature_name", "cm", "mb")]
gmax$sp = "gmax"
names(gmax) = c("chr", "marker", "cm", "mb", "sp")

#grai
grai.blast<-read.table("/Volumes/LaCie/Projects/Current/ne/final_data/maps/grai/marker_grai.out", header=F, sep="\t")
grai.blast = grai.blast[order(grai.blast$V1, grai.blast$V11),]
blast.ranks = as.data.frame(unlist(by(grai.blast$V11, grai.blast$V1, rank)))
names(blast.ranks) = c("rank")
grai.blast = cbind(grai.blast, blast.ranks)
grai.blast = subset(grai.blast, grai.blast$rank <= 8)


grai.blast$V2 = sub("Chr0", "", grai.blast$V2)
grai.blast$V2 = sub("Chr", "", grai.blast$V2)
grai.blast$mb = ceiling((grai.blast$V9+grai.blast$V10)/2)/1000000
grai.blast=grai.blast[,c(1,2,11,14)]
names(grai.blast)=c("marker", "seq.chr", "eval", "mb")

grai<-read.table("/Volumes/LaCie/Projects/Current/ne/final_data/maps/grai/cottonmap_initial.txt", header=F, fill=T)
grai = grai[,c(2,1,3)]
grai$V2 = sub("D", "", grai$V2)
names(grai) = c("chr", "marker", "cm")
map.info=rbind(map.info, data.frame(sp="grai", num.chr=length(unique(grai$chr)), map.length=sum(aggregate(grai$cm, by=list(chr=grai$chr), max)$x), markers=length(unique(grai$marker))))

grai = merge(grai, grai.blast, by="marker")
grai = subset(grai, grai$chr == grai$seq.chr)
grai$sp = "grai"
grai = grai[,c("chr", "marker", "cm", "mb", "sp")]

grai$cm[grai$chr==1] = abs(grai$cm[grai$chr==1] - max(grai$cm[grai$chr==1]))
grai$cm[grai$chr==2] = abs(grai$cm[grai$chr==2] - max(grai$cm[grai$chr==2]))
grai$cm[grai$chr==3] = abs(grai$cm[grai$chr==3] - max(grai$cm[grai$chr==3]))
grai$cm[grai$chr==4] = abs(grai$cm[grai$chr==4] - max(grai$cm[grai$chr==4]))
grai$cm[grai$chr==7] = abs(grai$cm[grai$chr==7] - max(grai$cm[grai$chr==7]))
grai$cm[grai$chr==8] = abs(grai$cm[grai$chr==8] - max(grai$cm[grai$chr==8]))
grai$cm[grai$chr==12] = abs(grai$cm[grai$chr==12] - max(grai$cm[grai$chr==12]))
grai$cm[grai$chr==13] = abs(grai$cm[grai$chr==13] - max(grai$cm[grai$chr==13]))

#hmel
hmel.map<-read.table("/Volumes/LaCie/Projects/Current/ne/final_data/maps/hmel/hmel_map.txt", header=F)
names(hmel.map) = c("chr", "marker", "cm")
hmel.map$cm = hmel.map$cm / 2 #correct for lack of recombination in males
map.info=rbind(map.info, data.frame(sp="hmel", num.chr=length(unique(hmel.map$chr)), map.length=sum(aggregate(hmel.map$cm, by=list(chr=hmel.map$chr), max)$x), markers=length(unique(hmel.map$marker))))
hmel.pos<-read.table("/Volumes/LaCie/Projects/Current/ne/final_data/maps/hmel/final_marker_positions.txt", header=F)
names(hmel.pos) = c("marker", "chr", "bp")
hmel.pos$chr = sub("chr", "", hmel.pos$chr, fixed=T)
hmel.max = aggregate(bp ~ marker + chr, hmel.pos, max)
names(hmel.max) = c("marker", "chr", "max")
hmel.min = aggregate(bp ~ marker + chr, hmel.pos, min)
names(hmel.min) = c("marker", "chr", "min")
hmel.marker = merge(hmel.max, hmel.min)
hmel.marker$mb = ((hmel.marker$max + hmel.marker$min)/2)/1000000
hmel.marker = hmel.marker[,c("marker", "chr", "mb")]
names(hmel.marker) = c("marker", "seq.chr", "mb")
hmel = merge(hmel.map, hmel.marker, by="marker", all=T)
hmel = hmel[,c("seq.chr", "marker", "cm", "mb")]
names(hmel)[1] = "chr"
hmel$sp = "hmel"
hmel = hmel[complete.cases(hmel),]
hmel$marker = as.factor(hmel$marker)

#human
hsap.m<-read.table("/Volumes/LaCie/Projects/Current/ne/final_data/maps/hsap/male.gmap.gmap", header=T)
hsap.f<-read.table("/Volumes/LaCie/Projects/Current/ne/final_data/maps/hsap/female.gmap.gmap", header=T)
hsap<-merge(hsap.m, hsap.f, by="snp")
hsap$cm.int = (hsap$cM.x + hsap$cM.y) / 2
hsap$chr = hsap$chr.x
hsap$mb = hsap$pos.x / 1000000
hsap = hsap[,c("chr", "snp", "cm.int", "mb")]
hsap = subset(hsap, !is.na(hsap$cm))
hsap$sp = "hsap"
hsap = hsap[order(hsap$chr, hsap$mb),]
hsap$cm = unlist(tapply(hsap$cm.int, hsap$chr, cumsum))
hsap$chr = sub("chr", "", hsap$chr, fixed=T)
hsap = hsap[,-3]
names(hsap) = c("chr", "marker", "mb", "sp", "cm")
hsap = hsap[,c("chr", "marker", "cm", "mb", "sp")]
hsap$marker = as.character(hsap$marker)
map.info=rbind(map.info, data.frame(sp="hsap", num.chr=length(unique(hsap$chr)), map.length=sum(aggregate(hsap$cm, by=list(chr=hsap$chr), max)$x), markers=length(unique(hsap$marker))))
#update positions
hsap.newpos<-read.table("/Volumes/LaCie/Projects/Current/ne/final_data/maps/hsap/hglft_genome_750e_55bf00.bed", head=F)
hsap.newpos$mb = hsap.newpos$V3/1000000
hsap.newpos$chr = sub("chr", "", hsap.newpos$V1, fixed=T)
hsap.newpos = hsap.newpos[,c(4,5,6)]
names(hsap.newpos)[1] = "marker"
hsap = merge(hsap, hsap.newpos, by="marker")
hsap = subset(hsap, hsap$chr.x==hsap$chr.y)
hsap = subset(hsap, !is.na(mb.y))
hsap = hsap[,c("chr.y", "marker", "cm", "mb.y", "sp")]
names(hsap) = c("chr", "marker", "cm", "mb", "sp")

#mmul
mmul.map<-read.table("/Volumes/LaCie/Projects/Current/ne/final_data/maps/mmul/mmul_map.txt", header=T)
map.info=rbind(map.info, data.frame(sp="mmul", num.chr=length(unique(mmul.map$chr)), map.length=sum(aggregate(mmul.map$cm, by=list(chr=mmul.map$chr), max)$x), markers=length(unique(mmul.map$marker))))

mmul.epcr<-read.table("/Volumes/LaCie/Projects/Current/ne/final_data/maps/mmul/mmul.primer.map", header=F)
names(mmul.epcr) = c("marker", "acc", "st", "start", "end", "mm", "gap", "len")
mmul.epcr$mb = ((mmul.epcr$start+mmul.epcr$end)/2)/1000000
mmul<-merge(mmul.map, mmul.epcr, by.x="marker", by.y="marker", all.x=T)
mmul.chrkey = read.table("/Volumes/LaCie/Projects/Current/ne/final_data/maps/mmul/chr_accessions_Mmul_051212")
mmul.chrkey=mmul.chrkey[,c(1,4)]
names(mmul.chrkey) = c("seq.chr", "acc")
mmul$acc = sub("gi\\|\\d+\\|gb\\|", "", mmul$acc, perl=T)
mmul$acc = sub("\\|", "", mmul$acc, perl=T)
mmul = merge(mmul, mmul.chrkey, all.x=T, by="acc")
mmul = subset(mmul, mmul$chr == mmul$seq.chr)
mmul = mmul[,c(3,2,4,11)]

#fix reverse chrs
mmul$cm[mmul$chr==2] = abs(mmul$cm[mmul$chr==2] - max(mmul$cm[mmul$chr==2]))
mmul$cm[mmul$chr==5] = abs(mmul$cm[mmul$chr==5] - max(mmul$cm[mmul$chr==5]))
mmul$cm[mmul$chr==7] = abs(mmul$cm[mmul$chr==7] - max(mmul$cm[mmul$chr==7]))
mmul$cm[mmul$chr==11] = abs(mmul$cm[mmul$chr==11] - max(mmul$cm[mmul$chr==11]))
mmul$cm[mmul$chr==14] = abs(mmul$cm[mmul$chr==14] - max(mmul$cm[mmul$chr==14]))
mmul$cm[mmul$chr==15] = abs(mmul$cm[mmul$chr==15] - max(mmul$cm[mmul$chr==15]))
mmul$cm[mmul$chr==17] = abs(mmul$cm[mmul$chr==17] - max(mmul$cm[mmul$chr==17]))

mmul=mmul[order(mmul$chr, mmul$cm),]
mmul$sp = "mmul"

#medicago
mtru.map<-read.table("/Volumes/LaCie/Projects/Current/ne/final_data/maps/mtru/MtYoungUMinn2006_cmap_LIS.txt", header=T, sep="\t", fill=F)
map.info=rbind(map.info, data.frame(sp="mtru", num.chr=length(unique(mtru.map$map_name)), map.length=sum(aggregate(mtru.map$feature_start, by=list(chr=mtru.map$map_name), max)$x), markers=length(unique(mtru.map$feature_acc))))
mtru.epcr = read.table("/Volumes/LaCie/Projects/Current/ne/final_data/maps/mtru//mtru.primer.map", header=F)

mtru = merge(mtru.map, mtru.epcr, by.x="feature_name", by.y="V1")
mtru$chr = sub("chr", "", mtru$V2, fixed=T)
mtru = subset(mtru, mtru$chr == mtru$map_name)
mtru$mb = ceiling((mtru$V4+mtru$V5)/2)/1000000
mtru = mtru[,c(19,1,8,20)]
names(mtru) = c("chr", "marker", "cm", "mb")
mtru$sp = "mtru"

#mgal
mgal<-read.table("/Volumes/LaCie/Projects/Current/ne/final_data/maps/mgal/turkey_map.txt", header=T, fill=F)
mgal$seqchr=sub("A0", "A", mgal$seqchr, fixed=T)
map.info=rbind(map.info, data.frame(sp="mgal", num.chr=length(unique(mgal$linkchr)), map.length=sum(aggregate(mgal$cm, by=list(chr=mgal$linkchr), max)$x), markers=length(unique(mgal$marker))))
mgal=subset(mgal, mgal$linkchr == mgal$seqchr)
mgal$chr = sub("MGA", "", mgal$seqchr, fixed=T)
mgal$mb = ceiling((mgal$start+mgal$end)/2)/1000000
mgal=mgal[,c("chr", "marker", "cm", "mb")]
mgal$sp = "mgal"
mgal = subset(mgal, !is.na(mgal$mb))

#mouse
mmus<-read.csv("/Volumes/LaCie/Projects/Current/ne/final_data/maps/mmus/Revised_HSmap_SNPs.csv", header=T)
mmus$cm = mmus$ave_cM
mmus = mmus[,c(2,1,7)]
mmus=subset(mmus, !is.na(mmus$cm))
map.info=rbind(map.info, data.frame(sp="mmus", num.chr=length(unique(mmus$chr)), map.length=sum(aggregate(mmus$cm, by=list(chr=mmus$chr), max)$x), markers=length(unique(mmus$snpID))))
mmus$sp = "mmus"
names(mmus)[2] = "marker"
#update positions
mmus.newpos<-read.table("/Volumes/LaCie/Projects/Current/ne/final_data/maps/mmus/hglft_genome_6cde_5886a0.bed", head=F)
mmus.newpos$mb = mmus.newpos$V3/1000000
mmus.newpos$chr = sub("chr", "", mmus.newpos$V1, fixed=T)
mmus.newpos = mmus.newpos[,c(4,5,6)]
names(mmus.newpos)[1] = "marker"
mmus = merge(mmus, mmus.newpos, by="marker")
mmus = subset(mmus, mmus$chr.x==mmus$chr.y)
mmus = mmus[,c("chr.y", "marker", "cm", "mb", "sp")]
names(mmus) = c("chr", "marker", "cm", "mb", "sp")
mmus$marker = as.character(mmus$marker)

#osat
osat.map<-read.table("/Volumes/LaCie/Projects/Current/ne/final_data/maps/osat/rice_map_fromRqtl.txt", header=T)
map.info=rbind(map.info, data.frame(sp="osat", num.chr=length(unique(osat.map$chr)), map.length=sum(aggregate(osat.map$pos, by=list(chr=osat.map$chr), max)$x), markers=length(unique(osat.map$marker))))
osat.pos1<-read.table("/Volumes/LaCie/Projects/Current/ne/final_data/maps/osat/rice_snps.bed", header=F)
osat.pos2<-read.table("/Volumes/LaCie/Projects/Current/ne/final_data/maps/osat/converted_rice_snps.gff", header=F, fill=T, colClasses=c("character", "character", "character", "numeric", "numeric", NA, NA), na.strings=c("NA", "."))
osat.pos = cbind(osat.pos1, osat.pos2)
osat.pos = osat.pos[,c(4,5,9)]
names(osat.pos) = c("marker", "chr.seq", "mb")
osat.pos$mb = osat.pos$mb/1000000
osat.pos = subset(osat.pos, !is.na(osat.pos$mb))

osat =  merge(osat.map, osat.pos, by.x="marker", by.y="marker")
osat = subset(osat, osat$chr == osat$chr.seq)
osat = osat[,c(2,1,3,5)]
osat$sp = "osat"
names(osat) = c("chr", "marker", "cm", "mb", "sp")

#oari
oari<-read.table("/Volumes/LaCie/Projects/Current/ne/final_data/maps/oari/Oari_final_map.txt", header=T)
oari$mb = oari$mb / 1000000
oari$sp = "oari"
oari = droplevels(subset(oari, chr != "X"))
map.info=rbind(map.info, data.frame(sp="oari", num.chr=length(unique(oari$chr)), map.length=sum(aggregate(oari$cm, by=list(chr=oari$chr), max)$x), markers=length(unique(oari$marker))))
oari = subset(oari, as.character(oari$chr) == as.character(oari$seq.chr))
oari = oari[,c(1,2,3,4,6)]

#medaka
olat.map<-read.table("/Volumes/LaCie/Projects/Current/ne/final_data/maps/olat/medaka_map.txt", header=F)
olat.map$map.chr=sub("CH0", "", olat.map$V1)
olat.map$map.chr=sub("CH", "", olat.map$map.chr)
names(olat.map) = c("oldchr", "marker", "cm", "map.chr")
map.info=rbind(map.info, data.frame(sp="olat", num.chr=length(unique(olat.map$map.chr)), map.length=sum(aggregate(olat.map$cm, by=list(chr=olat.map$map.chr), max)$x), markers=length(unique(olat.map$marker))))

olat.epcr<-read.table("/Volumes/LaCie/Projects/Current/ne/final_data/maps/olat/med.primer.map")
names(olat.epcr) = c("marker", "ncbi", "strand", "start", "end", "mm", "mm2", "len")
olat.epcr$mb = ceiling((olat.epcr$start + olat.epcr$end)/2)/1000000
olat.epcr = olat.epcr[,c(1,2,9)]
olat.chrkey = read.table("/Volumes/LaCie/Projects/Current/ne/final_data/maps/olat/chrkey.txt", header=F)
olat.chrkey$key = paste("gi|", olat.chrkey$V3, "|ref|", olat.chrkey$V2, "|", sep="")
olat.chrkey=olat.chrkey[,c(1,6)]
names(olat.chrkey) = c("chr", "key")
olat.pos = merge(olat.epcr, olat.chrkey, by.x="ncbi", by.y="key")
olat = merge(olat.pos, olat.map, by="marker")
olat = subset(olat, olat$chr==olat$map.chr)
olat$sp = "olat"
olat = droplevels(olat)
olat = olat[,c(4,1,6,3,8)]

#panu
panu.map<-read.table("/Volumes/LaCie/Projects/Current/ne/final_data/maps/panu/baboon_map_cm.txt", header=T)
map.info=rbind(map.info, data.frame(sp="panu", num.chr=length(unique(panu.map$chr)), map.length=sum(aggregate(panu.map$cm, by=list(chr=panu.map$chr), max)$x), markers=length(unique(panu.map$marker))))

panu.pos<-read.table("/Volumes/LaCie/Projects/Current/ne/final_data/maps/panu/panu.primer.run2.map", header=F)
names(panu.pos) = c("marker", "chr", "strand", "start", "end", "mm", "mm2", "len")
panu.pos$mb = ceiling((panu.pos$start + panu.pos$end)/2)/1000000
panu.pos$score = panu.pos$mm + panu.pos$mm2;
panu.pos = panu.pos[order(panu.pos$marker, panu.pos$score),]
blast.ranks = as.data.frame(unlist(by(panu.pos$score, panu.pos$marker, rank, ties="min")))
names(blast.ranks) = c("rank")
panu.pos = cbind(panu.pos, blast.ranks)
panu.pos = subset(panu.pos, panu.pos$rank <= 8)
panu.pos = panu.pos[,c(1,2,9)]

panu<-merge(panu.map, panu.pos, by.x="marker", by.y="marker")
panu = subset(panu, panu$chr.x == panu$chr.y)
panu = panu[,c(2,1,3,5)]
panu$sp = "panu"
names(panu) = c("chr", "marker", "cm", "mb", "sp")

panu$cm[panu$chr==10] = abs(panu$cm[panu$chr==10] - max(panu$cm[panu$chr==10]))
panu$cm[panu$chr==15] = abs(panu$cm[panu$chr==15] - max(panu$cm[panu$chr==15]))

#poplar
ptri.map<-read.table("/Volumes/LaCie/Projects/Current/ne/final_data/maps/ptri/poplar_map.txt", header=T)
map.info=rbind(map.info, data.frame(sp="ptri", num.chr=length(unique(ptri.map$chr)), map.length=sum(aggregate(ptri.map$cm, by=list(chr=ptri.map$chr), max)$x), markers=length(unique(ptri.map$Maker))))

ptri.blast<-read.table("/Volumes/LaCie/Projects/Current/ne/final_data/maps/ptri/marker_ptri.out")
ptri.blast = ptri.blast[order(ptri.blast$V1, ptri.blast$V11),]
blast.ranks = as.data.frame(unlist(by(ptri.blast$V11, ptri.blast$V1, rank)))
names(blast.ranks) = c("rank")
ptri.blast = cbind(ptri.blast, blast.ranks)
ptri.blast = subset(ptri.blast, ptri.blast$rank <= 8)

ptri.epcr<-read.table("/Volumes/LaCie/Projects/Current/ne/final_data/maps/ptri/ptri.primer.run2.map")
names(ptri.epcr) = c("marker", "chr", "strand", "start", "end", "mm", "mm2", "len")
ptri.epcr$mb = ceiling((ptri.epcr$start + ptri.epcr$end)/2)/1000000
ptri.epcr = ptri.epcr[,c(1,2,9)]
ptri.epcr$chr = sub("Chr0", "", ptri.epcr$chr)
ptri.epcr$chr = sub("Chr", "", ptri.epcr$chr)

ptri.blast$V2 = sub("Chr0", "", ptri.blast$V2)
ptri.blast$V2 = sub("Chr", "", ptri.blast$V2)
ptri.blast$mb = ceiling((ptri.blast$V9+ptri.blast$V10)/2)/1000000
ptri.blast=ptri.blast[,c(1,2,14)]
names(ptri.blast) = names(ptri.epcr)
ptri.pos = rbind(ptri.blast, ptri.epcr)

ptri=merge(ptri.map, ptri.pos, by.x="Maker", by.y="marker")
ptri=subset(ptri, ptri$chr.x==ptri$chr.y)
ptri=ptri[,c(2,1,4,6)]
ptri$sp = "ptri"
names(ptri) = c("chr", "marker", "cm", "mb", "sp")

ptri$cm[ptri$chr==3] = abs(ptri$cm[ptri$chr==3] - max(ptri$cm[ptri$chr==3]))
ptri$cm[ptri$chr==5] = abs(ptri$cm[ptri$chr==5] - max(ptri$cm[ptri$chr==5]))
ptri$cm[ptri$chr==7] = abs(ptri$cm[ptri$chr==7] - max(ptri$cm[ptri$chr==7]))
ptri$cm[ptri$chr==8] = abs(ptri$cm[ptri$chr==8] - max(ptri$cm[ptri$chr==8]))
ptri$cm[ptri$chr==9] = abs(ptri$cm[ptri$chr==9] - max(ptri$cm[ptri$chr==9]))

#millet
sita<-read.table("/Volumes/LaCie/Projects/Current/ne/final_data/maps/sita/millet_map.txt", header=T)
sita$chr = sub("Chr", "", sita$chr)
sita$chr = as.roman(sita$chr)
sita$chr = as.numeric(sita$chr)
sita$mb = as.numeric(sita$mb)
sita$mb = sita$mb / 1000000
sita$sp = "sita"
map.info=rbind(map.info, data.frame(sp="sita", num.chr=length(unique(sita$chr)), map.length=sum(aggregate(sita$cm, by=list(chr=sita$chr), max)$x), markers=length(unique(sita$marker))))

#pig
sscr<-read.table("/Volumes/LaCie/Projects/Current/ne/final_data/maps/sscr/pigmap.txt", sep="\t", header=T)
sscr=sscr[,c("ssc", "snp.start", "cm", "mb")]
names(sscr) = c("chr", "marker", "cm", "mb")
sscr$mb = sscr$mb / 1000000
sscr$sp = "sscr"
map.info=rbind(map.info, data.frame(sp="sscr", num.chr=length(unique(sscr$chr)), map.length=sum(aggregate(sscr$cm, by=list(chr=sscr$chr), max)$x), markers=length(unique(sscr$marker))))

#maize
nam.map<-read.table("/Volumes/LaCie/Projects/Current/ne/final_data/maps/zmay/nam_map.txt", header=T)
nam.pos<-read.table("/Volumes/LaCie/Projects/Current/ne/final_data/maps/zmay/nam_markers.txt", header=T)
map.info=rbind(map.info, data.frame(sp="zmay", num.chr=length(unique(nam.map$ch)), map.length=sum(aggregate(nam.map$cumulative, by=list(chr=nam.map$ch), max)$x), markers=length(unique(nam.map$marker))))

nam.merge<-merge(nam.map, nam.pos, all.x=T, all.y=F, by.x="marker", by.y="SNP_NAME")
nam.merge$chr.agree = ifelse(nam.merge$AGPv2_Chr == nam.merge$ch, T, F)
zmay = subset(nam.merge, nam.merge$chr.agree == T, select=c("ch", "marker", "cumulative", "AGPv2_pos"))
names(zmay) = c("chr", "marker", "cm", "mb")
zmay$mb = zmay$mb / 1000000
zmay$sp = "zmay"

#pper
pper.raw<-read.table("/Volumes/LaCie/Projects/Current/ne/final_data/maps/pper/pper_raw_map_data.txt", header=T, fill=T, sep="\t")
pper.raw<-subset(pper.raw, pper.raw$map.name=="TxE", select=c("marker.name", "scaffold", "start..bp.", "linkage.group", "start"))
names(pper.raw) = c("marker", "seq.chr", "mb", "map.chr", "cm")
pper.raw$seq.chr = sub("scaffold_", "", pper.raw$seq.chr)
pper.raw$map.chr = sub("G", "", pper.raw$map.chr)
map.info=rbind(map.info, data.frame(sp="pper", num.chr=length(unique(pper.raw$map.chr)), map.length=sum(aggregate(pper.raw$cm, by=list(chr=pper.raw$map.chr), max)$x), markers=length(unique(pper.raw$marker))))
pper = subset(pper.raw, pper.raw$seq.chr == pper.raw$map.chr, select=c("seq.chr", "marker", "cm", "mb"))
pper$sp = "pper"
pper$mb = pper$mb/1000000
names(pper)[1] = "chr"

#sorghun
sbic.map<-read.table("/Volumes/LaCie/Projects/Current/ne/final_data/maps/sbic/CIRAD-map.txt", header=T, sep="\t")
sbic.map$cm.start = as.numeric(as.character(sub("-\\d+\\.?\\d*", "", sbic.map$cm, perl=T)))
sbic.map$cm.end = as.numeric(as.character(sub("\\d+\\.?\\d*-", "", sbic.map$cm, perl=T)))
sbic.map$cm = (sbic.map$cm.start + sbic.map$cm.end) / 2
map.info=rbind(map.info, data.frame(sp="sbic", num.chr=length(unique(sbic.map$map.chr)), map.length=sum(aggregate(sbic.map$cm.end, by=list(chr=sbic.map$map.chr), max)$x), markers=length(unique(sbic.map$marker))))

sbic.pos<-read.table("/Volumes/LaCie/Projects/Current/ne/final_data/maps/sbic/markers_all.txt", header=T, sep="\t", colClasses=c("factor", "numeric", "numeric"))
sbic = merge(sbic.map, sbic.pos, by="marker")
sbic = subset(sbic, sbic$map.chr == sbic$chr)
sbic=sbic[,c(2,1,3,6)]
names(sbic)=c("chr", "marker", "cm", "mb")
sbic$mb = sbic$mb / 1000000
sbic$sp = "sbic"

sbic$cm[sbic$chr==9] = abs(sbic$cm[sbic$chr==9] - max(sbic$cm[sbic$chr==9]))

#locu
locu.map<-read.table("/Volumes/LaCie/Projects/Current/ne/final_data/maps/locu/locu_map.txt", header=T, sep="\t")
locu.map$chr = sub("LocLG", "", locu.map$Linkage.Group)
names(locu.map) = c("marker", "lg", "cm", "chr")
map.info=rbind(map.info, data.frame(sp="locu", num.chr=length(unique(locu.map$chr)), map.length=sum(aggregate(locu.map$cm, by=list(chr=locu.map$chr), max)$x), markers=length(unique(locu.map$marker))))

locu.blast<-read.table("/Volumes/LaCie/Projects/Current/ne/final_data/maps/locu//marker_locu.out", header=F, sep="\t")
locu.blast = locu.blast[order(locu.blast$V1, locu.blast$V11),]
blast.ranks = as.data.frame(unlist(by(locu.blast$V11, locu.blast$V1, rank)))
names(blast.ranks) = c("rank")
locu.blast = cbind(locu.blast, blast.ranks)
locu.blast = subset(locu.blast, locu.blast$rank <= 3)
locu.blast$chracc = sub("^gi\\|\\d+\\|gb\\|", "", locu.blast$V2, perl=T)
locu.blast$chracc = sub("\\|", "", locu.blast$chr, perl=T)

locu.chrkey = read.table("/Volumes/LaCie/Projects/Current/ne/final_data/maps/locu/chr2acc", header=F)
locu.blast = merge(locu.blast, locu.chrkey, all.x=T, by.x="chracc", by.y="V2")
locu.blast$chr = sub("LG", "", locu.blast$V1.y)
locu.blast$pos = ceiling((locu.blast$V9 + locu.blast$V10)/2)/1e6
locu.blast = locu.blast[,c("V1.x", "chr", "pos")]
names(locu.blast)[1] = "marker"

locu = merge(locu.map, locu.blast, by="marker")
locu = subset(locu, locu$chr.x == locu$chr.y)
locu = locu[,c("chr.x", "marker", "cm", "pos")]
locu$sp = "locu"
names(locu) = c("chr", "marker", "cm", "mb", "sp")

#fix chr
locu$cm[locu$chr==5] = abs(locu$cm[locu$chr==5] - max(locu$cm[locu$chr==5]))
locu$cm[locu$chr==9] = abs(locu$cm[locu$chr==9] - max(locu$cm[locu$chr==9]))
locu$cm[locu$chr==10] = abs(locu$cm[locu$chr==10] - max(locu$cm[locu$chr==10]))
locu$cm[locu$chr==12] = abs(locu$cm[locu$chr==12] - max(locu$cm[locu$chr==12]))
locu$cm[locu$chr==14] = abs(locu$cm[locu$chr==14] - max(locu$cm[locu$chr==14]))
locu$cm[locu$chr==16] = abs(locu$cm[locu$chr==16] - max(locu$cm[locu$chr==16]))
locu$cm[locu$chr==17] = abs(locu$cm[locu$chr==17] - max(locu$cm[locu$chr==17]))
locu$cm[locu$chr==18] = abs(locu$cm[locu$chr==18] - max(locu$cm[locu$chr==18]))
locu$cm[locu$chr==19] = abs(locu$cm[locu$chr==19] - max(locu$cm[locu$chr==19]))
locu$cm[locu$chr==21] = abs(locu$cm[locu$chr==21] - max(locu$cm[locu$chr==21]))
locu$cm[locu$chr==23] = abs(locu$cm[locu$chr==23] - max(locu$cm[locu$chr==23]))
locu$cm[locu$chr==25] = abs(locu$cm[locu$chr==25] - max(locu$cm[locu$chr==25]))

#ficAlb
falb.map<-read.table("/Volumes/LaCie/Projects/Current/ne/final_data/maps/falb/Chrompic_Results_ver20140121.txt", header=T, comment.char="", sep="\t")
falb.map<-falb.map[,c(1,4,13,14)]
falb.map$id = paste(falb.map$scaf_FAlb15, falb.map$pos_FAlb15, sep="-")
names(falb.map)=c("chr", "cm", "scaf", "pos", "id")
falb.map=falb.map[,c("chr", "cm", "id")]
falb.map$chr = tolower(falb.map$chr)
falb.map=subset(falb.map, chr != "chrz")
map.info=rbind(map.info, data.frame(sp="falb", num.chr=length(unique(falb.map$chr)), map.length=sum(aggregate(falb.map$cm, by=list(chr=falb.map$chr), max)$x), markers=length(unique(falb.map$id))))

falb.snp<-read.table("/Volumes/LaCie/Projects/Current/ne/final_data/maps/falb/snps.txt", header=F, sep="\t")
names(falb.snp)=c("id", "chr", "pos")
falb.snp=subset(falb.snp, !grepl("N\\d+", chr))
falb.snp=droplevels(falb.snp)

falb=merge(falb.map, falb.snp, by="id")
falb=falb[,c("chr.y", "id", "cm", "pos")]
names(falb)=c("chr", "marker", "cm", "mb")
falb$mb = falb$mb/1e6
falb$sp = "falb"

#csem
csem.acc<-read.table("/Volumes/LaCie/Projects/Current/ne/final_data/maps/csem/acc_key.txt", header=T, sep="\t")
names(csem.acc)=c("marker", "acc")
csem.map<-read.table("/Volumes/LaCie/Projects/Current/ne/final_data/maps/csem/csem_map_forR.txt", header=T, sep="\t")
names(csem.map)=c("chr", "cm", "marker")

map.info=rbind(map.info, data.frame(sp="csem", num.chr=length(unique(csem.map$chr)), map.length=sum(aggregate(csem.map$cm, by=list(chr=csem.map$chr), max)$x), markers=length(unique(csem.map$marker))))


csem.map<-merge(csem.map, csem.acc)

csem.blast<-read.table("/Volumes/LaCie/Projects/Current/ne/final_data/maps/csem/marker_csem.out", header=F, sep="\t")
csem.blast = csem.blast[order(csem.blast$V1, csem.blast$V11),]
blast.ranks = as.data.frame(unlist(by(csem.blast$V11, csem.blast$V1, rank)))
names(blast.ranks) = c("rank")
csem.blast = cbind(csem.blast, blast.ranks)
csem.blast = subset(csem.blast, csem.blast$rank <= 3)
csem.blast$acc = sub("^gi\\|\\d+\\|gb\\|", "", csem.blast$V1, perl=T)
csem.blast$acc = sub("\\.\\d+\\|", "", csem.blast$acc, perl=T)
csem.blast$chracc = sub("^gi\\|\\d+\\|gb\\|", "", csem.blast$V2, perl=T)
csem.blast$chracc = sub("\\.\\d+\\|", "", csem.blast$chracc, perl=T)

csem=merge(csem.map, csem.blast, by="acc")
csem.chr2acc<-read.table("/Volumes/LaCie/Projects/Current/ne/final_data/maps/csem/chr2acc.txt", header=T)
csem=merge(csem, csem.chr2acc, by.x="chracc", by.y="Accession.version")
csem=subset(csem, chr.x != "LG1")
csem$chr = sub("LG", "", csem$chr.x)
csem$chr = as.numeric(csem$chr) - 1
csem$chr[csem$chr==11]=112
csem$chr[csem$chr==12]=111
csem$chr[csem$chr==112]=12
csem$chr[csem$chr==111]=11
csem=subset(csem, csem$chr == csem$chr.y)
csem$mb = ceiling((csem$V9 + csem$V10)/2)/1e6
csem=csem[,c("chr", "marker", "cm", "mb")]

#fix chr
csem$cm[csem$chr==17] = abs(csem$cm[csem$chr==17] - max(csem$cm[csem$chr==17]))
csem$cm[csem$chr==16] = abs(csem$cm[csem$chr==16] - max(csem$cm[csem$chr==16]))
csem$cm[csem$chr==13] = abs(csem$cm[csem$chr==13] - max(csem$cm[csem$chr==13]))
csem$cm[csem$chr==11] = abs(csem$cm[csem$chr==11] - max(csem$cm[csem$chr==11]))
csem$cm[csem$chr==8] = abs(csem$cm[csem$chr==8] - max(csem$cm[csem$chr==8]))
csem$cm[csem$chr==5] = abs(csem$cm[csem$chr==5] - max(csem$cm[csem$chr==5]))
csem$cm[csem$chr==4] = abs(csem$cm[csem$chr==4] - max(csem$cm[csem$chr==4]))
csem$cm[csem$chr==2] = abs(csem$cm[csem$chr==2] - max(csem$cm[csem$chr==2]))
csem$sp="csem"

##write output##
#output is a data frame for each species (total 40) with the genetic and physical map information. All maps are filtered to include only autosomes, and in some cases only a subset of autosomes (where either there is no marker information for small chromosomes, or some chromsomes have inversions or other issues that prevents usability). Also produces a data frame, map.info, with information on each map

#clean up temporary objects
map.info.raw<-map.info

#concatencate all maps into a list
all_maps<-list(agam=agam, amel=amel, atha=atha, bdis=bdis, bmor=bmor, btau=btau, cbri=cbri, ccle=ccle, cele=cele, clan=clan, clup=clup, crub=crub, csat=csat, csem=csem, dmel=dmel, dpse=dpse, drer=drer, ecab=ecab, falb=falb, gacu=gacu, ggal=ggal, gmax=gmax, grai=grai, hmel=hmel, hsap=hsap, locu=locu, mgal=mgal, mmul=mmul, mmus=mmus, mtru=mtru, oari=oari, olat=olat, osat=osat, panu=panu, pper=pper, ptri=ptri, sbic=sbic, sita=sita, sscr=sscr, zmay=zmay)

#write raw map objects to file
save(list=ls(pattern="^\\w\\w\\w\\w$"), file="raw_genetic_maps.RData")

#write list to file
save(all_maps, file="all_maps_list.RData")

#write raw map info
write.table(map.info.raw, file="map_info_raw.txt", sep="\t", col.names=T, row.names=F, quote=F)
