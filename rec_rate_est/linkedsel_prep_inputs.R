#set up functions
source("../linkedsel_functions.R")

#load inputs
load(file="../prepare_input_files/all_maps_list.RData")
map.info.raw<-read.table("../prepare_input_files/map_info_raw.txt", header=T, sep="\t")

#########THIS CODE CLEANS MAPS AND ESTIMATES RECOMBINATION RATE###############

#clean up maps: make sure all chr are characters and not factors or numbers
#rename chr to add "chr" to beginning of string
all_maps_clean<-lapply(all_maps, cleanup_maps)

#get information about mapped markers for each map
map.info.mapped<-as.data.frame(sapply(all_maps_clean, function(x) length(unique(x$marker)), simplify=T))
names(map.info.mapped)<-c("mapped.markers")
map.info.mapped<-merge(map.info.raw, map.info.mapped, by.y="row.names", by.x="sp")

#now remove mismapped markers from all_maps_clean
all_maps_mm_removed<-lapply(all_maps_clean, remove_mismapped)

#plot all maps with duplicated/mismapped markers indicated
for (map in all_maps_mm_removed) {
  spe=unique(map$sp)
  num.chr = length(unique(map$chr));
  map$chr=sub("chr", "", map$chr)
  map$chr[grepl("^\\d\\D*$", map$chr, perl=T)]=paste0("0", map$chr[grepl("^\\d\\D*$", map$chr, perl=T)])
  map$chr[grepl("^\\d", map$chr, perl=T)]=paste0("chr", map$chr[grepl("^\\d", map$chr, perl=T)])
  chr.id = unique(map$chr)
  chr.id=sort(chr.id)
  #make pdf for each species
  pdf(paste("./plots/marey_maps/", spe, "_", "marey_map_all", ".pdf", sep=""))
  for (i in chr.id) {
    with(map[map$chr == i,], 
         plot(cm ~ mb, las=1, main = i, 
              pch=ifelse(marker.freq>1, 18, 16), 
              col=ifelse(good.marker, "black", "red")))
  }
  dev.off()
  rm(i)
  rm(map)
  rm(chr.id)
  rm(spe)
}

#now remove duplicates; this list now contains the final maps for each species
all_maps_final<-lapply(all_maps_mm_removed, remove_duplicate_markers)

#get info after removing duplicates
map.info.rmbad<-as.data.frame(sapply(all_maps_final, function(x) sum(x$good.marker)))
names(map.info.rmbad)=c("good.markers")
map.info.mb<-as.data.frame(sapply(all_maps_final, function(z) sum(aggregate(z$mb, list(chr=z$chr), function(x) max(x)-min(x))$x)))
names(map.info.mb)=c("mb.tot")
map.info.final<-merge(map.info.mapped, map.info.rmbad, by.y="row.names", by.x="sp")
map.info.final<-merge(map.info.final, map.info.mb, by.y="row.names", by.x="sp")
map.info.final$density.good = map.info.final$map.length/map.info.final$good.markers
map.info.final$prop.good = map.info.final$good.markers/map.info.final$markers
map.info.final$rec.rate = map.info.final$map.length/map.info.final$mb.tot
write.table(file="map_info_final.txt", map.info.final, row.names=F, sep="\t", quote=F)


#plot final marey maps
for (map in all_maps_final) {
  spe=unique(map$sp)
  num.chr = length(unique(map$chr));
  map$chr=sub("chr", "", map$chr)
  map$chr[grepl("^\\d\\D*$", map$chr, perl=T)]=paste0("0", map$chr[grepl("^\\d\\D*$", map$chr, perl=T)])
  map$chr[grepl("^\\d", map$chr, perl=T)]=paste0("chr", map$chr[grepl("^\\d", map$chr, perl=T)])
  chr.id = unique(map$chr)
  chr.id=sort(chr.id)
  #make pdf for each species
  pdf(paste("./plots/marey_maps/", spe, "_", "marey_map_final", ".pdf", sep=""))
  for (i in chr.id) {
    with(map[map$chr == i & map$good.marker==T,], 
         plot(cm ~ mb, las=1, main = i, 
              pch=16, 
              col="black"))
  }
  dev.off()
  rm(i)
  rm(map)
  rm(chr.id)
  rm(spe)
}

#make map mask - this is used later in the analysis to exclude windows in regions of the genome with poorly estimated recombination rate do to bad maps or lack of markers at chromosome ends
all_masks<-lapply(all_maps_final, make_map_mask)

#plot marey maps with masked regions indicated
for (map in all_maps_final) {
  spe=unique(map$sp)
  num.chr = length(unique(map$chr));
  map$chr=sub("chr", "", map$chr)
  map$chr[grepl("^\\d\\D*$", map$chr, perl=T)]=paste0("0", map$chr[grepl("^\\d\\D*$", map$chr, perl=T)])
  map$chr[grepl("^\\d", map$chr, perl=T)]=paste0("chr", map$chr[grepl("^\\d", map$chr, perl=T)])
  chr.id = unique(map$chr)
  chr.id=sort(chr.id)
  #make pdf for each species
  pdf(paste("./plots/marey_maps/", spe, "_", "marey_map_masked", ".pdf", sep=""))
  for (i in chr.id) {
    with(map[map$chr == i,], 
         plot(cm ~ mb, las=1, main = i, 
              pch=16, 
              col=ifelse(good.marker==T, "black", "red")))
    mask.chr=sub("chr0","chr", i)
    chr.mask<-subset(all_masks[[which(names(all_masks)==spe)]], chr==mask.chr)
    chr.mask.rows = nrow(chr.mask)
    if (chr.mask.rows>0) {
      for (mask.row in 1:chr.mask.rows) {
        min.cm=min(map[map$chr==i,"cm"])-1.5
        segments(x0=chr.mask[mask.row,"start"], y0=min.cm, x1 = chr.mask[mask.row,"end"], y1 = min.cm, lwd=2, col="black")	
      }
    }
    
  }
  dev.off()
  rm(i)
  rm(map)
  rm(chr.id)
  rm(spe)
  rm(num.chr)
  rm(mask.chr)
  rm(chr.mask)
  rm(chr.mask.rows)
  rm(min.cm)
  rm(mask.row)
}

#estimate recombination using both a piecewise regression and a high order polynomial
rec_piecewise<-lapply(all_maps_final, function(x) est_rec_piece(x[x$good.marker==T, c("chr", "cm", "mb")]))
rec_polynom<-lapply(all_maps_final, function(x) est_rec_poly(x[x$good.marker==T, c("chr", "cm", "mb")]))

#get r2
poly.r2<-unlist(lapply(rec_polynom, function(x) lapply(x, function(y) summary(y$func)$adj.r.squared)))
piece.r2<-unlist(lapply(rec_piecewise, function(x) lapply(x, function(y) summary(y$func)$adj.r.squared)))

#merged correlations
merged.r2=data.frame(piece=c(piece.r2), poly=c(poly.r2), species=names(poly.r2), stringsAsFactors=F)
merged.r2$species = sub(".\\w+$", "", merged.r2$species)
merged.r2$chr=sub("^\\w+.", "", row.names(merged.r2))
merged.r2$use.piece=F
merged.r2$use.piece[merged.r2$poly <= merged.r2$piece]=T

#plot recombination functions on marey maps
for (map in all_maps_final) {
  require(splines)
  require(polynom)
  map = subset(map, good.marker==T)
  spe=unique(map$sp)
  #clean up map and filter to remove chromosomes with fewer than 4 markers
  map=map[order(map$chr, map$mb),]
  map=merge(map, as.data.frame(table(map$chr)), all.x=T, by.x="chr", by.y="Var1")
  map=subset(map, map$Freq >= 3, select=c("chr", "cm", "mb"))
  num.chr = length(unique(map$chr));
  chr.id = unique(map$chr)
  rownames(map)=seq(1,nrow(map),1)
  
  map$chr=sub("chr", "", map$chr)
  map$chr[grepl("^\\d\\D*$", map$chr, perl=T)]=paste0("0", map$chr[grepl("^\\d\\D*$", map$chr, perl=T)])
  map$chr[grepl("^\\d", map$chr, perl=T)]=paste0("chr", map$chr[grepl("^\\d", map$chr, perl=T)])
  chr.id = unique(map$chr)
  chr.id=sort(chr.id)
  #make pdf for each species
  pdf(paste("./plots/marey_maps/", spe, "_", "marey_map_recrate", ".pdf", sep=""))
  for (i in chr.id) {
    func.chr=sub("chr0","chr", i)
    
    #get interpolation splines for plotting
    #check for monotonicity
    tmp.poly.pred = predict(rec_polynom[[spe]][[func.chr]][["func"]])
    tmp.piece.pred = predict(rec_piecewise[[spe]][[func.chr]][["func"]])
    
    tmp.poly.inc = (tmp.poly.pred==cummax(tmp.poly.pred))
    tmp.piece.inc = (tmp.piece.pred==cummax(tmp.piece.pred))
    
    tmp.mb.inc = map$mb[map$chr==i]
    
    tmp.poly.interp<-splinefun(x=tmp.mb.inc[tmp.poly.inc], y=tmp.poly.pred[tmp.poly.inc], method="monoH.FC")
    tmp.piece.interp<-splinefun(x=tmp.mb.inc[tmp.piece.inc], y=tmp.piece.pred[tmp.piece.inc], method="monoH.FC")
    tmp.p0<-polynomial(coef(rec_polynom[[spe]][[func.chr]][["func"]]))
    
    #make plot of actual vs fitting curve
    tmp.max.mb<-max(map$mb[map$chr==i])
    tmp.min.mb<-min(map$mb[map$chr==i])
    tmp.max.cm<-max(map$cm[map$chr==i])
    tmp.mb.pred = data.frame(mb=seq(tmp.min.mb,tmp.max.mb,0.01))
    tmp.mb.pred.poly = poly(tmp.mb.pred$mb, rec_polynom[[spe]][[func.chr]][["deg"]], raw=T)
    tmp.cm.piece = tmp.piece.interp(tmp.mb.pred$mb)
    tmp.cm.poly = tmp.poly.interp(tmp.mb.pred$mb)
    tmp.rate.piece = tmp.piece.interp(tmp.mb.pred$mb, deriv=1)
    tmp.rate.poly = tmp.poly.interp(tmp.mb.pred$mb, deriv=1)
    with(map[map$chr==i,], plot(cm ~ mb, main = i, xlim=c(0,tmp.max.mb+(tmp.max.mb*0.01)), ylim=c(0,tmp.max.cm+(tmp.max.cm*0.01)), cex=0.75))
    lines(tmp.p0, col="darkgreen", lwd=2)
    lines(tmp.mb.pred$mb, predict(rec_piecewise[[spe]][[func.chr]][["func"]], tmp.mb.pred), col="blue1", lwd=2)
    lines(tmp.cm.piece ~ tmp.mb.pred$mb, col="blue1", lwd=1, lty="dashed")
    lines(tmp.cm.poly ~ tmp.mb.pred$mb, col="darkgreen", lwd=1, lty="dashed")   
  }
  dev.off()
  
  #make rate plots
  pdf(paste("./plots/recrate/", spe, "_", "recrate", ".pdf", sep=""))
  for (i in chr.id) {
    func.chr=sub("chr0","chr", i)
    
    #get interpolation splines for plotting
    #check for monotonicity
    tmp.poly.pred = predict(rec_polynom[[spe]][[func.chr]][["func"]])
    tmp.piece.pred = predict(rec_piecewise[[spe]][[func.chr]][["func"]])
    
    tmp.poly.inc = (tmp.poly.pred==cummax(tmp.poly.pred))
    tmp.piece.inc = (tmp.piece.pred==cummax(tmp.piece.pred))
    
    tmp.mb.inc = map$mb[map$chr==i]
    
    tmp.poly.interp<-splinefun(x=tmp.mb.inc[tmp.poly.inc], y=tmp.poly.pred[tmp.poly.inc], method="monoH.FC")
    tmp.piece.interp<-splinefun(x=tmp.mb.inc[tmp.piece.inc], y=tmp.piece.pred[tmp.piece.inc], method="monoH.FC")
    
    tmp.max.mb<-max(map$mb[map$chr==i])
    tmp.min.mb<-min(map$mb[map$chr==i])
    tmp.max.cm<-max(map$cm[map$chr==i])
    tmp.mb.pred = data.frame(mb=seq(tmp.min.mb,tmp.max.mb,0.01))
    tmp.cm.piece = tmp.piece.interp(tmp.mb.pred$mb)
    tmp.cm.poly = tmp.poly.interp(tmp.mb.pred$mb)
    tmp.rate.piece = tmp.piece.interp(tmp.mb.pred$mb, deriv=1)
    tmp.rate.poly = tmp.poly.interp(tmp.mb.pred$mb, deriv=1)
    
    #piecewise deriv
    tmp.piece.deriv = predict(rec_piecewise[[spe]][[func.chr]][["func"]],tmp.mb.pred)
    tmp.dY<-diff(tmp.piece.deriv)/diff(tmp.mb.pred$mb)
    tmp.dY[tmp.dY<0] = 0
    tmp.dY[tmp.dY>quantile(tmp.dY,0.975)] = quantile(tmp.dY,0.975)
    tmp.dX<-rowMeans(embed(tmp.mb.pred$mb,2))
    
    tmp.max.rate = quantile(tmp.dY,0.975)*1.01
    
    #poly deriv
    tmp.p0<-polynomial(coef(rec_polynom[[spe]][[func.chr]][["func"]]))
    tmp.p1<-deriv(tmp.p0)
    
    #plot
    plot(tmp.p1, col="darkgreen", lwd=2, xlab="Mb", ylab="cM/Mb", type="l", main=i, ylim=c(0,tmp.max.rate), xlim=c(0,tmp.max.mb+(tmp.max.mb*0.01)))
    lines(tmp.dX,tmp.dY, col="blue1", lwd=2)
    lines(tmp.rate.piece ~ tmp.mb.pred$mb, col="blue1", lwd=2, lty="dashed")
    lines(tmp.rate.poly ~ tmp.mb.pred$mb, col="darkgreen", lwd=2, lty="dashed")
  }
  dev.off()
  
  rm(i)
  rm(map)
  rm(chr.id)
  rm(spe)
  rm(num.chr)
  rm(func.chr)
  rm(list=ls(pattern="^tmp.*"))
}

#read in gene density estimates 
gd.100<-read.table("../prepare_input_files/gd.100.txt", stringsAsFactors=F)
gd.100$window = 100
gd.500<-read.table("../prepare_input_files/gd.500.txt", stringsAsFactors=F)
gd.500$window = 500
gd.1000<-read.table("../prepare_input_files/gd.1000.txt", stringsAsFactors=F)
gd.1000$window = 1000
gd<-rbind(gd.100, gd.500, gd.1000)
#convert species names to be consistent with maps
names(gd)=c("species", "chr", "window.num", "window.pos", "exon.bp", "exon.frac", "window")
gd$species = tolower(gd$species)
gd$species = substring(gd$sp, 1, 4)
gd$species[gd$species=="bman"] = "bmor"
gd$species[gd$species=="cret"] = "ccle"
gd$species[gd$species=="mmca"] = "mmus"
gd$species[gd$species=="ocan"] = "oari"
gd$species[gd$species=="oruf"] = "osat"
gd$species[gd$species=="pdav"] = "pper"
#clean up temp objects
rm(list=ls(pattern="gd.\\d+"))
gd$chr = sub("^0", "", gd$chr, perl=T)
gd$chr[!grepl("LG", gd$chr)]=paste0("chr", gd$chr[!grepl("LG", gd$chr)])

#make window data frame to estimate rec rate for
windows<-gd
windows<-windows[,c("species", "chr", "window.pos", "window", "exon.bp")]
names(windows)[4]="window.size"
windows$window.id = paste(windows$window.size, windows$window.pos, sep="_")
windows=windows[order(windows$species, windows$chr, windows$window.pos),]

#compute recombination rate for windows by species and chromosome
require(plyr)

windows$poly.rr<-unlist(dlply(windows, .(species, chr), splat(est_rr_poly)))
windows$piece.rr<-unlist(dlply(windows, .(species, chr), splat(est_rr_piece)))
windows$poly.cm.i<-round(unlist(dlply(windows, .(species, chr), function(x) est_cm_poly(species=x$species, chr=x$chr, window.pos=x$window.pos-(x$window.size*500), window.size=x$window.size)))/100,5)
windows$poly.cm.k<-round(unlist(dlply(windows, .(species, chr), function(x) est_cm_poly(species=x$species, chr=x$chr, window.pos=x$window.pos, window.size=x$window.size)))/100,5)
windows$poly.cm.l<-round(unlist(dlply(windows, .(species, chr), function(x) est_cm_poly(species=x$species, chr=x$chr, window.pos=x$window.pos+(x$window.size*500), window.size=x$window.size)))/100,5)
windows$piece.cm.i<-round(unlist(dlply(windows, .(species, chr), function(x) est_cm_piece(species=x$species, chr=x$chr, window.pos=x$window.pos-(x$window.size*500), window.size=x$window.size)))/100,5)
windows$piece.cm.k<-round(unlist(dlply(windows, .(species, chr), function(x) est_cm_piece(species=x$species, chr=x$chr, window.pos=x$window.pos, window.size=x$window.size)))/100,5)
windows$piece.cm.l<-round(unlist(dlply(windows, .(species, chr), function(x) est_cm_piece(species=x$species, chr=x$chr, window.pos=x$window.pos+(x$window.size*500), window.size=x$window.size)))/100,5)

#export data for Gk calculations and popgen modeling
#Gk script expects: species, chromosome, window.id, mi_poly, mk_poly, ml_poly, {mi,mk,ml}_piece, fd 
#where mi, mk, ml are estimated positions of start, middle, and end of a window in Morgans
#and fd is the proportion of exonic sites in the genome that fall in a given window
#input to Gk should be sorted by species and chromosome (and position?) as Gk script does not santitize inputs
#Gk script also requires separate files for a given window size

all.exons<-read.table("../prepare_input_files/total.exons.txt")
names(all.exons)=c("species", "total.exons")
all.exons$species = tolower(all.exons$species)
all.exons$species = substring(all.exons$sp, 1, 4)
all.exons$species[all.exons$species=="bman"] = "bmor"
all.exons$species[all.exons$species=="cret"] = "ccle"
all.exons$species[all.exons$species=="mmca"] = "mmus"
all.exons$species[all.exons$species=="ocan"] = "oari"
all.exons$species[all.exons$species=="oruf"] = "osat"
all.exons$species[all.exons$species=="pdav"] = "pper"

windows<-merge(windows, all.exons, all.x=T, by="species")
windows$fd = windows$exon.bp/windows$total.exons
windows=windows[complete.cases(windows),] #remove regions w/o recomb estimates

wind100.poly.for.gk<-windows[windows$window.size==100,c("species", "chr", "window.id", "poly.cm.i", "poly.cm.k", "poly.cm.l","fd")]
wind500.poly.for.gk<-windows[windows$window.size==500,c("species", "chr", "window.id", "poly.cm.i", "poly.cm.k", "poly.cm.l", "fd")]
wind1000.poly.for.gk<-windows[windows$window.size==1000,c("species", "chr", "window.id", "poly.cm.i", "poly.cm.k", "poly.cm.l", "fd")]
wind100.piece.for.gk<-windows[windows$window.size==100,c("species", "chr", "window.id", "piece.cm.i", "piece.cm.k", "piece.cm.l", "fd")]
wind500.piece.for.gk<-windows[windows$window.size==500,c("species", "chr", "window.id", "piece.cm.i", "piece.cm.k", "piece.cm.l", "fd")]
wind1000.piece.for.gk<-windows[windows$window.size==1000,c("species", "chr", "window.id", "piece.cm.i", "piece.cm.k", "piece.cm.l", "fd")]

#add mask data to window file
windows$mb = windows$window.pos/1e6
mask.df=ldply(all_masks, data.frame)
names(mask.df)=c("species", "chr", "start", "end")
windows.filt=merge(windows, mask.df, all=T)
windows.filt$startP = windows.filt$mb - (windows.filt$window.size/1000/2)
windows.filt$endP = windows.filt$mb + (windows.filt$window.size/1000/2)
windows.filt = windows.filt[with(windows.filt, startP >= start & endP <= end),]
windows.filt = windows.filt[,c(1,2,5)]
windows.filt$use = F
windows.clean = merge(windows, windows.filt, all.x=T)
windows.clean$use[is.na(windows.clean$use)] = T

#windows clean now has recombination rate for each window

#make mutation rate files
#estimates for Arabidposis, Human, Drosophila, and C.elegans are used for plants, verts, insects, and nematodes respectively

mut.ests<-data.frame(class=c("plant", "vert", "ins", "nem"), est=c(0.000000007,0.0000000118,0.0000000028,0.0000000027))
u.df<-data.frame(species=names(all_maps_final), class=c("ins", "ins", "plant", "plant", "ins", "vert", "nem", "plant", "nem", "plant", "vert", "plant", "plant", "vert", "ins", "ins", "vert", "vert", "vert", "vert", "vert", "plant", "plant", "ins", "vert", "vert", "vert", "vert", "vert", "plant", "vert", "vert", "plant", "vert", "plant", "plant", "plant", "plant", "vert", "plant"))
genome.size<-read.table("../prepare_input_files/genome_size.txt", header=T)
u.df<-merge(u.df, mut.ests)
u.df<-merge(u.df, all.exons)
u.df<-merge(u.df, genome.size)
u.df$frac.conserved<-pmin((5*(u.df$total.exons/1e6)/u.df$size.mb),1)
u.df$u.min=u.df$est*u.df$total.exons*2 #diploid best est of mutation rate based on coding exons only
u.df$u.const=1
u.df$u.max=u.df$frac.conserved*u.df$size.mb*1e6*u.df$est*2

#add theta for each window in windows for theta
#read in theta estimates and clean up
theta.std.100<-read.table("../prepare_input_files/std.100.txt", stringsAsFactors=F)
theta.std.100$window = 100
theta.std.100$filt = "std"
theta.std.500<-read.table("../prepare_input_files/std.500.txt", stringsAsFactors=F)
theta.std.500$window = 500
theta.std.500$filt = "std"
theta.std.1000<-read.table("../prepare_input_files/std.1000.txt", stringsAsFactors=F)
theta.std.1000$window = 1000
theta.std.1000$filt = "std"
theta.q30.100<-read.table("../prepare_input_files/q30.100.txt", stringsAsFactors=F)
theta.q30.100$window = 100
theta.q30.100$filt = "q30"
theta.q30.500<-read.table("../prepare_input_files/q30.500.txt", stringsAsFactors=F)
theta.q30.500$window = 500
theta.q30.500$filt = "q30"
theta.q30.1000<-read.table("../prepare_input_files/q30.1000.txt", stringsAsFactors=F)
theta.q30.1000$window = 1000
theta.q30.1000$filt = "q30"
theta<-rbind(theta.std.100, theta.std.500, theta.std.1000, theta.q30.100, theta.q30.500, theta.q30.1000)
#convert species names to be consistent with maps
names(theta)=c("species", "chr", "window.num", "window.pos", "pi", "sites", "window", "filt")
theta$species = tolower(theta$species)
theta$species = substring(theta$sp, 1, 4)
theta$species[theta$species=="bman"] = "bmor"
theta$species[theta$species=="cret"] = "ccle"
theta$species[theta$species=="mmca"] = "mmus"
theta$species[theta$species=="ocan"] = "oari"
theta$species[theta$species=="oruf"] = "osat"
theta$species[theta$species=="pdav"] = "pper"

#fix chr names
theta$chr = sub("^0", "", theta$chr, perl=T)
theta$chr[!grepl("LG", theta$chr)]=paste0("chr", theta$chr[!grepl("LG", theta$chr)])
theta$window.id = paste(theta$window, theta$window.pos, sep="_")
windows.for.theta=windows.clean[,c("species", "chr", "window.id", "fd", "poly.rr", "piece.rr", "use")]
windows.for.theta=merge(windows.for.theta, theta, all=T, by=c("species", "chr", "window.id"))

#crude unstack of poly and piece
windows.poly=windows.for.theta
windows.poly$ratemtd="poly"
windows.poly=windows.poly[,c("species", "chr", "window.id", "poly.rr", "pi", "fd", "use", "filt", "ratemtd")]
windows.piece=windows.for.theta
windows.piece$ratemtd="piece"
windows.piece=windows.piece[,c("species", "chr", "window.id", "piece.rr", "pi", "fd", "use", "filt", "ratemtd")]
names(windows.piece)[4]="rr"
names(windows.poly)[4]="rr"
windows.for.theta=rbind(windows.piece, windows.poly)

####WRITE OUT 8 FILES FOR ANALYSIS#######
#windows for theta <- recombination rate data and masking data for each window
#windows for gk <- one file each per window size and rec estimator
#mutation rates <- mutation rate options to use for each species

write.table(windows.for.theta, file="windows_for_theta.out", row.names=F, quote=F, sep="\t")

write.table(wind100.poly.for.gk, file="wind100_poly_forgk.out", row.names=F, quote=F, sep="\t", col.names=F)
write.table(wind500.poly.for.gk, file="wind500_poly_forgk.out", row.names=F, quote=F, sep="\t", col.names=F)
write.table(wind1000.poly.for.gk, file="wind1000_poly_forgk.out", row.names=F, quote=F, sep="\t", col.names=F)
write.table(wind100.piece.for.gk, file="wind100_piece_forgk.out", row.names=F, quote=F, sep="\t", col.names=F)
write.table(wind500.piece.for.gk, file="wind500_piece_forgk.out", row.names=F, quote=F, sep="\t", col.names=F)
write.table(wind1000.piece.for.gk, file="wind1000_piece_forgk.out", row.names=F, quote=F, sep="\t", col.names=F)

write.table(u.df, file="mut_ests.txt", row.names=F, quote=F, sep="\t")
write.table(merged.r2, file="rec_fits.txt", row.names=F, col.names=T, quote=F, sep="\t")

############CODE ENDS####################

