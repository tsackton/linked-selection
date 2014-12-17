#best model selection
source("../linkedsel_functions.R")
require(plyr)

#load final table for analysis
ne<-read.table("../pogen_models/ne_final_table.out", header=T, stringsAsFactors=F)

#update qc params
qc.params<-read.table("../spec_props/qc_params.txt", header=T, sep="\t")
ne=ne[,-c(26,27,28,29,30)]
ne<-merge(ne, qc.params, by.x="spec", by.y="sp")
rm(qc.params)

#define some additional variables
ne$placed.frac = ne$placed.bp / ne$assembly.bp
ne$ungapped.frac = ne$placed.ungapped.bp / ne$placed.bp
ne$logsize = log10(ne$size.m)
ne$logrange = log10(ne$area)

#main subset described in manuscript
#500kb windows, standard filtering, U estimated as mutation rate * exonic bases, and variable hitchhiking model (accounting for differences in gene density across windows)

ne.main<-subset(ne, wind==500 & filt=="std" & mod.set=="set1" & U=="min")
####
##REQUIRES THE RAW MODEL FITS FROM THE SELECTION MODEL
##TO REPLICATE OUR ANALYSIS WITHOUT RERUNNING THE GK AND POPGEN MODEL CODE SKIP TO LINE 70 AND UNCOMMENT THE READ TABLE STATEMENT
####
select.results<-ddply(ne.main[,c("spec", "wind", "U", "rec", "selfing")], .(spec, wind, U, rec, selfing), splat(get_mod_results))
select.results=subset(select.results, filt=="std")

#reformat
require(reshape2)
select.modset.1<-droplevels(subset(select.results, model=="neutral" | model=="hh_only(var)" | model=="bgs_only(var)" | model=="full(hh_var)"))
select.modset.2<-droplevels(subset(select.results, model=="neutral" | model=="hh_only(const)" | model=="bgs_only(var)" | model=="full(hh_const)"))
select.modset.1$mod.set="set1" #set1 is variable hitchhiking, variable bgs
select.modset.2$mod.set="set2" #set2 is constant hitchhiking, variable bgs
select.modset.merge<-rbind(select.modset.1, select.modset.2)
select.modset.merge$model=sub("\\(\\w+\\)","", select.modset.merge$model)
select.modset.merge$model=sub("_only", "", select.modset.merge$model)

#get relative likelihoods of each 
rel.lik.models=ddply(.data=select.modset.merge, .variables=.(spec, mod.set), summarise, min.aic=min(best.aic, na.rm=T))
rel.lik.models=merge(rel.lik.models, select.modset.merge)
rel.lik.models=rel.lik.models[,c("spec", "mod.set", "model", "best.aic", "best.ll", "min.aic")]
rel.lik.models$rel.lik=exp((rel.lik.models$min.aic - rel.lik.models$best.aic)/2)
norm.factor=ddply(.data=rel.lik.models, .variables=.(spec,mod.set), summarise, norm.factor=sum(rel.lik))
rel.lik.models=merge(rel.lik.models, norm.factor)
rel.lik.models$norm.lik=rel.lik.models$rel.lik/rel.lik.models$norm.factor

#hh vs bgs
rel.lik.hhvsbgs=subset(rel.lik.models, mod.set=="set1" & (model=="bgs" | model=="full"))
hhvsbgs.final=ddply(.data=rel.lik.hhvsbgs, .variables=.(spec), summarise, min.aic=min(best.aic), full.aic=best.aic[model=="full"], bgs.aic=best.aic[model=="bgs"])
hhvsbgs.final$hh.rel=exp((hhvsbgs.final$min.aic - hhvsbgs.final$full.aic)/2)
hhvsbgs.final$bgs.rel=exp((hhvsbgs.final$min.aic - hhvsbgs.final$bgs.aic)/2)
hhvsbgs.final$hh.v.bgs = log2(hhvsbgs.final$hh.rel/hhvsbgs.final$bgs.rel)

#get best model for each spec, mod.set
best.model=ddply(.data=rel.lik.models, .variables=.(spec,mod.set), summarise, best.model=model[rel.lik==1], rel.lik.best=norm.lik[rel.lik==1], rel.lik.hh=norm.lik[model=="hh"]+norm.lik[model=="full"], rel.lik.bgs=norm.lik[model=="bgs"], rel.lik.neut=norm.lik[model=="neutral"])

#add best model info to ne.main
ne.best=merge(ne.main, best.model)
ne.best$best.model.plot=ne.best$best.model
ne.best$best.model.plot[ne.best$best.model=="full"]="hh"
ne.best$best.model.plot=factor(ne.best$best.model.plot, levels=c("neutral", "bgs", "hh"))

ne.best$hh.conf=cut(ne.best$rel.lik.hh, breaks=c(0,0.05,0.9,1), include.lowest=T, right=T, labels=c("low", "med", "high"))
ne.best$neut.conf=cut(ne.best$rel.lik.neut, breaks=c(0,0.05,0.9,1), include.lowest=T, right=T, labels=c("low", "med", "high"))
ne.best$bgs.conf=cut(ne.best$rel.lik.bgs, breaks=c(0,0.05,0.9,1), include.lowest=T, right=T, labels=c("low", "med", "high"))
##UNCOMMENT LINE BELOW TO LOAD NE.BEST DATA FRAME
#ne.best<-read.table("ne_best_table.out", header=T, sep="\t")

#significance tests, categorical
require(coin)
wilcox_test(ne.best$logrange ~ as.factor(ne.best$neut.conf != "low"))
wilcox_test(ne.best$logsize ~ as.factor(ne.best$neut.conf != "low"))
wilcox_test(ne.best$logsize[ne.best$neut.conf!="low"] ~ droplevels(ne.best$neut.conf[ne.best$neut.conf!="low"]))
wilcox_test(ne.best$logrange[ne.best$neut.conf!="low"] ~ droplevels(ne.best$neut.conf[ne.best$neut.conf!="low"]))
wilcox_test(ne.best$logsize[ne.best$neut.conf!="high"] ~ droplevels(ne.best$neut.conf[ne.best$neut.conf!="high"]))
wilcox_test(ne.best$logrange[ne.best$neut.conf!="high"] ~ droplevels(ne.best$neut.conf[ne.best$neut.conf!="high"]))
wilcox_test(ne.best$logsize[ne.best$neut.conf!="med"] ~ droplevels(ne.best$neut.conf[ne.best$neut.conf!="med"]))
wilcox_test(ne.best$logrange[ne.best$neut.conf!="med"] ~ droplevels(ne.best$neut.conf[ne.best$neut.conf!="med"]))

#hh
ne.best=merge(ne.best, hhvsbgs.final)
ne.best$rangeclass = cut(ne.best$logrange, c(-100,median(ne.best$logrange), 100), labels=c("low", "high"))
ne.best$sizeclass = cut(ne.best$logsize, c(-100,median(ne.best$logsize), 100), labels=c("low", "high"))
ne.best$nc.class = paste(ne.best$sizeclass, ne.best$rangeclass, sep=".")
ne.best$nc.class[ne.best$nc.class=="high.low"]="low"
ne.best$nc.class[ne.best$nc.class=="low.high"]="high"
ne.best$nc.class[ne.best$nc.class=="low.low"]="med"
ne.best$nc.class[ne.best$nc.class=="high.high"]="med"
ne.best$nc.class = factor(ne.best$nc.class, levels=c("low", "med", "high"))

ne.best$nc.class.2="med"
ne.best$nc.class.2[ne.best$logsize > quantile(ne.best$logsize, 0.60) & ne.best$logrange < quantile(ne.best$logrange, 0.40)]="low"
ne.best$nc.class.2[ne.best$logsize < quantile(ne.best$logsize, 0.40) & ne.best$logrange > quantile(ne.best$logrange, 0.60)]="high"


wilcox_test(ne.best$logrange ~ as.factor(ne.best$hh.conf != "low"))
wilcox_test(ne.best$logsize ~ as.factor(ne.best$hh.conf != "low"))

wilcox_test(ne.best$rel.lik.hh[ne.best$nc.class!="med"] ~ droplevels(ne.best$nc.class[ne.best$nc.class!="med"]))

##FIGURE 4##

#setup layout
pdf(file="Figure4.pdf", width=8, height=6)
par(mfrow=c(1,2))

#Figure 4A
par(mar=c(5,5,3,1), xpd=NA)
stripchart(ne.best$logrange ~ ne.best$neut.conf, vert=T, pch=16, at=c(0.2,0.7,1.2), xlim=c(0,1.4), ylim=c(3,8), las=1, ylab=expression('Log'[10]*' Range (sq km)'), xlab="Support for Neutral Model", group.names=c("Low", "Medium", "High"), frame.plot=F)
points(x=c(0.2, 0.7, 1.2), y=c(median(ne.best$logrange[ne.best$neut.conf=="low"]), median(ne.best$logrange[ne.best$neut.conf=="med"]), median(ne.best$logrange[ne.best$neut.conf=="high"])), pch=18, col="red", cex=2)
lines(x=c(0.2,0.7), y=c(8.1,8.1), col="black", lwd=1)
lines(x=c(0.2,(0.7+1.2)/2), y=c(8.25,8.25), col="black", lwd=1)
lines(x=c(0.2,1.2), y=c(8.4,8.4), col="black", lwd=1)
text(x=0.45, y=8.15, label="*", cex=1.5)
text(x=0.575, y=8.3, label="**", cex=1.5)
text(x=0.7, y=8.45, label="*", cex=1.5)

#Figure 4B
par(mar=c(5,5,3,1), xpd=NA)
stripchart(ne.best$logsize ~ ne.best$neut.conf, vert=T, pch=16, at=c(0.2,0.7,1.2), xlim=c(0,1.4), ylim=c(-3.5,2), las=1, ylab=expression('Log'[10]*' Size (m)'), xlab="Support for Neutral Model", group.names=c("Low", "Medium", "High"), frame.plot=F)
points(x=c(0.2, 0.7, 1.2), y=c(median(ne.best$logsize[ne.best$neut.conf=="low"]), median(ne.best$logsize[ne.best$neut.conf=="med"]), median(ne.best$logsize[ne.best$neut.conf=="high"])), pch=18, col="red", cex=2)
lines(x=c(0.2,0.7), y=c(1.9,1.9), col="black", lwd=1)
lines(x=c(0.2,(0.7+1.2)/2), y=c(2.05,2.05), col="black", lwd=1)
lines(x=c(0.2,1.2), y=c(2.2,2.2), col="black", lwd=1)
text(x=0.45, y=1.95, label="***", cex=1.5)
text(x=0.575, y=2.10, label="***", cex=1.5)
text(x=0.7, y=2.25, label="*", cex=1.5)

#labels
mtext("A", side=3, line=0.5, at=-2.5)
mtext("B", side=3, line=0.5, at=-.3)

dev.off()
