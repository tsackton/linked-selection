#load functions
source("../linkedsel_functions.R")
require(plyr)

#load in average pi calculations
avg.pi<-read.table("../prepare_input_files/std.average_pi.txt", header=F)
names(avg.pi)=c("spec", "pi.avg", "sites")
avg.pi$spec = tolower(avg.pi$spec)
avg.pi$spec = substring(avg.pi$spec, 1, 4)
avg.pi$spec[avg.pi$spec=="bman"] = "bmor"
avg.pi$spec[avg.pi$spec=="cret"] = "ccle"
avg.pi$spec[avg.pi$spec=="mmca"] = "mmus"
avg.pi$spec[avg.pi$spec=="ocan"] = "oari"
avg.pi$spec[avg.pi$spec=="oruf"] = "osat"
avg.pi$spec[avg.pi$spec=="pdav"] = "pper" 
avg.pi=avg.pi[,c("spec", "pi.avg")]

#set up plot
pdf(file="Figure1.pdf", width=8, height=8)
par(mfrow=c(2,1))

#PART A##
##ugly copy and pasted code starts below#
spec.to.plot="dmel"

plot.rec=read.table("../rec_rate_est/windows_for_theta.out", header=T, sep="\t")
plot.rec=droplevels(subset(plot.rec, species==spec.to.plot & grepl("500_", window.id, fixed=T) & use==T & filt=="std" & ratemtd=="piece"))
plot.rec$U="min"
plot.rec$wind=500
plot.rec$selfing="no"
plot.rec$rec="piece"
names(plot.rec)[1]="spec"

plot.sel=ddply(unique(plot.rec[,c("spec", "wind", "U", "rec", "selfing")]), .(spec, wind, U, rec, selfing), splat(get_mod_results_plot))
plot.sel=subset(plot.sel, model=="full(hh_var)" & filt=="std")

plot.gk=read.table("/Volumes/LaCie/Projects/Current/ne/final_data/gk/dmel.wind500_piece.min.gk.gz", header=F)
names(plot.gk)=c("spec", "chr", "window.id", "U", "sh", "P", "G")
plot.gk=droplevels(subset(plot.gk, P==1 & sh==plot.sel$sh))

plot.np=plot.sel$best.pi
plot.ap=avg.pi[avg.pi$spec==spec.to.plot, c("pi.avg")]

pred.plot = merge(plot.gk, plot.rec, by=c("spec", "chr", "window.id"))
pred.plot$rho = (1-exp(-1*(pred.plot$rr/1e6)*2/100))/2
pred.plot$alpha = plot.sel$param2
pred.plot$np = plot.np #neutral pi
pred.plot$ap = plot.ap #average pi
pred.plot$pp = (pred.plot$np*exp(-1*pred.plot$G)*pred.plot$rho)/(pred.plot$rho + (pred.plot$alpha*pred.plot$fd*exp(-1*pred.plot$G))) #predicted pi for each window
pred.plot$pi[pred.plot$pi>quantile(pred.plot$pi, 0.95)]=quantile(pred.plot$pi, 0.95) #for plotting truncate pi dist
pred.plot$reclass = cut(pred.plot$rr, breaks=quantile(pred.plot$rr, seq(0,1,0.02)), labels=F)
pred.plot$reclass[pred.plot$rr==0] = 0
pred.pi = by(pred.plot$pp, pred.plot$reclass, mean)
rec.rate = by(pred.plot$rr, pred.plot$reclass, mean)

par(mar=c(5,5,2,6), xpd=F)
plot(pi ~ rr, las=1, data=pred.plot, xlab="Recombination Rate (cM/Mb)", ylab="", bty="l", cex=0.5, col="gray70")
abline(h=unique(pred.plot$np), col="red", lwd=1.5)
abline(h=unique(pred.plot$ap), col="black", lty="dashed", lwd=1.5)
points(pred.pi ~ rec.rate, type="l", col="blue", lwd=2)
mtext("Neutral Diversity", 2, line=4)
mtext(text="Predicted Neutral Diversity", side=4, line=-3, at=unique(pred.plot$np)*1.05, cex=0.8, las=1, col="red")
mtext(text="Observed Neutral Diversity", side=4, line=-3, at=unique(pred.plot$ap)*1.05, cex=0.8, las=1, col="black")
mtext(text="Model Fit", side=4, line=-3, at=max(pred.plot$pp)*0.95, cex=0.8, las=1, col="blue")
mtext(text="A", side=3, line=0.5, at=-0.3)

##PART B##
spec.to.plot="ecab"

plot.rec=read.table("../rec_rate_est/windows_for_theta.out", header=T, sep="\t")
plot.rec=droplevels(subset(plot.rec, species==spec.to.plot & grepl("500_", window.id, fixed=T) & use==T & filt=="std" & ratemtd=="piece"))
plot.rec$U="min"
plot.rec$wind=500
plot.rec$selfing="no"
plot.rec$rec="piece"
names(plot.rec)[1]="spec"

plot.sel=ddply(unique(plot.rec[,c("spec", "wind", "U", "rec", "selfing")]), .(spec, wind, U, rec, selfing), splat(get_mod_results_plot))
plot.sel=subset(plot.sel, model=="full(hh_var)" & filt=="std")

plot.gk=read.table("/Volumes/LaCie/Projects/Current/ne/final_data/gk/ecab.wind500_piece.min.gk.gz", header=F)
names(plot.gk)=c("spec", "chr", "window.id", "U", "sh", "P", "G")
plot.gk=droplevels(subset(plot.gk, P==1 & sh==plot.sel$sh))

plot.np=plot.sel$best.pi
plot.ap=avg.pi[avg.pi$spec==spec.to.plot, c("pi.avg")]

pred.plot = merge(plot.gk, plot.rec, by=c("spec", "chr", "window.id"))
pred.plot$rho = (1-exp(-1*(pred.plot$rr/1e6)*2/100))/2
pred.plot$alpha = plot.sel$param2
pred.plot$np = plot.np #neutral pi
pred.plot$ap = plot.ap #average pi
pred.plot$pp = (pred.plot$np*exp(-1*pred.plot$G)*pred.plot$rho)/(pred.plot$rho + (pred.plot$alpha*pred.plot$fd*exp(-1*pred.plot$G))) #predicted pi for each window
pred.plot$pi[pred.plot$pi>quantile(pred.plot$pi, 0.95)]=quantile(pred.plot$pi, 0.95) #for plotting truncate pi dist
pred.plot$reclass = cut(pred.plot$rr, breaks=quantile(pred.plot$rr, seq(0,1,0.02)), labels=F)
pred.plot$reclass[pred.plot$rr==0] = 0
pred.pi = by(pred.plot$pp, pred.plot$reclass, mean)
rec.rate = by(pred.plot$rr, pred.plot$reclass, mean)

par(mar=c(5,5,2,6), xpd=F)
plot(pi ~ rr, las=1, data=pred.plot, xlab="Recombination Rate (cM/Mb)", ylab="", bty="l", cex=0.5, col="gray70")
abline(h=unique(pred.plot$np), col="red", lwd=1.5)
abline(h=unique(pred.plot$ap), col="black", lty="dashed", lwd=1.5)
points(pred.pi ~ rec.rate, type="l", col="blue", lwd=2)
mtext("Neutral Diversity", 2, line=4)
mtext(text="Predicted Neutral Diversity", side=4, line=-3.2, at=unique(pred.plot$np)*1.35, cex=0.8, las=1, col="red")
mtext(text="Observed Neutral Diversity", side=4, line=-2.7, at=unique(pred.plot$ap)*1.12, cex=0.8, las=1, col="black")
mtext(text="Model Fit", side=4, line=-2.9, at=max(pred.plot$pp)*0.85, cex=0.8, las=1, col="blue")
mtext(text="B", side=3, line=0.5, at=-0.3)

dev.off()
