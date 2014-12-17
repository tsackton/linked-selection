#packages
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

#basic default model
mod.bio.full<-lm(selstr.best ~ logrange+logsize+kingdom+logsize:kingdom+logrange:kingdom+logsize:logrange+logsize:logrange:kingdom, data=ne.main)
drop1(mod.bio.full, test="F") #test 3-way interaction
mod.bio.2nd<-lm(selstr.best ~ logrange+logsize+kingdom+logsize:kingdom+logrange:kingdom+logsize:logrange, data=ne.main)
drop1(mod.bio.2nd, test="F") #test all 2-way interactions
mod.bio.final<-lm(selstr.best ~ logrange+logsize+kingdom+logsize:kingdom, data=ne.main)

##ROBUSTNESS TESTS###

#nuisance + bio model
mod.nuis.all<-lm(selstr.best ~ genome.size.mb + prop.good + placed.frac + rec.rate + ungapped.frac + density.good, data=ne.main)
drop1(mod.nuis.all, test="F") # test all nuisance parameters
mod.nuis.bio<-lm(selstr.best ~ logrange+logsize+kingdom+logsize:kingdom+genome.size.mb + prop.good + placed.frac + rec.rate + ungapped.frac + density.good,data=ne.main)

#test whether bio parameters are significant after including nuisance parameters
anova(mod.nuis.all, mod.nuis.bio, test="F")

#residual model final
mod.res.final<-lm(residuals(lm(selstr.best ~ genome.size.mb + prop.good + placed.frac + rec.rate + ungapped.frac + density.good,data=ne.main))~ne.main$logrange+ne.main$logsize+ne.main$kingdom+ne.main$logsize:ne.main$kingdom)
mod.res.nogs<-lm(residuals(lm(selstr.best ~  prop.good + placed.frac + rec.rate + ungapped.frac + density.good,data=ne.main))~ne.main$logrange+ne.main$logsize+ne.main$kingdom+ne.main$logsize:ne.main$kingdom)

mod.final.nodom<-lm(selstr.best ~ logrange+logsize+kingdom+logsize:kingdom, data=ne.main[ne.main$domesticated=="no",])
mod.final.animal<-lm(selstr.best ~ logrange+logsize, data=ne.main[ne.main$kingdom=="animal",])
mod.final.plant<-lm(selstr.best ~ logrange+logsize, data=ne.main[ne.main$kingdom=="plant",])

#now we can getadjR2 from models for each variable set
robustness.test.main<-ddply(ne, .(wind, filt, U, mod.set), summarize, adjr2=summary(lm(selstr.best ~ logsize+logrange+kingdom+logsize:kingdom))$adj.r.squared, fstat=summary(lm(selstr.best ~ logsize+logrange+kingdom+logsize:kingdom))$fstatistic[1], numdf=summary(lm(selstr.best ~ logsize+logrange+kingdom+logsize:kingdom))$fstatistic[2], dendf=summary(lm(selstr.best ~ logsize+logrange+kingdom+logsize:kingdom))$fstatistic[3])
robustness.test.main$pvalue=apply(robustness.test.main[,c(6,7,8)], 1, function(x) pf(x[1], x[2], x[3], lower.tail=F))

#FIGURE 3#

#set up colors
ne.main$plotcol = NA
ne.main$plotcol[ne.main$kingdom=="animal"] = "blue1"
#ne.main$plotcol[ne.main$spec == "agam" | ne.main$spec=="amel" | ne.main$spec=="bmor" | ne.main$spec=="cbri" | ne.main$spec=="cele" | ne.main$spec == "dmel" | ne.main$spec=="dpse" | ne.main$spec == "hmel"] = "blue3"
ne.main$plotcol[ne.main$kingdom=="plant"] = "green4"
#ne.main$plotcol[ne.main$spec == "ccle" | ne.main$spec=="grai" | ne.main$spec=="pper" | ne.main$spec == "ptri"] = "green4"

#set up layout
fig3.out<-rbind(c(1,2),c(3,3))
pdf(file="Figure3.pdf")
layout(fig3.out)

#part A
par(mar=c(4,5,2,1))
par(xpd=F)
plot(ne.main$selstr.best ~ ne.main$logrange, col=ne.main$plotcol, pch=16, xlab=expression('Log'[10]*' Range (sq km)'), ylab="Impact of Selection", las=1, ylim=c(0,0.8), bty="l")
abline(lm(selstr.best ~ logrange, data=ne.main[ne.main$kingdom=="plant",]), col="green4")
abline(lm(selstr.best ~ logrange, data=ne.main[ne.main$kingdom=="animal",]), col="blue1")
mtext("A", 3, at=c(2.7), line=0.5)
legend("topleft", legend=c("Animals", "Plants"), pch=16, col=c("blue1", "green4"), bty="n")

#part B
par(mar=c(4,5,2,1))
plot(ne.main$selstr.best ~ ne.main$logsize, col=ne.main$plotcol, pch=16, xlab=expression('Log'[10]*' Size (meters)'), ylab="Impact of Selection", las=1, ylim=c(0,0.8), bty='l')
abline(lm(selstr.best ~ logsize, data=ne.main[ne.main$kingdom=="plant",]), col="darkgreen")
abline(lm(selstr.best ~ logsize, data=ne.main[ne.main$kingdom=="animal",]), col="blue")
mtext("B", 3, at=c(-3.5), line=0.5)

#part C
par(mar=c(8,5,3,1))
par(xpd=NA)
robustness.test.main$U=factor(robustness.test.main$U, levels=c("min", "const", "max"))
robustness.test.main$plotcol=NA
robustness.test.main$plotcol[robustness.test.main$pvalue<0.001 & is.na(robustness.test.main$plotcol)]="orangered3"
robustness.test.main$plotcol[robustness.test.main$pvalue<0.01 & is.na(robustness.test.main$plotcol)]="darkorange"
robustness.test.main$plotcol[robustness.test.main$pvalue<0.05 & is.na(robustness.test.main$plotcol)]="goldenrod1"
robustness.test.main$plotcol[is.na(robustness.test.main$plotcol)]="gray50"
robustness.test.main=robustness.test.main[order(robustness.test.main$U, robustness.test.main$wind, robustness.test.main$filt, robustness.test.main$mod.set),]
robustness.test.main$label=with(robustness.test.main, paste(U,filt,mod.set,wind,sep="."))
plot(robustness.test.main$adjr2, xaxt="n", ylab=expression("Adjusted R"^2), col="gray20", type="h", las=2, ylim=c(0,0.8), bty="l", xlab="")
points(robustness.test.main$adjr2, pch=16, type="p", col=robustness.test.main$plotcol)
axis(1, labels=c(as.character(robustness.test.main$label)), at=c(1:length(robustness.test.main$label)), las=2, cex.axis=0.65)
legend(x=29, y=0.95, col=c("orangered3", "darkorange", "goldenrod1", "gray50"), legend=c(expression("P "<=" 0.001"), expression("0.001 < P "<=" 0.01"), expression("0.01 < P "<=" 0.05"), "P > 0.05"), pch=16, bty="n")
mtext("Model Parameters (U, Filtering, Model Set, Window Size)", side=1, line=6, cex=0.8)
mtext("C", side=3, line=0.5, at=-2)

dev.off() #close PDF

