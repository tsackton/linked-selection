#correlations
require(plyr)
require(ppcor)
require(coin)

poly<-read.table("../rec_rate_est/windows_for_theta.out", header=T, stringsAsFactors=F)
map.params<-read.table("../rec_rate_est/map_info_final.txt", header=T, sep="\t")
poly<-merge(poly, map.params, by.x="species", by.y="sp")

#use poly fit if cM density good is >2, otherwise use piece
poly=subset(poly, (density.good > 2 & ratemtd == "poly") | (density.good <= 2 & ratemtd=="piece"))
poly$wind = sub("_\\d+", "", poly$window.id)

#subset for correlation tests
poly.for.cor<-poly[poly$use==T,c("filt", "wind", "species", "rr", "pi", "fd")]
poly.for.cor<-poly.for.cor[complete.cases(poly.for.cor),]
cor.res<-ddply(poly.for.cor, .(species, filt, wind), summarise, tau=pcor.test(rr,pi,fd,method="kendall")$estimate, pvalue=pcor.test(rr,pi,fd,method="kendall")$p.value)

class=read.table("../spec_props/ne_cors.txt", header=T, sep="\t")
class=class[,c("sp", "kingdom", "class.inter")]
cor.res=merge(cor.res, class, by.x="species", by.y="sp")

#main analysis - wilcox test
cor.res.main = subset(cor.res, filt=="std" & wind=="500")
wilcox_test(cor.res.main$tau[cor.res.main$kingdom=="animal"]~droplevels(cor.res.main$class.inter[cor.res.main$kingdom=="animal"]))
wilcox_test(cor.res.main$tau[cor.res.main$kingdom=="plant"]~droplevels(cor.res.main$class.inter[cor.res.main$kingdom=="plant"]))

#main analysis - permutation test
#permutation test
num.perms=100000
ne.animal = subset(cor.res.main, cor.res.main$kingdom == "animal")
animal.diff.mean = mean(cor.res.main$tau[cor.res.main$class.inter=="invert"]) - mean(cor.res.main$tau[cor.res.main$class.inter=="vert"])
animal.diff.rand = NULL
for (i in 1:num.perms) {
  invert.samp = sample(ne.animal$tau, sum(ne.animal$class.inter=="invert"), TRUE)
  vert.samp = sample(ne.animal$tau, sum(ne.animal$class.inter=="vert"), TRUE)
  animal.diff.rand[i]=mean(invert.samp) - mean(vert.samp)
}
sum(animal.diff.rand >= animal.diff.mean)/num.perms
sum(abs(animal.diff.rand) >= abs(animal.diff.mean))/num.perms

ne.plant = subset(cor.res.main, cor.res.main$kingdom == "plant")
plant.diff.mean = mean(cor.res.main$tau[cor.res.main$class.inter=="herb"]) - mean(cor.res.main$tau[cor.res.main$class.inter=="woody"])
plant.diff.rand = NULL
for (i in 1:num.perms) {
  herb.samp = sample(ne.plant$tau, sum(ne.plant$class.inter=="herb"), TRUE)
  woody.samp = sample(ne.plant$tau, sum(ne.plant$class.inter=="woody"), TRUE)
  plant.diff.rand[i]=mean(herb.samp) - mean(woody.samp)
}
sum(plant.diff.rand >= plant.diff.mean)/num.perms
sum(abs(plant.diff.rand) >= abs(plant.diff.mean))/num.perms


#robustness
cor.res.robust<-ddply(cor.res, .(kingdom, filt, wind), summarise, pvalue=pvalue(wilcox_test(tau ~ droplevels(class.inter))), effect.size=median(tau[class.inter=="invert" | class.inter=="herb"])-median(tau[class.inter=="vert" | class.inter=="woody"]))

#figure 2
pdf(file="Figure2.pdf")
cor.res.main$class.inter=factor(cor.res.main$class.inter, levels=c("invert", "vert", "herb", "woody"))
with(cor.res.main, stripchart(tau ~ class.inter, vert=T, ylim=c(-0.25,0.60), las=1, pch=16, col="gray40", xlab="Class", ylab="Partial Correlation (Rec Rate vs Neutral Diversity)", group.names=c("Invertebrates", "Vertebrates", "Herbaceous Plants", "Woody Plants"), frame.plot=F))
points(y=c(by(cor.res.main$tau, cor.res.main$class.inter, median)), x=c(1,2,3,4), pch=18, col="red", cex=2)
lines(y=c(0.57,0.57),x=c(1,2))
lines(y=c(0.57,0.57),x=c(3,4))
text(y=0.59, x=1.5, label="P=0.048", cex=1)
text(y=0.59, x=3.5, label="P=0.026", cex=1)
dev.off()
