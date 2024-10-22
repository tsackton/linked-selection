#script to summarize selection results
source("../linkedsel_functions.R")

#make data frame of all combinations to summarize
selfing.df<-read.table("selfing.txt", header=T, stringsAsFactors=F)
run.levels<-data.frame(wind=rep(c("100", "500", "1000"),6), U=c("min", "min", "min", "const", "const", "const", "max", "max", "max"), rec=c(rep("piece",9),rep("poly",9)))
runs.to.process<-merge(selfing.df, run.levels, by=NULL)

#process with plyr
require(plyr)
selection.results<-ddply(runs.to.process, .(spec, wind, U, rec, selfing), splat(get_mod_results))

#best pi and mean pi are very highly correlated:
plot(log10(selection.results$best.pi), log10(selection.results$mean.pi))
cor.test(selection.results$best.pi,selection.results$mean.pi, method="p")

#subset to only keep best pi for consistency with ll and aic
selection.results=selection.results[,c("spec", "selfing", "rec", "wind", "filt", "U", "model", "best.pi", "best.ll", "best.aic")]

#reformat
require(reshape2)
select.modset.1<-droplevels(subset(selection.results, model=="neutral" | model=="hh_only(var)" | model=="bgs_only(var)" | model=="full(hh_var)"))
select.modset.2<-droplevels(subset(selection.results, model=="neutral" | model=="hh_only(const)" | model=="bgs_only(var)" | model=="full(hh_const)"))
select.modset.1$mod.set="set1" #set1 is variable hitchhiking, variable bgs
select.modset.2$mod.set="set2" #set2 is constant hitchhiking, variable bgs
select.modset.merge<-rbind(select.modset.1, select.modset.2)
select.modset.merge$model=sub("\\(\\w+\\)","", select.modset.merge$model)
select.modset.merge$model=sub("_only", "", select.modset.merge$model)
select.modset.minaic=ddply(.data=select.modset.merge, .variables=.(spec, selfing, rec, wind, filt, U, mod.set), summarise, min.aic=min(best.aic, na.rm=T))
select.modset.merge=merge(select.modset.merge, select.modset.minaic)
select.modset.merge$rel.lik=exp((select.modset.merge$min.aic-select.modset.merge$best.aic)/2)

#three different pis for each model set: best by AIC, weighted mean, full model only
select.final.pi<-ddply(.data=select.modset.merge, .variables=.(spec, selfing, rec, wind, filt, U, mod.set), summarise, pi.best=best.pi[best.aic==min.aic], pi.wt=weighted.mean(best.pi, rel.lik), pi.full=best.pi[model=="full"])

#select.final.pi now has, for each species, 72 different data subsets with a best, full, and weighted pi for each

#now load in range and height data
ne.cors<-read.table("../spec_props/ne_cors.txt", header=T, sep="\t")
names(ne.cors)[1]="spec"
range<-read.table("../spec_props/range.final", header=T, sep="\t")
names(range)=c("spec", "area")
ne.cors<-merge(ne.cors, range)
ne<-merge(ne.cors, select.final.pi, by.x="spec", by.y="spec")

#nuisance parameters
map.params<-read.table("../rec_rate_est/map_info_final.txt", header=T, sep="\t")
qc.params<-read.table("../spec_props/qc_params.txt", header=T, sep="\t")

ne<-merge(ne, map.params, by.x="spec", by.y="sp")
ne<-merge(ne, qc.params, by.x="spec", by.y="sp")

#load in average pi calculations
std.pi<-read.table("../prepare_input_files/std.average_pi.txt", header=F)
q30.pi<-read.table("../prepare_input_files/q30.average_pi.txt", header=F)
std.pi$filt="std"
q30.pi$filt="q30"
avg.pi=rbind(std.pi, q30.pi)
names(avg.pi)=c("spec", "pi.avg", "sites", "filt")
avg.pi$spec = tolower(avg.pi$spec)
avg.pi$spec = substring(avg.pi$spec, 1, 4)
avg.pi$spec[avg.pi$spec=="bman"] = "bmor"
avg.pi$spec[avg.pi$spec=="cret"] = "ccle"
avg.pi$spec[avg.pi$spec=="mmca"] = "mmus"
avg.pi$spec[avg.pi$spec=="ocan"] = "oari"
avg.pi$spec[avg.pi$spec=="oruf"] = "osat"
avg.pi$spec[avg.pi$spec=="pdav"] = "pper" 
avg.pi=avg.pi[,c("spec", "filt", "pi.avg")]

ne<-merge(ne, avg.pi)

#set up 'intensity of selection'
ne$selstr.wt=pmax(0, 1-(ne$pi.avg/ne$pi.wt))
ne$selstr.full=pmax(0, 1-(ne$pi.avg/ne$pi.full))
ne$selstr.best=pmax(0, 1-(ne$pi.avg/ne$pi.best))

#again, weighted vs best makes very little difference; use best again for consistency
ne=subset(ne, select=-c(pi.wt, selstr.wt))

#use poly fit if cM density good is >5, otherwise use piece
ne=subset(ne, (density.good > 2 & rec == "poly") | (density.good <= 2 & rec=="piece"))

#ne is our final dataset, write out
write.table(ne, "ne_final_table.out", row.names=F, quote=F, sep="\t")

