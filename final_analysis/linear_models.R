

#now we can getadjR2 from models for each variable set
adj.r2.interact.best<-ddply(ne, .(rec, wind, filt, U, mod.set), summarize, adjr2=summary(lm(selstr.best ~ log10(size.m) + log10(area) + log10(size.m):kingdom + log10(area):kingdom))$adj.r.squared, fstat=summary(lm(selstr.best ~ log10(size.m) + log10(area) + log10(size.m):kingdom + log10(area):kingdom))$fstatistic[1], numdf=summary(lm(selstr.best ~ log10(size.m) + log10(area) + log10(size.m):kingdom + log10(area):kingdom))$fstatistic[2], dendf=summary(lm(selstr.best ~ log10(size.m) + log10(area) + log10(size.m):kingdom + log10(area):kingdom))$fstatistic[3])
adj.r2.interact.best$pvalue=apply(adj.r2.interact.best[,c(7,8,9)], 1, function(x) pf(x[1], x[2], x[3], lower.tail=F))

adj.r2.interact.full<-ddply(ne, .(rec, wind, filt, U, mod.set), summarize, adjr2=summary(lm(selstr.full ~ log10(size.m) + log10(area) + log10(size.m):kingdom + log10(area):kingdom))$adj.r.squared, fstat=summary(lm(selstr.full ~ log10(size.m) + log10(area) + log10(size.m):kingdom + log10(area):kingdom))$fstatistic[1], numdf=summary(lm(selstr.full ~ log10(size.m) + log10(area) + log10(size.m):kingdom + log10(area):kingdom))$fstatistic[2], dendf=summary(lm(selstr.full ~ log10(size.m) + log10(area) + log10(size.m):kingdom + log10(area):kingdom))$fstatistic[3])
adj.r2.interact.full$pvalue=apply(adj.r2.interact.full[,c(7,8,9)], 1, function(x) pf(x[1], x[2], x[3], lower.tail=F))


adj.r2.sep.best<-ddply(ne, .(rec, wind, filt, U, mod.set, kingdom), summarize, adjr2=summary(lm(selstr.best ~ log10(size.m) + log10(area)))$adj.r.squared, fstat=summary(lm(selstr.best ~ log10(size.m) + log10(area)))$fstatistic[1], numdf=summary(lm(selstr.best ~ log10(size.m) + log10(area)))$fstatistic[2], dendf=summary(lm(selstr.best ~ log10(size.m) + log10(area)))$fstatistic[3])
adj.r2.sep.best$pvalue=apply(adj.r2.sep.best[,c(8,9,10)], 1, function(x) pf(x[1], x[2], x[3], lower.tail=F))



