##FUNCTIONS FOR LINKED SELECTION PROJECT##

#clean up maps
cleanup_maps <- function(x) {
  #require data frame in correct format as input
  stopifnot(is.data.frame(x) && isTRUE(all.equal(colnames(x), c("chr", "marker", "cm", "mb", "sp"))));
  x$chr<-as.character(x$chr)
  x$marker<-as.character(x$marker)
  x$mb<-round(x$mb, 8)
  x$cm<-round(x$cm, 8)
  if (any(grepl("chr", x$chr, fixed=T))) {
    x<-unique(x)
    return(x)
  }
  else {
    x$chr=paste0("chr", x$chr)
    x<-unique(x)
    return(x)
  }
}

#remove mismapped
remove_mismapped <- function(map) {
  #require data frame in correct format as input
  stopifnot(is.data.frame(map) && isTRUE(all.equal(colnames(map), c("chr", "marker", "cm", "mb", "sp"))));
  #use LCS function from qualV package
  require(qualV)
  #remove species column
  spe=unique(map$sp)
  map=map[,-5]
  #set up good marker column
  map$good.marker = F;
  #get marker frequency (unique or multi-mapped)
  map=merge(map, as.data.frame(table(map$marker)), all.x=T, by.x="marker", by.y="Var1")
  names(map)[6] = "marker.freq"
  #order map by MB
  map=map[with(map, order(chr, mb)),]
  num.chr = length(unique(map$chr));
  chr.id = unique(map$chr)
  rownames(map)=seq(1,nrow(map),1)
  #process each chromosome
  for (i in chr.id) {
    #clean map - use longest common subsequence of ranks to define markers that are in the correct order
    r2<-rank(map[map$chr==i,4], ties="first")
    r1<-rank(map[map$chr==i,3], ties="first")
    lcs<-LCS(as.character(r1),as.character(r2))
    row.start=as.numeric(rownames(map[map$chr==i,])[1])-1
    good.rows = lcs$va + row.start
    map[good.rows,"good.marker"] = T
  }
  
  #remove mismapped markers and refilter
  num.chr = length(unique(map$chr));
  chr.id = unique(map$chr)
  rownames(map)=seq(1,nrow(map),1)
  map$sp = spe;
  map = map[,c("chr", "cm", "mb", "marker", "good.marker", "marker.freq", "sp")]
  return(map)
}

#remove duplicates
remove_duplicate_markers<-function(map) {
  stopifnot(is.data.frame(map) && isTRUE(all.equal(colnames(map), c("chr", "cm", "mb", "marker", "good.marker", "marker.freq", "sp"))));
  #remove duplicated markers
  #if there is a stretch of duplicated markers marked true with non-duplicated markers marked false within, set all the dups after the last false non-dup to false and set the non-dups to true, as long as the false markers are properly oriented
  map$spmark = paste(map$sp, map$marker, map$cm, sep=".")
  map = merge(map, as.data.frame(table(map$spmark)), all.x=T, by.x="spmark", by.y="Var1")
  map = map[order(map$chr, map$mb),]
  row.names(map)<-NULL
  
  dup.markers = unique(map$spmark[map$Freq > 1])
  rows.to.delete = numeric(0)
  rows.to.flip = numeric(0)
  
  for (marker in dup.markers) {
    num.true = sum(map$good.marker[map$spmark==marker])
    if (num.true == 0) {
      spec.sub = map[map$spmark==marker,]
      last = length(spec.sub$spmark)
      keep = floor((1+last)/2)
      rows.to.delete=c(rows.to.delete,rownames(spec.sub[-keep,]));		
    }
    else {
      rows.to.delete=c(rows.to.delete,rownames(map[map$spmark==marker & map$good.marker==F,]))
      row1 = min(as.numeric(rownames(map[map$spmark==marker & map$good.marker==T,])))
      row2 = max(as.numeric(rownames(map[map$spmark==marker & map$good.marker==T,])))
      spec.sub = map[as.character(seq(row1,row2)),]
      if (length(unique(spec.sub$spmark))==1) { #just one marker in the set
        last = length(spec.sub$spmark)
        keep = floor((1+last)/2)
        rows.to.delete=c(rows.to.delete,rownames(spec.sub[-keep,]));		
      }
      else {
        low.cm = spec.sub[1,"cm"]
        chr.id = unique(spec.sub$chr)
        stopifnot(length(chr.id)==1)
        high.cm.row = min(c(as.numeric(row2)+1),max(rownames(map[map$chr==chr.id,])))
        while(map[high.cm.row,"good.marker"]==F) {
          high.cm.row = as.character(as.numeric(high.cm.row)+1)
          high.cm.row = min(c(high.cm.row),max(as.numeric(rownames(map[map$chr==chr.id & map$good.marker==T,]))))
        }
        high.cm = map[high.cm.row, "cm"]
        spec.sub=spec.sub[-1,]
        rows.to.delete=c(rows.to.delete,rownames(spec.sub[spec.sub$spmark==marker,]))
        new.rows.to.flip = rownames(spec.sub[spec.sub$spmark!=marker & spec.sub$cm>=low.cm & spec.sub$cm <= high.cm & spec.sub$good.marker==F,]);
        if (length(unique(map[new.rows.to.flip,"cm"]))>1) {
          if (cor(map[new.rows.to.flip,"cm"],map[new.rows.to.flip,"mb"])>0) {
            rows.to.flip=c(rows.to.flip,new.rows.to.flip)			
          }
        }
      }	
    }
  }
  if (length(rows.to.flip) > 0) {
    map[rows.to.flip,"good.marker"] = T    
  }
  if (length(rows.to.delete > 0)) {
    map=map[-as.numeric(rows.to.delete),]    
  }
  map=map[,c("sp", "chr", "cm", "mb", "marker", "good.marker")]
  return(map)
}

#make mask of regions to ignore because of high density of low quality markers
make_map_mask<-function(map) {
  options(stringsAsFactors=F)
  stopifnot(is.data.frame(map) && isTRUE(all.equal(colnames(map), c("sp", "chr", "cm", "mb", "marker", "good.marker"))));
  map.mask = data.frame(chr=character(0), start=numeric(0), end=numeric(0))
  num.chr = length(unique(map$chr));
  chr.id = unique(map$chr)	
  for (i in chr.id) {
    chr.map = subset(map, map$chr==i)
    num.good.markers = sum(chr.map$good.marker)
    bad.tol = min(ceiling(num.good.markers*.2),5)
    #ignore whole chromosomes with fewer than 5 good markers
    if (num.good.markers < 5) {
      map.mask = rbind(map.mask, data.frame(chr=i, start=0, end=1000))
    }
    #otherwise look for runs
    else {
      runs = rle(chr.map$good.marker)
      runs$idx = head(cumsum(c(1,runs$lengths)),-1)
      bad.starts = runs$idx[which(runs$lengths>=bad.tol & runs$values==F)]
      bad.ends = bad.starts + runs$lengths[which(runs$lengths>=bad.tol & runs$values==F)] - 1
      if (length(bad.starts > 0)) {
        for (n in 1:length(bad.starts)) {
          bad.mb.start=ifelse(bad.starts[n]==1,0,chr.map[bad.starts[n]-1,"mb"])
          bad.mb.end=ifelse(bad.ends[n]==nrow(chr.map),1000,chr.map[bad.ends[n]+1,"mb"])
          map.mask=rbind(map.mask, data.frame(chr=i, start=bad.mb.start, end=bad.mb.end))		
        }
      }
      #add chr ends, if whole chromosome is not masked already
      chr.min = min(chr.map$mb[chr.map$good.marker])
      chr.max = max(chr.map$mb[chr.map$good.marker])
      map.mask=rbind(map.mask, data.frame(chr=i, start=0, end=chr.min))
      map.mask=rbind(map.mask, data.frame(chr=i, start=chr.max, end=1000))	
    }
  }
  map.mask=map.mask[order(map.mask$chr, map.mask$start, map.mask$end),]
  map.mask = subset(map.mask, !duplicated(map.mask))
  return(map.mask)
}

#estimate recombination rate using a polynomial
est_rec_poly<-function(map=NA, default.deg=NA) {
  #local function to optimize polynomial degree
  find.degree<-function(data) {
    x<-data$mb;
    y<-data$cm;
    polyfit<-function(i) x <- AIC(lm(y~poly(x,i,raw=T)))
    max.deg<-min(c(20,floor(length(unique(y))/3)))
    max.deg<-max(max.deg, 3)
    model.deg<-data.frame(deg=seq(1,max.deg,1))
    model.deg$aic=apply(model.deg, 1, polyfit)
    model.deg=model.deg[complete.cases(model.deg),]
    model.deg=model.deg[order(model.deg$aic, decreasing=F),]
    return(model.deg$deg[1])  
  }
  
  #requires a data frame with chromosome, cm, and mb positions
  stopifnot(is.data.frame(map) && isTRUE(all.equal(colnames(map), c("chr", "cm", "mb"))));
  
  #clean up map and filter to remove chromosomes with fewer than 4 markers
  map=map[order(map$chr, map$mb),]
  map=merge(map, as.data.frame(table(map$chr)), all.x=T, by.x="chr", by.y="Var1")
  map=subset(map, map$Freq >= 3, select=c("chr", "cm", "mb"))
  num.chr = length(unique(map$chr));
  chr.id = unique(map$chr)
  rownames(map)=seq(1,nrow(map),1)
  
  #list for functions for each chromosome
  rec.poly<-vector(num.chr, mode="list") #an optimized polynomial fit
  names(rec.poly)=chr.id
  for (j in chr.id) {
    chr.map = map[map$chr==j,]
    
    #optimize polynomial function
    if (is.na(default.deg)) {
        deg<-find.degree(chr.map); 
    }
    else {
        deg<-default.deg
    }
    poly.fit<-poly(chr.map$mb, deg, raw=T)
    poly.opt<-lm(chr.map$cm ~ poly.fit)
    while (sum(is.na(coef(poly.opt)))>0) {
      deg=deg-1
      poly.fit<-poly(chr.map$mb, deg, raw=T);
      poly.opt<-lm(chr.map$cm ~ poly.fit)
    }
    while (poly.opt$rank >= nrow(chr.map)) {
      deg=deg-1
      poly.fit<-poly(chr.map$mb, deg, raw=T);
      poly.opt<-lm(chr.map$cm ~ poly.fit)
    }
    rec.poly[[j]]$func=poly.opt;
    rec.poly[[j]]$deg=deg;
  }

  return(rec.poly)
}

#estimate recombination rate using a piecewise regression
est_rec_piece<-function(map=NA, default.df=NA) {
  #depends on splines package
  require(splines)
  
  #local function to optimize df in spline
  find.df<-function(data) {
    x<-data$mb;
    y<-data$cm;
    min.knot = 0;
    max.knot = max(x, na.rm=T)*1.05;
    bknots = c(min.knot, max.knot)
    piecefit<-function(i) x <- AIC(lm(y~0+bs(x,df=i,deg=1,Boundary.knots=bknots)))
    mn.df = floor(length(unique(y))/2)
    max.df=max(mn.df,2)
    max.df=min(100,max.df)
    stopifnot(max.df>0)
    model.df=data.frame(df=seq(1,max.df,1))
    model.df$aic=apply(model.df, 1, piecefit)
    model.df=model.df[complete.cases(model.df),]
    model.df=model.df[order(model.df$aic, decreasing=F),]
    return(model.df$df[1])  
  }
  
  
  #requires a data frame with chromosome, cm, and mb positions
  stopifnot(is.data.frame(map) && isTRUE(all.equal(colnames(map), c("chr", "cm", "mb"))));
  
  #clean up map and filter to remove chromosomes with fewer than 4 markers
  map=map[order(map$chr, map$mb),]
  map=merge(map, as.data.frame(table(map$chr)), all.x=T, by.x="chr", by.y="Var1")
  map=subset(map, map$Freq >= 3, select=c("chr", "cm", "mb"))
  num.chr = length(unique(map$chr));
  chr.id = unique(map$chr)
  rownames(map)=seq(1,nrow(map),1)
  
  #list for functions for each chromosome
  rec.piece<-vector(num.chr, mode="list") #an optimized polynomial fit
  names(rec.piece)=chr.id
  
  for (j in chr.id) {
    chr.map = map[map$chr==j,]
    
    #optimize df function
    if (is.na(default.df)) {
      df<-find.df(chr.map); 
    }
    else {
      df<-default.df
    }
    min.knot = 0;
    max.knot = max(chr.map$mb,na.rm=T)*1.05
    bknots = c(min.knot, max.knot)
    piece.opt<-with(chr.map, lm(cm ~ 0 + bs(mb, df=df, degree=1, Boundary.knots=bknots)))
    while (sum(is.na(coef(piece.opt)))>0) {
      df = df-1
      piece.opt<-with(chr.map, lm(cm ~ 0 + bs(mb, df=df, degree=1, Boundary.knots=bknots)))
    }

    while (piece.opt$rank >= nrow(chr.map)) {
      df=df-1
      piece.opt<-with(chr.map, lm(cm ~ 0 + bs(mb, df=df, degree=1, Boundary.knots=bknots)))
    }
    
    rec.piece[[j]]$func=piece.opt;
    rec.piece[[j]]$df=df;
  }
  
  return(rec.piece)
}

est_rr_poly<-function(species, chr, window.pos, ...) {
  require(polynom)  
  chr=unique(chr)
  species=unique(species)
  poly.fun=rec_polynom[[species]][[chr]][["func"]]
  if (!is.null(poly.fun)) {
    p0<-polynomial(coef(poly.fun))
    p1<-deriv(p0)
    mb<-as.numeric(window.pos)/1000000
    Y.poly<-as.function(p1)(mb)
    Y.poly[Y.poly < 0] = 0;
    Y.poly[Y.poly > quantile(Y.poly, probs=0.99, na.rm=T)]= quantile(Y.poly, probs=0.99, na.rm=T) 
  }
  else {
    Y.poly=rep(NA, length(window.pos))    
  }
  return(Y.poly)
}

est_cm_poly<-function(species, chr, window.pos, ...) {
  require(polynom)  
  require(splines)
  chr=unique(chr)
  species=unique(species)
  poly.fun=rec_polynom[[species]][[chr]][["func"]]
  if (!is.null(poly.fun)) {
    mb<-as.numeric(window.pos)/1000000
    poly.pred=predict(poly.fun)
    poly.inc=(poly.pred==cummax(poly.pred))
    mb.inc=model.frame(poly.fun)[,2][,1]
    poly.interp<-splinefun(x=mb.inc[poly.inc], y=poly.pred[poly.inc], method="monoH.FC")
    cm.poly = poly.interp(mb)
  }
  else {
    cm.poly=rep(NA, length(window.pos))    
  }
  return(cm.poly)
}


est_rr_piece<-function(species, chr, window.pos, window.size, ...) {
  require(splines)  
  chr=unique(chr)
  species=unique(species)
  piece.fun=rec_piecewise[[species]][[chr]][["func"]]
  if (!is.null(piece.fun)) {
    mb<-as.numeric(window.pos)/1000000
    Y.piece.data = data.frame(mb=mb-(window.size/2000))
    Y.piece.data = rbind(Y.piece.data, tail(mb, n=1)+(window.size/2000))
    Y.piece.raw = predict(piece.fun, newdata=Y.piece.data)
    Y.piece = diff(Y.piece.raw)/diff(Y.piece.data$mb)
    Y.piece[Y.piece < 0] = 0
    Y.piece[Y.piece>quantile(Y.piece, probs=0.99, na.rm=T)]= quantile(Y.piece, probs=0.99, na.rm=T)
  }
  else {
    Y.piece=rep(NA, length(window.pos)) 
  } 
  return(Y.piece)
}

est_cm_piece<-function(species, chr, window.pos, window.size, ...) {
  require(splines)  
  chr=unique(chr)
  species=unique(species)
  piece.fun=rec_piecewise[[species]][[chr]][["func"]]
  fitted.data=subset(all_maps_final[[species]], good.marker==T)
  fitted.data=fitted.data[fitted.data$chr==chr,]
  if (!is.null(piece.fun)) {
    mb<-as.numeric(window.pos)/1000000
    piece.pred=predict(piece.fun)
    piece.inc=(piece.pred==cummax(piece.pred))
    mb.inc=fitted.data$mb
    if (length(mb.inc) != length(piece.pred)) {
      stop("Fitted data not equal to predicted!")
    }
    piece.interp<-splinefun(x=mb.inc[piece.inc], y=piece.pred[piece.inc], method="monoH.FC")
    cm.piece = piece.interp(mb)
  }
  else {
    cm.piece=rep(NA, length(window.pos))
    
  }
    return(cm.piece)
}

get_mod_results<-function(spec, wind, U, rec, selfing, ...) {
  path="/Volumes/LaCie/Projects/Current/ne/final_data/modres/"
  file=paste(spec,wind,U,rec,"select.out", sep="_")
  print(file)
  mod.run<-read.table(paste0(path,file),header=T)
  #remove ill-conditioned fits (pi.neut = 1 is boundary of search space)
  mod.run=droplevels(subset(mod.run, model != "bgs_only(const)" & pi.neut < 1 & !is.na(pi.neut))) 
  if (is.na(selfing)) {
    mod.run=mod.run 
  } else if (selfing=="no") {
    mod.run=mod.run[!is.na(mod.run$P) & mod.run$P==1,]
  } else if (selfing=="partial") {
    mod.run=mod.run[!is.na(mod.run$P) & mod.run$P>=0.60 & mod.run$P<=0.72,]
  } else if (selfing=="yes"){
    mod.run=mod.run[!is.na(mod.run$P) & mod.run$P<=0.08,]  
  } 
  aic.min=aggregate(mod.run$aic, list(model=mod.run$model, filt=mod.run$filt), min, na.rm=T)
  names(aic.min)[3]="aic.min"
  mod.run=merge(mod.run, aic.min)
  mod.run=unique(mod.run[,c("model", "filt","aic", "ll", "pi.neut", "aic.min")])
  mod.run$rel.lik=exp((mod.run$aic.min-mod.run$aic)/2)
  results=ddply(.data=mod.run, .variables=.(model,filt), summarise, best.pi=pi.neut[aic==aic.min & !is.na(aic)], mean.pi=weighted.mean(pi.neut, rel.lik), best.aic=min(aic,na.rm=T), best.ll=ll[aic==aic.min & !is.na(aic)])
  results$U = U
  results$spec = spec
  results$wind = wind
  results$rec = rec
  return(results)
}

get_mod_results_plot<-function(spec, wind, U, rec, selfing, ...) {
  #same as above but also returns the sh parameter from the best model
  path="/Volumes/LaCie/Projects/Current/ne/final_data/modres/"
  file=paste(spec,wind,U,rec,"select.out", sep="_")
  print(file)
  mod.run<-read.table(paste0(path,file),header=T)
  #remove ill-conditioned fits (pi.neut = 1 is boundary of search space)
  mod.run=droplevels(subset(mod.run, model != "bgs_only(const)" & pi.neut < 1 & !is.na(pi.neut))) 
  if (is.na(selfing)) {
    mod.run=mod.run 
  } else if (selfing=="no") {
    mod.run=mod.run[!is.na(mod.run$P) & mod.run$P==1,]
  } else if (selfing=="partial") {
    mod.run=mod.run[!is.na(mod.run$P) & mod.run$P>=0.60 & mod.run$P<=0.72,]
  } else if (selfing=="yes"){
    mod.run=mod.run[!is.na(mod.run$P) & mod.run$P<=0.08,]  
  } 
  aic.min=aggregate(mod.run$aic, list(model=mod.run$model, filt=mod.run$filt), min, na.rm=T)
  names(aic.min)[3]="aic.min"
  mod.run=merge(mod.run, aic.min)
  mod.sh=mod.run[,c("model", "filt", "aic", "sh", "param2")]
  mod.run=unique(mod.run[,c("model", "filt","aic", "ll", "pi.neut", "aic.min")])
  mod.run$rel.lik=exp((mod.run$aic.min-mod.run$aic)/2)
  results=ddply(.data=mod.run, .variables=.(model,filt), summarise, best.pi=pi.neut[aic==aic.min & !is.na(aic)], mean.pi=weighted.mean(pi.neut, rel.lik), best.aic=min(aic,na.rm=T), best.ll=ll[aic==aic.min & !is.na(aic)])
  results$U = U
  results$spec = spec
  results$wind = wind
  results$rec = rec
  results = merge(results, mod.sh, by.x=c("model", "filt", "best.aic"), by.y=c("model", "filt", "aic"), all.x=T, all.y=F)
  return(results)
}



