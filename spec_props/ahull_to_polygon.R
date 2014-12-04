ah2sp<-function(x, rnd=6, proj4string=CRS(as.character(NA))) {
  require(rgeos)
  require(alphahull)
  
  if(class(x) != 'ahull') {
    stop('this function only works with `ahull` class objects')
  }
  
  # convert ashape edges to DF
  x.ah.df <- as.data.frame(x$arcs)
  
  # convert each arc to a line segment
  l.list <- list()
  
  for(i in 1:nrow(x.ah.df)) {
    
    # extract row i
    row_i <- x.ah.df[i,]
    
    # extract elements for arc()
    v <- c(row_i$v.x, row_i$v.y)
    theta <- row_i$theta
    r <- row_i$r
    cc <- c(row_i$c1, row_i$c2)
    
    # from arc()
    angles <- anglesArc(v, theta)
    seqang <- seq(angles[1], angles[2], length = 100)
    x <- round(cc[1] + r * cos(seqang),rnd)
    y <- round(cc[2] + r * sin(seqang),rnd)
    
    # convert to line segment
    l.list[[i]] <- Line(cbind(x,y))
  }
  
  # promote to Lines class, then to SpatialLines class
  l <- Lines(l.list, ID="ahull")
  
  # copy over CRS data from original point data
  l.spl <- SpatialLines(list(l), proj4string=proj4string)
  
  # merge lines that intersect
  l.merged<-gLineMerge(l.spl)
  
  #now convert to polygon
  l.segs<-slot(l.merged, "lines")[[1]]@Lines
  sppolys <- SpatialPolygons(list(Polygons(lapply(l.segs, function(x) { Polygon(slot(x, "coords")) }), ID = "1")), proj4string=proj4string)
  
  return(sppolys)
}