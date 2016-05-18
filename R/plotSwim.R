
plotSwim = function(object){
  
  requireNamespace("sp", quietly=TRUE)

  data(shpfiles)

  xlims = c(min(object$regData$lon)-2, max(object$regData$lon)+2) #longitude
  ylims = c(min(object$regData$lat)-2, max(object$regData$lat)+2) #latitude

  
  locs = data.frame(object$regData$lon, object$regData$lat)
  #make the spatial points data frame to get the location points
  spdf <- SpatialPointsDataFrame(coords = locs, data=locs, 
                                 proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
  #make the spatial lines to connect the dots (locations)
  line = Lines(list(Line(locs)), ID="path")
  sl = SpatialLines(list(line))
  
  #plot everything
  plot(world, col="mediumseagreen", xlim=xlims, ylim=ylims, axes=TRUE)
  plot(sl, add=TRUE, col="grey") 
  plot(spdf[which(object$states==1),], add=TRUE, pch=20, cex=0.7, type="o", 
       col="grey")
  plot(spdf[which(object$states==2),], add=TRUE, pch=20, cex=0.7, type="o", 
       col="blue")
  
}




