#plot the swim fits

plot.swim <- function(object,...){
  
  ##add shapefiles? Or define a whole new function to plot with shapefiles? 
  
  data = object$regData
  states = object$states
  
#   plot(data$lat~data$lon, type="o", pch=20, cex=0.7, 
#        col=ifelse(states==2, "blue", "grey"))
  
  
  plot(data[states==1,]$lat~data[states==1,]$lon, type="o", pch=20, cex=0.7, col="grey", 
       xlab="Longitude", ylab="Latitude")
  points(data[states==2,]$lat~data[states==2,]$lon, type="p", pch=20, cex=0.7, col="blue")
  
}
