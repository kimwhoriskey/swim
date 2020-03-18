#plot the swim fits

plot.swim <- function(x,...){
  
  data = x$regData
  states = x$states
  
  plot(data$lat~data$lon, type="o", pch=20, cex=0.7, 
       col=ifelse(states==2, "blue", "grey"), 
       xlab="Longitude", ylab="Latitude")
#   plot(data[states==1,]$lat~data[states==1,]$lon, type="o", pch=20, cex=0.7, col="grey", 
#        xlab="Longitude", ylab="Latitude")
#   points(data[states==2,]$lat~data[states==2,]$lon, type="p", pch=20, cex=0.7, col="blue")
  
}
