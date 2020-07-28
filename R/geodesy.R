#' @import dplyr
#' @import sf
#'
#' @useDynLib swim
dummy <- function() return(NULL)
############################################################
####################### step lengths #######################
############################################################


getStepLengths <- function(dat, units="m", convert_to_rad=TRUE, euclidean=FALSE) {
   
   # units can take one of "m" or "km"
   if(units=="m"){
      R=6371e3
   } else if(units=="km") {
      R=6371
   } else {
      stop("units not supported")
   }
   
   # have to convert lon/lat into radians first
   if(convert_to_rad) dat <- dat*pi/180
   
   # initialize
   n = nrow(dat)
   sl <- numeric(n)
   # sl[1] <- NA
   
   for(i in 2:n){
      a = sin((dat[i,2]-dat[i-1,2])/2)*sin((dat[i,2]-dat[i-1,2])/2) + 
          cos(dat[i-1,2])*cos(dat[i,2])*sin((dat[i,1]-dat[i-1,1])/2)*sin((dat[i,1]-dat[i-1,1])/2)
      b = 2*atan2(sqrt(a), sqrt(1-a))
      sl[i] = R*b
   }
   
   return(sl)
}






#############################################################
################ bearings and turning angles ################
#############################################################


getBearing <- function(loc1, loc2, euclidean=FALSE){
   
   # takes coordinates in radians, returns bearing in radians
   
   #lon,lat pairs
   if(euclidean){
      num <- loc2[2]-loc1[2]
      den <- loc2[1]-loc1[1]
   } else {
      num <- sin(loc2[1]-loc1[1])*cos(loc2[2])
      den <- cos(loc1[2])*sin(loc2[2]) - sin(loc1[2])*cos(loc2[2])*cos(loc2[1]-loc1[1])
      
   }
   
   bear <- atan2(num, den)
   
   return(bear)
   
}

getAngles <- function(dat, type="spherical", convert_to_rad=TRUE){
   
   # type can be one of "spherical", "spherical.adj", and "euclidean"
   
   # have to convert lon/lat into radians first
   if(convert_to_rad) dat <- dat*pi/180
   
   # dat <- obs$X
   n <- nrow(dat)
   ta <- numeric(n)
   # ta[1] <- ta[n] <- NA

   for(i in 2:(n-1)){ 
   if(type=="spherical.adj"){
         first.bearing <- getBearing(loc1=dat[i-1,], loc2=dat[i,], euclidean=FALSE)*180/pi
         second.bearing <- getBearing(loc1=dat[i,], loc2=dat[i+1,], euclidean=FALSE)*180/pi
         #this one is just the intial headings of both vectors
      } else {
      if(type=="spherical"){
         first.bearing <- 180 + getBearing(loc1=dat[i,], loc2=dat[i-1,], euclidean=FALSE)*180/pi
         second.bearing <- getBearing(loc1=dat[i,], loc2=dat[i+1,], euclidean=FALSE)*180/pi
         # this one is the initial heading of the second (in time) vector, and the final heading
         # of the first (in time) vector
      } else {
         if(type=="euclidean"){
            first.bearing <- getBearing(loc1=dat[i,], loc2=dat[i+1,], euclidean=TRUE)*180/pi
            second.bearing <- getBearing(loc1=dat[i-1,], loc2=dat[i,], euclidean=TRUE)*180/pi
            #reverse the order because the angles are counterclockwise for the euclidean bearings
         } else {
            stop("angle type not implemented")
         }
      }
      }
      turning.angle <- (second.bearing-first.bearing) 
      ta[i] <- turning.angle - 360*floor((turning.angle+180)/360)
         # ensure right turns are positive, euclidean angles are calculated counterclockwise from the x-axis, while
         # the angles on a sphere are calculated clockwise from North, so  we have to reverse the subtraction
         # relative to each other to be able to always call negative angles left turns
      }

   return(ta)  
   
}










#########################################################
##################### recreate path #####################
#########################################################



getNextLoc <-  function(start_loc, start_bearing, step_length){
   
   R=6371 # radius of the Earth
   
   # get the next latitude
   next.lat = asin( sin(start_loc[2])*cos(step_length/R)  + 
                       cos(start_loc[2])*sin(step_length/R)*cos(start_bearing) )
   
   # get the next longitude
   next.lon = start_loc[1] + atan2( sin(start_bearing)*sin(step_length/R)*cos(start_loc[2]),
                                    cos(step_length/R)-sin(start_loc[2])*sin(next.lat))
   return(c(next.lon, next.lat))
   
}



getPath  <- function(turning_angles, step_lengths, 
                      first_loc=NULL, first_bearing=NULL){
   
   # takes turning angles in degrees, but first.bearing and first.loc are in radians
   # spits out locations  in radians
   # turning angles and step lengths should be as long as the path that you want
   # this means that the first two turning angles and first step length won't be used

   # if first location and bearing aren't specified, set to c(0,0) and 0 radians 
   if(is.null(first_loc)) first_loc <- c(0,0)
   if(is.null(first_bearing)) first_bearing <- 0
   n=length(turning_angles)
   
   locations <- matrix(nrow=n, ncol=2)
   
   locations[1,] <- first_loc
   #call get next location
   locations[2,] <- getNextLoc(start_loc=first_loc, start_bearing=first_bearing,
                                 step_length=step_lengths[2])
   
   for(i in 3:n){
      # get the next bearing by adding the turning angle to the final bearing from x_{t-2} to x_{t-1} 
      next.bearing <- ((180 + getBearing(loc1=locations[i-1,], loc2=locations[i-2,], euclidean=FALSE)*180/pi) + turning_angles[i])*pi/180
      # get the next location
      locations[i,] <- getNextLoc(start_loc=locations[i-1,], start_bearing=next.bearing,
                                    step_length=step_lengths[i])
   }
      
   return(locations)
}


