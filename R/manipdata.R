
###############################     Data Functions     ###############################


##################### Split up a track based on data gaps #####################
#' split a track for a single individual based on temporal gaps
#' @export
#' @param dat the data, must have a date column in posix
#' @param the cutoff gap time in hours. Any gaps larger than this will force the track to be split. 
splitTrack <- function(dat, cutoff=24){
  timediff <- diff(as.numeric(dat$date))
  idx <- which(timediff >= cutoff*3600)
  suffixes <- rep(seq(length(idx)+1), times=diff(c(0, idx, nrow(dat))))
  dat$trackid <- paste("track", suffixes, sep='')
  return(dat)
}



#find foraging bouts
findBouts <- function(dat, numberofstatesinbout=5){
  rles <- data.frame(state=rle(dat$bhat)$values, run=rle(dat$bhat)$lengths) %>% 
    mutate(startidx=c(1,head(cumsum(run)+1, -1))) %>% 
    mutate(bout=case_when((run>=numberofstatesinbout & state==2) ~"bout", T~"notbout"))
  # idx <- which(rles$run>=3 & rles$state==2)
  return(dat %>% mutate(bout = rep(rles$bout, times=rles$run)) %>% 
           mutate(boutidx = paste(dat$id[1], dat$trackid[1], bout, rep(1:nrow(rles), times=rles$run), sep='')))
  }



##################### Regularize a track by linear interpolation #####################

#' regularize a track in time
#' @export
#' @param dat data frame, needs posix time and column names 'date', 'lon', 'lat' (even if eastings/northings)
#' @param ts time step in hours
regTrack <- function(dat, ts){
  #ts in hours
  delta = ts*60*60 #time step in seconds
  t0 = dat$date[1] #first time step
  tT = dat$date[nrow(dat)] #the final date time
  # dates = seq(t0, tT+delta, by=delta) #all the time values for interpolation
  dates = seq(t0, tT, by=delta) #all the time values for interpolation
  regx <- as.data.frame(cbind(lon=approx(dat$date, dat$lon, xout = dates)$y,
                              lat=approx(dat$date, dat$lat, xout = dates)$y))
  #best guess of final state is final observation, only needed for first step anyways
  regx <- rbind(regx, dat[nrow(dat),c("lon", "lat")])
  # return(list(regx=regx, xdates=c(dates, tail(dat,1)$date)))
  return(list(regx=regx, xdates=c(dates, tT+delta)))
}




##################### Calculate the idx's and jidx's for the SSM #####################

#' get the weighting factors for measurement error equation
#'  @export
#' @param dat a data frame with  date column in posix
#' @param ts time step in hours
getJidx <- function(dat, ts){
  
  # ts in hours
  delta = ts*60*60 #time step in seconds
  t0 = dat$date[1] #first time step
  tT = dat$date[nrow(dat)] #the final date time
  x_dates = seq(t0, tT+delta, by=delta) #all the time values for interpolation
  # x_dates = seq(t0, tT, by=delta) #all the time values for interpolation
  
  # interval of the y's in the xhats
  idx <- findInterval(dat$date, x_dates)
  
  # the proportion (how close y is to its floor xhat)
  jidx <- numeric()
  for(i in 1:nrow(dat)){
    jidx[i] <- as.numeric(difftime(dat$date[i], x_dates[idx[i]], units="hours")/ts)
  }
  
  return(list(idx=idx, jidx=jidx, x_dates=x_dates, nx=length(x_dates)))
  
}


####################### Get the argos errors for each t ##############################

#' get the Argos error terms for the measurement equation
#' @export
#' @param ac vector of argos classes for the observations
getAE <- function(ac, gpsdof=30, gpsdiv=10){
  
  all_classes <- c("GPS", "3", "2", "1", "0", "A", "B", "Z")
  sigmalon <- (c(0.289866/gpsdiv, 0.289866, 0.3119293, 0.9020423, 2.1625936, 0.507292, 4.2050261, 4.2050261)/6366.71 * 180)/pi#, 0.01)
  sigmalat <- (c(0.1220553/gpsdiv, 0.1220553, 0.2605126, 0.4603374, 1.607056, 0.5105468, 3.041276, 3.041276)/6366.71 * 180)/pi#, 0.01)
  nulon <- c(gpsdof, 3.070609, 2, 2.298819, 2, 2, 2, 2)#, 100000) #actual df values estimated
  nulat <- c(gpsdof, 2.075642, 6.314726, 3.896554, 2, 2, 2, 2)#, 100000)
  argos_error <- cbind(sigmalon, nulon, sigmalat, nulat)
  row.names(argos_error) <- all_classes
  
  ae <- as.data.frame(argos_error[as.character(ac),])
  return(as.matrix(ae))
}


