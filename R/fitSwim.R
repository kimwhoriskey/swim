#fit the SHMMM

fitSwim <- function(x, ...) UseMethod("swim")


fitSwim <- function(data, ts, regularize=TRUE){
  
  #data input, n X 2 vector
  #ts timestep (hours)

  if(any(!class(data$date) %in% c("POSIXt", "POSIXct", "POSIXlt"))){
    stop("Date-time vector is not a POSIX class")
  }
  
  #normally if one is NA, so is the other
  if(anyNA(data$lat)) {
    data = data[which(!is.na(data$lat)),]
    warning("Locations with NA value removed")
  }
  
  
  if(regularize==TRUE){
    #interpolate the data to regular time intervals
  delta = ts*60*60 #time step in seconds
  t0 = data$date[1] #first time step
  tT = data$date[length(data$date)] #the final date time
  t = seq(t0, tT, by=delta) #all the time values for interpolation
  iLoc <- cbind(approx(data$date, data$lon, xout = t)$y,
                approx(data$date, data$lat, xout = t)$y)
  } else {
    iLoc = cbind(data$lon, data$lat)
  }
  
  
  
  #load TMB
  requireNamespace("TMB", quietly=TRUE)
  
  
  #create a list of input data. The order must match the order in the c++ file!
  space = 0:1
  mu = c(0.5,0.5)
  dats <- list(x = t(iLoc),
               # b=as.integer(space),
               stateSpace=factor(space),
               initDist=as.double(mu))
  
  
  #create a list of input parameters. Again, order matters!
  parameters <- list(logitTheta1=0, logitTheta2=0, 
                     logitGamma1=0, logitGamma2=0, 
                     logSdlon=0, logSdlat=0, 
                     logA=matrix(log(1),ncol=1,nrow=2))
  
  
  #Then make an objective function, the negative log likelihood to minimize. 
  obj <- MakeADFun(dats,parameters,
                   DLL="swim")
  
  #Now we can pass the objective function and its derivatives into any regular R optimizer. 
  #We are using nlminb. 
  cat("Optimizing the objective function \n")
  time = system.time(opt <- nlminb(obj$par,obj$fn,obj$gr))

  #calculate the parameter results
  cat("Calculating the standard errors \n")
  srep <- summary(sdreport(obj))

  #calculate the latent behavioral states with the Viterbi algorithm
  states = obj$report()$states
  
  #calculate minimized nll
  nll = obj$fn()
  
  #return the object and the parameter results
  regData = data.frame(t, iLoc)
  names(regData) = c("date", "lon", "lat")
  rslts <- list(regData = regData, obj=obj, parameters=srep, states=states+1, time=time, nll=nll)

  class(rslts) <- "swim" #set class (summary, print, and plot functions)
  
  return(rslts)
  
}
