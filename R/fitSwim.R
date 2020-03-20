#fit the HMMM

fitSwim <- function(data, ts, regularize=TRUE, pars=list(working_theta = c(0, 0), 
                                                         working_gamma = c(0, 0), 
                                                         working_tau_lon=0, working_tau_lat=0, 
                                                         working_A=matrix(c(log(0.75), log(0.25)),ncol=1,nrow=2))){
  


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
  # requireNamespace("TMB", quietly=TRUE)
  
  #Then make an objective function, the negative log likelihood to minimize. 
  obj <- TMB::MakeADFun(data=list(x=t(iLoc)),
                   parameters = pars,
                   DLL="swim")
  
  #Now we can pass the objective function and its derivatives into any regular R optimizer. 
  #We are using nlminb. 
  cat("Optimizing the objective function \n")
  time = system.time(opt <- nlminb(obj$par,obj$fn,obj$gr))

  #calculate the parameter results
  cat("Calculating the standard errors \n")
  srep <- summary(TMB::sdreport(obj))

  #get the latent behavioral states with the Viterbi algorithm
  states = obj$report()$states
  
  # get the stationary distribution
  stationary <- obj$report()$delta
  
  # get the one step ahead (forecast pseudo) residuals
  pseudos <- as.data.frame(qnorm(obj$report()$pseudo[-1,]))
  names(pseudos) <- c("lon", "lat")
  
  #calculate minimized nll
  nll = obj$fn()
  
  #return the object and the parameter results
  regData = data.frame(t, iLoc)
  names(regData) = c("date", "lon", "lat")
  rslts <- list(regData = regData, 
                obj=obj, opt=opt, 
                parameters=srep, 
                states=states, stationary=stationary, 
                pseudos=pseudos,
                time=time, nll=nll)

  class(rslts) <- "swim" #set class (summary, print, and plot functions)
  
  return(rslts)
  
}
