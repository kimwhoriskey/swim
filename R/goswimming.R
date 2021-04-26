

############################### All of the fitting functions I will need ################################


#############################################################
#################### windsorize the data ####################
#############################################################

#' Winsorize some observations
#' @param obs dataframe of observations, 'id', date', 'lon', 'lat', 'lc'
#' @param inflate an inflation factor for the winsorizing cut off value
#' @export
winsorize <- function(obs, inflate=1){
  
  #get indices where the tracks start
  trackstartidx <- c(1, head(cumsum(rle(obs$trackid)$lengths), -1) +1)
  
  # generate new observation data frame
  newobs <- obs
  
  # step lengths first, units m
  if("sf" %in% class(obs)){
    empty <- st_as_sfc("POINT(EMPTY)", crs=st_crs(obs))
    d <- c(0, head(st_distance(obs$geometry, lead(obs$geometry, default=empty), by_element=TRUE), -1))
  } else {
    d <- c(0, argosfilter::distanceTrack(obs$lat, obs$lon))*1000
    # d <- getStepLengths(as.matrix(obs[,c("lon", "lat")]), convert_to_rad=TRUE, units="m")
  }
  # now divide by the time to get speed, units m/s
  secs <- c(1, as.numeric(diff(obs$date)))
  if(any(secs==0)) warning("there are multiple observations for some datetime stamps")
  sp <- d/secs 

  # get the quantile stats
  # remove the values corresponding to the temporal gaps
  # do this for each lc class
  if(length(trackstartidx)>1){
    speedstats <- tapply(sp[-trackstartidx[-1]], function(x)boxplot.stats(x)$stats, INDEX=obs[-trackstartidx[-1],]$lc)
  } else {
    speedstats <- tapply(sp, function(x)boxplot.stats(x)$stats, INDEX=obs$lc, simplify=FALSE)
  }
 
  # set an upperlimit to the step length value that we want to scale by
  # default is 1.5*IQR
  upperlim <- lapply(speedstats, function(x){x[4] + inflate*1.5*(x[4]-x[2])})
  
  outidx <- which(obs$outliers==TRUE)
  
  # remove any values that correspond to data gaps
  if(any(outidx %in% trackstartidx)){
    outidx <- outidx[-which(outidx %in% trackstartidx)]
  }
  
  # remove any outliers where the observed speed is less than 1.5xIQR
  # because that would pull the observed locations AWAY from the true 
  outidx <- outidx[as.numeric(upperlim[obs[outidx,]$lc]) < sp[outidx]]

  # throw a warning if the starting value of a track segment might be an outlier
  if(any(outidx %in% (trackstartidx+1))) print(paste(obs$id[1], obs$trackid[outidx[outidx %in% (trackstartidx+1)]], "might start with an outlying location"))
  if(any(which(sp>max(as.numeric(upperlim))) %in% (trackstartidx+1))) print(paste(obs$id[1], obs$trackid[which(sp>max(as.numeric(upperlim)))[which(sp>max(as.numeric(upperlim))) %in% (trackstartidx+1)]], "might start with an outlying location"))
  
  
  
  # calculate the speed ratio based on the upperlim
  # basically this is a proportion of how much we should scale an observation back, based on speed
  dr <- as.numeric(upperlim[obs[outidx,]$lc])/sp[outidx] 
  
  # replace outliers
  # basically its a projection of the vector onto itself, just scaling it back a bit
  # so this will technically affect the turning angles of the next two locations
  if(length(outidx) > 0 ){
      
    if("sf" %in% class(obs)){
        coos <- st_coordinates(newobs)
        newcoos <- coos
      for(i in 1:length(outidx)){
        newcoos[outidx[i], c("X", "Y")] <- coos[outidx[i]-1, c("X", "Y")] + 
          dr[i]*(coos[outidx[i], c("X", "Y")] - coos[outidx[i]-1, c("X", "Y")])
      }
      st_geometry(newobs) <- st_geometry(st_as_sf(data.frame(newcoos), 
                                                  coords=c("X", "Y"), 
                                                  crs=st_crs(obs)))
    } else {
      if(length(outidx) > 0 ){
        for(i in 1:length(outidx)){
          newobs[outidx[i], c("lon", "lat")] <- obs[outidx[i]-1, c("lon", "lat")] + 
            dr[i]*(obs[outidx[i], c("lon", "lat")] - obs[outidx[i]-1, c("lon", "lat")])
        }
      }
    }

    # print info on how many were scaled
    cat(" \n we winsorized", length(outidx), "outliers from seal", as.character(newobs$id[1]), "\n")
  } else {
    # or tell us that there aren't any outlying locations
    cat(" \n seal", as.character(newobs$id[1]), "does not have any outliers as you have defined them \n")
  }

  # add a column for the locations that were winsorized
  newobs$winsorized <- FALSE
  newobs[outidx, "winsorized"] <- TRUE
  # add the old locations in
  if("sf" %in% class(obs)){
    newobs$oldx <- st_coordinates(obs)[,1]
    newobs$oldy <- st_coordinates(obs)[,2]
  } else {
    newobs$oldlon <- obs$lon
    newobs$oldlat <- obs$lat
  }
  
  # return object
  return(newobs)
  
}




###############################################################
##################### First Step: Fit HMM #####################
###############################################################

#' fit an HMM 
#' @export
#' @param dat Named list: model denotes the model to use in integer format (check c++), x are the data, two rows, first one lon, second one lat
#' @param p_start Starting values for the parameters
#' @param silence if true, the tmb output isn't reported
fit_hmm <- function(dat,
                    p_start = list(working_theta = c(log(1), log(1)),
                                   working_gamma = c(log(1), log(1)),
                                   working_tau_lon=log(1), working_tau_lat=log(1),
                                   working_A=matrix(c(log(0.75/0.25), log(0.25/0.75)),ncol=1,nrow=2)),
                    silence=FALSE){
  
  # get track ids
  tracknames <- tail(names(dat), -2)
  
  # create TMB object
  obj <- MakeADFun(dat, p_start, DLL='swim', silent=silence)
  
  # optimize TMB object
  opt <- nlminb(obj$par, obj$fn, obj$gr)
  
  # likelihood value
  nll <- obj$fn(opt$par)
  
  # extract behavioural states
  # states <- obj$report()$states
  states <- list()
  for(i in 1:length(tracknames)) states[[tracknames[i]]] <- obj$report()[[tracknames[i]]]$states

  
  # extract parameter estimates
  params_hmm <- summary(sdreport(obj))
  
  # get the stationary distribution
  stationary <- obj$report()$delta
  
  # pseudoresiduals
  pseudos <- list()
  for(i in 1:length(tracknames)) pseudos[[tracknames[i]]] <- obj$report()[[tracknames[i]]]$pseudo
  
  # return results
  rslt <- list(mess = opt$mess, obj=obj, opt=opt, nll=nll, params = params_hmm, 
               bhat=states, pseudos=pseudos, 
               statdist = stationary)
  class(rslt) <- "swimhmm"
  return(rslt)
  
}











##############################################################################################
##################### Second Step: Fit SSM with fixed behavioural states #####################
##############################################################################################



#dat be a list of y, b_hat, idx, jidx, ae
#' fit an ssm
#' @export
#' @param dat list of data, model is the model integer (check cpp), y is the observations, b are the behavioural states, idx is the index linking the observations to the true locations, jidx is the proportion of time between x and y for the interpolation calculation
#' @param p_start the parameter list 
#' @param res the name of the ranndom effects, probs not needed anymore
#' @param mapping the map set of fixed parameters
#' @param optimizer the name of the outer optimizer
#' @param optimizer_arguments arguments for the outer optimizer
#' @param include_TMB_gr whether to include the tmb gradient in the optimization
#' @param inner_optim_method inner optimizer
#' @param inner_control inner optimizer options
#' @param allowErr whether to allow errors in the optimization
#' @param jiggle_err how much randomm noise to add to the parameter starting values i ncase of error
#' @param silence if true don't print the tmb output
fit_ssm <- function(dat,
                    p_start = list(working_theta = c(log(1), log(1)),
                                   working_gamma = c(log(1), log(1)),
                                   working_tau_lon=0, working_tau_lat=0,
                                   working_psi=log(1),
                                   x = matrix(log(1), nrow=2, ncol=n),
                                   working_A = matrix(c(log(0.75/0.25), log(0.25/0.75)))),
                    res = "x",
                    mapping=list(working_A=c(factor(NA), factor(NA))), 
                    optimizer="nlminb", optimizer_arguments = list(rel.tol = 1e-10, x.tol=1.5e-8), 
                    include_TMB_gr=TRUE, 
                    inner_optim_method = "newton", inner_control = list(maxit=1000), 
                    allowErr=FALSE, 
                    jiggle_err=0.02,  
                    silence=FALSE, 
                    maxtime=5*60){
  
  #extract the names of the fixed variables
  fixed_names <- names(mapping)
  # get the names of the tracks
  tracknames <- tail(names(dat), -2)
  # write the random effects
  res <- paste(tracknames, '.', res, sep='')
  
  
  # set up an error vector so that we can keep track of how many errors occur
  err_num = 0 # set initial error value at 0, base an iterative loop on this
  err = 0
  class(err) <- "error"
  
  while("error"%in%class(err) & allowErr==FALSE){
    
    if(err_num>10) stop("encountered errors > 10 times while attempting to fit the SSM") #cut optimization if errors > 10
    
    if(err_num>0){
      # print a message letting us know we're working on errors
      cat("\n -------------------- \n working on errors \n err: ", err$message, " \n -------------------- \n")
      
      # jiggle starting parameters - add a little bit of random error to the ones that aren't fixed
      # exclude the random effects as well
      idx <- which(!(names(p_start) %in% c("x", fixed_names)))
      for(i in idx) p_start[[i]]<- p_start[[i]] + rnorm(length(p_start[[i]]), 0, jiggle_err)
    }
    
    # re-create TMB object, remove previous just in case
    if(exists("obj")) rm(obj)

    obj <- MakeADFun(dat, p_start, map=mapping, random=res, DLL='swim',
                     inner.method = inner_optim_method,
                     inner.control = inner_control, silent=silence)
    
    # use trycatch loop again
    # optimizer="optim"
    # optimizer_arguments=list(control=list(reltol=1e-4))
     if(include_TMB_gr){
      argus <- append(list(obj$par, obj$fn, obj$gr), optimizer_arguments)
    } else {
      argus <- append(list(obj$par, obj$fn), optimizer_arguments)
    }
 
    err <- tryCatch({
      opt <- R.utils::withTimeout({do.call(optimizer, argus)}, timeout=maxtime, onTimeout="error") # optimize TMB object
      srep <- summary(sdreport(obj)) # extract parameter estimates
    },  error=function(e){e} #capture errors
    )
    
    err_num = err_num+1 #count the errors
    
  } # close while loop
  
    # extract the likelihood value
  nll <- obj$fn(opt$par)
  
  # parameter estimates, take the x's out of it
  params_ssm <- srep[! rownames(srep) %in% res,] 
  
  # extract location state estimates
  if(is.null(res)){
    xhat=NULL
  } else {
    xhat <- list()
    for(i in 1:length(tracknames)) xhat[[tracknames[i]]] <- obj$report()[[tracknames[i]]]$x
  }
    
  if(optimizer=='optim' & opt$convergence==0){
    mess <- "converged"
  } else {
    mess <- opt$mess
  }
  
  rslt <- list(mess = mess, obj=obj, opt=opt, nll=nll, 
               params=params_ssm, xhat=xhat)
  class(rslt) <- "swimssm"
  return(rslt)
  
}



######################################################################################
##################### Getting appropriate parameter values  #####################
######################################################################################

#' get the appropriate parameter values for each step, used in the issm model
#' @export
#' @param prev_hmm the previous hmm
#' @param prev_ssm the previous ssm 
#' @param curr_model what's next? an hmm or ssm
get_parms <- function(prev_hmm, prev_ssm, curr_model, ssm_type, init_psi){
  
  p <- list(working_theta = as.numeric(prev_hmm$params[row.names(prev_hmm$params) %in% 'working_theta', 'Estimate']), #have to check that these are the right thetas
            working_gamma = as.numeric(prev_hmm$params[row.names(prev_hmm$params) %in% 'working_gamma','Estimate']),
            working_tau_lon = prev_hmm$params['working_tau_lon', 'Estimate'], 
            working_tau_lat = prev_hmm$params['working_tau_lat', 'Estimate'])

  if(curr_model == "ssm"){
    
    # overwrite the previous process error 
    p$working_tau_lon <- p$working_tau_lat <- NULL
    p$working_sigma_lon <- prev_hmm$params['working_tau_lon', 'Estimate']
    p$working_sigma_lat <- prev_hmm$params['working_tau_lat', 'Estimate']
    
    if(ssm_type==1){
      p$working_psi <- log(init_psi)
    } else if (ssm_type == 2){
      p$working_gps_err <- c(log(1), log(1))
    } else if(ssm_type==3){
      p$working_psi <- log(init_psi)
      p$working_gps_err <- c(log(1), log(1))
    }
    
    tracknames <- names(prev_hmm$bhat)
    for(i in 1:length(tracknames)){
      p[[paste(tracknames[i], ".x", sep="")]] = matrix(log(1), nrow=2, ncol=length(prev_hmm$bhat[[tracknames[i]]]))
    }
    
    # if(!is.null(prev_ssm)){
    #   
    #   # measurement error
    #   p$working_psi <- as.numeric(prev_ssm$params[row.names(prev_ssm$params) %in% 'working_psi', 'Estimate'])
    #   # add in the locations
    #   p$x = matrix(log(1), nrow=2, ncol=length(prev_hmm$b_hat)) 
    #   # start at 0 regardless, but we can change this if we want
    # 
    # }
    
  } 
    
  p$working_A = matrix(c(prev_hmm$params[row.names(prev_hmm$params) %in% 'working_A', 'Estimate'][1:2]),ncol=1,nrow=2)

  return(p)
}

######################################################################################
##################### Iteratively fit multiple SSM and HMM steps #####################
######################################################################################

#' fit a switching ssm iteratively 
#' @export
#' @param obs the data: id, date (posix), lon, lat, ac
#' @param ts the time step
#' @param p_start_hmm starting values for the hmm step
#' @param init_psi starting values for the measurement error parameter
#' @param res random effect values
#' @param maxsteps max number of steps
#' @param fixsteps whether to fix the number of steps, should really be true always
#' @param allowfc whether to allow false convergence in the ssm step
#' @param jiggle_fc how much random noise to add to the starting values of the parameters if the ssm falsely covnerges
#' @param ssm_map mapping list for the ssm step
#' @param allsilent if true, don't print the tmb tape for the hmm or ssm
#' @param scaleobs value to scale the observations by (literally just divided), helpful for convergence if they are in eastings/northings (set to 1000 to put track in km)
#' @param ... further arguments to fit_ssm
fit_issm <- function(obs, 
                     ts, 
                     p_start_hmm,
                     update_pstarthmm = FALSE,
                     init_psi,
                     res="x",
                     maxsteps=10, fixsteps=FALSE, 
                     allowfc=FALSE, jiggle_fc=0.01, max_fc=3,
                     allsilent=TRUE,
                     scaleobs=1,
                     gpsdof=10, 
                     gpsdiv=10, 
                     ssm_args=NULL,
                     ssm_map=NULL
                     ){


  # obs=sim$obs
  # move="car"
  # ts=3
  # p_start_hmm = list(working_ta_pars = matrix(log(1), nrow=2, ncol=2),
  #                    working_sl_pars = matrix(log(1), nrow=2, ncol=2),
  #                    working_ac_pars = rep(log(1), 2),
  #                    working_A = matrix(c(log(0.75/0.25), log(0.25/0.75)), nrow=2, ncol=1))
  # ssm_map = list(working_ta_pars = factor(matrix(NA, nrow=2, ncol=2)),
  #                working_sl_pars = factor(matrix(NA, nrow=2, ncol=2)),
  #                working_ac_pars = rep(factor(NA),2))
  # res=c("working_ta", "working_sl")
  # init_psi=c(2,2)
  # maxsteps=3
  # fixsteps=TRUE
  # transform_random_effects=0
  # allowfc=TRUE
  # jiggle_fc=0.01
  # split=TRUE
  # allsilent=TRUE
  
  # obs = sim$obs
  # move = "dcrw"
  # ts = 3
  # p_start_hmm = list(working_theta = c(log(1), log(1)),
  #                    working_gamma = c(log(1), log(1)),
  #                    working_tau_lon=log(1), working_tau_lat=log(1),
  #                    working_A=matrix(c(log(0.75/0.25), log(0.25/0.75)),ncol=1,nrow=2))
  # init_psi = 1
  # res = c("x")
  # maxsteps = 5
  # fixsteps = TRUE
  # allowfc=FALSE
  # jiggle_fc=0.01
  # ssm_map=NULL
  # transform_random_effects=0
  # split=FALSE
  # allsilent=TRUE
  # ssm_map <- list(working_gamma = factor(c(NA, NA)),
  #                 working_theta = factor(c(NA, NA)),
  #                 working_sigma_lon=factor(NA), working_sigma_lat=factor(NA))
  # 
  
  fullargs <- match.call()
  
  #convert from sf
  if(inherits(obs, 'sf')){
    isprojected=TRUE
    obs <- obs %>% bind_cols(as_tibble(st_coordinates(obs))) %>% rename(lon=X, lat=Y) %>% as.data.frame() 
    firstloc <- obs[1, c("lon", "lat")]
  } else {
    isprojected=FALSE
    firstloc <- data.frame(lon=0, lat=0)
  }
  obs[,c("lon", "lat")] <- (obs[,c("lon", "lat")] - firstloc[rep(1, nrow(obs)),])/scaleobs
  # obs[,c("lon", "lat")] <- obs[,c("lon", "lat")]/scaleobs

  # set up accumulators
  # all of the results from each step
  hmm_results <- list()
  ssm_results <- list()
  # times it takes to fit each step
  hmm_time <- list() 
  ssm_time <- list()
  # convergence of the smm, cue a while loop on this if we don't allow false convergence
  ssm_conv <- list()
  fc_count <- numeric() # store the false convergence counts (how many at each step)
  
  # store the nll of the ssm, cue another while loop on this
  nll_ssm <- list()
  
  
  #####################################################################
  ############ initial step, get intial location estimates ############
  #####################################################################
  
  obs$trackid <- factor(obs$trackid, levels=unique(obs$trackid))
  splitobs <- split(obs, obs$trackid, )
  tracknames <- names(splitobs)
  
  kind <- split(obs$kindnum, obs$trackid)
  if(length(unique(obs$kind))==1 & "Argos" %in% unique(obs$kind)){
    ssm_type <- 1
  } else if (length(unique(obs$kind))==1 & "GPS" %in% unique(obs$kind)){
    ssm_type <- 2
  } else {
    ssm_type <- 3 # combo
  }
  
  
  idxs <- lapply(splitobs, getJidx, ts = ts)
  ae <- lapply(splitobs, function(x){as.matrix(getAE(x$lc, gpsdof=gpsdof, gpsdiv=gpsdiv))})

  regobsdat <- lapply(splitobs, function(x)regTrack(x, ts=ts))
  regobs <- lapply(regobsdat, function(x)t(x$regx))
  
  hmm_dat <- list(model=0, 
                  tracknames=as.list(tracknames))
  for(i in 1:length(tracknames)){
    hmm_dat[[tracknames[i]]]$x <- regobs[[tracknames[i]]]
  }
 


  
  ################################################################
  ############ initial step, fit HMM to observed data ############
  ################################################################
  
  
  # fit the hmm
  hmm_time[[1]] <- system.time(hmm_results[[1]] <- fit_hmm(dat = hmm_dat,
                                                           p_start = p_start_hmm,
                                                           silence=allsilent))
  
  # print out an update
  cat("\n ------------------------------------------------------------ \n finished HMM 1, convergence was:",
      hmm_results[[1]]$mess, 
      "\n ------------------------------------------------------------ \n ")
  
  ##############################################################################
  ########### second step, fit SSM to observed data with fixed states ###########
  ##############################################################################
  
  # put together the data list for the model
  ssm_dat <- list(model=ssm_type,
                  tracknames=as.list(tracknames))
  for(i in 1:length(tracknames)){
    ssm_dat[[tracknames[i]]] <- list(y = t(array(c(splitobs[[tracknames[i]]]$lon, 
                                                   splitobs[[tracknames[i]]]$lat), 
                                                 dim=c(nrow(splitobs[[tracknames[i]]]), 2))),
                                     b = hmm_results[[1]]$bhat[[tracknames[i]]]-1, 
                                     kind = kind[[tracknames[i]]],
                                     idx = idxs[[tracknames[i]]]$idx, 
                                     jidx = idxs[[tracknames[i]]]$jidx, 
                                     ae = ae[[tracknames[i]]])
  }
  # ssm_dat$ae = ae # options for either gaussian or t error? 

         
  # now get the parameters             
  p_start_ssm <- get_parms(prev_hmm = hmm_results[[1]], 
                           prev_ssm = NULL,
                           curr_model = "ssm", 
                           ssm_type = ssm_type, 
                           init_psi = init_psi)
 
  
  # for(i in 1:length(tracknames)){
  #   p_start_ssm[[paste(tracknames[i], ".x", sep="")]] = matrix(log(1), nrow=2, ncol=ncol(regobs[[tracknames[i]]]))
  # }
  
  # fit the ssm
  # block for false convergence
  # the idea is to jiggle the starting parameters by adding a teeny bit of random error
  # exclude the location states as well as anything that is fixed
  # if we fix all of the ssm parameters this may basically be moot
  # also count the numbers of false convergence so that we know, also so that we can cut off the optimizer
  ssm_fixed_names <- names(ssm_map) # record the names of the fixed parameters
  pidx <- which(!names(p_start_ssm) %in%  c("x", ssm_fixed_names)) #create an index to loop over, p for parameter
  # use a while loop acting on whether we allow fc and what the count of the false convergence is
  fc_count[1]=0
  ssm_conv[[1]]=0 #set it to 1 to initialize the loop
  runme = TRUE
  while(runme){
    # print a message letting us know we're working on false convergence
    
    if((sum(ssm_conv[[1]])>0) & (allowfc==FALSE) & (fc_count[1]>0)){
      cat("\n ---------------------------- \n working on false convergence  \n ---------------------------- \n")
      
      # jiggle the starting parameters, jiggle_fc is the st dev for random normal error
      for(i in pidx) p_start_ssm[[i]] <- p_start_ssm[[i]] + rnorm(length(p_start_ssm[[i]]), 0, jiggle_fc)
    }
 
    
    # fit SSM
    ssmargs <- ssm_args
    ssmargs$dat = ssm_dat
    ssmargs$p_start = p_start_ssm
    ssmargs$res=res
    ssmargs$silence=allsilent
    ssmargs$mapping = ssm_map
    ssm_time[[1]] <- system.time(
      ssm_results[[1]] <- do.call(fit_ssm, 
                                  ssmargs)
      )
    # ssm_time[[1]] <- system.time(ssm_results[[1]] <- fit_ssm(dat = ssm_dat, 
    #                                                          p_start = p_start_ssm,
    #                                                          res=res,
    #                                                          # mapping = ssm_map,
    #                                                          # inner_control = list(maxit=20000, step.tol=1e-4, grad.tol=1e-2),
    #                                                          # optimizer="optim", optimizer_arguments=list(control=list(reltol=1e-4)),
    #                                                          silence=allsilent,
    #                                                          # maxtime=1))#,
    #                                                          # ...))
      # save the nll 
      nll_ssm[[1]] = ssm_results[[1]]$nll
      # save the new convergence, cue on this
      ssm_conv[[1]] <- ssm_results[[1]]$opt$convergence
  
    if(ssm_conv[[1]]>0) fc_count[1] = fc_count[1] + 1 # add to the count every time we come back up to the beginning of this loop
    if(fc_count[1]>=max_fc) warning(paste("had false convergence", max_fc, "times while trying to get the SSM to fit"))
    # stop if greater than 10 times
    if(ssm_conv[[1]]==0 | allowfc==TRUE | fc_count[1]>=max_fc) runme=FALSE
  } # close while loop

  
  
  # update on where we are
  cat("\n ------------------------------------------------------------ \n finished SSM 1, convergence was:",
      ssm_results[[1]]$mess, 
      "\n ------------------------------------------------------------ \n ")
  

  
  
  ###########################################
  ############ rest of the steps ############
  ###########################################
  
  # initialize the indexer for the rest of the iteration
  i=2
  
  # use a while loop, and if fixsteps = false, then stop the optimization when nll_diff becomes positive
  # or i guess geq 0
  nll_diff = -1
  while(nll_diff<0){
    
    
    ### add extra data input
    for(j in 1:length(tracknames)) hmm_dat[[tracknames[j]]]$x = ssm_results[[i-1]]$xhat[[tracknames[j]]]
    # hmm_dat <- list(model=0, x=t(ssm_results[[i-1]]$xhat))
  
    # fit HMM again
    # update starting values of p_start_hmm if we want to
    if(update_pstarthmm){
      p_start_hmm <- get_parms(prev_hmm = hmm_results[[i-1]], 
                               prev_ssm = ssm_results[[i-1]],
                               curr_model = "hmm")
    }

    hmm_time[[i]] <- system.time(hmm_results[[i]] <- fit_hmm(dat=hmm_dat,
                                                             p_start = p_start_hmm,
                                                             silence=allsilent))

    # update on where we are (iteratively) in the optimization
    cat("\n ------------------------------------------------------------ \n finished HMM", 
        i,  
        ", convergence was:",
        hmm_results[[i]]$mess, 
        "\n ------------------------------------------------------------ \n ")
    
    
    # update the data for the ssm, the only thing that changes is the behavioural states
    # taken from the most recent ssm
    # ssm_dat$b = hmm_results[[i]]$b_hat[[tracknames[i]]]-1
    for(j in 1:length(tracknames)) ssm_dat[[tracknames[j]]]$b = hmm_results[[i]]$bhat[[tracknames[j]]]-1

    # get starting parameters for the ssm
    # starting values come from the most recent HMM, except for psi which comes from the last ssm
    # p_start_ssm <- get_parms(prev_hmm = hmm_results[[i]], 
    #                          prev_ssm = ssm_results[[i-1]],
    #                          curr_model = "ssm")
    # p_start_ssm$working_psi <- log(init_psi)
    # p_start_ssm$working_gps_err <- c(log(1), log(1))
    # for(j in 1:length(tracknames)) p_start_ssm[[paste(tracknames[j], ".x", sep="")]] = matrix(log(1), nrow=2, ncol=ncol(regobs[[tracknames[j]]]))
    # 
    
    p_start_ssm <- get_parms(prev_hmm = hmm_results[[i]], 
                             prev_ssm = NULL,
                             curr_model = "ssm", 
                             ssm_type = ssm_type, 
                             init_psi = init_psi)
    
    # block for false convergence, same deal as above
    fc_count[i]=0
    ssm_conv[[i]]=0 #set it to 1 to initialize the loop
    runme = TRUE
    while(runme){
      
      if(sum(ssm_conv[[i]])>0 & allowfc==FALSE & fc_count[i]>0){
        cat("\n --------------------------- \n working on false convergence \n --------------------------- \n")

        # jiggle parameters
        for(j in pidx){ 
          p_start_ssm[[j]] <- p_start_ssm[[j]] + rnorm(length(p_start_ssm[[j]]), 0, jiggle_fc)
        }
      }
      
      ssmargs$dat <- ssm_dat
      ssmargs$p_start <- p_start_ssm
      
      ssm_time[[i]] <- system.time(
        ssm_results[[i]] <- do.call(fit_ssm, 
                                    ssmargs)
      )
      # ssm_time[[i]] <- system.time(ssm_results[[i]] <- fit_ssm(ssm_dat, 
      #                                                          p_start = p_start_ssm,
      #                                                          res=res,
      #                                                          mapping = ssm_map,
      #                                                          # inner_control = list(maxit=20000, step.tol=1e-4, grad.tol=1e-2),
      #                                                          # optimizer="optim", optimizer_arguments=list(control=list(reltol=1e-4)),
      #                                                          silence=allsilent,
      #                                                          ...))
      
      # save the nll 
      nll_ssm[[i]] = ssm_results[[i]]$nll
      # save the new convergence, cue on this
      ssm_conv[[i]] <- ssm_results[[i]]$opt$convergence
      # calculate the difference between this and the last nll
      # should actually maybe make this a convergence criterion? 
      nll_diff <- nll_ssm[[i]]-nll_ssm[[i-1]]
    
      
      if(ssm_conv[[i]]>0) fc_count[i] = fc_count[i] + 1
      if(fc_count[i]>=max_fc) warning(paste("had false convergence", max_fc, "times while trying to get the SSM to fit"))
      
      # condition to close
      
      if(ssm_conv[[i]]==0 | allowfc==TRUE | fc_count[i]>=max_fc) runme=FALSE
      
    } # close fc while loop

    # if we've fixed the number of iterations then set this difference back to a negative number
    if(fixsteps) nll_diff <- -1
    # if we've reached the max number of steps allowed then also set this back to a negative value
    if(length(ssm_results)==maxsteps) nll_diff <- 1
    
    # update on where we are
    cat("\n ------------------------------------------------------------ \n finished SSM", i, ", convergence was:",
        ssm_results[[i]]$mess, 
        "\n ------------------------------------------------------------ \n ")
    
    
    # add to the iteration counter
    i=i+1
  } # close while loop

  
  #####################################
  ############ get results ############
  #####################################
  
  # nll values, convergence and messages, stuff that wasn't tracked in the for loops
  nll_hmm <- sapply(hmm_results, function(x)x$nll)
  hmm_conv <- sapply(hmm_results, function(x)x$opt$convergence)
  hmm_mess <- sapply(hmm_results, function(x)x$opt$message)
  ssm_mess <- sapply(ssm_results, function(x)x$mess)
  
  # put the results into a more convenient data frame
  ssm_nll <- data.frame(iter=1:length(ssm_results), mess = ssm_mess, conv = do.call(rbind, ssm_conv), nll = do.call(rbind, nll_ssm))
  hmm_nll <- data.frame(iter=1:length(ssm_results), mess = hmm_mess, conv = hmm_conv, nll = nll_hmm)
  
  # the best of the models based on proper convergence and minimum ssm nll
  winner = ssm_nll[ssm_nll$conv==0,'iter'][which.min(ssm_nll[ssm_nll$conv==0,'nll'])]
  if(length(winner)==0){
    winner=1
    warning("no winning iteration reached; setting winner to 1")
  }
  # new obs
  xhat <- data.frame(t(do.call(cbind, ssm_results[[winner]]$xhat)))
  names(xhat) <- c("lon", "lat")
  preds <- data.frame(id = rep(as.character(tracknames), times=sapply(regobsdat, function(x)nrow(x$regx))),
                      date = do.call(c, lapply(regobsdat, function(x)x$xdates)),
                      lon = t(do.call(cbind, ssm_results[[winner]]$xhat))[,1]*scaleobs+as.numeric(firstloc[1]),
                      lat = t(do.call(cbind, ssm_results[[winner]]$xhat))[,2]*scaleobs+as.numeric(firstloc[2]),
                      bhat = do.call(c, hmm_results[[winner]]$bhat))
  
  # final results list
  rslts <- list(fullargs=fullargs, obs=obs, idxs=idxs, maxsteps=maxsteps, ts=ts, 
                isprojected=isprojected, firstobs=firstloc, scaleobs=scaleobs,
                hmm_results = hmm_results, ssm_results = ssm_results, 
                hmm_time = hmm_time, ssm_time = ssm_time,
                hmm_nll = hmm_nll, ssm_nll = ssm_nll, 
                fc_count = fc_count, winner = winner,
                preds=preds)
  class(rslts) <- "issm" 
  return(rslts)
  
}


#



####### bootstrap the errors

#' Use a parametric bootstrap to calculate the standard errors of a model fit
#' @export
#' @param mod The model (issm object) to bootstrap errors for
#' @param startseed The starting seed to keep track of the random generation
#' @param nsamples The number of samples to take for the bootstrap
#' @param savepath The path to save the results in (full path, including .rda ext)
bootstrapCI <- function(mod, 
                        startseed=42, 
                        nsamples=50, 
                        nsteps=20, 
                        savepath=NULL, 
                        setREs=FALSE,
                        setjidx=TRUE,
                        setlc=TRUE,
                        usefirstlocs=FALSE,
                        issm_args=list(),
                        jumpsd=c(0.05, 0.05)){
  
  # get everything from the mod
  nx <- nrow(mod$preds)
  ny <- nrow(mod$obs)
  acprob <- table(mod$obs$lc)
  theta <- mod$hmm_results[[mod$winner]]$params[row.names(mod$hmm_results[[mod$winner]]$params) %in% 'theta', 'Estimate']
  gamma <- mod$hmm_results[[mod$winner]]$params[row.names(mod$hmm_results[[mod$winner]]$params) %in% 'gamma', 'Estimate']
  sdx <- mod$hmm_results[[mod$winner]]$params[c('tau_lon', 'tau_lat'), 'Estimate']
  alpha <- mod$hmm_results[[mod$winner]]$params[row.names(mod$hmm_results[[mod$winner]]$params) %in% 'A', 'Estimate']
  psi <- mod$ssm_results[[mod$winner]]$params['psi', 'Estimate']
  if("working_gps_err" %in% row.names(mod$ssm_results[[mod$winner]]$params)){
    gpserr = exp(mod$ssm_results[[mod$winner]]$params[row.names(mod$ssm_results[[mod$winner]]$params) %in% 'working_gps_err', 'Estimate'])
  } else {
    gpserr = NULL
  }
  
  if(setREs){
    res <- list(b=mod$preds$b,
                x=as.matrix(mod$preds[,c("lon", 'lat')]))
  } else {
    res <- list(b=NULL,
                x=NULL)
  }
  
  if(setjidx){
    tracknames <- unique(mod$obs$trackid)
    # trackstartidx <- c(1, head(cumsum(rle(as.character(mod$obs$trackid))$lengths),-1)+1)
    # jidx <- data.frame(idx=mod$idxs[[tracknames[1]]]$idx, ji=mod$idxs[[tracknames[1]]]$jidx)
    # if(length(tracknames)>1){
    #   for(i in 2:length(tracknames)){
    #     jidx <- rbind(jidx,data.frame(idx=jidx[trackstartidx[i-1], 1] + mod$idxs[[tracknames[i]]]$idx, ji=mod$idxs[[tracknames[i]]]$jidx))
    #   }
    # }
    idxtmp <- list()
    for(i in 1:length(mod$idxs)){
      idxtmp[[i]] <- cbind(mod$idxs[[i]]$idx, mod$idxs[[i]]$jidx)
    }
    jidx <- data.frame(do.call(rbind, idxtmp))
    names(jidx) <- c("idx", "ji")
    trackidy <- mod$obs$trackid
    trackidx <- mod$preds$id
    obsdates <- mod$obs$date
  } else {
    jidx <- NULL
    trackidy <- NULL
    trackidx <- NULL
    obsdates <- NULL
  }
  
  if(setlc){
    lc <- mod$obs$lc
  } else {
    lc <- NULL
  }
  
  if(usefirstlocs){
    firstlocs = mod$preds[1:2,c("lon", "lat")]
    firstlocs = as.matrix((firstlocs - firstlocs[c(1,1),])/1000000)
  } else {
    firstlocs = NULL
  }
  
  
  # initialize
  sims <- list()
  mods <- list()
  
  for(i in 1:nsamples){
    
    seed = startseed-1+i
    
    # need to generalize this
    sims[[i]] <- genTrack(nx=nx,
                          ny=ny,
                          process_pars = list(theta = theta,
                                              gamma = gamma,
                                              sdx = sdx),
                          start_seed=seed, 
                          alpha=matrix(alpha, nrow=sqrt(length(alpha))),
                          me="t", 
                          psi=psi,
                          gpserr=gpserr,
                          acprob=acprob,
                          res=res, 
                          jidx=jidx,
                          lc=lc,
                          trackidy=trackidy,
                          trackidx=trackidx,
                          obsdates=obsdates,
                          firstlocs=firstlocs,
                          startloc=mod$firstobs,
                          scaleobs=mod$scaleobs, 
                          jumpsd=jumpsd)


    issmargs = issm_args
    issmargs$obs = sims[[i]]$obs
    issmargs$ts = mod$ts
    issmargs$maxsteps = nsteps
    mods[[i]] <- tryCatch({
      do.call(fit_issm, issmargs)
      },
      error=function(e)e)
    # mods[[i]] <- tryCatch({fit_issm(obs = sims[[i]]$obs,
    #                                 ts = mod$ts,
    #                                 maxsteps = nsteps,
    #                                 ...)})
    # mods[[i]] <- tryCatch({fit_issm(obs = sims[[i]]$obs,
    #                                 ts = mod$ts,
    #                                 maxsteps = nsteps,
    #                                 p_start_hmm = list(working_theta = c(log(1), log(1)),
    #                                                    working_gamma = c(log(1), log(1)),
    #                                                    working_tau_lon=log(1), working_tau_lat=log(1),
    #                                                    working_A=matrix(c(log(0.75/0.25), log(0.25/0.75)),ncol=1,nrow=2)),
    #                                 init_psi = 1,
    #                                 res = c("x"),
    #                                 fixsteps = TRUE,
    #                                 allowfc=TRUE,
    #                                 jiggle_fc=0.01,
    #                                 allsilent=TRUE,
    #                                 ssm_map <- list(working_gamma = factor(c(NA, NA)),
    #                                                 working_theta = factor(c(NA, NA)),
    #                                                 working_sigma_lon=factor(NA), working_sigma_lat=factor(NA)))})
    
    # save it
    if(!is.null(savepath)) saveRDS(mods, file=paste(savepath, '_bootstrapmods.RDS', sep=''))
    
    # update on where we are in the iteration
    cat("\n ------------------------------------------------------------ \n finished sim study", i,
        "\n ------------------------------------------------------------ \n ")
    
  }
  
  nullidx <- which(sapply(mods, length)==1)  
  if(length(nullidx)>0) for(i in 1:length(nullidx)) mods[[nullidx[i]+1-i]] <- NULL
  
  m <- length(which(row.names(mod$hmm_results[[mod$winner]]$params) == "theta"))
  pars <- rbind(sapply(mods, function(x)x$hmm_results[[x$winner]]$params[row.names(x$hmm_results[[x$winner]]$params) %in% 'theta', 'Estimate']),
  sapply(mods, function(x)x$hmm_results[[x$winner]]$params[row.names(x$hmm_results[[x$winner]]$params) %in% 'gamma', 'Estimate']),
  sapply(mods, function(x)x$hmm_results[[x$winner]]$params['tau_lon', 'Estimate']),
  sapply(mods, function(x)x$hmm_results[[x$winner]]$params['tau_lat', 'Estimate']),
  sapply(mods, function(x)x$hmm_results[[x$winner]]$params[row.names(x$hmm_results[[x$winner]]$params) %in% 'A', 'Estimate']),
  sapply(mods, function(x)x$ssm_results[[x$winner]]$params['psi', 'Estimate'])
  ) %>% as.data.frame() 
  
  # parameters
  true.pars <- c(theta, gamma, sdx, alpha, psi)
  par.stats <- pars %>% 
    mutate(mean = rowMeans(.),
           median = apply(., 1, median), 
           sd = apply(., 1, sd),
           lower2.5 = apply(., 1, quantile, probs=0.025),
           upper97.5 = apply(., 1, quantile, probs=0.975),
           meanbias = (rowMeans(.)-true.pars),
           rmse = sqrt(rowMeans((. - matrix(rep(true.pars, length(mods)), ncol=length(mods)))^2)),
           par=c(rep("theta", m), rep("gamma", m), "sigma_lon", "sigma_lat", rep("A", m*m), "psi")) %>% 
    select(par, mean, median, sd, lower2.5, upper97.5, meanbias, rmse)
  
  # behavioural states
  err.rate = colMeans((sapply(mods, function(x)x$preds$bhat) - sapply(sims, function(x)x$b))^2) #only accurate for two states
  b.stats <- sapply(mods, function(x)x$preds$bhat) %>% as.data.frame %>%
    mutate(mean=rowMeans(.),
           mode=apply(., 1, function(x)unique(x)[which.max(table(x))]),
           err.rate= rowMeans((.-sapply(sims, function(x)x$b))^2)) %>% 
    select(mean, mode, err.rate) 

  #location states
  lon.rmse = sqrt(colMeans((sapply(mods, function(x)x$preds$lon) - sapply(sims, function(x)x$X[,1]))^2)) 
  x.stats.lon <- sapply(mods, function(x)x$preds$lon) %>% as.data.frame %>% 
    mutate(mean=rowMeans(.),
           median = apply(., 1, median), 
           sd=apply(., 1, sd),
           lower2.5 = apply(., 1, quantile, probs=0.025),
           upper97.5 = apply(., 1, quantile, probs=0.975),
           meanbias = (rowMeans(. - sapply(sims, function(x)x$X[,1]))),
           rmse = sqrt(rowMeans((. - sapply(sims, function(x)x$X[,1]))^2))) %>% 
    select(mean, median, sd, lower2.5, upper97.5, meanbias, rmse)
  lat.rmse = sqrt(colMeans((sapply(mods, function(x)x$preds$lat) - sapply(sims, function(x)x$X[,2]))^2)) 
  x.stats.lat <- sapply(mods, function(x)x$preds$lat) %>% as.data.frame %>% 
    mutate(mean=rowMeans(.),
           median = apply(., 1, median), 
           sd=apply(., 1, sd),
           lower2.5 = apply(., 1, quantile, probs=0.025),
           upper97.5 = apply(., 1, quantile, probs=0.975),
           meanbias = (rowMeans(. - sapply(sims, function(x)x$X[,2]))),
           rmse = sqrt(rowMeans((. - sapply(sims, function(x)x$X[,2]))^2))) %>% 
    select(mean, sd, lower2.5, upper97.5, meanbias, rmse)

  # also add a method for viewing
  boot <- list(sims=sims, 
               mods=mods,
               true.pars=true.pars,
               pars=pars, 
               par.stats=par.stats, 
               err.rate=err.rate, 
               b.stats=b.stats,
               lon.rmse=lon.rmse,
               x.stats.lon=x.stats.lon,
               lat.rmse=lat.rmse,
               x.stats.lat=x.stats.lat)
  class(boot) <- 'swim_bootstrap'
  if(!is.null(savepath)) saveRDS(boot, file=paste(savepath, '_bootstrap.RDS', sep=''))
  
  
  return(boot)
  
}


