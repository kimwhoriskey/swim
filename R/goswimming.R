

############################### All of the fitting functions I will need ################################


#############################################################
#################### windsorize the data ####################
#############################################################

#' Winsorize some observations
#' @param obs dataframe of observations, 'id', date', 'lon', 'lat', 'lc'
#' @param speedtol the speed tolerance for the animal. If "empirical" this is set to 
#' @param inflate the multiplicative inflation for the speedtol
#' @param upperlimscale the threshold value determining how far back the outliers should be scaled, = upperlimscale*IQR
#' @param proj are the data projected
#' @param returnfull return the index of the outliers as well as the winsorized data? 
#' @export
winsorize <- function(obs, speedtol=10, inflate=10, upperlimscale=1.5, proj=FALSE, returnfull=FALSE){
 
  # generate new observation data frame
  newobs <- obs
  
  # find the outliers based on the speed
  # step lengths first, units m
  if(proj){
    empty <- st_as_sfc("POINT(EMPTY)")
    sl <- c(0, head(st_distance(obs$geometry, lead(obs$geometry, default=empty), by_element=TRUE), -1))
  } else {
    sl <- getStepLengths(as.matrix(obs[,c("lon", "lat")]), convert_to_rad=FALSE, units="m")
  }
  # now divide by the time to get speed, units m/s
  if("date" %in% names(obs)){
    secs <- c(1, as.numeric(diff(obs$date)))
    if(any(secs==0)) warning("there are multiple observations for some datetime stamps")
    # secs[which(secs==0)] <- 1 # set to one second if it's 0 seconds
    sl <- sl/secs 
  } 
  
  # get the quantile stats
  bs <- boxplot.stats(sl)
  # set an upperlimit to the step length value that we want to scale by
  # default is 1.5*IQR
  upperlim <- bs$stats[4] + upperlimscale*(bs$stats[4]-bs$stats[2])
  # can set to an empirical value or a pre-determined value based on swim speed
  if(speedtol == "empirical"){
    tol = upperlim
  } else {
    tol = speedtol*inflate
  }
  outidx <- which(sl > tol) 
  
  # first figure out if the first location is an outlier
  # if so, remove it entirely
  # not sure if it is worth bringing it closer to the second location
  # also have to create an offset value for the indices for newobs
  offset <- 0 
  if(2 %in% outidx){ # need to have an and condition here in case the second one is really the outlier
    newobs <- newobs[-1,] 
    outidx <- outidx[-1] #should be the first element of outidx, hopefully this doesn't bite me in the butt
    cat("\n the first location of seal", as.character(newobs$id[1]), "was an outlier and was removed completely \n")
    offset <- 1
  }
 
  # this part is tricky, basically every single outlier based on step length will have a partner outlier
  # so basically what we are doing is finding all of the pairs of outliers
  # this pairs assumption is sometimes violated
  # in that case tho, you can run the windsorize code multiple times (jsut run it on output of the 
  # windsorize function), and it should work
  diffs <- diff(outidx)
  runs <- rle(diffs)
  runsidx <- which( (runs$values==1) & (runs$lengths>1) )
  runslen <- runs$lengths[runsidx]
  rmfromdups <- numeric()
  if(length(runslen)>0){ # condition on there being sets of outliers occuring one after the other
    for(i in 1:length(runsidx)){
      # keep the odd values, which are actual outliers, ditch the events which are pairs
      odds = seq(1, runslen[i], 2) 
      rmfromdups <- append(rmfromdups, cumsum(runs$lengths)[runsidx[i]-1] + odds)
    }
  }
 
  dupidx <- setdiff(which(diff(outidx)==1) + 1, rmfromdups) 
  # finally, remove the duplicated outliers from the original outlying index
  outidx <- outidx[-dupidx]
  
  # calculate the speed ratio based on the upperlim
  dr <- upperlim/sl[outidx] 
  
  # replace outliers
  # basically its a projection of the vector onto itself, just scaling it back a bit
  # so this will technically affect the turning angles of the next two locations
  # include the offset variable because the indices of newobs are one less than the indices
  # of obs if we have removed the first location
  # this should only kick in if you've removed the first location
  if(length(outidx) > 0 ){
      
    if(proj){
        coos <- st_coordinates(newobs)
        newcoos <- coos
      for(i in 1:length(outidx)){
        newcoos[outidx[i] - offset, c("X", "Y")] <- coos[outidx[i]-1, c("X", "Y")] + 
          dr[i]*(coos[outidx[i], c("X", "Y")] - coos[outidx[i]-1, c("X", "Y")])
      }
      st_geometry(newobs) <- st_geometry(st_as_sf(data.frame(newcoos), 
                                                  coords=c("X", "Y"), 
                                                  crs=st_crs(obs)))
    } else {
      if(length(outidx) > 0 ){
        for(i in 1:length(outidx)){
          newobs[outidx[i] - offset, c("lon", "lat")] <- obs[outidx[i]-1, c("lon", "lat")] + 
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

  # return object
  if(returnfull){
    return(list(outliers=outidx, obs=newobs))
  } else {
    return(newobs)
    
  }
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
  
  
  
  # create TMB object
  obj <- MakeADFun(dat, p_start, DLL='swim', silent=silence)
  
  # optimize TMB object
  opt <- nlminb(obj$par, obj$fn, obj$gr)
  
  # likelihood value
  nll <- obj$fn(opt$par)
  
  # extract behavioural states
  b_states <- obj$report()$states
  
  # extract parameter estimates
  params_hmm <- summary(sdreport(obj))
  
  # get the stationary distribution
  stationary <- obj$report()$delta
  
  # pseudoresiduals
  
  # return results
  rslt <- list(mess = opt$mess, obj=obj, opt=opt, nll=nll, params = params_hmm, 
               b_hat=b_states, statdist = stationary)
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
                    silence=FALSE){
  
  #extract the names of the fixed variables
  fixed_names <- names(mapping)
  
  # set up an error vector so that we can keep track of how many errors occur
  err_num = 0 # set initial error value at 0, base an iterative loop on this
  err = 0
  class(err) <- "error"
  
  while("error"%in%class(err) & allowErr==FALSE){
    
    if(err_num>10) stop("encountered errors > 10 times while attempting to fit the SSM") #cut optimization if errors > 10
    
    if(err_num>0){
      err_num = err_num+1 #count the errors
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
    if(include_TMB_gr){
      argus <- append(list(obj$par, obj$fn, obj$gr), optimizer_arguments)
    } else {
      argus <- append(list(obj$par, obj$fn), optimizer_arguments)
    }

    err <- tryCatch({
      opt <- do.call(optimizer, argus) # optimize TMB object
      srep <- summary(sdreport(obj)) # extract parameter estimates
    },  error=function(e){e} #capture errors
    )
    
  } # close while loop
  
    # extract the likelihood value
  nll <- obj$fn(opt$par)
  
  # parameter estimates, take the x's out of it
  params_ssm <- srep[! rownames(srep) %in% res,] 
  
  # extract location state estimates
  if(is.null(res)) res <- "nothing"
  if((!unique(res %in% fixed_names))){ # in case the states are fixed
    rehat = srep[rownames(srep) %in% res,]
    if(unique(res=="x")){
      xhat = data.frame(matrix(rehat[,1], nrow=nrow(rehat)/2, byrow=TRUE))
      names(xhat) = c("lon", "lat")
    } else {
      ta = c(0, 0, obj$report()$ta)
      sl = c(0, obj$report()$sl)
      xhat <- data.frame(ta=ta, sl=sl)
    }
  } else {
    xhat=NULL
  }

  rslt <- list(mess = opt$mess, obj=obj, opt=opt, nll=nll, 
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
get_parms <- function(prev_hmm, prev_ssm, curr_model){
  
  p <- list(working_theta = as.numeric(prev_hmm$params[row.names(prev_hmm$params) %in% 'working_theta', 'Estimate']), #have to check that these are the right thetas
            working_gamma = as.numeric(prev_hmm$params[row.names(prev_hmm$params) %in% 'working_gamma','Estimate']),
            working_tau_lon = prev_hmm$params['working_tau_lon', 'Estimate'], 
            working_tau_lat = prev_hmm$params['working_tau_lat', 'Estimate'])

  if(curr_model == "ssm"){
    
    # overwrite the previous process error 
    p$working_tau_lon <- p$working_tau_lat <- NULL
    p$working_sigma_lon <- prev_hmm$params['working_tau_lon', 'Estimate']
    p$working_sigma_lat <- prev_hmm$params['working_tau_lat', 'Estimate']
    
    if(!is.null(prev_ssm)){
      
      # measurement error
      p$working_psi <- as.numeric(prev_ssm$params[row.names(prev_ssm$params) %in% 'working_psi', 'Estimate'])
      # add in the locations
      p$x = matrix(log(1), nrow=2, ncol=length(prev_hmm$b_hat)) 
      # start at 0 regardless, but we can change this if we want

    }
    
  } 
    
  p$working_A = matrix(c(prev_hmm$params[row.names(prev_hmm$params) %in% 'working_A', 'Estimate'][1:2]),ncol=1,nrow=2)

  return(p)
}

######################################################################################
##################### Iteratively fit multiple SSM and HMM steps #####################
######################################################################################

#' fit a switching ssm iteratively 
#' @export
#' @param obs the data
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
#' @param scaleobs value to scale the observations by (literally just divided), helpful if they are in eastings/northings
#' @param ... further arguments to fit_ssm
fit_issm <- function(obs, 
                     ts, 
                     p_start_hmm,
                     init_psi,
                     res=c("x"),
                     maxsteps=10, fixsteps=FALSE, 
                     allowfc=FALSE, jiggle_fc=0.01, 
                     ssm_map=NULL,
                     allsilent=TRUE,
                     scaleobs=1, ...){


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
  
  
  #convert from sf
  if("sf" %in% class(obs)){
    obs <- obs %>% bind_cols(as_tibble(st_coordinates(obs))) %>% rename(lon=X, lat=Y) %>% as.data.frame() 
  }
  
  # if(scaleobs){
  #   obsscale <- obs[1, c("lon", "lat")]
  #   obs[,c("lon", "lat")] <- (obs[,c("lon", "lat")] - obsscale[rep(1, each=nrow(obs)),])/100000
  # } else {
  #   scaleobs <- 1
  # }
  obs[,c("lon", "lat")] <- obs[,c("lon", "lat")]/scaleobs
  # now it's in km
  
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
  
  # do everything else at the end using sapply instead
  
  #####################################################################
  ############ initial step, get intial location estimates ############
  #####################################################################
  
  idxs <- getJidx(obs, ts = ts)
  ae <- as.matrix(getAE(obs$ac))
  
  regobs <- regTrack(obs, ts=ts)
  nx <- nrow(regobs$regx)
  hmm_dat <- list(model=0, x=t(regobs$regx))


  
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
  ssm_dat <- list(model=1,
                  dat=list(y = t(array(c(obs$lon, obs$lat), dim=c(nrow(obs), 2))), 
                           b = hmm_results[[1]]$b_hat-1,
                           idx = idxs$idx, 
                           jidx = idxs$jidx,
                           ae=ae))
  # ssm_dat$ae = ae # options for either gaussian or t error? 

         
  # now get the parameters             
  p_start_ssm <- get_parms(prev_hmm = hmm_results[[1]], 
                           prev_ssm = NULL,
                           curr_model = "ssm")
  p_start_ssm$working_psi <- log(init_psi)
  p_start_ssm$x = matrix(log(1), nrow=2, ncol=nrow(regobs$regx)) #easy to start at 0

  # fit the ssm
  # block for false convergence
  # the idea is to jiggle the starting parameters by adding a teeny bit of random error
  # exclude the location states as well as anything that is fixed
  # if we fix all of the ssm parameters this may basically be moot
  # also count the numbers of false convergence so that we know, also so that we can cut off the optimizer
  ssm_fixed_names <- names(ssm_map) # record the names of the fixed parameters
  pidx <- which(!names(p_start_ssm) %in%  c("x", "sl", "ta", ssm_fixed_names)) #create an index to loop over, p for parameter
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
    ssm_time[[1]] <- system.time(ssm_results[[1]] <- fit_ssm(ssm_dat, 
                                                             p_start = p_start_ssm,
                                                             res=res,
                                                             mapping = ssm_map,
                                                             # inner_control = list(maxit=20000, step.tol=1e-4, grad.tol=1e-2),
                                                             silence=allsilent,
                                                             ...))
      # save the nll 
      nll_ssm[[1]] = ssm_results[[1]]$nll
      # save the new convergence, cue on this
      ssm_conv[[1]] <- ssm_results[[1]]$opt$convergence
  
    if(ssm_conv[[1]]>0) fc_count[1] = fc_count[1] + 1 # add to the count every time we come back up to the beginning of this loop
    if(fc_count[1]>=10) warning("had false convergence 10 times while trying to get the SSM to fit")
    # stop if greater than 10 times
    if(ssm_conv[[1]]==0 | allowfc==TRUE | fc_count[1]>=10) runme=FALSE
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
    hmm_dat <- list(model=0, x=t(ssm_results[[i-1]]$xhat))
  
    # fit HMM again
    p_start_hmm <- get_parms(prev_hmm = hmm_results[[i-1]], 
                             prev_ssm = ssm_results[[i-1]],
                             curr_model = "hmm")
    
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
    ssm_dat$dat$b = hmm_results[[i]]$b_hat-1
    
    # get starting parameters for the ssm
    # starting values come from the most recent HMM, except for psi which comes from the last ssm
    p_start_ssm <- get_parms(prev_hmm = hmm_results[[i]], 
                             prev_ssm = ssm_results[[i-1]],
                             curr_model = "ssm")
    
    
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
      
      ssm_time[[i]] <- system.time(ssm_results[[i]] <- fit_ssm(ssm_dat, 
                                                               p_start = p_start_ssm,
                                                               res=res,
                                                               mapping = ssm_map,
                                                               # inner_control = list(maxit=20000, step.tol=1e-4, grad.tol=1e-2),
                                                               silence=allsilent,
                                                               ...))
        # save the nll 
      nll_ssm[[i]] = ssm_results[[i]]$nll
      # save the new convergence, cue on this
      ssm_conv[[i]] <- ssm_results[[i]]$opt$convergence
      # calculate the difference between this and the last nll
      # should actually maybe make this a convergence criterion? 
      nll_diff <- nll_ssm[[i]]-nll_ssm[[i-1]]
    
      
      if(ssm_conv[[i]]>0) fc_count[i] = fc_count[i] + 1
      if(fc_count[i]>=10) warning("had false convergence >10 times while trying to get the SSM to fit")
      
      # condition to close
      
      if(ssm_conv[[i]]==0 | allowfc==TRUE | fc_count[i]>=10) runme=FALSE
      
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
  ssm_mess <- sapply(ssm_results, function(x)x$opt$message)
  
  # put the results into a more convenient data frame
  ssm_nll <- data.frame(iter=1:length(ssm_results), mess = ssm_mess, conv = do.call(rbind, ssm_conv), nll = do.call(rbind, nll_ssm))
  hmm_nll <- data.frame(iter=1:length(ssm_results), mess = hmm_mess, conv = hmm_conv, nll = nll_hmm)
  
  # the best of the models based on proper convergence and minimum ssm nll
  winner = ssm_nll[ssm_nll$conv==0,'iter'][which.min(ssm_nll[ssm_nll$conv==0,'nll'])]
  

  
  # final results list
  rslts <- list(obs = obs, scaleobs=scaleobs, maxsteps = maxsteps, ts = ts, 
                hmm_results = hmm_results, ssm_results = ssm_results, 
                hmm_time = hmm_time, ssm_time = ssm_time,
                hmm_nll = hmm_nll, ssm_nll = ssm_nll, 
                fc_count = fc_count, winner = winner)
  class(rslts) <- "issm" 
  return(rslts)
  
}


#' mmm not sure about this one yet
runSimStudy <- function(startseed, nsims, nsteps, path, ...){
  
  sims <- list()
  mods <- list()
  
  for(i in 1:nsims){
    
    seed = startseed-1+i
    
    # need to generalize this
    sims[[i]] <- genTrack(n_x=644,
                    n_y=1314,
                    process_pars = list(theta = c(0.0132, 0.9345),
                                        gamma = c(0.9433,0.2045),
                                        sdx = c(0.0272, 0.0235)),
                    move="dcrw", 
                    start_seed=seed, 
                    alpha=matrix(c(0.6952, 0.3048, 0.3322, 0.6678), nrow=2, byrow=TRUE),
                    me="t", 
                    psi=0.5531)

    
    #maybe try 10 and 13
    
    mods[[i]] <- tryCatch({fit_issm(obs = sims[[i]]$obs,
                                    move = "dcrw",
                                    ts = 3,
                                    p_start_hmm = list(working_theta = c(log(1), log(1)),
                                                       working_gamma = c(log(1), log(1)),
                                                       working_tau_lon=log(1), working_tau_lat=log(1),
                                                       working_A=matrix(c(log(0.75/0.25), log(0.25/0.75)),ncol=1,nrow=2)),
                                    init_psi = 1,
                                    res = c("x"),
                                    maxsteps = 20,
                                    fixsteps = TRUE,
                                    allowfc=FALSE,
                                    jiggle_fc=0.01,
                                    transform_random_effects=0,
                                    split=FALSE,
                                    allsilent=FALSE,
                                    ssm_map <- list(working_gamma = factor(c(NA, NA)),
                                                    working_theta = factor(c(NA, NA)),
                                                    working_sigma_lon=factor(NA), working_sigma_lat=factor(NA)))})
    # save it
    save(mods, file=path)
    # 
    # update on where we are in the iteration
    cat("\n ------------------------------------------------------------ \n finished sim study", i, 
        "\n ------------------------------------------------------------ \n ")    
    
  }
  
  return(list(sims=sims, mods=mods))
  
}


