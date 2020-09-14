#' Simulate Correlated Random Walk Data
#'
#'  Simulate data from a DCRWS movement process with or without measurement error. 
#' @export 
#' @param nx Number of simulated true locations. 
#' @param ny Number of simulated observatioons, ny does not necessarily = nx.
#' @param date1 Starting date.
#' @param ts Time step.
#' @param start_seed Starting seed.
#' @param process_pars List of relevant parameters. 
#' @param alpha Matrix of switching probabilities. Determines the number of states you wish to have. 
#' @param start_loc Starting location in c(lon, lat). If NULL assumed to be (0,0).
#' @param me Measurement error distribution, Gaussian or t.
#' @param acprob Vector of relative proportions for drawing the numbers of observed locations from each Argos class   
#' @return Simulated data set.
genTrack <- function(nx=500, ny=1500, 
                     date1=Sys.time(), ts=3, start_seed=42,
                     process_pars = list(theta = c(0, pi),
                                         gamma = c(0.8,0.3),
                                         sdx = rep(0.1, 2)), 
                     alpha = matrix(c(0.95, 0.05, 0.05, 0.95), nrow=2, ncol=2),
                     me='t', acprob=NULL, psi = rep(0.3, 2),
                     res=NULL, lc=NULL, jidx=NULL, trackid=NULL, obsdates=NULL){
  
  set.seed(start_seed)
  
  
  ##########################
  ######  behaviours  ######
  ##########################
  
  if(is.null(res$b)){
    m=nrow(alpha)
    B <-  1:m
    b <- numeric(nx)
    b[1] <- sample(B, size=1, prob=rep(1/m, m))
    for(i in 2:nx) b[i] = sample(B, size=1, prob=alpha[b[i-1],])
  } else {
    b = res$b
  }

  
  
  ##############################
  ######  true locations  ######
  ##############################
  
  if(is.null(res$x)){
    X <- matrix(nrow=nx, ncol=2)
    X[1,] <- rep(0,2)
    X[2,] = mvtnorm::rmvnorm(1, X[1,], diag(process_pars$sdx, nrow=2)) 
    for(i in 3:nx){
      Tr <- matrix(c(cos(process_pars$theta[b[i]]), -sin(process_pars$theta[b[i]]),
                     sin(process_pars$theta[b[i]]), cos(process_pars$theta[b[i]])), 
                   nrow=2, byrow=TRUE)
      X[i,] <- mvtnorm::rmvnorm(1, X[i-1,] + process_pars$gamma[b[i]] * Tr %*% (X[i-1,] - X[i-2,]), diag(process_pars$sdx)) 
    } 
  } else {
    X <- res$x
  }


  
  #observed locations
  
  # regularization 
  if(is.null(jidx)){
    idx <- sample(1:(nx-1), ny, replace=TRUE)
    ji <- runif(ny)
    jidx <- data.frame(cbind(idx,ji))
    jidx <- jidx[do.call("order",jidx[c("idx","ji")]),]
    if(jidx[nrow(jidx),1] != (nx-1)) jidx[nrow(jidx),1] = nx-1
    # uniform dists for the idx and jidx
    # set the last one to the last data point to avoid simulation problems
    # sometimes there are no simulated observations in the final x interval
    # to avoid model fitting problems I need to have data in the last interval;
    # set seed to 1 or 7 before running gentrack , then check tail idx to see example
  } else {
    jidx <- jidx
  }

  
  
  Y <- matrix(nrow=ny, ncol=2)
  if(me == 'gau'){
    for(i in 1:ny){
      y <- (1 - jidx[i,2]) * X[jidx[i,1],] + jidx[i,2] * X[jidx[i,1]+1,]
      Y[i,] = y  + mvtnorm::rmvnorm(1, rep(0,2), diag(psi^2, nrow=2))
    }
  } else if(me == 't'){
      all_classes <- c("Z", "B", "A", "0", "1", "2", "3")
      if(is.null(lc)){
        lc <- sample(all_classes, ny, replace=TRUE, prob=acprob)
      } else {
        lc <- lc
      }

      #from marie's code
      ############# parameter values for the argos classes
      # Make Z class = B class
      sigmalon <- (c(4.2050261, 4.2050261, 0.507292, 2.1625936, 0.9020423, 0.3119293, 0.289866)/6366.71 * 180)/pi
      sigmalat <- (c(3.041276, 3.041276, 0.5105468, 1.607056, 0.4603374, 0.2605126, 0.1220553)/6366.71 * 180)/pi
      # If nu < 2 bsam sets it to 2
      nulon <- c(2, 2, 2, 2, 2.298819, 2, 3.070609)
      nulat <- c(2, 2, 2, 2, 3.896554, 6.314726, 2.075642)
      
      argos_error <- cbind(sigmalon, nulon, sigmalat, nulat)
      row.names(argos_error) <- all_classes
      data_error <- as.data.frame(argos_error[lc,])

      for(i in 1:ny){
        y <- (1 - jidx[i,2]) * X[jidx[i,1],] + jidx[i,2] * X[jidx[i,1]+1,]
        Y[i,1] <- y[1] + rt(1, data_error$nulon[i])*data_error$sigmalon[i]/psi    # longitude
        Y[i,2] <- y[2] + rt(1, data_error$nulat[i])*data_error$sigmalat[i]/psi    # latitude
      }
  }


  # get the simulated dates
  # check the indexing on this
  date_x <- seq(date1, date1+ts*3600*nx, ts*3600)
  if(is.null(obsdates)){
    date_y <- date1
    for(i in 2:ny) date_y = append(date_y, (jidx[i,2]) * 3600 * 3 + date_x[jidx[i,1]-1])
  } else {
    date_y <- obsdates
  }

  obs <- data.frame(date=date_y, lon=Y[,1], lat=Y[,2])
  if(is.null(trackid)){
    obs$trackid = 'track1'
  } else {
    obs$trackid = trackid
  }
  if(me=='t') obs$lc = lc
  
  rslt <- list(b=b, X=X, Y=Y, jidx=jidx,
                 date_x=date_x, date_y=date_y, obs=obs)
  return(rslt)
  
}

