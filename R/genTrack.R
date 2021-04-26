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
                     date1=Sys.time(),
                     ts=3, 
                     start_seed=42,
                     process_pars = list(theta = c(0, pi),
                                         gamma = c(0.8,0.3),
                                         sdx = rep(0.1, 2)), 
                     alpha = matrix(c(0.95, 0.05, 0.05, 0.95), nrow=2, ncol=2),
                     me='t', 
                     gpserr=NULL, 
                     acprob=NULL, 
                     psi = rep(0.3, 2),
                     res=NULL, 
                     lc=NULL, 
                     jidx=NULL, 
                     trackidx=NULL,
                     trackidy=NULL, 
                     obsdates=NULL,
                     firstlocs=NULL,
                     jumpsd=c(1,1),
                     startloc=matrix(c(0,0), nrow=1), 
                     scaleobs=1){
  
  set.seed(start_seed)
  
  if(is.null(trackidx)){
    tracknames="track1"
    trackidx <- rep("track1", nx)
  } else {
    tracknames <- unique(trackidx)
  }
  if(length(tracknames)==1){
    trackstartidx <- c(1, nx+1)
  } else {
    trackstartidx <- c(1, head(cumsum(rle(as.character(trackidx))$lengths),-1)+1, length(trackidx)+1)
  }

  # tmp <- numeric()
  # for(i in 1:length(mod$idxs)) tmp <- append(tmp, mod$idxs[[i]]$nx)
  # tmp2 <- numeric()
  # for(i in 1:length(mod$idxs)) tmp2 <- append(tmp2, max(mod$idxs[[i]]$idx))

  
  ##########################
  ######  behaviours  ######
  ##########################
  
  if(is.null(res$b)){
    m=nrow(alpha)
    B <-  1:m
    b <- numeric(nx)
    for(i in 1:(length(trackstartidx)-1)){
      b[trackstartidx[i]] <- sample(B, size=1, prob=rep(1/m, m))
      for(j in (trackstartidx[i]+1):(trackstartidx[i+1]-1)) b[j] = sample(B, size=1, prob=alpha[b[j-1],])
    }
    # b[1] <- sample(B, size=1, prob=rep(1/m, m))
    # for(i in 2:nx) b[i] = sample(B, size=1, prob=alpha[b[i-1],])
  } else {
    b = res$b
  }

  
  
  ##############################
  ######  true locations  ######
  ##############################
  
  if(is.null(res$x)){
    X <- matrix(nrow=nx, ncol=2)
    if(is.null(firstlocs)){
      X[1,] <- rep(0,2) + mvtnorm::rmvnorm(1, rep(0,2), diag(process_pars$sdx^2, nrow=2)) 
      X[2,] <- mvtnorm::rmvnorm(1, X[1,], diag(process_pars$sdx^2, nrow=2)) 
    } else {
      X[1:2,] <- firstlocs
    }
    # X[1,] <- rep(0,2) + mvtnorm::rmvnorm(1, rep(0,2), diag(process_pars$sdx, nrow=2)) 
    # X[2,] = mvtnorm::rmvnorm(1, X[1,], diag(process_pars$sdx, nrow=2)) 
    for(i in 3:(trackstartidx[2]-1)){
      Tr <- matrix(c(cos(process_pars$theta[b[i]]), -sin(process_pars$theta[b[i]]),
                     sin(process_pars$theta[b[i]]), cos(process_pars$theta[b[i]])), 
                   nrow=2, byrow=TRUE)
      X[i,] <- mvtnorm::rmvnorm(1, X[i-1,] + process_pars$gamma[b[i]] * Tr %*% (X[i-1,] - X[i-2,]), diag(process_pars$sdx^2)) 
    } 
    if(length(trackstartidx)>2){
      for(j in 2:(length(trackstartidx)-1)){
        X[trackstartidx[j],] <- X[trackstartidx[j]-1,] + mvtnorm::rmvnorm(1, rep(0,2), diag(jumpsd, nrow=2)) 
        X[trackstartidx[j]+1,] <- mvtnorm::rmvnorm(1, X[trackstartidx[j],], diag(process_pars$sdx^2, nrow=2)) 
        for(i in (trackstartidx[j]+2):(trackstartidx[j+1]-1)){
          Tr <- matrix(c(cos(process_pars$theta[b[i]]), -sin(process_pars$theta[b[i]]),
                         sin(process_pars$theta[b[i]]), cos(process_pars$theta[b[i]])), 
                       nrow=2, byrow=TRUE)
          X[i,] <- mvtnorm::rmvnorm(1, X[i-1,] + process_pars$gamma[b[i]] * Tr %*% (X[i-1,] - X[i-2,]), diag(process_pars$sdx^2)) 
        }
      }
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

  
  if(length(tracknames)==1){
    trackstartjidx <- c(1, ny+1)
  } else {
    trackstartjidx <- c(1, head(cumsum(rle(as.character(trackidy))$lengths),-1)+1, length(trackidy)+1)
  }
  
  Y <- matrix(nrow=ny, ncol=2)
  if(me == 'gau'){
    for(j in 1:(length(trackstartidx)-1)){
      for(i in trackstartjidx[j]:(trackstartjidx[j+1]-1)){
        y <- (1 - jidx[i,2]) * X[trackstartidx[j]+jidx[i,1]-1,] + jidx[i,2] * X[trackstartidx[j]+jidx[i,1],]
        Y[i,] = y  + mvtnorm::rmvnorm(1, rep(0,2), diag(psi^2, nrow=2))
      }
    }
  } else if(me == 't'){
      all_classes <- c("Z", "B", "A", "0", "1", "2", "3", "GPS")
      if(is.null(lc)){
        lc <- sample(all_classes, ny, replace=TRUE, prob=acprob)
      } else {
        lc <- lc
      }

      data_error <- as.data.frame(getAE(lc))

      if(is.null(gpserr)){
        for(j in 1:(length(trackstartidx)-1)){
          for(i in trackstartjidx[j]:(trackstartjidx[j+1]-1)){
            xhat <- (1 - jidx[i,2]) * X[trackstartidx[j]+jidx[i,1]-1,] + jidx[i,2] * X[trackstartidx[j]+jidx[i,1],]
            Y[i,1] <- xhat[1] + rt(1, data_error$nulon[i])*data_error$sigmalon[i]/sqrt(psi)    # longitude
            Y[i,2] <- xhat[2] + rt(1, data_error$nulat[i])*data_error$sigmalat[i]/sqrt(psi)    # latitude
          }   
        }
      } else {
        for(j in 1:(length(trackstartidx)-1)){
          for(i in trackstartjidx[j]:(trackstartjidx[j+1]-1)){
            xhat <- (1 - jidx[i,2]) * X[trackstartidx[j]+jidx[i,1]-1,] + jidx[i,2] * X[trackstartidx[j]+jidx[i,1],]
            if(lc[i]=='GPS'){
              Y[i,1] <- xhat[1] + rnorm(1, 0, gpserr[1])    # longitude
              Y[i,2] <- xhat[2] + rnorm(1, 0, gpserr[2])    # latitude
            } else {
              Y[i,1] <- xhat[1] + rt(1, data_error$nulon[i])*data_error$sigmalon[i]/sqrt(psi)    # longitude
              Y[i,2] <- xhat[2] + rt(1, data_error$nulat[i])*data_error$sigmalat[i]/sqrt(psi)    # latitude
            }
          }          
        }
      }
 
  }

  X <- X*scaleobs + startloc[rep(1, nrow(X)),]
  Y <- Y*scaleobs + startloc[rep(1, nrow(Y)),]

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
  if(is.null(trackidy)){
    obs$trackid = 'track1'
  } else {
    obs$trackid = trackidy
  }
  if(me=='t') obs$lc = lc
  obs <- obs %>% mutate(kind=case_when(lc=="GPS"~"GPS", T~"Argos")) %>% 
    mutate(kindnum=case_when(kind=="GPS"~0, T~1))
  

  rslt <- list(b=b, X=X, trackidx=trackidx, Y=Y, jidx=jidx,
                 date_x=date_x, date_y=date_y, obs=obs)
  class(rslt) <- 'swimsim'
  return(rslt)
  
}

