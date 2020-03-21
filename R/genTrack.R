#' Simulate Correlated Random Walk Data
#'
#'  Simulate data from either a DCRWS movement process or a conditional autoregressive (step length and turning angle) process. Can simulate with or without measurement error. 
#' @export 
#' @param n_x Number of simulated true locations. 
#' @param start_seed Starting seed
#' @param movement_process Type of process model to simulate from. "DCRWS" implies the first-difference correlated random walk, while "carHMM" is the conditional autoregressive correlated random walk. 
#' @param process_pars List of relevant parameters depending on the movement process. 
#' @param alpha Matrix of switching probabilities. Determines the number of states you wish to have. 
#' @param start_loc Starting location in c(lon, lat). If NULL assumed to be (0,0).  
#' @return Simulated data set.
genTrack <- function(n_x=500, 
                     start_seed=42,
                     movement_process="DCRWS",
                     process_pars = list(theta = c(0, pi),
                                         gamma = c(0.8,0.3),
                                         sdx = rep(0.1, 2)), 
                     alpha = matrix(c(0.95, 0.05, 0.05, 0.95), nrow=2, ncol=2),
                     start_loc=NULL){

  
  set.seed(start_seed)
  
  
  ##########################
  ######  behaviours  ######
  ##########################
  
  m=nrow(alpha)
  B <-  1:m
  b <- numeric(n_x)
  b[1] <- sample(B, size=1, prob=rep(1/m, m))
  for(i in 2:n_x) b[i] = sample(B, size=1, prob=alpha[b[i-1],])
  
  
  
  ##############################
  ######  true locations  ######
  ##############################
  
  if(movement_process=="DCRWS"){
    
    X <- matrix(nrow=n_x, ncol=2)
    if(is.null(start_loc)){
      X[1,] <- rep(0,2)
    } else {
      X[1,] <- start_loc
    }
    X[2,] = mvtnorm::rmvnorm(1, X[1,], diag(process_pars$sdx, nrow=2)) 
    for(i in 3:n_x){
      Tr <- matrix(c(cos(process_pars$theta[b[i]]), -sin(process_pars$theta[b[i]]),
                     sin(process_pars$theta[b[i]]), cos(process_pars$theta[b[i]])), 
                   nrow=2, byrow=TRUE)
      X[i,] <- mvtnorm::rmvnorm(1, X[i-1,] + process_pars$gamma[b[i]] * Tr %*% (X[i-1,] - X[i-2,]), diag(process_pars$sdx)) 
    } 
    
  } else if (movement_process=="carHMM") {
    
    # step lengths is gamma distribution
    sl <- numeric(n_x)
    mu <- process_pars$mrl[b[i]]
    sigma  <- process_pars$sigma[b[i]]
    sl[2] <- rgamma(1, (mu/sigma)^2, scale=sigma^2/mu)
    for(i in 3:n_x){
      mu <- process_pars$acf[b[i]]*sl[i-1] + (1-process_pars$acf[b[i]])*process_pars$mrl[b[i]]
      sigma  <- process_pars$sigma[b[i]]
      sl[i] <- rgamma(1, shape=(mu/sigma)^2, scale=sigma^2/mu)
    }
    
    # turning angles are wrapped cauchy
    # parameters should be a list of vectors of the form list(location=c(pi,0,0), rho=c(0.2, 0.6, 0.8))
    ta <- numeric(n_x)
    for(i in 2:(n_x-1)) ta[i] <- CircStats::rwrpcauchy(1, location=process_pars$location[b[i]], rho=process_pars$rho[b[i]])
    # ta[i] <- ta[i] - (2*pi)*floor((ta[i]+pi)/(2*pi))
    
    
    X.geo <- data.frame(ta=ta, sl=sl)
    
    # turn the turning angles and step lengths into locations
    # assume the first location is (0,0), and the first bearing is 0 degrees
    
    X  <- getPath(angles=ta*180/pi, type='turning', step_lengths=sl, 
                  first_loc=start_loc, first_bearing=NULL)#*180/pi
    
  } else {
    
    stop("model not implemented...")
    
  }
  

  
  if(movement_process=="DCRWS"){
    rslt <- list(b=b, X=X)
  } else { 
    rslt <- list(b=b, X=X, X.geo=X.geo)
  }
  
  return(rslt)
  
}