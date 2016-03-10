#the simulation function
simTrack = function(T=500, theta=c(0,pi), gamma=c(0.7,0.2), 
                    Sigma=matrix(c(0.2*0.2,0,0,0.1*0.1),nrow=2), 
                    alphas=c(0.8,0.9)){

  #accumulators
  x = matrix(NA, T, 2)
  emp = matrix(NA, T-1, 2)
  b = c()
  alpha=c(1-alphas[1], alphas[2])
  
  #Initialize locations randomly
  x[1, ] = mvtnorm::rmvnorm(1, c(0,0), Sigma) 
  x[2, ] = mvtnorm::rmvnorm(1, x[1,], Sigma)
  
  #Initialize behavioral state
  b[1] = rbinom(1, 1, 0.5) + 1
  
  
  #Process equation
  for(i in 2:(T-1)){
    #evolve states
    b[i] = rbinom(1, 1, alpha[b[i-1]]) + 1
    #the idea is that a success is defined as being in state 2, and a failure being
    #in state 1. So we want the first alpha probability here to be low, since it 
    #is the probability of not switching to state 2 given you were in state 1. 
    
    #evolve movement process
    emp[i,1] = cos(theta[b[i]]) * (x[i,1] - x[i-1,1]) + sin(theta[b[i]]) * (x[i,2] - x[i-1,2])
    emp[i,2] = -sin(theta[b[i]]) * (x[i,1] - x[i-1,1]) + cos(theta[b[i]]) * (x[i,2] - x[i-1,2])
    
    #add random error
    x[i+1,] = mvtnorm::rmvnorm(1, x[i,] + emp[i,] * gamma[b[i]], Sigma)
  }
  
  #Get the last behavioral state
  b[T] = rbinom(1, 1, alpha[b[T-1]]) + 1
  
  dat = data.frame(lon=x[,1], lat = x[,2], b)
  
  return(dat)
  
}