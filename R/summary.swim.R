#Get the Wald tests and the Confidence intervals for the parameter values

summary.swim <- function(object, ...){
  
  object = test
  
  #Summary (Z statistics)
  Ztest = as.data.frame(object$parms, 
                       row.names= c("logitTheta1", "logitTheta2", 
                                    "logitGamma1", "logitGamma2", 
                                    "logSdlat", "logSdlon", 
                                    "logA11", "logA21", 
                                    "theta1", "theta2", 
                                    "gamma1", "gamma2", 
                                    "sdLat", "sdLon", 
                                    "a11", "a21", "a12", "a22"))
  Ztest$z.val = Ztest[,1]/Ztest[,2]
  Ztest$p.val = 2*pnorm(-abs(Ztest$z.val))
  
  #Confidence Intervals
  CIs = as.data.frame(object$parms, 
                               row.names= c("logitTheta1", "logitTheta2", 
                                            "logitGamma1", "logitGamma2", 
                                            "logSdlat", "logSdlon", 
                                            "logA11", "logA21", 
                                            "theta1", "theta2", 
                                            "gamma1", "gamma2", 
                                            "sdLat", "sdLon", 
                                            "a11", "a21", "a12", "a22"))
  CIs$lower95 = CIs[,1] - 2*CIs[,2]
  CIs$upper95 = CIs[,1] + 2*CIs[,2]
  CIs
  
  sum <- list(Ztest, CIs)
  names(sum) <- c("Ztest", "CIs")
  class(sum) <- "summary.swim"
  
  sum
}