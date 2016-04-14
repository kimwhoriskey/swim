#Get the Wald tests and the Confidence intervals for the parameter values

summary.swim <- function(object, ...){
  
  time = object$time
  
  #object = test
  
  #Summary (Z statistics)
  Ztest = as.data.frame(object$parms, 
                       row.names= c("logitTheta1", "logitTheta2", 
                                    "logitGamma1", "logitGamma2", 
                                    "logSdlon", "logSdlat",  
                                    "logA11", "logA21", 
                                    "theta1", "theta2", 
                                    "gamma1", "gamma2", 
                                    "sdLon", "sdLat",  
                                    "a11", "a21", "a12", "a22"))
  Ztest$z.val = Ztest[,1]/Ztest[,2]
  Ztest$p.val = 2*pnorm(-abs(Ztest$z.val))
  
  #Confidence Intervals
  CIs = as.data.frame(object$parms, 
                               row.names= c("logitTheta1", "logitTheta2", 
                                            "logitGamma1", "logitGamma2", 
                                            "logSdlon", "logSdlat", 
                                            "logA11", "logA21", 
                                            "theta1", "theta2", 
                                            "gamma1", "gamma2", 
                                            "sdLon", "sdLat",  
                                            "a11", "a21", "a12", "a22"))
  
  lower95 = CIs[1:8,1] - 1.96*CIs[1:8,2]
  upper95 = CIs[1:8,1] + 1.96*CIs[1:8,2]
  lower95 = append(lower95, 2*pi/(1.0+exp(-lower95[1])) - pi) #theta1
  upper95 = append(upper95, 2*pi/(1.0+exp(-upper95[1])) - pi)
  lower95 = append(lower95, 2*pi/(1.0+exp(-lower95[2]))) #theta2
  upper95 = append(upper95, 2*pi/(1.0+exp(-upper95[2])))
  lower95 = append(lower95, 1/(1.0+exp(-lower95[3]))) #gamma1
  upper95 = append(upper95, 1/(1.0+exp(-upper95[3])))
  lower95 = append(lower95, 1/(1.0+exp(-lower95[4]))) #gamma2
  upper95 = append(upper95, 1/(1.0+exp(-upper95[4])))
  lower95 = append(lower95, exp(lower95[5])) #sigmalon
  upper95 = append(upper95, exp(upper95[5])) 
  lower95 = append(lower95, exp(lower95[6])) #sigmalat
  upper95 = append(upper95, exp(upper95[6])) 
  lower95 = append(lower95, 1/(1+exp(-lower95[7]))) #a11, logit
  upper95 = append(upper95, 1/(1+exp(-upper95[7]))) 
  lower95 = append(lower95, 1/(1+exp(-lower95[8]))) #a21, note have to logit these ones but I'm not sure if that's correct
  upper95 = append(upper95, 1/(1+exp(-upper95[8]))) 
  
  lower95 = append(lower95, 1-upper95[15])
  upper95 = append(upper95, 1-lower95[15])
  
  lower95 = append(lower95, 1-upper95[16])
  upper95 = append(upper95, 1-lower95[16])

  CIs$lower95 = lower95
  CIs$upper95 = upper95
  
  
  CIs
  
  sum <- list(time, Ztest, CIs)
  names(sum) <- c("Time", "Ztest", "CIs")
  class(sum) <- "summary.swim"
  
  sum
}