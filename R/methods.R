#' @import ggplot2

############################################
############ set a print method ############
############################################

# self explanatory
print.issm <- function(mod){
  
  cat("----------------------------------------------------------- \n the best iteration is ... \n\n")
  print(mod$winner)
  cat("\n\n\n")
  
  cat("----------------------------------------------------------- \n the associated HMM parameters are ... \n\n")
  print(round(mod$hmm_results[[mod$winner]]$params,4))
  cat("\n\n\n")
  
  cat("----------------------------------------------------------- \n the associated SSM parameters are ... \n\n")
  print(round(mod$ssm_results[[mod$winner]]$params,4))
  cat("\n\n\n")
  
  cat("----------------------------------------------------------- \n the convergence information for the SSM is ... \n\n")
  print(mod$ssm_nll)
  cat("\n\n\n")
  
  
  
}





############################################
############ set a plot method #############
############################################

# this function gets the slope and intercept of a qqline for the residuals below
get_qqline <- function(vec){
  y <- quantile(vec[!is.na(vec)], c(0.25, 0.75))
  x <- qnorm(c(0.25, 0.75))
  m <- diff(y)/diff(x)
  b <- y[1] - m * x[1]
  return(c(m,b))
}


# plot the predicted states along with the pseudoresiduals
plot.issm <- function(mod, sim=NULL, main="seal"){
  
  fmf <- wesanderson::wes_palette("FantasticFox1", type = "discrete")
  
  if(!is.null(sim)){
    lplot <- ggplot() + 
      ylab("Northing") +
      xlab("Easting") +
      geom_path(data=sim$obs, aes(x=lon, y=lat), col='navy') +
      geom_path(data=mod$preds, aes(x=lon, y=lat), col='cyan3') +
      ggtitle("True and Predicted Path") +
      theme_bw()
  }
  
  # plot the predicted locations and their predicted behavioural states
  if(mod$isprojected){
    obs <- as.data.frame(st_coordinates(st_as_sf(mod$obs)))
    names(obs) <- c("lon", "lat")
  } else {
    obs <- mod$obs
  }
  bplot <- ggplot() + 
    ylab("Northing") +
    xlab("Easting") +
    geom_path(data=obs, aes(x=lon, y=lat), col='grey') +
    geom_path(data=mod$preds, aes(x=lon, y=lat), col='black') +
    geom_point(data=mod$preds, aes(x=lon, y=lat, col=as.factor(bhat))) +
    scale_color_manual(values=c(fmf[2], fmf[3])) +
    guides(colour=guide_legend(title="b state")) +
    ggtitle("Predicted Path and Predicted States") +
    theme_bw()
  
  # take the residuals, gotta get rid of the rows with Inf
  qqplots <- plotqq(mod)
  
  # put it all together
  if(is.null(sim)){
    gridExtra::grid.arrange(bplot,
                 qqplots$xqqplot, qqplots$yqqplot,
                 widths=c(2,1),
                 layout_matrix=rbind(c(1,2),
                                     c(1,3)),
                 top = paste(main, "; winning iteration =", mod$winner))
  } else {
    gridExtra::grid.arrange(lplot, bplot,
                 qqplots$xqqplot, qqplots$yqqplot,
                 widths=c(2,1),
                 layout_matrix=rbind(c(1,3),
                                     c(2,4)),
                 top = paste(main, "; winning iteration =", mod$winner))
  }
  
}

# retrieve just the qqplots
plotqq <- function(mod){
  
  fmf <- wesanderson::wes_palette("FantasticFox1", type = "discrete")
  
  pseudo <- as.data.frame(do.call(rbind, mod$hmm_results[[mod$winner]]$pseudos))
  names(pseudo) <- c("lon", "lat")
  if(length(which(is.infinite(rowSums(pseudo))))>0){
    warning("infinite pseudo values detected at ", 
            paste(which(is.infinite(rowSums(pseudo))), collapse=","),
            " and were removed")
    pseudo <- pseudo[is.finite(rowSums(pseudo)),]
  } 
  
  linetmp <- get_qqline(pseudo[,1]) # get qqline parameters
  xresidplot <- ggplot(data=pseudo, aes(sample=lon)) +
    geom_qq() +
    geom_abline(slope=linetmp[1], intercept=linetmp[2], col=fmf[4]) + 
    ggtitle("X Coordinate Residuals") +
    theme_bw() #+ theme(axis.title.x=element_blank(), axis.title.y=element_blank()) 

  linetmp <- get_qqline(pseudo[,2])
  yresidplot <- ggplot(data=pseudo, aes(sample=lat)) +
    geom_qq() +
    geom_abline(slope=linetmp[1], intercept=linetmp[2], col=fmf[4]) +
    ggtitle("Y Coordinate Residuals") +
    theme_bw() #+ theme(axis.title.x=element_blank(), axis.title.y=element_blank()) 
  
  rslt <- list(xqqplot=xresidplot, yqqplot=yresidplot)
  return(rslt)
}



# retrieve just the qqplots
checkqq <- function(mod){
  
  fmf <- wesanderson::wes_palette("FantasticFox1", type = "discrete")
  
  pseudo <- as.data.frame(do.call(rbind, mod$hmm_results[[mod$winner]]$pseudos))
  names(pseudo) <- c("lon", "lat")
  if(length(which(is.infinite(rowSums(pseudo))))>0){
    warning("infinite pseudo values detected at ", 
            paste(which(is.infinite(rowSums(pseudo))), collapse=","),
            " and were removed")
    pseudo <- pseudo[is.finite(rowSums(pseudo)),]
  } 
  
  linetmp <- get_qqline(pseudo[,1]) # get qqline parameters
  ecdf(pseudo[,1]) %>% str
  ks.test(pseeudo[,1], pnorm)
  ecdf <- data.frame(y=sort(unique(pseudo[,1])), ecdf=knots(ecdf(pseudo[,1])))
  lonqqpred <- linetmp[2] + linetmp[1]*   
  
  linetmp <- get_qqline(pseudo[,2])
  latqqpred
  
  rslt <- list(xqqplot=xresidplot, yqqplot=yresidplot)
  return(rslt)
}




##########################################################
############ set a plot method for bootstrap #############
##########################################################

plot.swim_bootstrap <- function(mod, paramvec=c('theta1', 'theta2', 'gamma1', 'gamma2', 'sigmalon', 'sigmalat', 'a11', 'a21', 'a12', 'a22', 'psi')){
  t(mod$pars) %>% 
    as.data.frame() %>% 
    rename(theta1=theta, theta2=theta.1, gamma1=gamma, gamma2=gamma.1, 
           sigmalon=X, sigmalat=X.1, a11=A, a21=A.1, a12=A.2, a22=A.3, psi=X.2) %>% 
    pivot_longer(cols=c(theta1, theta2, gamma1, gamma2, a11, a21, a12, a22, sigmalon, sigmalat, psi), values_to='val') %>% 
    ggplot() + 
    geom_violin(aes(x=name, y=val), fill=fmf[2]) + 
    facet_wrap(~name, scales='free') +
    theme_bw() + 
    geom_hline(data=data.frame("name"=paramvec,
                               "val"=mod$true.pars), aes(yintercept=val), col=fmf[5])
  
}



###########################################################
############ set a plot method for simulation #############
###########################################################


plot.swimsim <- function(sim){
  
  fmf <- wesanderson::wes_palette("FantasticFox1", type = "discrete")

  devAskNewPage(TRUE)
  
  print(
  sim$X %>% ggplot + 
    geom_path(aes(x=lon, y=lat), col='black')+
    geom_point(aes(x=lon, y=lat, col=factor(sim$b))) + 
    theme_bw()+
    scale_color_manual(values=c(fmf[2], fmf[5]))+
    coord_fixed()
  )  
  
  print(
  sim$X %>% ggplot + 
    geom_path(data=sim$X, aes(x=lon, y=lat, col=sim$trackidx))+
    theme_bw()+
    coord_fixed()
  )

  print(
  sim$Y %>% ggplot + geom_path(aes(x=lon, y=lat), col='grey')+
    geom_path(data=sim$X, aes(x=lon, y=lat), col='black')+
    theme_bw()+
    coord_fixed()
  )
  
  devAskNewPage(FALSE)
  
  
}



