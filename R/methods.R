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
plot.issm <- function(mod, sim=NULL){
  
  fmf <- wesanderson::wes_palette("FantasticFox1", type = "discrete")
  
  if(!is.null(sim)){
    lplot <- ggplot() + 
      ylab("Northing") +
      xlab("Easting") +
      geom_path(data=as.data.frame(sim$locs), aes(x=coords.x1, y=coords.x2), col='navy', lwd=1.1) +
      geom_path(data=as.data.frame(mod$ssm_results[[mod$winner]]$l_hat), aes(x=x, y=y), col='cyan3', lwd=1.1) +
      ggtitle("True and Predicted Path") +
      theme_bw()
  }
  
  # plot the predicted locations and their predicted behavioural states
  bplot <- ggplot() + 
    ylab("Northing") +
    xlab("Easting") +
    geom_path(data=mod$obs, aes(x=lon, y=lat), col='grey') +
    geom_path(data=as.data.frame(mod$ssm_results[[mod$winner]]$xhat), aes(x=lon, y=lat), col='black') +
    geom_point(data=as.data.frame(mod$ssm_results[[mod$winner]]$xhat), aes(x=lon, y=lat, col=as.factor(mod$hmm_results[[mod$winner]]$b_hat))) +
    scale_color_manual(values=c(fmf[2], fmf[3])) +
    guides(colour=guide_legend(title="b state")) +
    ggtitle("Predicted Path and Predicted States") +
    theme_bw()
  
  # take the residuals, gotta get rid of the rows with Inf
  # get rid of the first one, CHECK THIS LATER
  pseudo <- as.data.frame(mod$hmm_results[[mod$winner]]$obj$report()$pseudo)[-1,]
  names(pseudo) <- c("lon", "lat")
  pseudo <- pseudo[is.finite(rowSums(pseudo)),]
  
  # plot the residuals
  linetmp <- get_qqline(pseudo[,1]) # get qqline parameters
  xresidplot <- ggplot(data=pseudo, aes(sample=lon)) +
    geom_qq() +
    geom_abline(slope=linetmp[1], intercept=linetmp[2], col=fmf[4]) + 
    ggtitle("X Coordinate Residuals") +
    theme_bw()
  
  linetmp <- get_qqline(pseudo[,2])
  yresidplot <- ggplot(data=pseudo, aes(sample=lat)) +
    geom_qq() +
    geom_abline(slope=linetmp[1], intercept=linetmp[2], col=fmf[4]) +
    ggtitle("Y Coordinate Residuals") +
    theme_bw()
  
  # put it all together
  if(is.null(sim)){
    grid.arrange(bplot,
                 xresidplot, yresidplot,
                 widths=c(2,1),
                 layout_matrix=rbind(c(1,2),
                                     c(1,3)),
                 top = paste("winning iteration =", mod$winner))
  } else {
    grid.arrange(lplot, bplot,
                 xresidplot, yresidplot,
                 widths=c(2,1),
                 layout_matrix=rbind(c(1,3),
                                     c(2,4)),
                 top = paste("winning iteration =", mod$winner))
  }
  
}

# retrieve just the qqplots
plotqq <- function(mod){
  
  fmf <- wesanderson::wes_palette("FantasticFox1", type = "discrete")
  
  pseudo <- as.data.frame(mod$hmm_results[[mod$winner]]$obj$report()$pseudo)[-1,]
  names(pseudo) <- c("lon", "lat")
  pseudo <- pseudo[is.finite(rowSums(pseudo)),]
  
  linetmp <- get_qqline(pseudo[,1]) # get qqline parameters
  xresidplot <- ggplot(data=pseudo, aes(sample=lon)) +
    geom_qq() +
    geom_abline(slope=linetmp[1], intercept=linetmp[2], col=fmf[4]) + 
    # ggtitle("X Coordinate Residuals") +
    theme_bw() #+ theme(axis.title.x=element_blank(), axis.title.y=element_blank()) 
  
  linetmp <- get_qqline(pseudo[,2])
  yresidplot <- ggplot(data=pseudo, aes(sample=lat)) +
    geom_qq() +
    geom_abline(slope=linetmp[1], intercept=linetmp[2], col=fmf[4]) +
    # ggtitle("Y Coordinate Residuals") +
    theme_bw() #+ theme(axis.title.x=element_blank(), axis.title.y=element_blank()) 
  
  rslt <- list(xqqplot=xresidplot, yqqplot=yresidplot)
  return(rslt)
}
