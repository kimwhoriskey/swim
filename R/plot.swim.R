#plot the swim fits

# this function gets the slope and intercept of a qqline for the residuals below
get_qqline <- function(vec){
  y <- quantile(vec[!is.na(vec)], c(0.25, 0.75))
  x <- qnorm(c(0.25, 0.75))
  m <- diff(y)/diff(x)
  b <- y[1] - m * x[1]
  return(c(m,b))
}

plot.swim <- function(x,...){
  
  # locations
  bplot <- ggplot(data=x$regData, aes(x=lon, y=lat)) + 
    geom_path() +
    geom_point(aes(x=lon, y=lat, col=as.factor(x$states)))+
    guides(colour=guide_legend(title="b state")) +
    ggtitle("Interpolated Path and Predicted States")
  
  # take the residuals, get rid of the rows with Inf
  if(any(is.finite(rowSums(x$pseudos)))) warning("Inf occurred at least once in the pseudoresiduals - all instances were removed.")
  pseudos <- x$pseudos[is.finite(rowSums(x$pseudos)),]
  
  # plot the residuals
  linetmp <- get_qqline(pseudos[,1]) # get qqline parameters
  xresidplot <- ggplot(data=pseudos, aes(sample=lon)) +
    geom_qq() +
    geom_abline(slope=linetmp[1], intercept=linetmp[2]) + 
    ggtitle("X Coordinate Residuals")
  
  linetmp <- get_qqline(pseudos[,2])
  yresidplot <- ggplot2::ggplot(data=pseudos, aes(sample=lat)) +
    geom_qq() +
    geom_abline(slope=linetmp[1], intercept=linetmp[2]) +
    ggtitle("Y Coordinate Residuals")
  
  # put it all together
  gridExtra::grid.arrange(bplot,
                 xresidplot, yresidplot,
                 widths=c(2,1),
                 layout_matrix=rbind(c(1,2),
                                     c(1,3)))
  

}
