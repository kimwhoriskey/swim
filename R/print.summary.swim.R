#prints the output from the summary of the swim object

print.summary.swim <- function(object, ...){
  
  cat("Switching Hidden Markov Movement Model:\n")
  cat("\nWald Tests: \n")
  
  printCoefmat(object[[1]][9:18,], P.values=TRUE, has.Pvalue=TRUE)
  
  cat("\n Confidence Intervals: \n")
  print(round(summary(test)[[2]][9:18,-2], 7))

}
