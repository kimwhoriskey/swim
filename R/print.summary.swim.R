#prints the output from the summary of the swim object

print.summary.swim <- function(sum, ...){
  
  cat("Switching Hidden Markov Movement Model:\n")
  
  cat("\n Time to fit: \n")
  
  print(sum[[1]])
  
  cat("\n Confidence Intervals: \n")
  print(round(sum[[2]][9:18,-2], 7))
  
#   cat("\nWald Tests: \n")
#   
#   printCoefmat(sum[[2]][9:18,], P.values=TRUE, has.Pvalue=TRUE)

}
