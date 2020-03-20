#prints the output from the summary of the swim object

print.summary.swim <- function(x){
  
  cat("Hidden Markov Movement Model:\n")
  
  cat("\n Time to fit: \n")
  
  print(x[[1]])
  
  cat("\n Confidence Intervals: \n")
  print(round(x[[2]][9:18,-2], 7))
  
#   cat("\nWald Tests: \n")
#   
#   printCoefmat(x[[2]][9:18,], P.values=TRUE, has.Pvalue=TRUE)

}
