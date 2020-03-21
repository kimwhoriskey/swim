#' @export
print.swim <- function(object){
  
  cat("----------------------------------------------------------- \n the convergence was ... \n\n")
  print(object$opt$message)
  cat("\n\n\n")
  
  cat("----------------------------------------------------------- \n the associated movement parameters and their standard errors are ... \n\n")
  print(round(object$parameters[9:14, ], 4))
  cat("\n\n\n")
  
  cat("----------------------------------------------------------- \n the associated TPM is ... \n\n")
  print(matrix(round(object$parameters[15:18, 1], 4), nrow=2, byrow=FALSE))
  cat("\n\n\n")
  
  cat("----------------------------------------------------------- \n the stationary distribution is ... \n\n")
  print(object$stationary)
  cat("\n\n\n")
  
}
