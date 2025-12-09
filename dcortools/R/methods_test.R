#' @export
print.dctest <- function(x, ...) {
  dctest <- x
  if (dctest$method == "none")
    testout <- "None"
  else if (dctest$method == "permutation")
    testout <- paste("Permutation test with ",dctest$b," permutations")
  else if (dctest$method == "gamma")
    testout <- "Simple gamma approximation"               
  else if (dctest$method == "conservative")
    testout <- "Conservative two-moment approximation"  
  else if (dctest$method == "bb3")
    testout <- "Three moment approximation by Berschneider and Boettcher"
  else if (dctest$method == "wildbs1")
    testout <- "Wild bootstrap by Chwialkowksi, et al., Method 1"
  else if (dctest$method == "wildbs2")
    testout <- "Wild bootstrap by Chwialkowksi, et al., Method 2"
  
  #if (dctest$pvalue<1e-6)
 #   pvalout <- "< 1E-6"
 # else
#pvalout <- round(dctest$pvalue,6)
  
  pvalout <- dctest$pvalue
  
  cat(paste("pvalue: ",pvalout,"   dcov: ",round(dctest$dcov,6),"   dcor: ",round(dctest$dcor,6)))
  cat("\n")
  cat("\n")
  cat(paste("Method:",testout, "\n") )
  cat("\n")
  cat("Call:", "\n")
  cat(paste(deparse(dctest$call), "\n"))
  cat("\n")
  cat(paste("Bias correction: ", ifelse(dctest$bias.corr,"Yes \n", "No \n")) )
  cat(paste("Affinely invariant: ", ifelse(dctest$affine,"Yes \n", "No \n")) )
  cat(paste("Standardization applied: ", ifelse(dctest$standardize,"Yes \n", "No \n")) )
  cat(paste("Metric for X: ", dctest$metr.X, "\n") )
  cat(paste("Metric for Y: ", dctest$metr.Y, "\n") )
}


is.dctest <- function(x) inherits(x, "dctest")