#' Performs Distance correlation Sure Independence Screening (DC-SIS)
#'
#' @param X A dataframe or matrix with n rows and p columns.
#' @param Y A vector-valued response of length n
#' @param k Number of variables that are selected (only used when threshold is not provided).
#' @param threshold If provided, variables with a distance correlation larger than threshold are selected.
#' @param calc.cor If set as "pearson", "spearman" or "kendall", a corresponding correlation matrix is addionally calculated.
#' @param calc.pval.cor IF TRUE, a p-value based on the Pearson or Spearman correlation matrix is calculated (not implemented for calc.cor ="kendall") using Hmisc::rcorr
#' @param return.data IF TRUE, X and Y are contained in the resulting dcmatrix object.
#' @param test specifies the type of test that is performed, "permutation" performs a Monte Carlo Permutation test. "gamma" performs a test based on a gamma approximation of the test statistic under the null. "conservative" performs a conservative two-moment approximation. "bb3" performs a quite precise three-moment approximation and is recommended when computation time is not an issue.
#' @param adjustp If setting this parameter to "holm", "hochberg", "hommel", "bonferroni", "BH", "BY" or "fdr", corresponding adjusted p-values are returned.
#' @param b specifies the number of random permutations used for the permutation test. Ignored when test="gamma"
#' @param affine logical; indicates if the affinely transformed distance covariance should be calculated or not.
#' @param bias.corr logical; indicates if the bias corrected version of the sample distance covariance should be calculated
#' @param use : "all" uses all observations, "complete.obs" excludes NA's, "pairwise.complete.obs" uses pairwise complete observations for each comparison.
#' @param algorithm: One of "auto", "fast", "memsave" and "standard". "memsave" is typically very inefficient for dcmatrix and should only be applied in exceptional cases.
#' @return dcmatrix object
#' @export


dcsis <- function(X, 
                  Y,
                  k = floor(nrow(X)/log(nrow(X))),
                  threshold = NULL,
                  calc.cor = "spearman",
                  calc.pval.cor = FALSE,
                  return.data = FALSE,
                  test = "none",
                  adjustp = "none",
                  b = 499,
                  bias.corr = FALSE,
                  use="everything",
                  algorithm ="auto") {
  
  p <- ncol(X)
  q <- ncol(Y)
  
  if (is.null(q))
    q <- 1
  
  if (k  < 0) {
    stop("k must be positive")
  }
  
  if (q != 1) 
    stop("Dimension of response must be 1 for screening algorithm!")
  
  output <- dcmatrix(X = X,
                     Y = Y,
                     calc.cor = calc.cor,
                     calc.pval.cor = calc.pval.cor,
                     return.data = return.data,
                     test = test,
                     adjustp = adjustp,
                     b = b,
                     bias.corr = bias.corr,
                     use = use,
                     algorithm = algorithm)
  
  if (!is.null(threshold)) {
    output$selected <- which(calc.cor > threshold)
    output$dcor.selected <- output$dcor[output$selected]
  } else {
    k <- min(round(k),p)
    
    output$selected <- order(output$dcor, decreasing = TRUE)[1:k]
    output$dcor.selected <- output$dcor[output$selected]
  }
  
  return(output)  

  
}
