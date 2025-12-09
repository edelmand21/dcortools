#' Performs distance correlation sure independence screening \insertCite{li2012feature}{dcortools} with some additional options (such as calculating corresponding tests).
#'
#' @param X A dataframe or matrix.
#' @param Y A vector-valued response having the same length as the number of rows of X.
#' @param k Number of variables that are selected (only used when threshold is not provided).
#' @param threshold If provided, variables with a distance correlation larger than threshold are selected.
#' @param calc.cor If set as "pearson", "spearman" or "kendall", a corresponding correlation matrix is additionally calculated.
#' @param calc.pvalue.cor logical; IF TRUE, a p-value based on the Pearson or Spearman correlation matrix is calculated (not implemented for calc.cor = "kendall") using Hmisc::rcorr.
#' @param return.data logical; specifies if the dcmatrix object should contain the original data.
#' @param test Allows for additionally calculating a test based on distance Covariance. Specifies the type of test that is performed, "permutation" performs a Monte Carlo Permutation test. "gamma" performs a test based on a gamma approximation of the test statistic under the null. "conservative" performs a conservative two-moment approximation. "bb3" performs a quite precise three-moment approximation and is recommended when computation time is not an issue.
#' @param adjustp If setting this parameter to "holm", "hochberg", "hommel", "bonferroni", "BH", "BY" or "fdr", corresponding adjusted p-values are additionally returned for the distance covariance test.
#' @param b specifies the number of random permutations used for the permutation test. Ignored for all other tests.
#' @param bias.corr logical; specifies if the bias corrected version of the sample distance covariance \insertCite{huo2016fast}{dcortools} should be calculated.
#' @param use "all" uses all observations, "complete.obs" excludes NAs, "pairwise.complete.obs" uses pairwise complete observations for each comparison.
#' @param algorithm specifies the algorithm used for calculating the distance covariance. 
#' 
#' "fast" uses an O(n log n) algorithm if the observations are one-dimensional and metr.X and metr.Y are either "euclidean" or "discrete", see also \insertCite{huo2016fast;textual}{dcortools}. 
#' 
#' "memsave" uses a memory saving version of the standard algorithm with computational complexity O(n^2) but requiring only O(n) memory. 
#' 
#' "standard" uses the classical algorithm. User-specified metrics always use the classical algorithm.
#' 
#' "auto" chooses the best algorithm for the specific setting using a rule of thumb.
#' 
#' "memsave" is typically very inefficient for dcsis and should only be applied in exceptional cases.
#' @return dcmatrix object with the following two additional slots:
#' \item{name selected}{description indices of selected variables.} 
#' \item{name dcor.selected}{distance correlation of the selected variables and the response Y.}
#' @export
#' @references
#' \insertRef{berschneider2018complex}{dcortools}

#' \insertRef{dueck2014affinely}{dcortools}
#' 
#' \insertRef{huang2017statistically}{dcortools}
#' 
#' \insertRef{huo2016fast}{dcortools}
#' 
#' \insertRef{li2012feature}{dcortools}
#' 
#' \insertRef{szekely2007}{dcortools}
#' 
#' \insertRef{szekely2009brownian}{dcortools}
#' @examples
#' X <- matrix(rnorm(1e5), ncol = 1000)
#' Y <- sapply(1:100, function(u) sum(X[u, 1:50])) + rnorm(100)
#' a <- dcsis(X, Y)

dcsis <- function(X, 
                  Y,
                  k = floor(nrow(X)/log(nrow(X))),
                  threshold = NULL,
                  calc.cor = "spearman",
                  calc.pvalue.cor = FALSE,
                  return.data = FALSE,
                  test = "none",
                  adjustp = "none",
                  b = 499,
                  bias.corr = FALSE,
                  use="all",
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
                     calc.pvalue.cor = calc.pvalue.cor,
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
