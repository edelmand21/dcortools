#' @useDynLib dcortools, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom Rfast standardise
NULL

#' @importFrom Rdpack reprompt
NULL

#' Calculates the distance covariance \insertCite{szekely2007,szekely2009brownian}{dcortools}.
#'
#' @param X contains either the first  sample or its corresponding distance matrix.
#'
#' In the first case, X can be provided either as a vector (if one-dimensional), a matrix or a data.frame (if two-dimensional or higher). 
#'
#' In the second case, the input must be a distance matrix corresponding to the sample of interest.
#'
#' If X is a sample, type.X must be specified as "sample". If X is a distance matrix, type.X must be specified as "distance".
#' @param Y see X.
#' @param affine logical; specifies if the affinely invariant distance covariance \insertCite{dueck2014affinely}{dcortools} should be calculated or not.
#' @param standardize logical; specifies if X and Y should be standardized dividing each component by its standard deviations. No effect when affine = TRUE.
#' @param bias.corr logical; specifies if the bias corrected version of the sample distance covariance \insertCite{huo2016fast}{dcortools} should be calculated.
#' @param type.X For "distance", X is interpreted as a distance matrix. For "sample", X is interpreted as a sample.
#' @param type.Y see type.X.
#' @param metr.X specifies the metric which should be used to compute the distance matrix for X (ignored when type.X = "distance").
#' 
#'  Options are "euclidean", "discrete", "alpha", "minkowski", "gaussian", "gaussauto", "boundsq" or user-specified metrics (see examples).
#'  
#'  For "alpha", "minkowski", "gaussian", "gaussauto" and "boundsq", the corresponding parameters are specified via "c(metric, parameter)", e.g. c("gaussian", 3) for a Gaussian metric with bandwidth parameter 3; the default parameter is 2 for "minkowski" and "1" for all other metrics.
#'  
#'  See \insertCite{lyons2013distance,sejdinovic2013equivalence,bottcher2017detecting;textual}{dcortools} for details.
#' @param metr.Y see metr.X.
#' @param use specifies how to treat missing values. "complete.obs" excludes observations containing NAs, "all" uses all observations.
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
#' @return numeric; the distance covariance between samples X and Y.
#' @references
#' \insertRef{bottcher2017detecting}{dcortools}
#' 
#' \insertRef{dueck2014affinely}{dcortools}
#' 
#' \insertRef{huo2016fast}{dcortools}
#' 
#' \insertRef{lyons2013distance}{dcortools}
#' 
#' \insertRef{sejdinovic2013equivalence}{dcortools}
#' 
#' \insertRef{szekely2007}{dcortools}
#' 
#' \insertRef{szekely2009brownian}{dcortools}
#' @examples 
#' X <- rnorm(100)
#' Y <- X + 3 * rnorm(100)
#' distcov(X, Y) # standard distance covariance
#' 
#' distcov(X, Y, metr.X = "gaussauto", metr.Y = "gaussauto") 
#' # Gaussian distance with bandwidth choice based on median heuristic
#' 
#' distcov(X, Y, metr.X = c("alpha", 0.5), metr.Y = c("alpha", 0.5)) 
#' # alpha distance covariance with alpha = 0.5.
#' 
#' 
#' #Define a user-specified (slow) version of the alpha metric
#' 
#' alpha_user <- function(X, prm = 1, kernel = FALSE) {
#'     as.matrix(dist(X)) ^ prm
#' }
#' 
#' distcov(X, Y, metr.X = c("alpha", 0.5), metr.Y = c("alpha", 0.5)) 
#' # Gives the same result as before.
#'    
#'
#' #User-specified Gaussian kernel function  
#'      
#' gauss_kernel <- function(X, prm = 1, kernel = TRUE)  {
#'     exp(as.matrix(dist(X)) ^ 2 / 2 / prm ^ 2)
#' }  
#' 
#' distcov(X, Y, metr.X = c("gauss_kernel", 2), metr.Y = c("gauss_kernel", 2)) 
#' # calculates the distance covariance using the corresponding kernel-induced metric
#' 
#' distcov(X, Y, metr.X = c("gaussian", 2), metr.Y = c("gaussian", 2)) 
#' # same result
#' 
#' Y <- matrix(nrow = 100, ncol = 2)
#' X <- rnorm(300)
#' dim(X) <- c(100, 3)
#' Z <- rnorm(100)
#' Y <- matrix(nrow = 100, ncol = 2)
#' Y[, 1] <- X[, 1] + Z
#' Y[, 2] <- 3 * Z
#' 
#' distcov(X, Y) 
#' 
#' distcov(X, Y, affine = TRUE) 
#' # affinely invariant distance covariance
#' 
#' distcov(X, Y, standardize = TRUE) 
#' ## distance covariance standardizing the components of X and Y
#' 
#' @export
distcov <-
  function(X,
           Y,
           affine = FALSE,
           standardize=FALSE,
           bias.corr = FALSE,
           type.X = "sample",
           type.Y = "sample",
           metr.X = "euclidean",
           metr.Y = "euclidean",
           use = "all",
           algorithm = "auto") {
   
    #extract dimensions and sample sizes
    ss.dimX <- extract_np(X, type.X)
    ss.dimY <- extract_np(Y, type.Y)
    
    n <- ss.dimX$Sample.Size
    p <- ss.dimX$Dimension
    
    m <- ss.dimY$Sample.Size
    q <- ss.dimY$Dimension
    
    if (n != m) {
      stop("Samples X and Y must have the same sizes!")
    }
    
    if  (use == "complete.obs") {
      ccX <- ccY <- cc <- 1:n
      ccX <- which(stats::complete.cases(X))
      ccY <- which(stats::complete.cases(Y))
      cc <- intersect(ccX, ccY)
      n <- m <- length(cc)
      if (type.X == "sample" && p == 1) {
        X <- X[cc]
      } else if (type.X == "sample" && p > 1) {
        X <- X[cc, ]
      }
      if (type.Y == "sample" && q == 1) {
        Y <- Y[cc]
      } else if (type.X == "sample" && q > 1) {
        Y <- Y[cc, ]
      }
 
      
      if (type.X == "distance") {
        X <- X[cc,cc]
      }
      if (type.Y == "distance") {
        Y <- Y[cc,cc]
      }
    }
    
    ## normalize samples if calculation of affinely invariant distance covariance is desired
    if (affine) {
      if (type.X == "distance" | type.Y == "distance") {
        stop("Affinely invariant distance covariance cannot be calculated for type distance")
      }
      if (p > n | q > n) {
        stop("Affinely invariant distance covariance cannot be calculated for p>n")
      }
      if (p > 1) {
        X <- X %*% Rfast::spdinv(mroot(stats::var(X)))
      } else {
        X <- X / stats::sd(X)
      }
      if (q > 1) {
        Y <- Y %*% Rfast::spdinv(mroot(stats::var(Y)))
      } else {
        Y <- Y / stats::sd(Y)
      }
    } else if (standardize) {
      if (type.X == "distance" | type.Y == "distance") {
        stop("Standardization cannot be applied for type distance.")
      }
      if (p > 1) {
        X <- standardise(X, center = FALSE)
      } else {
        X <- X / stats::sd(X)
      }
      if (q > 1) {
        Y <- standardise(Y, center = FALSE)
      } else {
        Y <- Y / stats::sd(Y)
      }
    }
    
    
    alg.fast <- alg.standard <- alg.memsave <- FALSE
    
    if (algorithm == "auto") {
      if (p == 1 & q == 1 & metr.X[1] %in% c("euclidean", "discrete") 
          & metr.Y[1] %in% c("euclidean", "discrete") & type.X == "sample" & type.Y == "sample" & n > 100) {
        alg.fast <- TRUE
      } else if (metr.X[1] %in% c("euclidean", "alpha", "gaussian", "boundsq", "minkowski", "discrete") &
                 metr.Y[1] %in% c("euclidean", "alpha", "gaussian", "boundsq", "minkowski", "discrete") & type.X == "sample" & type.Y == "sample") {
        alg.memsave <- TRUE
      } else {
        alg.standard <- TRUE
      }
    } else if (algorithm == "fast") {
      alg.fast <- TRUE 
    } else if (algorithm == "standard") {
      alg.standard <- TRUE
    }  else if (algorithm == "memsave") {
      alg.memsave <- TRUE
    } else
      stop ("Algorithm must be one of \"fast\", \"standard\", \"memsave\" or \"auto\"")
    

    terms <- dcovterms(X, Y, n, calc.dcor = FALSE, doperm = FALSE, dobb3 = FALSE, alg.fast = alg.fast, alg.memsave = alg.memsave, alg.standard = alg.standard, p = p, q = q, metr.X = metr.X, metr.Y =metr.Y, type.X = type.X, type.Y = type.Y)

   
    if (bias.corr) {
      dcov2 <- terms$aijbij / n / (n - 3) - 2 * terms$Sab / n / (n - 2) / (n - 3) + terms$Tab / n / (n - 1) / (n - 2) / (n - 3)
      dcov <- sqrt(abs(dcov2)) * sign(dcov2) 
    } else {
      dcov2 <- terms$aijbij / n / n - 2 * terms$Sab / n / n / n  + terms$Tab / n / n / n / n
      dcov <- sqrt(abs(dcov2))
    }
    return(dcov)
}




#' Calculates the distance standard deviation \insertCite{edelmann2017distance}{dcortools}.
#'
#' @param X contains either the sample or its corresponding distance matrix.
#'
#' In the first case, X can be provided either as a vector (if one-dimensional), a matrix or a data.frame (if two-dimensional or higher). 
#'
#' In the second case, the input must be a distance matrix corresponding to the sample of interest.
#'
#' If X is a sample, type.X must be specified as "sample". If X is a distance matrix, type.X must be specified as "distance".
#' @param affine logical; specifies if the affinely invariant distance standard deviation \insertCite{dueck2014affinely}{dcortools} should be calculated or not.
#' @param standardize logical; specifies if X and Y should be standardized dividing each component by its standard deviations. No effect when affine = TRUE.
#' @param bias.corr logical; specifies if the bias corrected version of the sample distance standard deviation \insertCite{huo2016fast}{dcortools} should be calculated.
#' @param type.X For "distance", X is interpreted as a distance matrix. For "sample", X is interpreted as a sample.
#' @param metr.X specifies the metric which should be used to compute the distance matrix for X (ignored when type.X = "distance").
#' 
#'  Options are "euclidean", "discrete", "alpha", "minkowski", "gaussian", "gaussauto", "boundsq" or user-specified metrics (see examples).
#'  
#'  For "alpha", "minkowski", "gaussian", "gaussauto" and "boundsq", the corresponding parameters are specified via "c(metric, parameter)", e.g. c("gaussian", 3) for a Gaussian metric with bandwidth parameter 3; the default parameter is 2 for "minkowski" and "1" for all other metrics.
#'  
#'  See \insertCite{lyons2013distance,sejdinovic2013equivalence,bottcher2017detecting;textual}{dcortools} for details.
#' @param use specifies how to treat missing values. "complete.obs" excludes observations containing NAs, "all" uses all observations.
#' @param algorithm specifies the algorithm used for calculating the distance standard deviation. 
#' 
#' "fast" uses an O(n log n) algorithm if the observations are one-dimensional and metr.X and metr.Y are either "euclidean" or "discrete", see also \insertCite{huo2016fast;textual}{dcortools}. 
#' 
#' "memsave" uses a memory saving version of the standard algorithm with computational complexity O(n^2) but requiring only O(n) memory. 
#' 
#' "standard" uses the classical algorithm. User-specified metrics always use the classical algorithm.
#' 
#' "auto" chooses the best algorithm for the specific setting using a rule of thumb.
#' 
#' @return numeric; the distance standard deviation of X.
#' @references
#' \insertRef{bottcher2017detecting}{dcortools}
#' 
#' \insertRef{dueck2014affinely}{dcortools}
#' 
#' \insertRef{edelmann2017distance}{dcortools}
#' 
#' \insertRef{huo2016fast}{dcortools}
#' 
#' \insertRef{lyons2013distance}{dcortools}
#' 
#' \insertRef{sejdinovic2013equivalence}{dcortools}
#' 
#' \insertRef{szekely2007}{dcortools}
#' 
#' \insertRef{szekely2009brownian}{dcortools}
#' @export
#' @examples 
#' X <- rnorm(100)
#' distsd(X) # for more examples on the options see the documentation of distcov.
distsd <-
  function(X,
           affine = FALSE,
           standardize=FALSE,
           bias.corr = FALSE,
           type.X = "sample",
           metr.X = "euclidean",
           use = "all",
           algorithm = "auto") {
    #extract dimensions and sample sizes
    ss.dimX <- extract_np(X, type.X)
    
    n <- ss.dimX$Sample.Size
    p <- ss.dimX$Dimension
    
    if (use == "complete.obs") {
      ccX <-  1:n
      n <- length(ccX)
      if (type.X == "sample") {
        ccX <- which(stats::complete.cases(X))
      }
      if (type.X == "sample" && p == 1) {
        X <- X[ccX]
      } else if (type.X == "sample" && p > 1) {
        X <- X[ccX, ]
      }
      if (type.X == "distance") {
        X <- X[ccX,ccX]
      }
    }
    
    
    ## normalize samples if calculation of affinely invariant distance covariance is desired
    if (affine == TRUE) {
      if (type.X == "distance") {
        stop("Affinely invariant distance variance cannot be calculated for type distance")
      }
      if (p > n) {
        stop("Affinely invariant distance variance cannot be calculated for p>n")
      }
      if (p > 1) {
        X <- X %*% Rfast::spdinv(mroot(stats::var(X)))
      } else {
        X <- X / stats::sd(X)
      }
    } else if (standardize) {
      if (type.X == "distance") {
        stop("Standardization cannot be applied for type distance.")
      }
      if (p > 1) {
        X <- standardise(X, center = FALSE)
      } else {
        X <- X / stats::sd(X)
      }
    }
    
    
    alg.fast <- alg.standard <- alg.memsave <- FALSE
    
    if (algorithm == "auto") {
      if (p == 1 & metr.X[1] %in% c("euclidean", "discrete") 
            & type.X == "sample" & n > 100) {
        alg.fast <- TRUE
      } else if (metr.X[1] %in% c("euclidean", "alpha", "gaussian", "boundsq", "minkowski", "discrete") & type.X == "sample") {
        alg.memsave <- TRUE
      } else {
        alg.standard <- TRUE
      }
    } else if (algorithm == "fast") {
      alg.fast <- TRUE 
    } else if (algorithm == "standard") {
      alg.standard <- TRUE
    }  else if (algorithm == "memsave") {
      alg.memsave <- TRUE
    } else
      stop ("Algorithm must be one of \"fast\", \"standard\", \"memsave\" or \"auto\"")
    
    if (alg.fast) {
      if (p == 1) {
        if (metr.X[1] == "euclidean") {
          terms <- dcovterms.fast(X = X, n = n, calc.dvar = TRUE)
        } else if (metr.X[1] == "discrete") {
          X <- as.factor(X)
          terms <- dcovterms.fast.discrete(X = X, n = n, calc.dvar =TRUE)
        } else {
          stop("metr.X has to be \"euclidean\" or \"discrete\" for fast algorithms")
        }
      } else {
        stop("Dimensions of X must be 1 for fast algorithms.")
      }
    } else if (alg.memsave) {
      if (metr.X[1] %in% c("euclidean", "alpha", "gaussian", "boundsq", "minkowski", "discrete")) {
        terms <- dvarterms.memsave(X, metr.X, p)
      } else {
        stop("Memory efficient algorithms cannot be run with user-defined metrics")
      }
    } else if (alg.standard) {
      terms <- dcovterms.standard(X = X, type.X = type.X, metr.X = metr.X, p = p, calc.dvar = TRUE)
    }
    
    
    if (bias.corr) {
      dvar <- terms$aijaij / n / (n - 3) - 2 * terms$Saa / n / (n - 2) / (n - 3) + terms$Taa / n / (n - 1) / (n - 2) / (n - 3)
      dsd <- sqrt(abs(dvar))
    } else {
      dvar <- terms$aijaij / n / n - 2 * terms$Saa / n / n / n  + terms$Taa / n / n / n / n
      dsd <-  sqrt(abs(dvar))
    }
    return(dsd) 

  }





#' Calculates the distance correlation \insertCite{szekely2007,szekely2009brownian}{dcortools}.
#'
#' @param X contains either the first  sample or its corresponding distance matrix.
#'
#' In the first case, X can be provided either as a vector (if one-dimensional), a matrix or a data.frame (if two-dimensional or higher). 
#'
#' In the second case, the input must be a distance matrix corresponding to the sample of interest.
#'
#' If X is a sample, type.X must be specified as "sample". If X is a distance matrix, type.X must be specified as "distance".
#' @param Y see X.
#' @param affine logical; specifies if the affinely invariant distance correlation \insertCite{dueck2014affinely}{dcortools} should be calculated or not.
#' @param standardize logical; specifies if X and Y should be standardized dividing each component by its standard deviations. No effect when affine = TRUE.
#' @param bias.corr logical; specifies if the bias corrected version of the sample distance correlation \insertCite{huo2016fast}{dcortools} should be calculated.
#' @param type.X For "distance", X is interpreted as a distance matrix. For "sample", X is interpreted as a sample.
#' @param type.Y see type.X.
#' @param metr.X specifies the metric which should be used to compute the distance matrix for X (ignored when type.X = "distance").
#' 
#'  Options are "euclidean", "discrete", "alpha", "minkowski", "gaussian", "gaussauto", "boundsq" or user-specified metrics (see examples).
#'  
#'  For "alpha", "minkowski", "gaussian", "gaussauto" and "boundsq", the corresponding parameters are specified via "c(metric, parameter)", e.g. c("gaussian", 3) for a Gaussian metric with bandwidth parameter 3; the default parameter is 2 for "minkowski" and "1" for all other metrics.
#'  
#'  See \insertCite{lyons2013distance,sejdinovic2013equivalence,bottcher2017detecting;textual}{dcortools} for details.
#' @param metr.Y see metr.X.
#' @param use specifies how to treat missing values. "complete.obs" excludes observations containing NAs, "all" uses all observations.
#' @param algorithm specifies the algorithm used for calculating the distance correlation. 
#' 
#' "fast" uses an O(n log n) algorithm if the observations are one-dimensional and metr.X and metr.Y are either "euclidean" or "discrete", see also \insertCite{huo2016fast;textual}{dcortools}. 
#' 
#' "memsave" uses a memory saving version of the standard algorithm with computational complexity O(n^2) but requiring only O(n) memory. 
#' 
#' "standard" uses the classical algorithm. User-specified metrics always use the classical algorithm.
#' 
#' "auto" chooses the best algorithm for the specific setting using a rule of thumb.
#' 
#' @return numeric; the distance correlation between samples X and Y.
#' @references
#' \insertRef{bottcher2017detecting}{dcortools}
#' 
#' \insertRef{dueck2014affinely}{dcortools}
#' 
#' \insertRef{huo2016fast}{dcortools}
#' 
#' \insertRef{lyons2013distance}{dcortools}
#' 
#' \insertRef{sejdinovic2013equivalence}{dcortools}
#' 
#' \insertRef{szekely2007}{dcortools}
#' 
#' \insertRef{szekely2009brownian}{dcortools}
#' @export
#' @examples 
#' X <- rnorm(200)
#' Y <- rnorm(200)
#' Z <- X + rnorm(200)
#' dim(X) <- dim(Y) <- dim(Z) <- c(20, 10)
#' 
#' #Demonstration that biased-corrected distance correlation is
#' #often more meaningful than without using bias-correction
#' 
#' distcor(X, Y) 
#' distcor(X, Z) 
#' distcor(X, Y, bias.corr = TRUE)
#' distcor(X, Z, bias.corr = TRUE)
#' 
#' #For more examples of the different options,
#' #see the documentation of distcov.
distcor <-
  function(X,
           Y,
           affine = FALSE,
           standardize=FALSE,
           bias.corr = FALSE,
           type.X = "sample",
           type.Y = "sample",
           metr.X = "euclidean",
           metr.Y = "euclidean",
           use = "all",
           algorithm = "auto") {
    #extract dimensions and sample sizes
    ss.dimX <- extract_np(X, type.X)
    ss.dimY <- extract_np(Y, type.Y)
    
    n <- ss.dimX$Sample.Size
    p <- ss.dimX$Dimension
    
    m <- ss.dimY$Sample.Size
    q <- ss.dimY$Dimension
    
    if  (use == "complete.obs") {
      ccX <- ccY <- cc <- 1:n
      ccX <- which(stats::complete.cases(X))
      ccY <- which(stats::complete.cases(Y))
      cc <- intersect(ccX, ccY)
      n <- m <- length(cc)
      if (type.X == "sample" && p == 1) {
        X <- X[cc]
      } else if (type.X == "sample" && p > 1) {
        X <- X[cc, ]
      }
      if (type.Y == "sample" && q == 1) {
        Y <- Y[cc]
      } else if (type.X == "sample" && q > 1) {
        Y <- Y[cc, ]
      }
      if (type.X == "distance") {
        X <- X[cc,cc]
      }
      if (type.Y == "distance") {
        Y <- Y[cc,cc]
      }
    }
    
    
    if (n != m) {
      stop("Samples X and Y must have the same sizes!")
    }
    
    
    ## normalize samples if calculation of affinely invariant distance correlation is desired
    if (affine == TRUE) {
      if (type.X == "distance" | type.Y == "distance") {
        stop("Affinely invariant distance correlation cannot be calculated for type distance")
      }
      if (p > n | q > n) {
        stop("Affinely invariant distance correlation cannot be calculated for p>n")
      }
      if (p > 1) {
        X <- X %*% Rfast::spdinv(mroot(stats::var(X)))
      } else {
        X <- X / stats::sd(X)
      }
      if (q > 1) {
        Y <- Y %*% Rfast::spdinv(mroot(stats::var(Y)))
      } else {
        Y <- Y / stats::sd(Y)
      }
    } else if (standardize) {
      if (type.X == "distance" | type.Y == "distance") {
        stop("Standardization cannot be applied for type distance.")
      }
      if (p > 1) {
        X <- standardise(X, center = FALSE)
      } else {
        X <- X / stats::sd(X)
      }
      if (q > 1) {
        Y <- standardise(Y, center = FALSE)
      } else {
        Y <- Y / stats::sd(Y)
      }
    }
    
    
    alg.fast <- alg.standard <- alg.memsave <- FALSE
    
    if (algorithm == "auto") {
      if (p == 1 & q == 1 & metr.X[1] %in% c("euclidean", "discrete") 
          & metr.Y[1] %in% c("euclidean", "discrete") & type.X == "sample" & type.Y == "sample" & n > 100) {
        alg.fast <- TRUE
      } else if (metr.X[1] %in% c("euclidean", "alpha", "gaussian", "boundsq", "minkowski", "discrete") &
                 metr.Y[1] %in% c("euclidean", "alpha", "gaussian", "boundsq", "minkowski", "discrete") & type.X == "sample" & type.Y == "sample") {
        alg.memsave <- TRUE
      } else {
        alg.standard <- TRUE
      }
    } else if (algorithm == "fast") {
      alg.fast <- TRUE 
    } else if (algorithm == "standard") {
      alg.standard <- TRUE
    }  else if (algorithm == "memsave") {
      alg.memsave <- TRUE
    } else
      stop ("Algorithm must be one of \"fast\", \"standard\", \"memsave\" or \"auto\"")
    
    
    terms <- dcovterms(X, Y, n, calc.dcor = TRUE, doperm = FALSE, dobb3 = FALSE, alg.fast = alg.fast, alg.memsave = alg.memsave, alg.standard = alg.standard, p = p, q = q, metr.X = metr.X, metr.Y =metr.Y, type.X = type.X, type.Y = type.Y)
    
    
    if (bias.corr) {
      dvarX <- terms$aijaij / n / (n - 3) - 2 * terms$Saa / n / (n - 2) / (n - 3) + (terms$adotdot * terms$adotdot) / n / (n - 1) / (n - 2) / (n - 3)
      dvarY <- terms$bijbij / n / (n - 3) - 2 * terms$Sbb / n / (n - 2) / (n - 3) + (terms$bdotdot * terms$bdotdot) / n / (n - 1) / (n - 2) / (n - 3)
      dcov2 <- terms$aijbij / n / (n - 3) - 2 * terms$Sab / n / (n - 2) / (n - 3) + (terms$adotdot * terms$bdotdot) / n / (n - 1) / (n - 2) / (n - 3)
      dcov <-  sqrt(abs(dcov2)) * sign(dcov2)
      dcor <-  dcov / sqrt(sqrt(dvarX * dvarY))
    } else {
      dvarX <- terms$aijaij / n / n - 2 * terms$Saa / n / n / n + (terms$adotdot * terms$adotdot) / n / n / n / n
      dvarY <- terms$bijbij / n / n - 2 * terms$Sbb / n / n / n + (terms$bdotdot * terms$bdotdot) / n / n / n / n
      dcov2 <- terms$aijbij / n / n - 2 * terms$Sab / n / n / n  + (terms$adotdot * terms$bdotdot) / n / n / n / n
      dcov <- sqrt(abs(dcov2))
      dcor <- dcov / sqrt(sqrt(dvarX * dvarY))
    }
    
 
    return(dcor)
  }
