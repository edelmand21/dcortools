#' Calculate the distance covariance
#'
#' @param X contains either the first sample or its corresponding distance matrix. In the first case, this input can be either a vector of positive length, a matrix with one column or a data.frame with one column. In this case, type.X must be specified as "sample". In the second case, the input must be a distance matrix corresponding to the sample of interest. In this second case, type.X must be "distance".
#' @param Y see X.
#' @param affine logical; indicates if the affinely transformed distance covariance should be calculated or not.
#' @param bias.corr logical; indicates if the bias corrected version of the sample distance covariance should be calculated.
#' @param type.X For "distance", X is interpreted as a distance matrix. For "sample" (or any other value), X is intepreted as a sample
#' @param type.Y see type.X.
#' @param metr.X specifies the metric which should be used for X to analyse the distance covariance. TO DO: Provide details for this.
#' @param metr.Y see metr.X.
#' @param use specifies how to treat missing values. "complete.obs" excludes NA's, "all" (or any other value) uses all observations.
#' @param algorithm : "fast" uses an O(n log n) algorithm if the observations are one-dimensional and metr.X and metr.Y are either "euclidean" or discrete. "memsave" uses a memory saving algorithm, which does not save the distance matrices. "standard" uses the classical algorithm. "auto" chooses the best algorithm for the specific setting using a rule of thumb.
#' 
#' @return numeric giving the distance covariance between samples X and Y.
#' @export
distcov <-
  function(X,
           Y,
           affine = FALSE,
           bias.corr = TRUE,
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
      if (type.X == "sample") {
        ccX <- which(complete.cases(X))
      }
      if (type.Y == "sample") {
        ccY <- which(complete.cases(Y))
      }
      cc <- intersect(ccX, ccY)
      if (type.X == "sample" && p == 1) {
        X <- X[cc]
      } else if (type.X == "sample" && p > 1) {
        X <- X[cc, ]
      }
      if (type.Y == "sample" && p == 1) {
        Y <- Y[cc]
      } else if (type.X == "sample" && p > 1) {
        Y <- Y[cc, ]
      }
      n <- m <- length(cc)
      
      if (type.X == "distance") {
        X <- X[cc,cc]
      }
      if (type.Y == "distance") {
        Y <- Y[cc,cc]
      }
    }
    
    ## normalize samples if calculation of affinely invariant distance covariance is desired
    if (affine) {
      if (p > n | q > n) {
        stop("Affinely invariant distance covariance cannot be calculated for p>n")
      }
      if (type.X == "distance" | type.Y == "distance") {
        stop("Affinely invariant distance covariance cannot be calculated for type distance")
      }
      if (p > 1) {
        X <- X %*% Rfast::spdinv(mroot(var(X)))
      } else {
        X <- X / sd(X)
      }
      if (q > 1) {
        Y <- Y %*% Rfast::spdinv(mroot(var(Y)))
      } else {
        Y <- Y / sd(Y)
      }
    }
    
    if (algorithm == "auto") {
      if (p == 1 & q == 1 & metr.X[1] %in% c("euclidean", "discrete") 
          & metr.Y[1] %in% c("euclidean", "discrete") & n > 100) {
        algorithm <- "fast"
      } else if (metr.X[1] %in% c("euclidean", "alpha", "gaussian", "boundsq", "minkowski", "discrete") &
                 metr.Y[1] %in% c("euclidean", "alpha", "gaussian", "boundsq", "minkowski", "discrete")) {
        algorithm <- "memsave"
      } else {
        algorithm <- "standard"
      }
    }
    
    if (algorithm == "fast") {
      if (p == 1 & q == 1) {
        if (metr.X == "euclidean" & metr.Y == "euclidean") {
          terms <- dcovterms.fast(X, Y, n)
        } else if (metr.X == "euclidean" & metr.Y == "discrete") {
          Y <- as.factor(Y)
          terms <- dcovterms.fast.numdisc(X, Y, n)
        } else if (metr.X == "discrete" & metr.Y == "euclidean") {
          X <- as.factor(X)
          terms <- dcovterms.fast.numdisc(Y, X, n)
        } else if (metr.X == "discrete" & metr.Y == "discrete") {
          X <- as.factor(X)
          Y <- as.factor(Y)
          terms <- dcovterms.fast.discrete(X ,Y, n)
        } else {
          stop("metr.X and metr.Y have to be \"euclidean\" or \"discrete\" for fast algorithms")
        }
      } else {
          stop("Dimensions of X and Y must be 1 for fast algorithms.")
        }
    } else if (algorithm == "memsave") {
      if (metr.X[1] %in% c("euclidean", "alpha", "gaussian", "boundsq", "minkowski", "discrete") &
          metr.Y[1] %in% c("euclidean", "alpha", "gaussian", "boundsq", "minkowski", "discrete")) {
        terms <- dcovterms.memsave(X, Y, metr.X, metr.Y, p, q)
      } else {
        stop("Memory efficient algorithms cannot be run with user-defined metrics")
      }
    } else if (algorithm == "standard") {
      terms <- dcovterms.standard(X, Y, type.X, type.Y, metr.X, metr.Y, p, q)
    }
    
   
    if (bias.corr) {
      dcov2 <- terms$aijbij / n / (n - 3) - 2 * terms$Sab / n / (n - 2) / (n - 3) + terms$Tab / n / (n - 1) / (n - 2) / (n - 3)
      dcov <- sqrt(abs(dcov2)) * sign(dcov2) 
    } else {
      dcov2 <- terms$aijbij / n / n - 2 * terms$Sab / n / n / n  + terms$Tab / n / n / n / n
      dcov <- sqrt(abs(dcov2))
    }
    return(dcov)
}




#' Calculates the distance standard deviation
#'
#' @param X contains either the first sample or its corresponding distance matrix.
#'
#' In the first case, this input can be either a vector of positive length,
#'
#' a matrix with one column or a data.frame with one column.
#'
#' In this case, type.X must be specified as "sample".
#'
#' In the second case, the input must be a distance matrix corresponding to the sample of interest.
#'
#' In this second case, type.X must be "distance".
#' @param affine logical; indicates if the affinely transformed distance standard deviation should be calculated or not.
#' @param bias.corr logical; indicates if the bias corrected version of the sample distance variance should be calculated.
#' @param type.X either "sample" or "distance"; specifies the type of input for X.
#' @param metr.X specifies the metric which should be used for X to analyse the distance variance TO DO: Provide details for this.
#' @param use : "all" uses all observations, "complete.obs" excludes NA's
#' @return numeric giving the distance variance of the sample X..
#' @export
distsd <-
  function(X,
           affine = FALSE,
           bias.corr = TRUE,
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
      if (type.X == "sample") {
        ccX <- which(complete.cases(X))
      }
      if (type.X == "sample" && p == 1) {
        X <- X[ccX]
      } else if (type.X == "sample" && p > 1) {
        X <- X[ccX, ]
      }
      n <- length(ccX)
      if (type.X == "distance") {
        X <- X[cc,cc]
      }
    }
    
    
    ## normalize samples if calculation of affinely invariant distance covariance is desired
    if (affine == TRUE) {
      if (p > n) {
        stop("Affinely invariant distance variance cannot be calculated for p>n")
      }
      if (type.X == "distance") {
        stop("Affinely invariant distance variance cannot be calculated for type distance")
      }
      if (p > 1) {
        X <- X %*% Rfast::spdinv(mroot(var(X)))
      } else {
        X <- X / sd(X)
      }
    }
    
    if (algorithm == "auto") {
      if (p == 1 & metr.X[1] %in% c("euclidean", "discrete") & n > 100) {
        algorithm <- "fast"
      } else if (metr.X[1] %in% c("euclidean", "alpha", "gaussian", "boundsq", "minkowski", "discrete")) {
        algorithm <- "memsave"
      } else {
        algorithm <- "standard"
      }
    }
    
    if (algorithm == "fast") {
      if (p == 1) {
        if (metr.X == "euclidean") {
          terms <- dcovterms.fast(X = X, n = n, calc.dvar = TRUE)
        } else if (metr.X == "discrete") {
          X <- as.factor(X)
          terms <- dcovterms.fast.discrete(X = X, n = n, calc.dvar =TRUE)
        } else {
          stop("metr.X has to be \"euclidean\" or \"discrete\" for fast algorithms")
        }
      } else {
        stop("Dimensions of X must be 1 for fast algorithms.")
      }
    } else if (algorithm == "memsave") {
      if (metr.X[1] %in% c("euclidean", "alpha", "gaussian", "boundsq", "minkowski", "discrete")) {
        terms <- dvarterms.memsave(X, metr.X, p)
      } else {
        stop("Memory efficient algorithms cannot be run with user-defined metrics")
      }
    } else if (algorithm == "standard") {
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





#' Calculates the distance correlation
#'
#' @param X contains either the first sample or its corresponding distance matrix. In the first case, this input can be either a vector of positive length, a matrix with one column or a data.frame with one column. In this case, type.X must be specified as "sample". In the second case, the input must be a distance matrix corresponding to the sample of interest. In this second case, type.X must be "distance".
#' @param Y see X.
#' @param affine logical; indicates if the affinely transformed distance correlation should be calculated or not.
#' @param bias.corr logical; indicates if the bias corrected version of the sample distance correlation should be calculated.
#' @param type.X either "sample" or "distance"; specifies the type of input for X.
#' @param type.Y see type.X.
#' @param metr.X specifies the metric which should be used for X to analyse the distance correlation. TO DO: Provide details for this.
#' @param metr.Y see metr.X.
#' @param use : "all" uses all observations, "complete.obs" excludes NA's
#' @return numeric giving the distance correlation between samples X and Y.
#' @export

distcor <-
  function(X,
           Y,
           affine = FALSE,
           bias.corr = TRUE,
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
    
    # if (type.X == "distance") {
    #   if (is.null(dim.X)) {
    #     stop("When X is a distance matrix the dimension of X must be specified as input.")
    #   }
    #   p <- dim.X
    # }
    # 
    # if (type.Y == "distance") {
    #   if (is.null(dim.X)) {
    #     stop("When Y is a distance matrix the dimension of X must be specified as input.")
    #   }
    #   q <- dim.Y
    # }
    
    if (use == "complete.obs") {
      ccX <- ccY <- cc <- 1:n
      if (type.X == "sample") {
        ccX <- which(complete.cases(X))
      }
      if (type.Y == "sample") {
        ccY <- which(complete.cases(Y))
      }
      cc <- intersect(ccX, ccY)
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
      n <- m <- length(cc)
      
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
      if (p > n | q > n) {
        stop("Affinely invariant distance correlation cannot be calculated for p>n")
      }
      if (type.X == "distance" | type.Y == "distance") {
        stop("Affinely invariant distance correlation cannot be calculated for type distance")
      }
      if (p > 1) {
        X <- X %*% Rfast::spdinv(mroot(var(X)))
      } else {
        X <- X / sd(X)
      }
      if (q > 1) {
        Y <- Y %*% Rfast::spdinv(mroot(var(Y)))
      } else {
        Y <- Y / sd(Y)
      }
    }
    
    
    if (algorithm == "auto") {
      if (p == 1 & q == 1 & metr.X[1] %in% c("euclidean", "discrete") 
          & metr.Y[1] %in% c("euclidean", "discrete") & n > 100) {
        algorithm <- "fast"
      } else if (metr.X[1] %in% c("euclidean", "alpha", "gaussian", "boundsq", "minkowski", "discrete") &
                 metr.Y[1] %in% c("euclidean", "alpha", "gaussian", "boundsq", "minkowski", "discrete")) {
        algorithm <- "memsave"
      } else {
        algorithm <- "standard"
      }
    }
    
    if (algorithm == "fast") {
      if (p == 1 & q == 1) {
        if (metr.X == "euclidean" & metr.Y == "euclidean") {
          terms <- dcovterms.fast(X, Y, n, calc.dcor = TRUE)
        } else if (metr.X == "euclidean" & metr.Y == "discrete") {
          Y <- as.factor(Y)
          terms <- dcovterms.fast.numdisc(X, Y, n, calc.dcor = TRUE)
        } else if (metr.X == "discrete" & metr.Y == "euclidean") {
          X <- as.factor(X)
          terms <- dcovterms.fast.numdisc(Y, X, n, calc.dcor = TRUE)
        } else if (metr.X == "discrete" & metr.Y == "discrete") {
          X <- as.factor(X)
          Y <- as.factor(Y)
          terms <- dcovterms.fast.discrete(X ,Y, n, calc.dcor = TRUE)
        } else {
          stop("metr.X and metr.Y have to be \"euclidean\" or \"discrete\" for fast algorithms")
        }
      } else {
        stop("Dimensions of X and Y must be 1 for fast algorithms.")
      }
    } else if (algorithm == "memsave") {
      if (metr.X[1] %in% c("euclidean", "alpha", "gaussian", "boundsq", "minkowski", "discrete") &
          metr.Y[1] %in% c("euclidean", "alpha", "gaussian", "boundsq", "minkowski", "discrete")) {
        terms <- dcovterms.memsave(X, Y, metr.X, metr.Y, p, q, calc.dcor = TRUE)
      } else {
        stop("Memory efficient algorithms cannot be run with user-defined metrics")
      }
    } else if (algorithm == "standard") {
      terms <- dcovterms.standard(X, Y, type.X, type.Y, metr.X, metr.Y, p, q, calc.dcor = TRUE)
    }
    
    
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
