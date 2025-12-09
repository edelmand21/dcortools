

#' Performs a distance covariance test.
#'
#' @param X contains either the first sample or its corresponding distance matrix.
#'
#' In the first case, X can be provided either as a vector (if one-dimensional), a matrix or a data.frame (if two-dimensional or higher). 
#'
#' In the second case, the input must be a distance matrix corresponding to the sample of interest.
#'
#' If X is a sample, type.X must be specified as "sample". If X is a distance matrix, type.X must be specified as "distance".
#' @param Y see X.
#' @param method specifies the type of test that is performed.
#'  
#' "permutation" performs a Monte Carlo Permutation test. 
#' 
#' "gamma" performs a test based on a gamma approximation of the test statistic under the null \insertCite{huang2017statistically}{dcortools}. This test tends to be anti-conservative, if the ``real'' p-value is small
#' 
#' "conservative" performs a conservative two-moment approximation \insertCite{berschneider2018complex}{dcortools}.
#'   
#' "bb3" performs a  three-moment approximation  \insertCite{berschneider2018complex}{dcortools}. This is the most precise parametric option, but only available with the standard algorithm.
#' 
#' "wildbs1" and "wilbs2" perform wild bootstrap tests \insertCite{chwialkowski2014wild}{dcortools}; experimental at the moment.
#' @param b integer; specifies the number of random permutations/bootstrap samples used for the permutation or wild bootstraps tests. Ignored for other tests.
#' @param ln numeric; block size parameter for wild bootstrap tests. Ignored for other tests.
#' @param affine logical; specifies if the affinely invariant distance covariance \insertCite{dueck2014affinely}{dcortools} should be calculated or not.
#' @param standardize logical; specifies if X and Y should be standardized dividing each component by its standard deviations. No effect when affine = TRUE.
#' @param bias.corr logical; specifies if the bias corrected version of the sample distance covariance \insertCite{huo2016fast}{dcortools} should be calculated.
#' @param type.X For "distance", X is interpreted as a distance matrix. For "sample", X is interpreted as a sample.
#' @param type.Y see type.X.
#' @param metr.X specifies the metric which should be used to compute the distance matrix for X (ignored when type.X = "distance").
#' 
#'  Options are "euclidean", "discrete", "alpha", "minkowski", "gaussian", "gaussauto", "boundsq" or user-specified metrics (see examples).
#'  
#'  For "alpha", "minkowski", "gauss", "gaussauto" and "boundsq", the corresponding parameters are specified via "c(metric, parameter)", c("gaussian", 3) for example uses a Gaussian metric with bandwidth parameter 3; the default parameter is 2 for "minkowski" and "1" for all other metrics.
#'  
#'  See \insertCite{lyons2013distance,sejdinovic2013equivalence,bottcher2017detecting;textual}{dcortools} for details.
#' @param metr.Y see metr.X.
#' @param use specifies how to treat missing values. "complete.obs" excludes NAs, "all" uses all observations.
#' @param return.data logical; specifies if the test object should contain the original data.
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
#' @return distcov.test object
#' @export
#' @references
#' \insertRef{berschneider2018complex}{dcortools}
#' 
#' \insertRef{bottcher2017detecting}{dcortools}
#' 
#' \insertRef{chwialkowski2014wild}{dcortools}
#' 
#' \insertRef{dueck2014affinely}{dcortools}
#' 
#' \insertRef{huang2017statistically}{dcortools}
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
#'
distcov.test <- function(X,
                         Y,
                         method = "permutation",
                         b = 499L,
                         ln = 20,
                         affine = FALSE,
                         standardize = FALSE,
                         bias.corr = FALSE,
                         type.X = "sample",
                         type.Y = "sample",
                         metr.X = "euclidean",
                         metr.Y = "euclidean",
                         use = "all",
                         return.data = FALSE,
                         algorithm = "auto") {
  
  smp <- NULL
  #extract dimensions and sample sizes
  ss.dimX <- extract_np(X, type.X)
  ss.dimY <- extract_np(Y, type.Y)
  
  n <- ss.dimX$Sample.Size
  p <- ss.dimX$Dimension
  
  m <- ss.dimY$Sample.Size
  q <- ss.dimY$Dimension
  
  dogamma <- docons  <- dobb3 <-  doperm <- dowild1 <- dowild2 <-  FALSE
  
  output <- list()
  
  if (method == "gamma")
    dogamma <- TRUE 
    else if (method == "conservative")
    docons <- TRUE
    else if (method == "bb3")
    dobb3 <- TRUE
    else if (method == "permutation")
    doperm <- TRUE
    else if (method == "wildbs1")
    dowild1 <- TRUE
    else if (method == "wildbs2")
    dowild2 <- TRUE
    else
    stop ("Test must be one of \"permutation\", \"gamma\", \"bb3\", \"conservative\", \"wildbs1\" or \"wildbs2\" ")
    
  
  
  
  if (n != m) {
    stop("Samples X and Y must have the same sizes!")
  }
  
  if(return.data) {
    output$X <- X
    output$Y <- Y
  } else {
    output <- NULL
  }
    
  
  if (bias.corr == TRUE) {
    termstodcov2 <- function(aijbij, Sab, Tab, n) {
      aijbij / n / (n - 3) - 2 * Sab / n / (n - 2) / (n - 3) + Tab / n / (n - 1) / (n - 2) / (n - 3) 
    }
    dcov2todcov <- function(dcov2) {
        sqrt(abs(dcov2)) * sign(dcov2)
    }
    dcov2todcor <- function(dcov2, dvarX, dvarY) {
        (sqrt(abs(dcov2)) * sign(dcov2)) / sqrt(sqrt(dvarX * dvarY))
    }
  } else  {
      termstodcov2 <- function(aijbij, Sab, Tab, n) {
        aijbij / n / n - 2 * Sab / n / n / n + Tab / n / n / n / n
      }
      dcov2todcov <- function(dcov2) {
        sqrt(dcov2)
      }
      dcov2todcor <- function(dcov2, dvarX, dvarY) {
        sqrt(dcov2) / sqrt(sqrt(dvarX * dvarY))
      }
  }
  
  
  
  if (dogamma) {
    testfunc <- function(dcov2, smp, terms, n, moms, ...) {
      dvarX <- terms$aijaij / n / (n - 3) - 2 * terms$Saa / n / (n - 2) / (n - 3) + terms$adotdot * terms$adotdot / n / (n - 1) / (n - 2) / (n - 3) 
      dvarY <- terms$bijbij / n / (n - 3) - 2 * terms$Sbb / n / (n - 2) / (n - 3) + terms$bdotdot * terms$bdotdot / n / (n - 1) / (n - 2) / (n - 3) 
      dcov2 <- terms$aijbij / n / (n - 3) - 2 * terms$Sab / n / (n - 2) / (n - 3) + terms$adotdot * terms$bdotdot / n / (n - 1) / (n - 2) / (n - 3) 
      U1 <- dvarX  * dvarY
      U2 <- terms$adotdot / n / (n - 1)
      U3 <- terms$bdotdot / n / (n - 1)
      alph <- 1/2*(U2 ^ 2 * U3 ^ 2) / U1
      beta <- 1/2*(U2 * U3) / U1
      stat <- n *  dcov2 + U2 * U3
      pval <- stats::pgamma(stat, alph, beta, lower.tail = FALSE) 
      return(pval)
    }
  } else if (doperm) {
    testfunc <- function(dcov2, smp, terms, n, moms, ...) {
      if (is.na(dcov2))
        return(NA)
      Tab <- terms$adotdot * terms$bdotdot
      reps <- lapply(1:b, function(t) {
        terms.sample <- terms.smp(terms = terms, smp = smp[[t]])
        return(termstodcov2(terms.sample$aijbij, terms.sample$Sab, Tab, n))
      })
      pval <- (1 + length(which(reps >= dcov2))) / (1 + b)
      return(pval)
    }
  } else if (dowild1 | dowild2) {
    testfunc <- function(dcov2, smp, terms, n, moms, ...) {
      cX <- doublecent(terms$distX)
      cY <- doublecent(terms$distY)
      dXY <- cX*cY
      stat <- sum(dXY)
      reps <- sapply(1:b, function(t) {
        return(t(smp[[t]]) %*% dXY %*% smp[[t]])
      })
      pval <- (1 + length(which(reps >= stat))) / (1 + b)
      return(pval)
    }
  } else if (docons) {
    testfunc <- function(dcov2, smp, terms, n, moms, ...) {
      est.m2 <- sum(moms) / n ^ 10
      est.m1 <- terms$adotdot * terms$bdotdot / n ^ 3 / (n - 1)
      est.var <- (est.m2 - est.m1 ^ 2)
      alpha <- sqrt(est.var / 2 / est.m1 ^ 2)
      stat <- terms$aijbij / n - 2 * terms$Sab / n ^ 2 + terms$adotdot * terms$bdotdot / n ^ 3
      pval <- stats::pchisq(stat * sqrt(2) / sqrt(est.var), df = 1 / alpha, lower.tail = FALSE) 
    }
  } else if (dobb3) {
    testfunc <- function(dcov2, smp, terms, n, moms,...) {
      est.m2 <- sum(moms$vc) / n ^ 10
      est.m1 <- terms$adotdot * terms$bdotdot / n ^ 3 / (n - 1)
      est.var <- (est.m2 - est.m1 ^ 2)
      est.skw <- moms$skw
      beta <- est.skw / sqrt(8)
      stat <- terms$aijbij / n - 2 * terms$Sab / n ^ 2 + terms$adotdot * terms$bdotdot / n ^ 3
      centstat <- (stat - est.m1) /  sqrt(est.var)
      pval <- stats::pchisq((centstat * sqrt(2) + 1 / beta) / beta , df = 1 / beta ^ 2, lower.tail = FALSE)  
      return(pval)
    } 
  }
  
  
  
  

  if (use == "complete.obs") {
    ccX <- ccY <- cc <- 1:n
    if (type.X == "sample") {
      ccX <- which(stats::complete.cases(X))
    }
    if (type.Y == "sample") {
      ccY <- which(stats::complete.cases(Y))
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
  
  
  ## normalize samples if calculation of affinely invariant distance covariance is desired
  if (affine == TRUE) {
    if (type.X == "distance" | type.Y == "distance") {
      stop("Affinely invariant distance covariance cannot be calculated for type 'distance'")
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
  
  if (algorithm == "auto") {
    if (p == 1 & q == 1 & metr.X[1] %in% c("euclidean", "discrete") 
        & metr.Y %in% c("euclidean", "discrete") & n > 100 &  type.X == "sample" & type.Y == "sample" & !(dobb3|dowild1|dowild2)) {
      algorithm <- "fast"
    } else if (metr.X[1] %in% c("euclidean", "alpha", "gaussian", "boundsq", "minkowski", "discrete") &
               metr.Y[1] %in% c("euclidean", "alpha", "gaussian", "boundsq", "minkowski", "discrete") & type.X == "sample" & type.Y == "sample" & !(dobb3|dowild1|dowild2)) {
      algorithm <- "memsave"
    } else {
      algorithm <- "standard"
    }
  }
  
  alg.fast <- alg.standard <- alg.memsave <- FALSE
  
  if (algorithm == "fast") {
    alg.fast <- TRUE 
    if (doperm) 
      terms.smp <- function(terms, smp) {sampleterms.fast.discrete(terms, smp)}
  } else if (algorithm == "standard") {
    alg.standard <- TRUE
    if (doperm) 
      terms.smp <- function(terms, smp) {sampleterms.standard(terms, smp)}
  }  else if (algorithm == "memsave") {
    alg.memsave <- TRUE
    if (doperm) 
      terms.smp <- function(terms, smp) {sampleterms.memsave(terms, smp)}
  } 
  else
    stop ("Algorithm must be one of \"fast\", \"standard\", \"memsave\" or \"auto\"")
  
  
  if (!alg.standard & (dobb3|dowild1|dowild2))
    stop("bb3 and wild bootstrap p-value calculation is only possible with algorithm=standard!")
  
  
  if (dowild1 | dowild2) {
   perms <- lapply(1:b, function(t) {
      epsilon <- stats::rnorm(n+1)
      W <- rep(NA,n)
      W[1] <- exp(-1/ln) * epsilon[1]+sqrt(1-exp(-2/ln)) *epsilon[2]
      for (r in 2:n) {
        W[r] <- exp(-1/ln) * W[r-1] + sqrt(1-exp(-2/ln))*epsilon[r+1]
      }
      if (dowild2)
        W <- W - mean(W)
      return(W)
    })
  } else if (doperm) {
    perms <- lapply(1:b, function(t) sample(1:n))
  } else {
    perms <- NULL
  }
  
  if (type.X == "distance")
    metr.X <- "distance matrix provided"
  
  if (type.Y == "distance")
    metr.Y <- "distance matrix provided"
  
  terms <- dcovterms(X,Y,n, calc.dcor = TRUE, doperm = doperm, dobb3 = dobb3|dowild1|dowild2, alg.fast = alg.fast, alg.memsave = alg.memsave, alg.standard = alg.standard, p = p, q = q, metr.X = metr.X, metr.Y =metr.Y, type.X = type.X, type.Y = type.Y)
  
  
  
  

  
  dcov2 <- termstodcov2(terms$aijbij, terms$Sab, terms$adotdot * terms$bdotdot, n)
  output$dcov <- dcov2todcov(dcov2)
  dvarX <- termstodcov2(terms$aijaij, terms$Saa, terms$adotdot * terms$adotdot, n)
  dvarY <- termstodcov2(terms$bijbij, terms$Sbb, terms$bdotdot * terms$bdotdot, n)
  output$dsdX <- sqrt(dvarX)
  output$dsdY <- sqrt(dvarY)
  output$dcor <- dcov2todcor(dcov2, dvarX, dvarY)
  
  if (docons) {
    moms <- calcmom(terms$aijaij, terms$Saa, terms$adotdot, terms$bijbij, terms$Sbb, terms$bdotdot, n = n, dobb3 = FALSE)
  } else if (dobb3) {
    moms <- calcmom(terms$aijaij, terms$Saa, terms$adotdot, terms$bijbij, terms$Sbb, terms$bdotdot, terms$distX, terms$distY, terms$aidot, terms$bidot,  n, dobb3 = TRUE)
  } else {
    moms <- NULL
  }
    
  
  output$pvalue <-  testfunc(dcov2 = dcov2, terms = terms, moms = moms, n = n, smp = perms)
  
  class(output) <- "dctest"
  
  output$call <- match.call()
  output$method <- method
  output$affine <- affine
  output$bias.corr <- bias.corr
  output$standardize <- standardize
  output$metr.X <- metr.X
  output$metr.Y <- metr.Y
  output$b <- b
  
 # if (dogamma) 
  #  warning("The simple gamma approximation can be anticonservative, in particular for small p-values.")
  
  return(output)
}




calcmom  <- function(aijaij, Saa, adotdot, bijbij = NULL, Sbb = NULL, bdotdot = NULL, distX = NULL, distY =NULL,  aidot = NULL, bidot =NULL,  n, dobb3) {
  
  n2 <- n ^ 2
  n3 <- n ^ 3
  n4 <- n ^ 4
  
  
  vc.C <- c(n * (n - 1) * (n - 2) * (n - 3),
            2 * n * (n - 1),
            4 * n * (n - 1) * (n - 2),
            n * (n - 1),
            4 * n * (n - 1),
            2 * n * (n - 1) * (n - 2),
            n)
  
  
  vc.b <- c(6 * n + 2 * n2,
            6 * n - 2 * n2 - 2 * n3 + n4,
            6 * n - n3,
            6 * n - 2 * n2,
            6 * n - 4 * n2,
            6 * n,
            6 * n - 10 * n2 + 4 * n3)
  
  vc.c <- c(-24 * n - 4 * n2,
            -24 * n + 12 * n2 + 4 * n3 - 2 * n4,
            -24 * n + 4 * n2 + 2 * n3,
            -24 * n + 12 * n2,
            -24 * n + 20 * n2 - 4 * n3,
            -24 * n + 4 * n2,
            -24 * n + 44 * n2 - 24 * n3 + 4 * n4)
  
  
  vc.d <- c(18 * n + 3 * n2,
            18 * n - 9 * n2 - 2 * n3 + n4,
            18 * n - 3 * n2 -  n3,
            18 * n - 9 * n2 - 2 * n3 + n4,
            18 * n - 15 * n2 + 3 * n3,
            18 * n - 3 * n2 -  n3,
            18 * n - 33 * n2 + 18 * n3 - 3 * n4)
  
  biX <- aijaij / n / (n - 1)
  
  ciX <- (Saa - aijaij) / n / (n - 1) / (n - 2)
  
  diX <- (adotdot ^ 2 + 2 * aijaij - 4 * Saa) / n / (n - 1) / (n - 2) / (n - 3)
  
  vc.X <- vc.b * biX + vc.c * ciX + vc.d * diX
  
  if (!dobb3) {
    if (is.null(bijbij)) {
      return(sqrt(vc.C) * vc.X)
    }  else {
      biY <- bijbij / n / (n - 1)
      ciY <- (Sbb - bijbij) / n / (n - 1) / (n - 2)
      diY <- (bdotdot ^ 2 + 2 * bijbij - 4 * Sbb) / n / (n - 1) / (n - 2) / (n - 3)
      vc.Y <- vc.b * biY + vc.c * ciY + vc.d * diY
      return(vc.C * vc.X * vc.Y)
    }
  } else {
    vc.X.lim <- biX - 2 * ciX + diX
    est.var.X.lim <- sqrt(2) * vc.X.lim
    distXmatp2 <- squaresym(distX)
    amatp2dot <- Rfast::colsums(distXmatp2)
    distX2 <- distX * distX
    a2dot <- Rfast::colsums(distX2)
    
    BiX.C.C <- sum(amatp2dot * aidot)
    BiX.C.H <- matrix_prod_sum(distXmatp2, distX)
    BiX.H.C <- sum(a2dot * aidot)
    BiX.H.H <- sum_hadamard_power3(distX)
    BiX.H <- aijaij
    BiX.C <- Saa
    CS3X <- sum(aidot ^ 3)
    
    eiX <- BiX.C.H / n / (n - 1) / (n - 2)
    fiX <- (BiX.C.C - BiX.C.H - 2 * BiX.H.C + BiX.H.H) / n / (n - 1) / (n - 2) / (n - 3)
    yiX <- (BiX.C * adotdot - BiX.H * adotdot - 2 * CS3X - 4 * BiX.H.H -
              4 * BiX.C.C + 2 * BiX.C.H + 10 * BiX.H.C) /
      n / (n - 1) / (n - 2) / (n - 3) / (n - 4)
    uiX <- (adotdot ^ 3 + 16 * BiX.H.H - 48 * BiX.H.C - 8 * BiX.C.H +
              6 * adotdot * BiX.H + 24 * BiX.C.C + 16 * CS3X - 12 * BiX.C * adotdot) /
      n / (n - 1) / (n - 2) / (n - 3) / (n - 4) / (n - 5)
    
    
    est.m3cent.X <- - 1 * eiX + 3 * fiX - 3 * yiX + uiX
    
    est.m3cent.X <- sqrt(8) * est.m3cent.X
    
    est.skw.X <- est.m3cent.X / est.var.X.lim ^ (3 / 2)
    
    if (!is.na(est.skw.X)) {
      if (est.skw.X < 0)
        est.skw.X <- 1e-3
    }
    
    
    if (is.null(bijbij)) { 
      return(list("vc" = sqrt(vc.C) * vc.X, "skw" = est.skw.X))
    } else {
      
      biY <- bijbij / n / (n - 1)
      ciY <- (Sbb - bijbij) / n / (n - 1) / (n - 2)
      diY <- (bdotdot ^ 2 + 2 * bijbij - 4 * Sbb) / n / (n - 1) / (n - 2) / (n - 3)
      vc.Y <- vc.b * biY + vc.c * ciY + vc.d * diY
      
      vc.Y.lim <- biY - 2 * ciY + diY
      est.var.Y.lim <- sqrt(2) * vc.Y.lim
      distYmatp2 <- squaresym(distY)
      bmatp2dot <- Rfast::colsums(distYmatp2)
      distY2 <- distY * distY
      b2dot <- Rfast::colsums(distY2)
      
      BiY.C.C <- sum(bmatp2dot * bidot)
      BiY.C.H <- matrix_prod_sum(distYmatp2, distY)
      BiY.H.C <- sum(b2dot * bidot)
      BiY.H.H <- sum_hadamard_power3(distY)
      BiY.H <- bijbij
      BiY.C <- Sbb
      CS3Y <- sum(bidot ^ 3)
      
      eiY <- BiY.C.H / n / (n - 1) / (n - 2)
      fiY <- (BiY.C.C - BiY.C.H - 2 * BiY.H.C + BiY.H.H) / n / (n - 1) / (n - 2) / (n - 3)
      yiY <- (BiY.C * bdotdot - BiY.H * bdotdot - 2 * CS3Y - 4 * BiY.H.H -
                4 * BiY.C.C + 2 * BiY.C.H + 10 * BiY.H.C) /
        n / (n - 1) / (n - 2) / (n - 3) / (n - 4)
      uiY <- (bdotdot ^ 3 + 16 * BiY.H.H - 48 * BiY.H.C - 8 * BiY.C.H +
                6 * bdotdot * BiY.H + 24 * BiY.C.C + 16 * CS3Y - 12 * BiY.C * bdotdot) /
        n / (n - 1) / (n - 2) / (n - 3) / (n - 4) / (n - 5)
      
      
      est.m3cent.Y <- - 1 * eiY + 3 * fiY - 3 * yiY + uiY
      
      est.m3cent.Y <- sqrt(8) * est.m3cent.Y
      
      est.skw.Y <- est.m3cent.Y / est.var.Y.lim ^ (3 / 2)
      
      if (!is.na(est.skw.Y)) {
        if (est.skw.Y < 0)
          est.skw.Y <- 1e-3
      }
      
      return(list("vc" = vc.C * vc.X * vc.Y, "skw" = est.skw.X * est.skw.Y))
    }
  }
}
