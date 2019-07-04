#' Performs a distance covariance test
#'
#' @param X contains either the first sample or its corresponding distance matrix. In the first case, this input can be either a vector of positive length, a matrix with one column or a data.frame with one column. In this case, type.X must be specified as "sample". In the second case, the input must be a distance matrix corresponding to the sample of interest. In this second case, type.X must be "distance".
#' @param Y see X.
#' @param test specifies the type of test that is performed, "permutation" performs a Monte Carlo Permutation test. "gamma" performs a test based on a gamma approximation of the test statistic under the null. "conservative" performs a conservative two-moment approximation. "bb3" performs a quite precise three-moment approximation and is recommended when computation time is not an issue.
#' @param b specifies the number of random permutations used for the permutation test. Ignored when test="gamma"
#' @param affine logical; indicates if the affinely transformed distance covariance should be calculated or not.
#' @param bias.corr logical; indicates if the bias corrected version of the sample distance covariance should be calculated.
#' @param type.X either "sample" or "distance"; specifies the type of input for X.
#' @param type.Y see type.X.
#' @param metr.X specifies the metric which should be used for X to analyse the distance covariance. One of "euclidean", "gaussian", "boundsq", "alpha", "discrete" or a user-defined metric or kernel (see examples). It is also possible to specify parameters for "gaussian", "boundsq" and "alpha" via c(metr,param), e.g. c("gaussian",2) will result in a Gaussian metric with bandwith parameter 2.
#' @param metr.Y see metr.X.
#' @param use : "all" uses all observations, "complete.obs" excludes NA's
#' @param return.data : IF TRUE, X and Y are contained in the resulting distcov.test object.
#' @param algorithm: One of "auto", "fast", "memsave" and "standard". "fast" has a lower computational complexity than "memsave" and "standard", but cannot be used with vector-valued data, other metrics than "discrete" or "euclidean" or the "bb3" test. "memsave" can be used in most settings and is slightly faster than the "standard" algorithm. Moreover, different to "standard", it can be used for very large sample size. "standard" can be used with all metrics and tests, however is the slowest algorithm and breaks down with large sample sizes. "auto" chooses an algorithm using feasibility and rules of thumb.
#' @return distcov.test object
#' @export

distcov.test <- function(X,
                         Y,
                         method = "permutation",
                         b = 499L,
                         affine = FALSE,
                         bias.corr = TRUE,
                         type.X = "sample",
                         type.Y = "sample",
                         metr.X = "euclidean",
                         metr.Y = "euclidean",
                         use = "all",
                         return.data = FALSE,
                         algorithm = "auto") {
  
  #extract dimensions and sample sizes
  ss.dimX <- extract_np(X, type.X)
  ss.dimY <- extract_np(Y, type.Y)
  
  n <- ss.dimX$Sample.Size
  p <- ss.dimX$Dimension
  
  m <- ss.dimY$Sample.Size
  q <- ss.dimY$Dimension
  
  dogamma <- docons  <- dobb3 <-  doperm <- FALSE
  
  output <- list()
  
  if (method == "gamma")
    dogamma <- TRUE 
    else if (method == "conservative")
    docons <- TRUE
    else if (method == "bb3")
    dobb3 <- TRUE
    else if (method == "permutation")
    doperm <- TRUE
    else
    stop ("Test must be one of \"permutation\", \"gamma\", \"bb3\" or \"conservative\"")
    
  
  
  
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
    testfunc <- function(terms, n, ...) {
      dvarX <- terms$aijaij / n / (n - 3) - 2 * terms$Saa / n / (n - 2) / (n - 3) + terms$adotdot * terms$adotdot / n / (n - 1) / (n - 2) / (n - 3) 
      dvarY <- terms$bijbij / n / (n - 3) - 2 * terms$Sbb / n / (n - 2) / (n - 3) + terms$bdotdot * terms$bdotdot / n / (n - 1) / (n - 2) / (n - 3) 
      dcov2 <- terms$aijbij / n / (n - 3) - 2 * terms$Sab / n / (n - 2) / (n - 3) + terms$adotdot * terms$bdotdot / n / (n - 1) / (n - 2) / (n - 3) 
      U1 <- dvarX  * dvarY
      U2 <- terms$adotdot / n / (n - 1)
      U3 <- terms$bdotdot / n / (n - 1)
      alph <- 1 / 2 * (U2 ^ 2 * U3 ^ 2) / U1
      beta <- 1 / 2 * (U2 * U3) / U1
      stat <- n *  dcov2 + U2 * U3
      pval <- pgamma(stat, alph, beta, lower.tail = FALSE) 
      return(pval)
    }
  } else if (doperm) {
    testfunc <- function(dcov2, smp, terms, n, ...) {
      if (is.na(dcov2))
        return(NA)
      Tab <- terms$adotdot * terms$bdotdot
      reps <- lapply(1:b, function(t) {
        terms.sample <- terms.smp(terms,smp[[t]])
        return(termstodcov2(terms.sample$aijbij, terms.sample$Sab, Tab, n))
      })
      pval <- (1 + length(which(reps > dcov2))) / (1 + b)
      return(pval)
    }
  } else if (docons) {
    testfunc <- function(terms, moms, n, ...) {
      est.m2 <- sum(moms) / n ^ 10
      est.m1 <- terms$adotdot * terms$bdotdot / n ^ 3 / (n - 1)
      est.var <- (est.m2 - est.m1 ^ 2)
      alpha <- sqrt(est.var / 2 / est.m1 ^ 2)
      stat <- terms$aijbij / n - 2 * terms$Sab / n ^ 2 + terms$adotdot * terms$bdotdot / n ^ 3
      pval <- pchisq(stat * sqrt(2) / sqrt(est.var), df = 1 / alpha, lower.tail = FALSE)  
      return(pval)
    }
  } else if (dobb3) {
    testfunc <- function(terms, moms, n,...) {
      est.m2 <- sum(moms$vc) / n ^ 10
      est.m1 <- terms$adotdot * terms$bdotdot / n ^ 3 / (n - 1)
      est.var <- (est.m2 - est.m1 ^ 2)
      est.skw <- moms$skw
      beta <- est.skw / sqrt(8)
      stat <- terms$aijbij / n - 2 * terms$Sab / n ^ 2 + terms$adotdot * terms$bdotdot / n ^ 3
      centstat <- (stat - est.m1) /  sqrt(est.var)
      pval <- pchisq((centstat * sqrt(2) + 1 / beta) / beta , df = 1 / beta ^ 2, lower.tail = FALSE)  
      return(pval)
    } 
  }
  
  
  
  

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
  
  
  ## normalize samples if calculation of affinely invariant distance covariance is desired
  if (affine == TRUE) {
    if (p > n | q > n) {
      stop("Affinely invariant distance covariance cannot be calculated for p>n")
    }
    if (type.X == "distance" | type.Y == "distance") {
      stop("Affinely invariant distance covariance cannot be calculated for type 'distance'")
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
    if (p == 1 & q == 1 & metr.X %in% c("euclidean", "discrete") 
        & metr.Y %in% c("euclidean", "discrete") & n > 100 & dobb3 == FALSE) {
      algorithm <- "fast"
    } else if (metr.X[1] %in% c("euclidean", "alpha", "gaussian", "boundsq", "minkowski", "discrete") &
               metr.Y[1] %in% c("euclidean", "alpha", "gaussian", "boundsq", "minkowski", "discrete") & dobb3 == FALSE) {
      algorithm <- "memsave"
    } else {
      algorithm <- "standard"
    }
  }
  
  alg.fast <- alg.standard <- alg.memsave <- FALSE
  
  if (algorithm == "fast") {
    alg.fast <- TRUE 
    if (doperm) 
      terms.smp <- function(terms, smp) {sampleterms.fast.matr(terms, smp)}
  } else if (algorithm == "standard") {
    alg.standard <- TRUE
    if (doperm) 
      terms.smp <- function(terms, smp, ndisc = NULL) {sampleterms.standard(terms, smp)}
  }  else if (algorithm == "memsave") {
    alg.memsave <- TRUE
    if (doperm) 
      terms.smp <- function(terms, smp, ndisc = NULL) {sampleterms.memsave(terms, smp)}
  } 
  else
    stop ("Algorithm must be one of \"fast\", \"standard\", \"memsave\" or \"auto\"")
  
  
  if (!alg.standard & dobb3) 
    stop("bb3 p-value calculation is only possible with algorithm=standard!")
  
  
  if (doperm) {
    perms <- lapply(1:b, function(t) sample(1:n))
  } else {
    perms <- NULL
  }
  
  if (type.X == "distance")
    metr.X <- "distance matrix provided"
  
  if (type.Y == "distance")
    metr.Y <- "distance matrix provided"
  
  
  
  if (alg.fast) {
    if (p == 1 & q == 1) {
      if (metr.X == "euclidean" & metr.Y == "euclidean") {
        terms <- dcovterms.fast(X, Y, n, calc.dcor = TRUE, calc.perm = doperm)
        if (doperm) 
          terms.smp <- function(terms, smp) {sampleterms.fast(terms, smp)}
      } else if (metr.X == "euclidean" & metr.Y == "discrete") {
        Y <- as.factor(Y)
        terms <- dcovterms.fast.numdisc(X, Y, n, calc.dcor = TRUE, calc.perm = doperm)
        if (doperm) 
          terms.smp <- function(terms, smp) {sampleterms.fast.numdisc(terms, smp)}
      } else if (metr.X == "discrete" & metr.Y == "euclidean") {
        X <- as.factor(X)
        terms <- dcovterms.fast.numdisc(Y, X, n, calc.dcor = TRUE, calc.perm = doperm)
        if (doperm) 
          terms.smp <- function(terms, smp) {sampleterms.fast.numdisc(terms, smp)}
      } else if (metr.X == "discrete" & metr.Y == "discrete") {
        X <- as.factor(X)
        Y <- as.factor(Y)
        terms <- dcovterms.fast.discrete(X ,Y, n, calc.dcor = TRUE, calc.perm = doperm)
        if (doperm) 
          terms.smp <- function(terms, smp) {sampleterms.fast.discrete(terms, smp)}
      } else {
        stop("metr.X and metr.Y have to be \"euclidean\" or \"discrete\" for fast algorithms")
      }
    } else {
      stop("Dimensions of X and Y must be 1 for fast algorithms.")
    }
  } else if (alg.memsave) {
    if (metr.X[1] %in% c("euclidean", "alpha", "gaussian", "boundsq", "minkowski", "discrete") &
        metr.Y[1] %in% c("euclidean", "alpha", "gaussian", "boundsq", "minkowski", "discrete")) {
      terms <- dcovterms.memsave(X, Y, metr.X, metr.Y, p, q, calc.dcor = TRUE, calc.perm = doperm)
      if (doperm) 
        terms.smp <- function(terms, smp) {sampleterms.memsave(terms, smp)}
    } else {
      stop("Memory efficient algorithms cannot be run with user-defined metrics")
    }
  } else if (alg.standard) {
    terms <- dcovterms.standard(X, Y, type.X, type.Y, metr.X, metr.Y, p, q, calc.dcor = TRUE, calc.bb3 = dobb3 | doperm)
    if (doperm) 
      terms.smp <- function(terms, smp) {sampleterms.standard(terms, smp)}
  }  
  
  

  
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
  output$metr.X <- metr.X
  output$metr.Y <- metr.Y
  output$b <- b
  
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
