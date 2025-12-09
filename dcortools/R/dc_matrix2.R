#' Calculates distance covariance and distance correlation matrices
#'
#' @param X A data.frame or matrix.
#' @param Y Either NULL or a data.frame or a matrix with the same number of rows as X. If only X is provided, distance covariances/correlations are calculated between all groups in X. If X and Y are provided, distance covariances/correlations are calculated between all groups in X and all groups of Y.
#' @param calc.dcov logical; specifies if the distance covariance matrix is calculated.
#' @param calc.dcor logical; specifies if the distance correlation matrix is calculated.
#' @param calc.cor If set as "pearson", "spearman" or "kendall", a corresponding correlation matrix is additionally calculated.
#' @param calc.pvalue.cor logical; IF TRUE, a p-value based on the Pearson or Spearman correlation matrix is calculated (not implemented for calc.cor ="kendall") using Hmisc::rcorr.
#' @param return.data logical; specifies if the dcmatrix object should contain the original data.
#' @param test specifies the type of test that is performed, "permutation" performs a Monte Carlo Permutation test. "gamma" performs a test based on a gamma approximation of the test statistic under the null. "conservative" performs a conservative two-moment approximation. "bb3" performs a quite precise three-moment approximation and is recommended when computation time is not an issue.
#' @param adjustp If setting this parameter to "holm", "hochberg", "hommel", "bonferroni", "BH", "BY" or "fdr", corresponding adjusted p-values are additionally returned for the distance covariance test.
#' @param b specifies the number of random permutations used for the permutation test. Ignored for all other tests.
#' @param affine logical; indicates if the affinely transformed distance covariance should be calculated or not.
#' @param standardize specifies if data should be standardized dividing each component by its standard deviations. No effect when affine = TRUE.
#' @param bias.corr logical; specifies if the bias corrected version of the sample distance covariance \insertCite{huo2016fast}{dcortools} should be calculated.
#' @param group.X A vector, each entry specifying the group membership of the respective column in X. Each group is handled as one sample for calculating the distance covariance/correlation matrices. If NULL, every sample is handled as an individual group.
#' @param group.Y A vector, each entry specifying the group membership of the respective column in Y. Each group is handled as one sample for calculating the distance covariance/correlation matrices. If NULL, every sample is handled as an individual group.
#' @param metr.X Either a single metric or a list providing a metric for each group in X (see examples).
#' @param metr.Y see metr.X.
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
#' "memsave" is typically very inefficient for dcmatrix and should only be applied in exceptional cases.
#' 
#' @param fc.discrete logical; If TRUE, "discrete" metric is applied automatically on samples of type "factor" or "character".
#' @param calc.dcov.pw logical; If TRUE, a distance covariance matrix between the univariate observations/columns is additionally calculated. Not meaningful if group.X and group.Y are not specified.
#' @param calc.dcor.pw logical; If TRUE, a distance correlation matrix between the univariate observations/columns is additionally calculated. Not meaningful if group.X and group.Y are not specified.
#' @param test.pw specifies a test (see argument "test") that is performed between all single observations.
#' @param metr.pw.X Either a single metric or a list providing a metric for each single observation/column in X (see metr.X).
#' @param metr.pw.Y See metr.pw.Y.
#' @return S3 object of class "dcmatrix" with the following components
#' \item{name X, Y}{description original data (if return.data = TRUE).} 
#' \item{name dcov, dcor}{distance covariance/correlation matrices between the groups specified in group.X/group.Y (if calc.dcov/calc.dcor = TRUE).} 
#' \item{name corr}{correlation matrix between the univariate observations/columns (if cal.cor is "pearson", "spearman" or "kendall").}
#' \item{name pvalue}{matrix of p-values based on a corresponding distance covariance test based on the entries in dcov (if argument test is not "none").} 
#' \item{name pvalue.adj}{matrix of p-values adjusted for multiple comparisons using the method specified in argument adjustp.} 
#' \item{name pvalue.cor}{matrix of pvalues based on "pearson"/"spearman" correlation (if calc.cor is "pearson" or "spearman" and calc.pvalue.cor = TRUE).}
#' \item{name dcov.pw,dcor.pw}{distance covariance/correlation matrices between the univariate observations (if calc.dcov.pw/calc.dcor.pw = TRUE.)} 
#' \item{name pvalue.pw}{matrix of p-values based on a corresponding distance covariance test based on the entries in dcov.pw (if argument test is not "none").} 
#' @export
#' @references
#' \insertRef{berschneider2018complex}{dcortools}
#' 
#' \insertRef{bottcher2017detecting}{dcortools}
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
#'@examples
#' X <- matrix(rnorm(1000), ncol = 10)
#' 
#' dcm <- dcmatrix(X, test="bb3",calc.cor = "pearson",
#'  calc.pvalue.cor = TRUE, adjustp = "BH") 
#'  
#' dcm <- dcmatrix(X, test="bb3",calc.cor = "pearson", 
#'  calc.pvalue.cor = TRUE, adjustp = "BH", 
#'  group.X = c(rep(1, 5), rep(2, 5)), 
#'  calc.dcor.pw = TRUE, test.pw = "bb3")
#' 
#' 
#' Y <- matrix(rnorm(600), ncol = 6)
#' 
#' Y[,6] <- rbinom(100, 4, 0.3)
#' 
#' dcm <- dcmatrix(X, Y, test="bb3",calc.cor = "pearson",
#'  calc.pvalue.cor = TRUE, adjustp = "BH")
#'   
#' dcm <- dcmatrix(X, Y, test="bb3",calc.cor = "pearson",
#'  calc.pvalue.cor = TRUE, adjustp = "BH",
#'  group.X = c(rep("group1", 5), rep("group2", 5)),
#'  group.Y = c(rep("group1", 5), "group2"), 
#'  metr.X = "gaussauto",
#'  metr.Y = list("group1" = "gaussauto", "group2" = "discrete"))
dcmatrix <- function (X,
                      Y = NULL,
                      calc.dcov = TRUE,
                      calc.dcor = TRUE,
                      calc.cor = "none",
                      calc.pvalue.cor = FALSE,
                      return.data = TRUE,
                      test = "none",
                      adjustp = "none",
                      b = 499,
                      affine = FALSE,
                      standardize = FALSE,
                      bias.corr = FALSE,
                      group.X = NULL,
                      group.Y = NULL,
                      metr.X = "euclidean",
                      metr.Y = "euclidean",
                      use="all",
                      algorithm ="auto",
                      fc.discrete = FALSE,
                      calc.dcor.pw = FALSE,
                      calc.dcov.pw = FALSE,
                      test.pw = "none",
                      metr.pw.X = "euclidean",
                      metr.pw.Y = "euclidean"
                      ) {
  
  output <- .dcmatrixmain(X, Y, calc.dcov, calc.dcor, calc.cor, calc.pvalue.cor, return.data, test, adjustp, b, affine, standardize, bias.corr, group.X, group.Y, metr.X, metr.Y, use, algorithm, fc.discrete)
  output$call <- match.call()
  
  if (calc.dcor.pw | calc.dcov.pw | test.pw != "none") {
    output2 <- .dcmatrixmain(X, Y, calc.dcov = calc.dcov.pw, calc.cor = "no", calc.dcor = calc.dcor.pw, calc.pvalue.cor = FALSE, return.data, test = test.pw, adjustp, b, affine, standardize, bias.corr, group.X = NULL, group.Y = NULL, metr.X = metr.pw.X, metr.Y = metr.pw.Y, use, algorithm, fc.discrete)
    output$dcor.pw <- output2$dcor
    output$dcov.pw <- output2$dcov
    output$pvalue.pw <- output2$pvalue
    output$metr.pw.X <- metr.pw.X
    output$metr.pw.Y <- metr.pw.Y
    output$dX.pw <-  output2$dX
    output$dY.pw <-  output2$dY
    names.X <- names.Y <- colnames(X)
    if (!is.null(Y))
      names.Y <- colnames(Y)
    if (!is.null(output$dcor.pw)) {
      rownames(output$dcor.pw) <- names.X
      colnames(output$dcor.pw) <- names.Y
    }
    if (!is.null(output$dcov.pw)) {
      rownames(output$dcov.pw) <- names.X
      colnames(output$dcov.pw) <- names.Y
    }
    
    if (!is.null(output$pvalue.pw)) {
      rownames(output$pvalue.pw) <- names.X
      colnames(output$pvalue.pw) <- names.Y
    }
    
    
  }
  
  output$test.pw <- test.pw
  output$calc.dcor.pw <- calc.dcor.pw
  output$calc.dcov.pw <- calc.dcov.pw
  
  return(output)
  
}


.dcmatrixmain <- function (X,
                           Y = NULL,
                           calc.dcov = TRUE,
                           calc.dcor = TRUE,
                           calc.cor = "none",
                           calc.pvalue.cor = FALSE,
                           return.data = TRUE,
                           test = "none",
                           adjustp = "none",
                           b = 499,
                           affine = FALSE,
                           standardize=FALSE,
                           bias.corr = FALSE,
                           group.X = NULL,
                           group.Y = NULL,
                           metr.X = "euclidean",
                           metr.Y = "euclidean",
                           use="all",
                           algorithm ="auto",
                           fc.discrete = FALSE) {
  output <- list()
  
  
  if(return.data) {
    output$X <- X
    output$Y <- Y
  } else {
    output <- NULL
  }
  
  
  withY <- ifelse(is.null(Y), FALSE, TRUE)
  
  
  dogamma <- docons <- dobb3 <- doperm <- donotest <- FALSE
  
  if (test == "none")
    donotest <- TRUE 
  else if (test == "gamma")
    dogamma <- TRUE 
  else if (test == "conservative")
    docons <- TRUE
  else if (test == "bb3")
    dobb3 <- TRUE
  else if (test == "permutation")
    doperm <- TRUE
  else
    stop ("Test must be one of \"none\", \"permutation\", \"gamma\", \"bb3\" or \"conservative\"")
  
  
  use.all <- use.pw <- FALSE
  
  if (use == "complete.obs") {
    cc <- which(stats::complete.cases(X))
    if (withY) {
      cc <- intersect(cc,which(stats::complete.cases(Y)))
      Y <- Y[cc,]
    }
    X <- X[cc,]
    use.all <- TRUE
  } else if (use == "all") {
    use.all <- TRUE
  } else if (use == "pairwise.complete.obs") {
    use.pw <- TRUE
  } else {
    stop("use must be one of \"all\", \"complete.obs\" or \"pairwise.complete.obs\"")
  }
  
  
  
  
  if (is.vector(X)) {
    X <- as.matrix(X)
  }
  
  p <- ncol(X)
  n <- nrow(X)
  
  if (is.null(group.X)) {
    if (is.null(colnames(X)))
     group.X <- 1:p
     else
     group.X <- colnames(X)
  } 
  
  names.X <- names.Y <-  unique(group.X)
  dX <- dY <- length(unique(group.X))
  ms.grpX <- ms.grpY <- NULL
  groupslistX <- lapply(1:dX, function(t) which(group.X == names.X[t]))
  pX <- sapply(1:dX, function(t) length(groupslistX[[t]]))
  prepX <- as.list(rep(NA,dX))
  dvarX <- rep(NA,dX)
    
    
  # 
  # tblX <- table(group.X)
  # labelsX <- names.X <- names.Y <- as.factor(names(tblX))
  # #if (sum(group.X - 1:p) == 0)
  #  # names.X <- names.Y <- colnames(X)
  # pX <- as.numeric(tblX)
  # dX <- dY <- length(pX)
  # ms.grpX <- ms.grpY <- NULL
  # groupslistX <- lapply(1:dX, function(t) which(group.X == labelsX[t]))
  #prepX <- as.list(rep(NA,dX))
  #dvarX <- rep(NA,dX)
  

  
  lmX <- length(metr.X)
  if (lmX ==1) {
    metr.X <- as.list(replicate(dX, metr.X))
  } else if (lmX == 2) {
    ischar <- suppressWarnings(is.na(as.numeric(metr.X[2])))
    if (!ischar)
      metr.X <- lapply(1:dX, function(u) metr.X)
  }
  
  if (is.character(metr.X) & lmX == dX) {
    ischar <- suppressWarnings(is.na(as.numeric(metr.X[2])))
    if (ischar)
      metr.X <- as.list(metr.X)
  }
  
  if (use.all) {
    ms.X <- sapply(1:dX, function(t) any(!stats::complete.cases(X[,t])))
    ms.grpX <- which(ms.X)
  } else {
    ms.X <- rep(FALSE,dX)
    ms.grpX <- numeric(0)
  }
  
  ## normalize samples if calculation of affinely invariant distance covariance is desired
  if (affine) {
    for (j in 1 : dX) {
      if (use.all) {
        X[,groupslistX[[j]]] <- normalize.sample(X[,groupslistX[[j]]], n, pX[j])
      } else {
        cc <- stats::complete.cases(X[,groupslistX[[j]]])
        ncc <- length(cc)
        X[cc,groupslistX[[j]]] <- normalize.sample(X[cc,groupslistX[[j]]], n, pX[j])
      }
    }
  } else if (standardize) {
    if (use.all) {
      X[,groupslistX[[j]]] <- scale.sample(X[,groupslistX[[j]]], n, pX[j])
    } else {
      cc <- stats::complete.cases(X[,groupslistX[[j]]])
      ncc <- length(cc)
      X[cc,groupslistX[[j]]] <- scale.sample(X[cc,groupslistX[[j]]], n, pX[j])
    }
  }
  
  
  
  
  if (withY) {
    
    lmY <- length(metr.Y)
    
    
    if (is.vector(Y)) {
      Y <- as.matrix(Y)
    }
    
    q <- ncol(Y)
    m <- nrow(Y)
    
    if (is.null(group.Y)) {
      if (is.null(colnames(Y)))
        group.Y <- 1 : q
        else
        group.Y <- colnames(Y)
    }
   #   names.Y <- colnames(Y)
   # } else 
   #   names.Y <- unique(group.Y)
    
    names.Y <- unique(group.Y)
    dY <- length(unique(group.Y))
    groupslistY <- lapply(1:dY, function(t) which(group.Y == names.Y[t]))
    pY <- sapply(1:dY, function(t) length(groupslistY[[t]]))
    prepY <- as.list(rep(NA,dY))
    dvarY <- rep(NA,dY)
    
    #tblY <- table(group.Y)
    #labelsY <- as.factor(names(tblY))
    #pY <- as.numeric(tblY)
    #dY <- length(pY)
    #groupslistY <- lapply(1:dY, function(t) which(group.Y == labelsY[t]))
    #prepY <- as.list(rep(NA,dY))
    #dvarY <- rep(NA,dY)
    if (use.all) {
      ms.Y <- sapply(1:dY, function(t) any(!stats::complete.cases(Y[,t])))
      ms.grpY <- which(sapply(1:dY, function(t) any(!stats::complete.cases(Y[,t]))))
    } else {
      ms.Y <- rep(FALSE,dY)
      ms.grpY <- numeric(0)
    }
    
    if (lmY ==1) {
      metr.Y <- as.list(replicate(dY, metr.Y))
    } else if (lmX == 2) {
      ischar <- suppressWarnings(is.na(as.numeric(metr.Y[2])))
      if (!ischar)
        metr.Y <- lapply(1:dY, function(u) metr.Y)
    }
    
    if (is.character(metr.Y) & lmY == dY) {
      ischar <- suppressWarnings(is.na(as.numeric(metr.Y[2])))
      if (ischar)
        metr.Y <- as.list(metr.Y)
    }
    
    if (m != n) 
      stop("X and Y must have same number of rows (samples)")
    
    if (affine) {
      for (j in 1 : dY) {
        if (use.all) {
          Y[,groupslistY[[j]]] <- normalize.sample(Y[,groupslistY[[j]]], n, pY[j])
        } else {
          cc <- stats::complete.cases(Y[,groupslistY[[j]]])
          ncc <- length(cc)
          Y[cc,groupslistX[[j]]] <- normalize.sample(Y[cc,groupslistX[[j]]], ncc, pY[j])
        }
      }
    } else if (standardize) {
      if (use.all) {
        Y[,groupslistY[[j]]] <- scale.sample(Y[,groupslistY[[j]]], n, pY[j])
      } else {
        cc <- stats::complete.cases(X[,groupslistX[[j]]])
        ncc <- length(cc)
        Y[cc,groupslistY[[j]]] <- scale.sample(Y[cc,groupslistY[[j]]], n, pY[j])
      }
    }
    
  }
  
  
  
  
  if (algorithm == "auto") {
    gofast <- (((p == length(names.X))) * (n>200)) & (!dobb3) * all(metr.X %in% c("euclidean", "discrete"))
    if (withY) 
      gofast <- gofast * (q == length(names.Y)) * all(metr.Y %in% c("euclidean", "discrete"))
    if (gofast) {
      algorithm <- "fast"
    } else {
      algorithm <- "standard"
    }
  }
  
  alg.fast <- alg.standard <- alg.memsave <- FALSE
  
  if (algorithm == "fast") {
    alg.fast <- TRUE 
    if (doperm) 
      terms.smp <- function(terms, smp, ndisc = NULL) {sampleterms.fast.matr(terms, smp)}
  } else if (algorithm == "standard") {
    alg.standard <- TRUE
    if (doperm) 
      terms.smp <- function(terms, smp, ndisc = NULL) {sampleterms.standard(terms, smp)}
  }  else if (algorithm == "memsave") {
    alg.memsave <- TRUE
    if (doperm) 
      terms.smp <- function(terms, smp, ndisc = NULL) {sampleterms.memsave(terms, smp)}
  } else
    stop ("Algorithm must be one of \"fast\", \"standard\", \"memsave\" or \"auto\"")
  
  
  if (!alg.standard & dobb3) 
    stop("bb3 p-value calculation is only possible with algorithm=\"standard\"!")
  
  
  if (bias.corr == TRUE) {
    termstodcov2 <- function(aijbij,Sab,Tab,n) {
      aijbij/ n / (n - 3) - 2 * Sab / n / (n - 2) / (n - 3) + Tab / n / (n - 1) / (n - 2) / (n - 3) 
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
    testfunc <- function(terms, moms.X, moms.Y, dcov2, smp, ...) {
      n <- terms$ncc
      Saa <- vector_prod_sum(terms$aidot,terms$aidot)
      Sbb <- vector_prod_sum(terms$bidot,terms$bidot)
      Sab <- vector_prod_sum(terms$aidot,terms$bidot)
      dvarX <- terms$aijaij / n / (n - 3) - 2 * Saa/ n / (n - 2) / (n - 3) + terms$adotdot * terms$adotdot / n / (n - 1) / (n - 2) / (n - 3) 
      dvarY <- terms$bijbij / n / (n - 3) - 2 * Sbb / n / (n - 2) / (n - 3) + terms$bdotdot * terms$bdotdot / n / (n - 1) / (n - 2) / (n - 3) 
      dcov2 <- terms$aijbij / n / (n - 3) - 2 * Sab / n / (n - 2) / (n - 3) + terms$adotdot * terms$bdotdot / n / (n - 1) / (n - 2) / (n - 3) 
      U1 <- dvarX  * dvarY
      U2 <- terms$adotdot / n / (n - 1)
      U3 <- terms$bdotdot / n / (n - 1)
      alph <- 1 / 2 * (U2 ^ 2 * U3 ^ 2) / U1
      beta <- 1 / 2 * (U2 * U3) / U1
      stat <- n *  dcov2 + U2 * U3
      pval <- stats::pgamma(stat, alph, beta, lower.tail = FALSE) 
      return(pval)
    }
  } else if (doperm) {
    testfunc <- function(terms, moms.X, moms.Y, dcov2, smp,  ...) {
      n <- terms$ncc
      if (is.na(dcov2))
        return(NA)
      Tab <- terms$adotdot * terms$bdotdot
      reps <- lapply(1:b, function(t) {
        terms.sample <- terms.smp(terms,smp[[t]])
        return(termstodcov2(terms.sample$aijbij, terms.sample$Sab, Tab, n))
      })
      pval <- (1 + length(which(reps >= dcov2))) / (1 + b)
      return(pval)
    }
  } else if (docons) {
    testfunc <- function(terms, moms.X, moms.Y, dcov2, smp,...) {
      n <- terms$ncc
      est.m2 <- sum((moms.X * moms.Y)) / n ^ 10
      est.m1 <- terms$adotdot * terms$bdotdot / n ^ 3 / (n - 1)
      est.var <- (est.m2 - est.m1 ^ 2)
      alpha <- sqrt(est.var / 2 / est.m1 ^ 2)
      stat <- terms$aijbij / n - 2 * vector_prod_sum(terms$aidot,terms$bidot) / n ^ 2 + terms$adotdot * terms$bdotdot / n ^ 3
      pval <- stats::pchisq(stat * sqrt(2) / sqrt(est.var), df = 1 / alpha, lower.tail = FALSE)  
      return(pval)
    }
  } else if (dobb3) {
    testfunc <- function(terms, moms.X, moms.Y, dcov2, smp, ...) {
      n <- terms$ncc
      est.m2 <- sum((moms.X$vc * moms.Y$vc)) / n ^ 10
      est.m1 <- terms$adotdot * terms$bdotdot / n ^ 3 / (n - 1)
      est.var <- (est.m2 - est.m1 ^ 2)
      est.skw <- moms.X$skw * moms.Y$skw
      beta <- est.skw / sqrt(8)
      stat <- terms$aijbij / n - 2 * vector_prod_sum(terms$aidot,terms$bidot) / n ^ 2 + terms$adotdot * terms$bdotdot / n ^ 3
      centstat <- (stat - est.m1) /  sqrt(est.var)
      pval <- stats::pchisq((centstat * sqrt(2) + 1 / beta) / beta , df = 1 / beta ^ 2, lower.tail = FALSE)  
      return(pval)
    } 
  } else if (donotest) {
    testfunc <- function(terms, moms.X, moms.Y, dcov2, smp, ...) {}
  }
  
  if (!calc.dcov) {
    dcov2todcov <- function(dcov2) {}
  }
  
  if (!calc.dcor) {
    dcov2todcor <- function(dcov2, dvarX, dvarY) {}
  }
  
  if (doperm & use.all) {
    perms <- lapply(1:b, function(t) sample(1:n))
  } else {
    perms <- NULL
  }
  
  extendoutput <- doperm| ((dobb3|docons)*use.pw)
  
  
  
  if (fc.discrete) {
    for (j in 1:dX) {
      if (is.factor(X[,groupslistX[[j]]]) | is.character(X[,groupslistX[[j]]]))
        metr.X[[j]] <- "discrete"
    }
    if (withY) {
      for (j in 1:dY) {
        if (is.factor(Y[,groupslistY[[j]]]) | is.character(Y[,groupslistY[[j]]]))
          metr.Y[[j]] <- "discrete"
      }
    }
  }
  
  
  
  if (calc.cor %in% c("spearman","kendall", "pearson")) {
    output$corr <- stats::cor(X,Y, use = use, method = calc.cor)
    if (calc.pvalue.cor) {
      if (calc.cor %in% c("spearman", "pearson")) {
        if (!withY) {
          corrp <- Hmisc::rcorr(X, type = calc.cor)
          output$pvalue.cor <- corrp$P
          diag(output$pvalue.cor) <- 0
          if (use.all)
            output$pvalue.cor[which(corrp$n<n,arr.ind=TRUE)] <- NA
        } else {
          corrp <- Hmisc::rcorr(X,Y)
          output$pvalue.cor <- corrp$P[1:dX,(dX+1):(dX+dY)]
          if (use.all)
            output$pvalue.cor[which(corrp$n[1:dX,(dX+1):(dX+dY)]<n,arr.ind=TRUE)] <- NA
        }  
      } else 
        warning("P-Value calculation for Kendall correlation not implemented")
    }
  }
  
  
  
  if (alg.fast) {
    discrete.X <- (metr.X == "discrete")
    if (withY)
      discrete.Y <- (metr.Y == "discrete")
  }
  
  
  
  if (calc.dcov) {
    output$dcov <- matrix(nrow = dX, ncol = dY)
    rownames(output$dcov) <- names.X
    colnames(output$dcov) <- names.Y
  }
  
  if (calc.dcor) {
    output$dcor <- matrix(nrow = dX, ncol = dY)
    rownames(output$dcor) <- names.X
    colnames(output$dcor) <- names.Y
    if (!withY)
      diag(output$dcor) <- 1
  }
  
  if (!donotest) {
    output$pvalue <- matrix(nrow = dX, ncol = dY)
    rownames(output$pvalue) <- names.X
    colnames(output$pvalue) <- names.Y
    if (!withY)
      diag(output$pvalue) <- 0
  }
  
  momsX <- momsY <- NULL
  
  if ((docons | dobb3) & !use.pw) {
    momsX <- as.list(rep(NA,dX))
    if (withY)
      momsY <- as.list(rep(NA,dY))
  }
  
  
  
  for (j in setdiff(1:dX,ms.grpX)) {
    if (alg.fast) {
      prepX[[j]] <- prep.fast(X[,j], n, discrete = discrete.X[j], pairwise = use.pw)
    } else if (alg.memsave) {
      prepX[[j]] <- prep.memsave(X[,groupslistX[[j]]], n, pX[j],  metr.X = metr.X[[j]], pairwise = use.pw)
    } else if (alg.standard) {
      prepX[[j]] <- prep.standard(X[,groupslistX[[j]]], n, pX[j],  metr.X = metr.X[[j]], pairwise = use.pw)
    }  
    
    Saa <- vector_prod_sum(prepX[[j]]$aidot,prepX[[j]]$aidot)
    
    if ((docons | dobb3) & !use.pw) {
      momsX[[j]] <- calcmom(aijaij = prepX[[j]]$aijaij, Saa = Saa, adotdot = prepX[[j]]$adotdot, aidot = prepX[[j]]$aidot, distX = prepX[[j]]$distX, n = n, dobb3 = dobb3)
    }
    dvarX[j] <- termstodcov2(prepX[[j]]$aijaij, Saa, prepX[[j]]$adotdot*prepX[[j]]$adotdot, prepX[[j]]$ncc) 
  }
  
  if (!withY & calc.dcov) 
    diag(output$dcov) <- sqrt(dvarX)
  
  
  if (withY) {
    for (j in setdiff(1:dY,ms.grpY)) {
      if (alg.fast) {
        prepY[[j]] <- prep.fast(Y[,j], n, discrete = discrete.Y[j], pairwise = use.pw)
      } else if (alg.memsave) {
        prepY[[j]] <- prep.memsave(Y[,groupslistY[[j]]], n, pY[j], metr.X = metr.Y[[j]], pairwise = use.pw)
      } else if (alg.standard) {
        prepY[[j]] <- prep.standard(Y[,groupslistY[[j]]], n, pY[j], metr.X = metr.Y[[j]], pairwise = use.pw)
      }  
      Sbb <- vector_prod_sum(prepY[[j]]$aidot, prepY[[j]]$aidot)
      
      if ((docons | dobb3) & !use.pw) {
        momsY[[j]] <- calcmom(aijaij = prepY[[j]]$aijaij, Saa = Sbb, adotdot = prepY[[j]]$adotdot, aidot = prepY[[j]]$aidot, distX = prepY[[j]]$distX, n = n, dobb3 = dobb3)
      }
      dvarY[j] <- termstodcov2(prepY[[j]]$aijaij, Sbb, prepY[[j]]$adotdot*prepY[[j]]$adotdot, prepY[[j]]$ncc) 
    }
  }
  
  
  if (!withY) {
    if (dX > 1) {
      for (i in setdiff(1:(dX-1),ms.grpX)) {
        for (j in setdiff((i+1):dX,ms.grpX)) {
          if (alg.fast) {
            terms <- preptoterms.fast(prepX[[i]], prepX[[j]], n, pairwise = use.pw, discrete.X[[i]], discrete.X[[j]], perm = extendoutput)
          } else if (alg.memsave) {
            terms <- preptoterms.memsave(prepX[[i]], prepX[[j]], metr.X[[i]], metr.X[[j]], n, pairwise = use.pw, perm = extendoutput) 
          } else if (alg.standard) {
            terms <- preptoterms.standard(prepX[[i]], prepX[[j]], n, pairwise = use.pw, perm = extendoutput)
          }  
          dcov2XY <- termstodcov2(terms$aijbij, vector_prod_sum(terms$aidot, terms$bidot), terms$adotdot * terms$bdotdot, terms$ncc)
          output$dcov[i,j] <- output$dcov[j,i] <- dcov2todcov(dcov2 = dcov2XY)
          if (use.pw) {
            Saa <- vector_prod_sum(terms$aidot, terms$aidot)
            Sbb <- vector_prod_sum(terms$bidot, terms$bidot)
            dvX <-  termstodcov2(terms$aijaij, Saa, terms$adotdot * terms$adotdot, terms$ncc)
            dvY <-  termstodcov2(terms$bijbij, Sbb, terms$bdotdot * terms$bdotdot, terms$ncc)
            if (docons | dobb3) {
              moms.X <- calcmom(aijaij = terms$aijaij, Saa = Saa, adotdot = terms$adotdot, distX = terms$distX, n = terms$ncc, aidot = terms$aidot, dobb3 = dobb3)
              moms.Y <- calcmom(aijaij = terms$bijbij, Saa = Sbb, adotdot = terms$bdotdot, distX = terms$distY, n = terms$ncc, aidot = terms$bidot, dobb3 = dobb3)
            }
            if (doperm) {
              perms <- lapply(1:b, function(t) sample(1:terms$ncc))
            }
            
          } else {
            dvX <- dvarX[i]
            dvY <- dvarX[j]
            if (docons | dobb3) {
              moms.X <- momsX[[i]]
              moms.Y <- momsX[[j]]
            }
          }
          
          
          
          output$dcor[i,j] <- output$dcor[j,i] <- dcov2todcor(dcov2 = dcov2XY, dvX, dvY)
          
          
          
          output$pvalue[i,j] <- output$pvalue[j,i] <- testfunc(dcov2 = dcov2XY, terms = terms, moms.X = moms.X, moms.Y = moms.Y, n = n, smp = perms, prepX[[i]], prepX[[j]])
        }
      }
    }
  } else {
    for (i in setdiff(1:dX,ms.grpX)) {
      for (j in setdiff(1:dY,ms.grpY)) {
        if (alg.fast) {
          terms <- preptoterms.fast(prepX[[i]], prepY[[j]], n, pairwise = use.pw, discrete.X[[i]], discrete.Y[[j]], perm = extendoutput)
        } else if (alg.memsave) {
          terms <- preptoterms.memsave(prepX[[i]], prepY[[j]], metr.X[[i]], metr.Y[[j]], n, pairwise = use.pw, perm = extendoutput) 
        } else if (alg.standard) {
          terms <- preptoterms.standard(prepX[[i]], prepY[[j]], n, pairwise = use.pw, perm = extendoutput)
        }  
        dcov2XY <- termstodcov2(terms$aijbij, vector_prod_sum(terms$aidot, terms$bidot), terms$adotdot * terms$bdotdot, terms$ncc)
        output$dcov[i,j] <- dcov2todcov(dcov2 = dcov2XY)
        if (use.pw) {
          Saa <- vector_prod_sum(terms$aidot, terms$aidot)
          Sbb <- vector_prod_sum(terms$bidot, terms$bidot)
          dvX <-  termstodcov2(terms$aijaij, Saa, terms$adotdot * terms$adotdot, terms$ncc)
          dvY <-  termstodcov2(terms$bijbij, Sbb, terms$bdotdot * terms$bdotdot, terms$ncc)
          if (docons | dobb3) {
            moms.X <- calcmom(aijaij = terms$aijaij, Saa = Saa, adotdot = terms$adotdot, distX = terms$distX, aidot = terms$aidot, n = terms$ncc, dobb3 = dobb3)
            moms.Y <- calcmom(aijaij = terms$bijbij, Saa = Sbb, adotdot = terms$bdotdot, distX = terms$distY, aidot = terms$bidot, n = terms$ncc, dobb3 = dobb3)
          }
          if (doperm) {
            perms <- lapply(1:b, function(t) sample(1:terms$ncc))
          }
        } else {
          dvX <- dvarX[i]
          dvY <- dvarY[j]
          if (docons | dobb3) {
            moms.X <- momsX[[i]]
            moms.Y <- momsY[[j]]
          }
        }
        output$dcor[i,j] <- dcov2todcor(dcov2 = dcov2XY, dvX, dvY)
        
        
        output$pvalue[i,j] <- testfunc(dcov2 = dcov2XY, terms = terms, moms.X = moms.X, moms.Y = moms.Y, smp = perms, prepX[[i]], prepY[[j]])
      }
    }
  }
  
  
  
  
  
  if (adjustp %in% c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr")) {
    if (withY) {
      output$pvalue.adj <- matrix(stats::p.adjust(output$pvalue,method = adjustp), ncol = dY)
    } else {
      ind <- which(lower.tri(output$pvalue), arr.ind=TRUE)
      pvec <- as.vector(output$pvalue[ind])
      pvec <- stats::p.adjust(pvec, method = adjustp)
      output$pvalue.adj <- diag(0,dX)
      ind2 <- ind[,2:1]
      output$pvalue.adj[ind] <- output$pvalue.adj[ind2] <- pvec
    }  
  } else if (adjustp != "none")
    warning ("adjustp should be one of \"holm\", \"hochberg\", \"hommel\", \"bonferroni\", \"BH\", \"BY\", \"fdr\" \n
               No p-value correction performed")
  
  class(output) <- "dcmatrix"
  output$withY <- withY
  output$dX <- dX
  output$n <- n
  output$b <- b
  output$test <- test
  output$calc.dcov <- calc.dcov
  output$calc.dcor <- calc.dcor
  output$bias.corr <- bias.corr
  output$affine <- affine
  output$calc.cor <- calc.cor
  output$group.X  <-  group.X
  output$names.X <-  names.X
  output$groupslistX <-  groupslistX
  
  if (withY) {
    output$group.Y <- group.Y
    output$dY <- dY
    output$names.Y <- names.Y
    output$groupslistY <- groupslistY
  } 
  
  return(output)

}