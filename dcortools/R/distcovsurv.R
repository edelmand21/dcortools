#' Calculates an inverse-probability-of-censoring weighted (IPCW) distance correlation based on IPCW U-statistics \insertCite{datta2010inverse}{dcortools}.
#'
#' @param Y A matrix with two columns, where the first column contains the survival times and the second column the status indicators (a survival object will work).
#' @param X A vector or matrix containing the covariate information.
#' @param affine logical; specifies if X should be transformed such that the result is invariant under affine transformations of X
#' @param standardize logical; should X be standardized using the standard deviations of single observations?. No effect when affine = TRUE.
#' @param timetrafo specifies a transformation applied on the follow-up times. Can be "none", "log" or a user-specified function.
#' @param type.X For "distance", X is interpreted as a distance matrix. For "sample", X is interpreted as a sample.
#' @param metr.X specifies the metric which should be used to compute the distance matrix for X (ignored when type.X = "distance").
#' 
#'  Options are "euclidean", "discrete", "alpha", "minkowski", "gaussian", "gaussauto", "boundsq" or user-specified metrics (see examples).
#'  
#'  For "alpha", "minkowski", "gaussian", "gaussauto" and "boundsq", the corresponding parameters are specified via "c(metric,parameter)", c("gaussian",3) for example uses a Gaussian metric with bandwidth parameter 3; the default parameter is 2 for "minkowski" and "1" for all other metrics.
#' @param use specifies how to treat missing values. "complete.obs" excludes observations containing NAs, "all" uses all observations.
#' @param cutoff If provided, all survival times larger than cutoff are set to the cutoff and all corresponding status indicators are set to one. Under most circumstances, choosing a cutoff is highly recommended.
#' @return An inverse-probability of censoring weighted estimate for the distance correlation between X and the survival times.
#' @export
#' @references
#' \insertRef{bottcher2017detecting}{dcortools}
#' 
#' \insertRef{datta2010inverse}{dcortools}
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
#' survtime <- rgamma(100, abs(X))
#' cens <- rexp(100)
#' status <- as.numeric(survtime < cens)
#' time <- sapply(1:100, function(u) min(survtime[u], cens[u]))
#' surv <- cbind(time, status)
#' ipcw.dcor(surv, X)

ipcw.dcor <- function(Y, X, affine = FALSE, standardize = FALSE, timetrafo = "none", type.X = "sample", metr.X = "euclidean", use = "all", cutoff = NULL) {
  
  dcor <- ipcw.dcov(Y, X, affine, standardize, timetrafo, type.X, metr.X, use, cutoff) / sqrt(ipcw.dcov(Y, Y[,1], affine, standardize, timetrafo, type.X = "sample", metr.X, use, cutoff)) / sqrt(distsd(X =  X, affine = affine, standardize = standardize, type.X = type.X, metr.X = metr.X, use = use, bias.corr=TRUE))
  
  return(dcor)
}



#' Calculates an inverse-probability-of-censoring weighted (IPCW) distance covariance based on IPCW U-statistics \insertCite{datta2010inverse}{dcortools}.
#'
#' @param Y A column with two rows, where the first row contains the survival times and the second row the status indicators (a survival object will work).
#' @param X A vector or matrix containing the covariate information.
#' @param affine logical; indicates if X should be transformed such that the result is invariant under affine transformations of X
#' @param standardize logical; should X be standardized using the standard deviations of single observations?. No effect when affine = TRUE.
#' @param timetrafo specifies a transformation applied on the follow-up times. Can be "none", "log" or a user-specified function.
#' @param type.X For "distance", X is interpreted as a distance matrix. For "sample" (or any other value), X is interpreted as a sample
#' @param metr.X metr.X specifies the metric which should be used for X to analyze the distance covariance. Options are "euclidean", "discrete", "alpha", "minkowski", "gaussian", "gaussauto" and "boundsq". For "alpha", "minkowski", "gauss", "gaussauto" and "boundsq", the corresponding parameters are specified via "c(metric,parameter)" (see examples); the standard parameter is 2 for "minkowski" and "1" for all other metrics.
#' @param use specifies how to treat missing values. "complete.obs" excludes observations containing NAs, "all" uses all observations.
#' @param cutoff If provided, all survival times larger than cutoff are set to the cutoff and all corresponding status indicators are set to one. Under most circumstances, choosing a cutoff is highly recommended.
#' @return An inverse-probability of censoring weighted estimate for the distance covariance between X and the survival times.
#' @export
#' @references
#' \insertRef{bottcher2017detecting}{dcortools}
#' 
#' \insertRef{datta2010inverse}{dcortools}
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
#' survtime <- rgamma(100, abs(X))
#' cens <- rexp(100)
#' status <- as.numeric(survtime < cens)
#' time <- sapply(1:100, function(u) min(survtime[u], cens[u]))
#' surv <- cbind(time, status)
#' ipcw.dcov(surv, X)

ipcw.dcov <- function(Y, X, affine = FALSE, standardize = FALSE,  timetrafo = "none", type.X = "sample", metr.X = "euclidean", use = "all", cutoff = NULL) {

    # extract sample size and dimension

    ss.dimX <- extract_np(X, "sample")

    n <- ss.dimX$Sample.Size
    p <- ss.dimX$Dimension
    


    # sort data in order of follow-up times
    
    if (ncol(Y) != 2)
      stop("Y must be a matrix with two columns with time in the first column and status indicator in the second.")
    
    
    if (use == "complete.obs") {
      ccY <- stats::complete.cases(Y)
      Y <- Y[ccY,]
      n <- length(ccY)
      if (type.X == "sample" && p == 1) {
        X <- X[ccY]
      } else if (type.X == "sample" && p > 1) {
        X <- X[ccY, ]
      }
      if (type.X == "distance") {
        X <- X[ccY,ccY]
      }
    }
    
    time <- Y[,1]
    status <- Y[,2]

    IX <- Rfast::Order(time)

    time <- time[IX]
    status <- status[IX]
    
    if (timetrafo != "none") {
      if (timetrafo == "log")
        time <- log(time)
        else {
        trafo <- match.fun(timetrafo)
        time <- trafo(time)
        }
    }

    
    if (type.X == "sample" && p == 1) {
      X <- X[IX]
    } else if (type.X == "sample" && p > 1) {
      X <- X[IX, ]
    }
    if (type.X == "distance") {
      X <- X[IX,IX]
    }
    

    
    if (use == "complete.obs") {
      ccX <- which(stats::complete.cases(X))
      n <- length(ccX)
      if (type.X == "sample" && p == 1) {
        X <- X[ccX]
      } else if (type.X == "sample" && p > 1) {
        X <- X[ccX, ]
      }
      if (type.X == "distance") {
        X <- X[ccX,ccX]
      }
      time <- time[ccX]
      status <- status[ccX]
    }
    
    
    
    
    if (!is.null(cutoff)) {
      time[time>=cutoff] <- cutoff
      status[time>=cutoff] <- 1
    }
    
    events <- which(status == 1)
    
    # calculate IPC weights
    
    
    ipcw <- rep(0, n)
    ipcw[events[1]] <- 1 / (n - events[1] + 1)
    
    help <- sapply(1:n, function(x) ((n - x) / (n - x + 1)) ^ status[x])
    
    for (i in (events[1] + 1):n) {
      ipcw[i]<- prod(help[1:(i - 1)]) * (status[i] / (n - i + 1))
    }
    
    ipcw <- ipcw * n
    
    ipcw <- ipcw[events]
    time <- time[events]
    status <- status[events]
    
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
    


    

    #distance matrices

    if (type.X == "distance") {
      distXall <- X
    } else {
      distXall <- distmat(X, metr.X, p)
    }
    distX <- distXall[events, events]
    
    distT <- Rfast::Dist(time)
    k <- sum(status)

    # auxiliary terms

    m <- sum(ipcw)
    M <- sum(ipcw^2)
    W2mat <- ipcw %*% t(ipcw)
    ipcwm <- matrix(rep(ipcw, k), ncol=k)
    Wmatsum <- ipcwm + t(ipcwm)

    rowsumsX_w <- Rfast::colsums(ipcw * distX)
    rowsumsT_w <- Rfast::colsums(ipcw * distT)
    rowsumsX_w2 <- Rfast::colsums(ipcw^2 * distX)
    rowsumsT_w2 <- Rfast::colsums(ipcw^2 * distT)

    #summands of squared dcov

    term1 <- sum(W2mat * distX * distT * (m ^ 2 -m * Wmatsum - M + 2 * W2mat))
    term2 <- sum(ipcw * (m + ipcw) * rowsumsX_w * rowsumsT_w)
    term3 <- sum(ipcw * rowsumsX_w2 * rowsumsT_w)
    term4 <- sum(ipcw * rowsumsX_w * rowsumsT_w2)
    term5 <- sum(W2mat * distX) * sum(W2mat * distT)

    #calculate dcov and return

    dcov2 <- 1 / (n * (n-1) * (n-2) * (n-3)) * (term1 - 2 * (term2 - term3 - term4) + term5)

    dcov <- sign(dcov2) * sqrt(abs(dcov2))
    return(dcov)
}



#' Performs a permutation test based on the IPCW distance covariance.
#'
#' @param Y A column with two rows, where the first row contains the survival times and the second row the status indicators (a survival object will work).
#' @param X A vector or matrix containing the covariate information.
#' @param affine logical; indicates if X should be transformed such that the result is invariant under affine transformations of X.
#' @param standardize logical; should X be standardized using the standard deviations of single observations. No effect when affine = TRUE.
#' @param timetrafo specifies a transformation applied on the follow-up times. Can be "none", "log" or a user-specified function.
#' @param type.X For "distance", X is interpreted as a distance matrix. For "sample" (or any other value), X is interpreted as a sample.
#' @param metr.X metr.X specifies the metric which should be used for X to analyze the distance covariance. Options are "euclidean", "discrete", "alpha", "minkowski", "gaussian", "gaussauto" and "boundsq". For "alpha", "minkowski", "gauss", "gaussauto" and "boundsq", the corresponding parameters are specified via "c(metric,parameter)" (see examples); the standard parameter is 2 for "minkowski" and "1" for all other metrics.
#' @param use specifies how to treat missing values. "complete.obs" excludes observations containing NAs, "all" uses all observations.
#' @param cutoff If provided, all survival times larger than cutoff are set to the cutoff and all corresponding status indicators are set to one. Under most circumstances, choosing a cutoff is highly recommended.
#' @param B The number of permutations used for the permutation test
#' @return An list with two arguments, $dcov contains the IPCW distance covariance, $pvalue the corresponding p-value
#' @export
#' @references
#' \insertRef{bottcher2017detecting}{dcortools}
#' 
#' \insertRef{datta2010inverse}{dcortools}
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
#' survtime <- rgamma(100, abs(X))
#' cens <- rexp(100)
#' status <- as.numeric(survtime < cens)
#' time <- sapply(1:100, function(u) min(survtime[u], cens[u]))
#' surv <- cbind(time, status)
#' ipcw.dcov.test(surv, X)
#' ipcw.dcov.test(surv, X, cutoff = quantile(time, 0.8)) 
#' # often better performance when using a cutoff time
#' 
ipcw.dcov.test <- function(Y, X, affine = FALSE, standardize = FALSE, timetrafo = "none", type.X = "sample", metr.X = "euclidean", use = "all", cutoff = NULL, B=499)
{
    ss.dimX <- extract_np(X, "sample")

    n <- ss.dimX$Sample.Size
    p <- ss.dimX$Dimension
    
    if (ncol(Y) != 2)
      stop("Y must be a matrix with two columns with time in the first column and status indicator in the second.")
    
    
    if (use == "complete.obs") {
      ccY <- stats::complete.cases(Y)
      Y <- Y[ccY,]
      n <- length(ccY)
      if (type.X == "sample" && p == 1) {
        X <- X[ccY]
      } else if (type.X == "sample" && p > 1) {
        X <- X[ccY, ]
      }
      if (type.X == "distance") {
        X <- X[ccY,ccY]
      }
    }
    
    time <- Y[,1]
    status <- Y[,2]
    
   
      

    IX <- Rfast::Order(time)

    time <- time[IX]
    status <- status[IX]
    
    
    if (timetrafo != "none") {
      if (timetrafo == "log")
        time <- log(time)
      else {
        trafo <- match.fun(timetrafo)
        time <- trafo(time)
      }
    }

    if (type.X == "sample" && p == 1) {
      X <- X[IX]
    } else if (type.X == "sample" && p > 1) {
      X <- X[IX, ]
    }
    if (type.X == "distance") {
      X <- X[IX,IX]
    }
    


 

    if (use == "complete.obs") {
      ccX <- which(stats::complete.cases(X))
      n <- length(ccX)
    if (type.X == "sample" && p == 1) {
        X <- X[ccX]
      } else if (type.X == "sample" && p > 1) {
        X <- X[ccX, ]
      }
      if (type.X == "distance") {
        X <- X[ccX,ccX]
      }
      time <- time[ccX]
      status <- status[ccX]
    }
    
    
    
    if (!is.null(cutoff)) {
      time[time>=cutoff] <- cutoff
      status[time>=cutoff] <- 1
    }
    
    events <- which(status == 1)
    # calculate IPC weights
    
    
    ipcw <- rep(0, n)
    ipcw[events[1]] <- 1 / (n - events[1] + 1)
    
    help <- sapply(1:n, function(x) ((n - x) / (n - x + 1)) ^ status[x])
    
    for (i in (events[1] + 1):n) {
      ipcw[i]<- prod(help[1:(i - 1)]) * (status[i] / (n - i + 1))
    }
    
    ipcw <- ipcw * n
    
    ipcw <- ipcw[events]
    time <- time[events]
    status <- status[events]
    
    
    
    
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

    if (type.X == "distance") {
      distXall <- X
    } else {
      distXall <- distmat(X, metr.X, p)
    }
    distX <- distXall[events, events]
    distT <- Rfast::Dist(time)
    k <- sum(status)

    m <- sum(ipcw)
    M <- sum(ipcw^2)
    W2mat <- ipcw%*%t(ipcw)
    Wmatsum <- matrix(rep(ipcw, k), ncol=k) + matrix(rep(ipcw, k), ncol=k, byrow=TRUE)

    term1 <- sum(W2mat * distX * distT * (m ^ 2 -m * Wmatsum - M + 2 * W2mat))

    rowsumsX_w <- Rfast::colsums(ipcw*distX)
    rowsumsT_w <- Rfast::colsums(ipcw*distT)
    rowsumsX_w2 <- Rfast::colsums(ipcw^2*distX)
    rowsumsT_w2 <- Rfast::colsums(ipcw^2*distT)
    term2 <- sum(ipcw * (m + ipcw) * rowsumsX_w * rowsumsT_w)
    term3 <- sum(ipcw * rowsumsX_w2 * rowsumsT_w)
    term4 <- sum(ipcw * rowsumsX_w * rowsumsT_w2)
    term5 <- sum(W2mat * distX) * sum(W2mat * distT)

    dcov2 <- 1 / (n * (n-1) * (n-2) * (n-3)) * (term1 - 2 * (term2 - term3 - term4) + term5)

    dcov <- sign(dcov2) * sqrt(abs(dcov2))




    samp <- lapply(1:B, function(x) sample(1:n))

    reps <- sapply(1:B,
                 function(x) {
                   events_samp <- samp[[x]][events]
                   distX_samp <- distXall[events_samp,events_samp]

                   rowsumsX_w <- Rfast::colsums(ipcw*distX_samp)
                   rowsumsX_w2 <- Rfast::colsums(ipcw^2*distX_samp)


                   term1 <- sum(W2mat * distX_samp * distT * (m ^ 2 -m * Wmatsum - M + 2 * W2mat))


                   term2 <- sum(ipcw * (m + ipcw) * rowsumsX_w * rowsumsT_w)
                   term3 <- sum(ipcw * rowsumsX_w2 * rowsumsT_w)
                   term4 <- sum(ipcw * rowsumsX_w * rowsumsT_w2)
                   term5 <- sum(W2mat * distX_samp) * sum(W2mat * distT)

                   res2 <- 1 / (n * (n-1) * (n-2) * (n-3)) * (term1 - 2 * (term2 - term3 - term4) + term5)
                   res <- sign(res2) * sqrt(abs(res2))
                   return(res)
                 })

  pval <- (1 + length(which(reps>dcov))) / (1 + B)


  return(list("dcov"=dcov, "pvalue"=pval))
}




