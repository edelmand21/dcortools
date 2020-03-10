#' Calculates the inverse-probability-of-censoring weighted (IPCW) distance correlation
#'
#' @param Y A column with two rows, where the first row contains the survival times and the second row the status indicators (a survival object will work).
#' @param X A vector or matrix containing the covariate information.
#' @param cutoff If provided, all survival times larger than cutoff are set to the cutoff and all corresponding status indicators are set to one. Under most circumstances, choosing a cutoff is highly recommended.
#' @return An inverse-probability of censoring weighted estimate for the distance correlation between X and the survival times.
#' @export


ipcw.dcor <- function(Y, X, cutoff = NULL) {
  
  dcor <- ipcw.dcov(Y, X, cutoff) / sqrt(ipcw.dcov(Y, Y[,1], cutoff)) / sqrt(distsd(X,bias.corr=TRUE))
  
  return(dcor)
}



#' Calculates the inverse-probability-of-censoring weighted (IPCW) distance covariance
#'
#' @param Y A column with two rows, where the first row contains the survival times and the second row the status indicators (a survival object will work).
#' @param X A vector or matrix containing the covariate information.
#' @param cutoff If provided, all survival times larger than cutoff are set to the cutoff and all corresponding status indicators are set to one. Under most circumstances, choosing a cutoff is highly recommended.
#' @return An inverse-probability of censoring weighted estimate for the distance covariance between X and the survival times.
#' @export


ipcw.dcov <- function(Y, X, cutoff = NULL) {

    # extract sample size and dimension

    ss.dimX <- dcortools:::extract_np(X, "sample")

    n <- ss.dimX$Sample.Size
    p <- ss.dimX$Dimension
    


    # sort data in order of follow-up times
    
    if (ncol(Y) != 2)
      stop("Y must be a matrix with two columns with time in the first column and status indicator in the second.")
    
    time <- Y[,1]
    status <- Y[,2]

    IX <- Rfast::Order(time)

    time <- time[IX]
    status <- status[IX]


    if (p==1)
        X <- X[IX]
    else
        X <- X[IX,]


    events <- which(status == 1)
    
    if (!is.null(cutoff)) {
      time[time>=cutoff] <- cutoff
      status[time>=cutoff] <- 1
    }
      
  
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


    if (p==1)
        X <- X[events]
    else
        X <- X[events,]


    #distance matrices

    distX <- Rfast:::Dist(X)
    distT <- Rfast:::Dist(time)
    k <- sum(status)

    # auxiliary terms

    m <- sum(ipcw)
    M <- sum(ipcw^2)
    W2mat <- ipcw %*% t(ipcw)
    ipcwm <- matrix(rep(ipcw, k), ncol=k)
    Wmatsum <- ipcwm + t(ipcwm)

    rowsumsX_w <- Rfast:::colsums(ipcw * distX)
    rowsumsT_w <- Rfast:::colsums(ipcw * distT)
    rowsumsX_w2 <- Rfast:::colsums(ipcw^2 * distX)
    rowsumsT_w2 <- Rfast:::colsums(ipcw^2 * distT)

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



#' Performs a permutation test based on the IPCW distance covariance
#'
#' @param Y A column with two rows, where the first row contains the survival times and the second row the status indicators (a survival object will work).
#' @param X A vector or matrix containing the covariate information.
#' @param cutoff If provided, all survival times larger than cutoff are set to the cutoff and all corresponding status indicators are set to one. Under most circumstances, choosing a cutoff is highly recommended.
#' @param B The number of permutations used for the permutation test
#' @return An list with two arguments, $dcov contains the IPCW distance covariance, $pvalue the corresponding p-value
#' @export


ipcw.dcov.test <- function(Y, X, cutoff = NULL, B=499)
{
    ss.dimX <- dcortools:::extract_np(X, "sample")

    n <- ss.dimX$Sample.Size
    p <- ss.dimX$Dimension
    
    if (ncol(Y) != 2)
      stop("Y must be a matrix with two columns with time in the first column and status indicator in the second.")
    
    time <- Y[,1]
    status <- Y[,2]

    IX <- Rfast::Order(time)

    time <- time[IX]
    status <- status[IX]

    if (p==1)
        X <- X[IX]
    else
        X <- X[IX,]


    events <- which(status == 1)
    
    if (!is.null(cutoff)) {
      time[time>=cutoff] <- cutoff
      status[time>=cutoff] <- 1
    }


    ipcw <- rep(0, n)
    ipcw[events[1]] <- 1 / (n - events[1]+1)

    help <- sapply(1:n, function(x) ((n - x) / (n - x +1)) ^ status[x])

    for (i in (events[1] + 1):n)
    {
        ipcw[i]<- prod(help[1:(i - 1)]) * (status[i] / (n - i + 1))
    }

    ipcw <- ipcw * n

    ipcw <- ipcw[events]
    time <- time[events]
    status <- status[events]



    distXall <- Rfast:::Dist(X)
    distX <- distXall[events, events]
    distT <- Rfast:::Dist(time)
    k <- sum(status)

    m <- sum(ipcw)
    M <- sum(ipcw^2)
    W2mat <- ipcw%*%t(ipcw)
    Wmatsum <- matrix(rep(ipcw, k), ncol=k) + matrix(rep(ipcw, k), ncol=k, byrow=TRUE)

    term1 <- sum(W2mat * distX * distT * (m ^ 2 -m * Wmatsum - M + 2 * W2mat))

    rowsumsX_w <- Rfast:::colsums(ipcw*distX)
    rowsumsT_w <- Rfast:::colsums(ipcw*distT)
    rowsumsX_w2 <- Rfast:::colsums(ipcw^2*distX)
    rowsumsT_w2 <- Rfast:::colsums(ipcw^2*distT)
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

                   rowsumsX_w <- Rfast:::colsums(ipcw*distX_samp)
                   rowsumsX_w2 <- Rfast:::colsums(ipcw^2*distX_samp)


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




