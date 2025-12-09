centmat <- function(X,
                    metr.X = "euclidean",
                    p) {
  
  distX <- distmat(X, metr.X, p)
  return(doublecent(distX))
}



doublecent <- function(distX) {
  n <- nrow(distX)
  rm <- Rfast::rowmeans(distX)
  totm <- mean(rm)
  rmmat <- matrix(rep(rm,n), ncol=n)
  centmat <- distX - rmmat - t(rmmat) + totm
  return(centmat)
}




distmat <- function(X,
                    metr.X = "euclidean",
                    p) {
  len <- length(metr.X) 
  
  if (len == 1) {
    if (p == 1) {
      X <- as.numeric(X)
      if (metr.X == "euclidean" | metr.X == "alpha" | metr.X == "minkowski") {
        distX <- Rfast::vecdist(X)
      } else if (metr.X == "gaussian") {
        distX <- 1 - exp(-0.5 * Rfast::vecdist(X) ^ 2)
      } else if (metr.X == "gaussauto") {
        preX <- Rfast::vecdist(X)^2
        bw <- stats::median(stats::as.dist(preX), na.rm = TRUE)
        distX <- 1 - exp(- preX / bw)
      } else if (metr.X == "boundsq") {
        distX <- Rfast::vecdist(X) 
        distX <- (distX ^ 2) / (1 + distX ^ 2)
      } else if (metr.X == "discrete") {
        distX <- 1 * (Rfast::vecdist(X)>0)
      } else {
        distfunc <- match.fun(metr.X)
        distX <- distfunc(X = X)
        if (formals(distfunc)$kernel) {
          distX <- kerntodist(distX)
        }
      } 
    } else {
      if (metr.X == "euclidean" | metr.X == "alpha" | metr.X == "minkowski") {
        distX <- Rfast::Dist(X)
      } else if (metr.X == "gaussian") {
        distX <- 1 - exp(-0.5 * Rfast::Dist(X) ^ 2)
      } else if (metr.X == "gaussauto") {
        preX <- Rfast::Dist(X)^2
        bw <- stats::median(stats::as.dist(preX), na.rm= TRUE)
        distX <- 1 - exp(- preX / bw)
      } else if (metr.X == "boundsq") {
        distX <- Rfast::Dist(X) 
        distX <- (distX ^ 2) / (1 + distX ^ 2)
      } else if (metr.X == "discrete") {
        distX <- 1*(Rfast::Dist(X)>0)
      } else {
        distfunc <- match.fun(metr.X)
        distX <- distfunc(X = X)
        if (formals(distfunc)$kernel) {
          distX <- kerntodist(distX)
        }
      } 
    }
  } else {
    prm <- as.numeric(metr.X[2:len])
    metr.X <- metr.X[1]
    if (p == 1) {
      if (metr.X == "alpha") {
        distX <- Rfast::vecdist(X) ^ prm
      } else if (metr.X == "gaussian") {
        distX <- 1 - exp(-0.5 * Rfast::vecdist(X) ^ 2 / prm ^ 2)
      } else if (metr.X == "gaussauto") {
        preX <-Rfast::vecdist(X)^2
        bw <- stats::median(stats::as.dist(preX), na.rm= TRUE)
        distX <- 1 - exp(- preX / bw / prm ^ 2)
        } else if (metr.X == "boundsq") {
        distX <- Rfast::vecdist(X) 
        distX <- (distX ^ 2) / (prm ^ 2 + distX ^ 2)
      } else if (metr.X == "minkowski") {
        distX <- Rfast::vecdist(X) 
      } else {
        distfunc <- match.fun(metr.X)
        distX <- distfunc(X = X, prm = prm)
        if (formals(distfunc)$kernel) {
          distX <- kerntodist(distX)
        }
      } 
    } else {
      if (metr.X == "alpha") {
        distX <- Rfast::Dist(X) ^ prm
      } else if (metr.X == "gaussian") {
        distX <- 1 - exp(-0.5 * Rfast::Dist(X) ^ 2 / prm ^ 2)
      } else if (metr.X == "gaussauto") {
        preX <-Rfast::Dist(X)^2
        bw <- stats::median(stats::as.dist(preX), na.rm= TRUE)
        distX <- 1 - exp(- preX / bw / prm ^ 2)
      } else if (metr.X == "boundsq") {
        distX <- Rfast::Dist(X) 
        distX <- (distX ^ 2) / (prm ^ 2 + distX ^ 2)
      } else if (metr.X == "minkowski") {
        if (prm > 2 | prm <= 1)
          warning("Distance covariance with Minkowski distance outside the interval (1,2] does not define independence")
        distX <- as.matrix(stats::dist(X, method="minkowski", p=prm)) 
      } else {
        distfunc <- match.fun(metr.X)
        distX <- distfunc(X = X, prm = prm)
        if (formals(distfunc)$kernel) {
          distX <- kerntodist(distX)
        }
      } 
    }
  }
  return(distX)
}


kerntodist <- function(kernmat) {
  n <- ncol(kernmat)
  dmat <- matrix(rep(diag(kernmat), n), ncol = n)
  return((dmat + t(dmat) - 2 * kernmat) / 2)
}


#' Extract the dimensions of X.
#'
#' @param X a numeric vector or a numeric matrix.
#' @param type.X either "sample" or "distance". If type.X = "sample", X must be
#' a numeric vector or numeric matrix with the corresponding observations. If metr.X = "distance",
#' X must be a distance matrix.
#'
#' @return The centralized distance matrix corresponding to X.
extract_np <- function(X, type.X) {
  if (type.X == "sample") {
    if (is.vector(X) | is.factor(X))
    {
      n <- length(X)
      p <- 1L
    } else if (is.matrix(X) | is.data.frame(X)) {
      n <- nrow(X)
      p <- ncol(X)
    } else {
      stop("X must be a vector, matrix or dataframe for type 'sample'!")
    }
  } else if (type.X == "distance") {
    if (is.matrix(X)) {
      n <- nrow(X)
      p <- NA
    } else {
      stop("X must be a matrix for type 'distance'!")
    }
  } else {
    stop("type.X must be either 'sample' or 'distance'.")
  }
  return(list("Sample.Size" = n, "Dimension" = p))
}


decom.metr<- function(metr.X) {
  if (length(metr.X) == 1) {
    prmX <- 0
  } else {
    prmX <- as.numeric(metr.X[2])
    metr.X <- as.character(metr.X[1])
    if (is.na(prmX)) 
      stop("Parameters of standard metrics must be numeric")
  }
  return(c(metr.X,prmX))
}

normalize.sample <- function(X, n, p) {
  if (p > n) {
    stop("Affinely invariant distance covariance cannot be calculated for p>n")
  }
  if (p > 1) {
    X <- X %*% Rfast::spdinv(mroot(stats::var(X)))
  } else {
    X <- X / stats::sd(X)
  }
}


scale.sample <- function(X, n, p) {
  if (p > 1) {
    X <- standardise(X, center = FALSE)
  } else {
    X <- X / stats::sd(X)
  }
}

mroot <- function(A) {
  e <- eigen(A)
  V <- e$vectors
  V %*% diag(e$values) %*% t(V)
  
  
  B <- V %*% diag(sqrt(e$values)) %*% t(V)
  return(B)
}

