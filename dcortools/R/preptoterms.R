preptoterms.fast <- function(prepX, prepY, n, pairwise, discreteX, discreteY, perm) {
  terms <- list()
  
  if (!pairwise) {
    terms$aidot <- prepX$aidot
    terms$bidot <- prepY$aidot
    terms$adotdot <- prepX$adotdot
    terms$bdotdot <- prepY$adotdot
    terms$aijaij <- prepX$aijaij
    terms$bijbij <- prepY$aijaij
    terms$ncc <- n
  } else {
    cc <- intersect(prepX$cc, prepY$cc)
    ncc <- length(cc)
    terms$ncc <- ncc
  }
  
  if (!discreteX & !discreteY) {
    if (pairwise) {
  
      
      IX <- IY <- temp <- IY.s <-  1:ncc
      
      IX0here <- rep(NA,n)
      IXhere <- prepX$IX
      IXhere[prepY$mc] <- NA
      IX0here[IXhere[cc]] <- temp 
      cc0 <- !is.na(IX0here)
      IX0 <- IX0here[cc0]
      IX[IX0]<- temp
      
      IY0here <- rep(NA,n)
      IYhere <- prepY$IX
      IY0here[IYhere[cc]] <- temp 
      IYhere[prepX$mc] <- NA
      cc0 <- !is.na(IY0here)
      IY0 <- IY0here[cc0]
      IY[IY0]<- temp
      
      Xcc <- prepX$X[cc]
      Ycc <- prepY$X[cc]
      
      X.sort.X <- Xcc[IX0]
      Y.sort.Y <- Ycc[IY0]
      
      
      sX.X <- cumsum(X.sort.X)
      sY.Y <- cumsum(Y.sort.Y)
      
      alphaX <- IX - 1
      alphaY <- IY - 1
      betaX <- sX.X[IX] - X.sort.X[IX]
      betaY <- sY.Y[IY] - Y.sort.Y[IY]
      
      Xdot <- sX.X[ncc]
      Ydot <- sY.Y[ncc]
      
      terms$aidot <- Xdot + (2 * alphaX - ncc) * Xcc - 2 * betaX
      terms$bidot <- Ydot + (2 * alphaY - ncc) * Ycc - 2 * betaY
  
      terms$adotdot <- sum(terms$aidot)
      terms$bdotdot <- sum(terms$bidot)
      
      XXdot <- sum(Xcc^2)
      terms$aijaij <- 2 * ncc * XXdot - 2* Xdot^2
      
      YYdot <- sum(Ycc^2)
      terms$bijbij <- 2 * ncc * YYdot - 2* Ydot^2
      
  
      Y.sortX <- Ycc[IX0]
      XY.sortX <- Y.sortX * X.sort.X
      
      IY0.s <- IX[IY0]
      IY.s[IY0.s] <- temp
      
      
      
      
      } else {
      IY.s <- temp <- 1:n
   
      X.sort.X <- prepX$X.sort.X
      IX0 <- prepX$IX0
      
      Y.sortX <- prepY$X[IX0]
      XY.sortX <- Y.sortX *  X.sort.X
      
      IY0.s <- prepX$IX[prepY$IX0]
      IY.s[IY0.s] <- temp
 
      
      sX.X <- prepX$sX.X
      sY.Y <- prepY$sX.X
      
      Ycc <- prepY$X
      IY <- prepY$IX
   
      }
    
   sY.X <- cumsum(X.sort.X[IY0.s])
   sX.Y <- cumsum(Y.sortX)
   sY.XY <- cumsum(XY.sortX[IY0.s]) 
   sX.XY <- cumsum(XY.sortX)
      
   terms$aijbij <- SUMAIJBIJ(IY.s, X.sort.X, Y.sortX, XY.sortX, sX.X, sX.Y, sX.XY, sY.X[IY.s], sY.Y[IY.s], sY.XY[IY.s])
   
   if (perm) {
     terms$n <- terms$ncc
     terms$Y <- Ycc  
     terms$IX0 <- IX0
     terms$X.sort.X <- X.sort.X
     terms$IY <- IY
     terms$sX.X <- sX.X
     terms$sY.Y <- sY.Y
     terms$ndisc <- 0
   }

  
 } else if (discreteX + discreteY == 1) {
   if (discreteX) {
     prX <- prepY
     prY <- prepX
     prepX <- prX
     prepY <- prY
   }
   
   if (pairwise) {

     IX <- IY <- temp <- IY.s <-  1:ncc
     
     IX0here <- rep(NA,n)
     IXhere <- prepX$IX
     IXhere[prepY$mc] <- NA
     IX0here[IXhere[cc]] <- temp 
     cc0 <- !is.na(IX0here)
     IX0 <- IX0here[cc0]
     IX[IX0]<- temp
     
   
     Y <- prepY$X[cc]
     Y.sortX <- Y[IX0]
     X <- prepX$X[cc]
     X.sort.X <- X[IX0]
     
     sX.X <- cumsum(X.sort.X)
     alphaX <- 0:(ncc-1)
     betaX <- sX.X - X.sort.X
     Xdot <-sX.X[ncc]
     Y.sortX <- Y[IX0]
     
     terms$aidot <- Xdot + (2 * alphaX - ncc) * X.sort.X - 2 * betaX
     terms$adotdot <- sum(terms$aidot)
     XXdot <- sum(X^2)
     terms$aijaij <- 2 * ncc * XXdot - 2* Xdot^2

     nY <- table(Y)
     terms$bidot <- as.numeric(ncc - nY[Y.sortX])
     terms$bdotdot <- sum(terms$bidot)
     terms$bijbij <- ncc * (ncc - 1) - sum(nY * (nY - 1))
     
   } else {
     Y <- prepY$X
     IX0 <- prepX$IX0
     Y.sortX <- Y[IX0]
     X.sort.X <- prepX$X.sort.X
     nY <- table(Y)
     terms$aidot <- terms$aidot[IX0]
     terms$bidot <- terms$bidot[IX0]
   }
     
   levY <-levels(Y)
   aijbij <- terms$adotdot
     
   for (lvl in levY) {
     set0 <- which(Y.sortX == lvl)
     n0 <- length(set0)
     alphaX0 <- 0:(n0-1)
     X0 <- X.sort.X[set0]
     sX.X0 <- cumsum(X0)
     betaX0 <- sX.X0 - X0
     X0dot <- sX.X0[n0]
     aidot0 <- X0dot + (2 * alphaX0 - n0) * X0 - 2 * betaX0
     adotdot0 <- sum(aidot0)
     aijbij <- aijbij - adotdot0
   }
     
   terms$aijbij <- aijbij
   
   if (perm) {
     terms$X.sort.X <- X.sort.X
     terms$Y.sortX <- as.factor(Y.sortX)
     terms$levY <- levY
     terms$ndisc <- 1
   }
   
 } else if (discreteX & discreteY) {
      if (pairwise) {
        X <- prepX$X[cc]
        Y <- prepY$X[cc]
        nXY <- table(X,Y)
        nX <- rowSums(nXY)
        nY <- colSums(nXY)
        
        aidot <-ncc - nX
        terms$adotdot <- sum(aidot * nX)
        bidot <- ncc - nY
        terms$bdotdot <- sum(bidot * nY)
        
        terms$aidot <- aidot[X]
        terms$bidot <- bidot[Y] 
        terms$aijaij <- ncc * (ncc - 1) - sum(nX * (nX - 1)) 
        terms$bijbij <- ncc * (ncc - 1) - sum(nY * (nY - 1))
        
        terms$aijbij <- ncc * ncc - sum(nX * nX) - sum(nY * nY) + sum(nXY * nXY)
      
        if (perm) {
          terms$aidotshort <- aidot
          terms$bidotshort <- bidot
          terms$nX <- nX
          terms$nY <- nY
          terms$X <- X
          terms$Y <- Y
          terms$n <- terms$ncc
          terms$ndisc <- 2
        }   
      
    } else {
        nXY <- table(prepX$X,prepY$X)
        nX <- prepX$nX
        nY <- prepY$nX
        
        terms$aijbij <- n * n - sum(nX * nX) - sum(nY * nY) + sum(nXY * nXY)
        
        if (perm) {
          terms$aidotshort <- prepX$aidotshort
          terms$bidotshort <- prepY$aidotshort
          terms$nX <- prepX$nX
          terms$nY <- prepY$nX
          terms$X <- prepX$X
          terms$Y <- prepY$X
          terms$n <- terms$ncc
          terms$ndisc <- 2
        } 
        
    }
    
 }
 
  return(terms)
  
}




preptoterms.standard <- function(prepX, prepY, n, pairwise, perm) {
  distX <- prepX$distX
  distY <- prepY$distX

  if (!pairwise) {
    aijbij <- matrix_prod_sum(distX, distY)
    
    if (perm) {
      return(list("aijbij" = aijbij,"aijaij" = prepX$aijaij, "aidot" = prepX$aidot,"bijbij" = prepY$aijaij, "bidot" = prepY$aidot, "adotdot" = prepX$adotdot, "bdotdot" = prepY$adotdot, "ncc" = n, "distX" = distX, "distY" = distY))
    } else {
      return(list("aijbij" = aijbij,"aijaij" = prepX$aijaij, "aidot" = prepX$aidot,"bijbij" = prepY$aijaij, "bidot" = prepY$aidot, "adotdot" = prepX$adotdot, "bdotdot" = prepY$adotdot, "ncc" = n))
    }
  } else {
    cc <- intersect(prepX$cc, prepY$cc)
    ncc <- length(cc)
    aidot <- colsums_subset(distX,cc)
    adotdot <- sum(aidot)
    bidot <- colsums_subset(distY,cc)
    bdotdot <- sum(bidot)
    aijaij <- matrix_prod_sum_subset(distX, distX, cc)
    bijbij <- matrix_prod_sum_subset(distY, distY, cc)
    aijbij <- matrix_prod_sum_subset(distX, distY, cc)
    
    
    if (perm) {
      return(list("aijbij" = aijbij,"aijaij" = aijaij, "aidot" = aidot,"bijbij" = bijbij, "bidot" = bidot, "adotdot" = adotdot, "bdotdot" = bdotdot, "ncc" = ncc, "distX" = distX[cc,cc], "distY" = distY[cc,cc]))
    } else {
    return(list("aijbij" = aijbij, "aijaij" = aijaij, "bijbij" = bijbij, "aidot" = aidot, "bidot" = bidot, "adotdot" = adotdot, "bdotdot" = bdotdot, "ncc" = ncc))
    }
  }
}


preptoterms.memsave <- function(prepX, prepY, metr.X, metr.Y, n, pairwise, perm) {

  
   if (pairwise) {
    cc <- intersect(prepX$cc, prepY$cc)
    X <- prepX$X[cc]
    Y <- prepY$X[cc]
    ncc <- length(cc)
    dvartermsX <- dvarterms.memsave(X, metr.X, prepX$p)
    dvartermsY <- dvarterms.memsave(Y, metr.Y, prepY$p)
    aijbij <- aijbij.memsave(X, Y, metr.X, metr.Y, prepX$p, prepY$p)
 
    if (perm) { 
      dec.metr.X <- decom.metr(metr.X)
      metr.X <- dec.metr.X[1]
      prmX <- dec.metr.X[2]
      
      dec.metr.Y <- decom.metr(metr.Y)
      metr.Y <- dec.metr.Y[1]
      prmY <- dec.metr.Y[2]
      
      return(list("aijbij" = aijbij, "aijaij" = dvartermsX$aijaij, "bijbij" = dvartermsY$aijaij, "aidot" = dvartermsX$aidot, "bidot" = dvartermsY$aidot, "adotdot" = dvartermsX$adotdot, "bdotdot" = dvartermsY$adotdot, "ncc" = ncc, "X" = X, "Y" = Y, "p" = prepX$p, "q" = prepY$p, "metr.X" = metr.X, "metr.Y" = metr.X, "prmX" = prmX, "prmY" = prmY))
   }  else 
      return(list("aijbij" = aijbij, "aijaij" = dvartermsX$aijaij, "bijbij" = dvartermsY$aijaij, "aidot" = dvartermsX$aidot, "bidot" = dvartermsY$aidot, "adotdot" = dvartermsX$adotdot, "bdotdot" = dvartermsY$adotdot, "ncc" = ncc))
  } else {
    ncc <- n
    X <- prepX$X
    Y <- prepY$X

  aijbij <- aijbij.memsave(X, Y, metr.X, metr.Y, prepX$p, prepY$p)
  if (perm) {
    dec.metr.X <- decom.metr(metr.X)
    metr.X <- dec.metr.X[1]
    prmX <- dec.metr.X[2]
    
    dec.metr.Y <- decom.metr(metr.Y)
    metr.Y <- dec.metr.Y[1]
    prmY <- dec.metr.Y[2]
    return(list("aijbij" = aijbij, "aijaij" = prepX$aijaij, "aidot" = prepX$aidot, "bidot" = prepY$aidot, "bijbij" = prepY$aijaij, "adotdot" = prepX$adotdot, "bdotdot" = prepY$adotdot, "ncc" = ncc, "X" = X, "Y" = Y, "p" = prepX$p, "q" = prepY$p, "metr.X" = metr.X, "metr.Y" = metr.X, "prmX" = prmX, "prmY" = prmY))
  } else
    return(list("aijbij" = aijbij, "aijaij" = prepX$aijaij, "aidot" = prepX$aidot, "bidot" = prepY$aidot, "bijbij" = prepY$aijaij, "adotdot" = prepX$adotdot, "bdotdot" = prepY$adotdot, "ncc" = n, "X" = X, "Y" = Y, "p" = prepX$p, "q" = prepY$p, "metr.X" = metr.X, "metr.Y" = metr.X))
  }
}