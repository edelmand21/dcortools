#' Calculates distance correlation from Pearson correlation under assumption of a bivariate normal distribution
#' @param rho Pearson correlation.
#' @return Distance correlation assuming a bivariate normal distribution
#' @export
dcorgaussianbiv <- function(rho)
{
  nenner <-   rho * asin(rho) + sqrt(1-rho^2) - rho * asin(rho/2) - sqrt(4-rho^2) +1
  zaehler <- 1+ pi/3 - sqrt(3)
  return(sqrt(nenner/zaehler))
}


#' @export
print.dcmatrix <- function(x, ...) {
  dcmat <- x
  if (dcmat$test == "none")
    testout <- "None"
    else if (dcmat$test == "permutation")
    testout <- paste("Permutation test with ",dcmat$b," permutations")
    else if (dcmat$test == "gamma")
    testout <- "Simple gamma approximation"               
    else if (dcmat$test == "conservative")
    testout <- "Conservative two-moment approximation"  
    else if (dcmat$test == "bb3")
    testout <- "Three moment approximation by Berschneider and Boettcher"
    
    
    if (dcmat$test.pw == "none")
      testout.pw <- "None"
    else if (dcmat$test.pw == "permutation")
      testout.pw  <- paste("Permutation test with ",dcmat$b," permutations")
    else if (dcmat$test.pw == "gamma")
      testout.pw  <- "Simple gamma approximation"               
    else if (dcmat$test.pw == "conservative")
      testout.pw  <- "Conservative two-moment approximation"  
    else if (dcmat$test.pw == "bb3")
      testout.pw  <- "Three moment approximation by Berschneider and Boettcher"
    
    pw <- (dcmat$test.pw != "none") | dcmat$calc.dcor.pw | dcmat$calc.dcov.pw
    
    if (dcmat$calc.cor == "pearson")
      corout <- "Pearson"
    else if (dcmat$calc.cor == "spearman")
    corout <- "Spearman"
    else if (dcmat$calc.cor == "kendall")
    corout <- "Kendall"               
    else 
    corout <- "None"     
    
  dX <- dY <- dcmat$dX
  dX.pw <- dY.pw <- dcmat$dX.pw
    
  if (dcmat$withY) {  
    dY.pw <- dcmat$dY.pw
    dY <- dcmat$dY
  }
  cat(paste("dcmatrix object of size", dX, "x", dY, "\n"))
  cat("\n")
  cat("Call:", "\n")
  cat(paste(deparse(dcmat$call), "\n"))
  cat("\n")
  cat(paste("Distance correlation matrix: ", ifelse(dcmat$calc.dcor,"Available \n", "Not available \n")) )
  cat(paste("Distance covariance matrix: ", ifelse(dcmat$calc.dcov,"Available \n", "Not available \n")) )
  cat(paste("Bias correction: ", ifelse(dcmat$bias.corr,"Yes \n", "No \n")) )
  cat(paste("Affinely invariant: ", ifelse(dcmat$affine,"Yes \n", "No \n")) )
  cat(paste("Test: ", testout, "\n") )
  cat(paste("Standard correlation matrix: ", corout, "\n"))
  if (pw) {
    cat("\n")
    cat("Pairwise results:\n")
    cat("\n")
    cat(paste("Size:", dX.pw, "x", dY.pw, "\n"))
    cat("\n")
    cat(paste("Distance correlation matrix: ", ifelse(dcmat$calc.dcor.pw,"Available \n", "Not available \n")) )
    cat(paste("Distance covariance matrix: ", ifelse(dcmat$calc.dcov.pw,"Available \n", "Not available \n")) )
    cat(paste("Test: ", testout.pw, "\n") )
  }
}

#' Plots a heatmap from a dcmatrix object using the function "pheatmap" from the package "pheatmap".
#'
#' @param x a dcmatrix object.
#' @param type specifies what should be displayed in the heatmap. One of "dcor", "dcov", "logp" (-log10 of corresponding p-values), "cor", "abscor" (absolute correlation), "logp.cor", "dcor.pw", "dcov.pw" or "logp.pw".
#' @param trunc.up truncates the values to be plotted; if set to numeric, all values larger than trunc.up are set to trunc.up.
#' @param trunc.low truncates the values to be plotted; if set to numeric, all values smaller than trunc.low are set to trunc.low.
#' @param cluster_rows,cluster_cols,display_numbers passed to pheatmap().
#' @param ... passed to pheatmap().
#' @return a heatmap plotting the entries of the slot specified in type of the object specified in dcmat.
#' @export
plot.dcmatrix <- function(x, type = "dcor", trunc.up = NULL, trunc.low = NULL, cluster_rows = FALSE, cluster_cols = FALSE, display_numbers=TRUE,  ...) {
  dcmat <- x
  if (type =="dcor") {
    mat <- dcmat$dcor 
  } else if (type == "dcov") {
    mat <- dcmat$dcov 
  } else if (type == "logp") {
    mat <- -log10(dcmat$pvalue)
  } else if (type == "cor") {
    mat <- dcmat$cor
  } else if (type == "abscor") {
    mat <- abs(dcmat$cor)
  } else if (type == "logp.cor") {
    mat <- -log10(dcmat$pval.cor)
  } else if (type == "dcor.pw") {
    mat <- dcmat$dcor.pw 
  } else if (type == "dcov.pw") {
    mat <- dcmat$dcov.pw 
  } else if (type == "logp.pw") {
    mat <- -log10(dcmat$pvalue.pw)
  } else {
    stop(" \"type\" must be one of \"dcor\",\"logp\",\"cor\",\"abscor\", \"logp.cor\",\"dcor.pw\",\"dcov.pw\" or \"logp.pw\" ")
  }
  
  if (!is.null(trunc.up) | !is.null(trunc.low)) {
    mat[mat > trunc.up] <- trunc.up
    mat[mat < trunc.low] <- trunc.low
  }  
  
  if (anyNA(mat)|any(mat==Inf)) {
    warning("NAs and Infs were set to 0 for plotting.")
    mat[is.na(mat)] <- 0
    mat[mat == Inf] <- 0
  }
  

  pheatmap::pheatmap(mat, cluster_rows = cluster_rows, cluster_cols = cluster_cols, display_numbers= display_numbers, ...)

}



#' Plots Pearson/Spearman/Kendall correlation against distance correlation (often resembling a horseshoe(hs)).
#'
#' @param dcmat A dcmatrix object.
#' @param maxcomp Maximum number of associations, for which distance correlation is plotted against correlation. If the number of associations in the dcmat object is larger, only the maxcomp associations with the largest difference between distance correlation and absolute (Pearson/Spearman/Kendall) correlation are plotted.
#' @param col color of the plot.
#' @param alpha alpha parameter of the plot.
#' @param cortrafo Either "none" or "gaussiandcor". If "gaussiandcor", the distance correlation under assumption of normality is calculated and plotted against the actual distance correlation. 
#' 
#' Note that this is only sensible for Pearson correlation!
#' @return Plot of (possibly transformed) Pearson/Spearman/Kendall correlation against distance correlation.
#' @export
hsplot <- function(dcmat, maxcomp = 1e5, col = "blue", alpha = 1, cortrafo = "none") {
  pw <- FALSE
  
  if (is.null(dcmat$corr) | is.null(dcmat$dcor))
    stop("Correlation and (pairwise) distance correlation matrix needed for horseshoe plot!")

  if (anyDuplicated(dcmat$group.X)) {
   pw <- TRUE
   if (is.null(dcmat$dcor.pw))
     stop("Pairwise distance correlation matrix needed for horseshoe plot!")
  }#
  
  if (dcmat$withY) {
    if (anyDuplicated(dcmat$group.Y)) {
      pw <- TRUE
      if (is.null(dcmat$dcor.pw))
        stop("Pairwise distance correlation matrix needed for horseshoe plot!")
    }
  }
    
  
  

  if (!dcmat$withY) {
    corr <- c(stats::as.dist(dcmat$corr))
    if (pw)
      dcor <- c(stats::as.dist(dcmat$dcor.pw))
    else
      dcor <- c(stats::as.dist(dcmat$dcor))
  } else {
    
    corr <- c(dcmat$corr)
    if (pw)
      dcor <- c(dcmat$dcor.pw)
    else 
      dcor <- c(dcmat$dcor)
  }
  
  
  
  xl <- "Correlation"

  if (cortrafo == "gaussiandcor") {
    corr <- dcorgaussianbiv(corr)
    xl <- "Distance Correlation under asssumption of bivariate normal"
  }

  if (length(corr) > maxcomp) {
    sel <- order(dcor-abs(corr), decreasing=TRUE)[1:maxcomp]
    warning("Not all comparisons plotted. Increase argument maxcomp to view all comparisons")
    dcor <- dcor[sel]
    corr <- corr[sel]
  }
  print(ggplot2::ggplot() + ggplot2::geom_point(ggplot2::aes(dcor,corr), alpha = alpha, color =col) + ggplot2::ylab(xl) + ggplot2::xlab("Distance Correlation"))
}

#' #' @export
#' plotrelations <- function(dcmat, candidates = NULL, diff, smooth = "loess", cortrafo = "gaussiandcor", switchXY = FALSE, names = "nonlin", width = 480, height = 480, criterion = "dcor") {
#'   if (is.null(dcmat$X)) {
#'     stop("Data needed to plot interesting nonlinear associations! Run dcmatrix again or add data manually to dcmat object.")
#'   } 
#'   
#'   X <- as.matrix(dcmat$X)
#'   
#'   if (dcmat$withY) 
#'     Y <- as.matrix(dcmat$Y)
#'   
#'   if (criterion == "logp") {
#'     if (is.null(dcmat$pvalue) | is.null(dcmat$pval.cor)) {
#'       stop("P-values missing in dcmatrix object.")
#'     } 
#'     corr <- -log10(dcmat$pval.cor)
#'     dcor <- -log10(dcmat$pvalue)
#'   } else if (criterion == "dcor") {
#'     if (is.null(dcmat$dcor) | is.null(dcmat$corr)) {
#'       stop("Distance correlation or correlation missing in dcmatrix object.")
#'     }
#'     corr <- abs(dcmat$corr)
#'     dcor <- dcmat$dcor
#'   }
#' 
#'   
#'   if (!dcmat$withY) {  
#'     lowtr <- lower.tri(corr)
#'     corr[lowtr] <- NA
#'     dcor[lowtr] <- NA
#'     Y <- X
#'   }
#'   
#'   if (cortrafo == "gaussiandcor") {
#'     corr <- dcorgaussianbiv(corr)
#'     xl <- "Distance Correlation under asssumption of bivariate normal"
#'   }  
#'   
#'   if (is.null(candidates)) {
#'     candidates <- which((dcor - corr) > diff, arr.ind=TRUE)
#'   }
#'   ncand <- nrow(candidates)
#'   
#'   
#'   if (ncand<0.5) {
#'     stop("No associations with specified difference between two tests/association measures found.")
#'   }
#'   
#'   names.Y <- "Y"
#'   
#'   if (is.null(colnames(X)))
#'     colnames(X) <- paste("X",1:ncol(X))
#'   
#'   
#'   
#'   
#'   
#'   if (!dcmat$withY) {  
#'     Y <- X
#'     names.Y <- "X"
#'     if (switchXY) {
#'       candidates <- cbind(candidates[,2],candidates[,1])
#'     }
#'   } else {
#'     if (dcmat$withY) {
#'       if (is.null(colnames(Y)))
#'         colnames(Y) <- paste("Y",1:ncol(Y), sep = "")
#'     }
#'     if (switchXY) {
#'       Xnew <- Y
#'       Ynew <- X
#'       X <- Xnew
#'       Y <- Ynew
#'       candidates <- cbind(candidates[,2],candidates[,1])
#'     }
#'   }
#'   
#'   
#'   for (i in 1:ncand) {
#'     name.X <- colnames(X)[candidates[i,1]]
#'     name.Y <- colnames(Y)[candidates[i,2]]
#'     
#'     title <- paste("Association of ", name.X," with ", name.Y,sep="")
#'     filename <- paste(names,"X-",name.X,"_Y-",name.Y,".png",sep="")
#'     
#'     png(filename, width = width, height = height)
#'     plotobj <- ggplot2::ggplot()+ggtitle(title)+ ggplot2::xlab(name.X) + ggplot2::ylab(name.Y)
#'     plotobj <- plotobj + ggplot2::geom_point(ggplot2::aes(X[,candidates[i,1]], Y[,candidates[i,2]]))
#'     if (smooth != "none")
#'       plotobj <- plotobj + geom_smooth(ggplot2::aes(X[,candidates[i,1]], Y[,candidates[i,2]]), method=smooth, se=F)  
#'     print(plotobj)
#'     dev.off()
#'   }
#' }
#' 
#' #' @export
#' extract.nonlins <- function(dcmat, diff, criterion = "dcor") {
#'   if (is.null(dcmat$X)) {
#'     stop("Data needed to plot interesting nonlinear associations! Run dcmatrix again or add data manually to dcmat object.")
#'   } 
#'   
#'   if (criterion == "logp") {
#'     if (is.null(dcmat$pvalue) | is.null(dcmat$pval.cor)) {
#'       stop("P-values missing in dcmatrix object.")
#'     } 
#'     corr <- -log10(dcmat$pval.cor)
#'     dcor <- -log10(dcmat$pvalue)
#'   } else if (criterion == "dcor") {
#'     if (is.null(dcmat$dcor) | is.null(dcmat$corr)) {
#'       stop("Distance correlation or correlation missing in dcmatrix object.")
#'     }
#'     corr <- abs(dcmat$corr)
#'     dcor <- dcmat$dcor
#'   }
#'   
#'   
#'   if (!dcmat$withY) {  
#'     lowtr <- lower.tri(corr)
#'     corr[lowtr] <- NA
#'     dcor[lowtr] <- NA
#'     Y <- X
#'   }
#'   
#'   selint <- which((dcor - corr) > diff, arr.ind=TRUE)
#'   
#'   return(selint)
#'  
#' }
#' 
#' #' @export
#' plotsubgr <- function(dcmat, candidates, data, cov.subgroup = names(data), splitnum = 3, method = "color", smooth = "lm", switchXY = FALSE, names = "plot", width = 480, height = 480) {
#' 
#'   if (is.null(dcmat$X)) {
#'     stop("Data needed to plot subgroup associations! Run dcmatrix again or add data manually to dcmat object.")
#'   } 
#'   
#'   if ("X" %in% cov.subgroup | "Y" %in% cov.subgroup)  {
#'     stop("Variable name X and Y not allowed in cov.subgroup. Please rename variables.")
#'   } 
#'   
#'   
#'   X <- as.matrix(dcmat$X)
#'   if (dcmat$withY) 
#'     Y <- as.matrix(dcmat$Y)
#'   
#'   if (length(candidates)==2)
#'     candidates <- matrix(candidates,ncol=2)
#'   
#'   names.Y <- "Y"
#'   
#'   if (is.null(colnames(X)))
#'     colnames(X) <- paste("X",1:ncol(X))
#'   
#' 
#' 
#'   
#'   
#'   if (!dcmat$withY) {  
#'     Y <- X
#'     names.Y <- "X"
#'     if (switchXY) {
#'       candidates <- cbind(candidates[,2],candidates[,1])
#'     }
#'   } else {
#'     if (dcmat$withY) {
#'       if (is.null(colnames(Y)))
#'         colnames(Y) <- paste("Y",1:ncol(Y), sep = "")
#'     }
#'    if (switchXY) {
#'       Xnew <- Y
#'       Ynew <- X
#'       X <- Xnew
#'       Y <- Ynew
#'       candidates <- cbind(candidates[,2],candidates[,1])
#'     }
#'   }
#'   
#'   nsubgr <- length(cov.subgroup)
#'   
#'   for (j in 1:nsubgr) {
#'     if (is.numeric(data[,cov.subgroup[j]])) {
#'       data[,cov.subgroup[j]] <- gtools::quantcut(data[,cov.subgroup[j]], q= splitnum)
#'     }
#'   }
#' 
#'   
#' 
#'   ncand <- nrow(candidates)
#'   for (i in 1:ncand) {
#'     for (j in 1:nsubgr) {
#'       name.X <- colnames(X)[candidates[i,1]]
#'       name.Y <- colnames(Y)[candidates[i,2]]
#'       
#'       title <- paste("Association of ", name.X," with ", name.Y, " in ",cov.subgroup[j]," subgroups",sep="")
#'       filename <- paste(names,"X-",name.X,"_Y-",name.Y,"_subgroup",cov.subgroup[j],".png",sep="")
#'       png(filename, width = width, height = height)
#'       plotobj <- ggplot2::ggplot(data)+ggtitle(title)+ ggplot2::xlab(name.X) + ggplot2::ylab(name.Y)
#'       
#'       plotobj <- plotobj + ggplot2::geom_point(ggplot2::aes(X[,candidates[i,1]], Y[,candidates[i,2]], color = data[,cov.subgroup[j]]))
#'       
#'       if (method == "facet")
#'         plotobj <- plotobj + ggplot2::geom_point(ggplot2::aes(X[,candidates[i,1]], Y[,candidates[i,2]])) + facet_grid(rows = vars(data[,cov.subgroup[j]]))
#'       
#'       if (smooth != "none")
#'         plotobj <- plotobj + geom_smooth(ggplot2::aes(X[,candidates[i,1]], Y[,candidates[i,2]], group=data[,cov.subgroup[j]], color = data[,cov.subgroup[j]]), method=smooth, se=F)  
#'       
#'       print(plotobj)
#'       
#'       dev.off()
#'     }
#'   }
#'   
#' }
#' 
#' is.dcmatrix <- function(x) inherits(x, "dcmatrix")