#' Plots a linked chart for exploratory data analysis from a dcmatrix object using the functionality provided by the "rlc" package. 
#'
#' @param dcmatrix A dcmatrix object
#' @param X should only be specified if the dcmatrix object does not contain the original data (see argument "return.data" of dcmatrix); should contain the same data as used in the calculation of the dcmatrix object.
#' @param Y see X.
#' @param heatmap specifies the dependence measure to use for the heatmap. Options are "dcor", "dcov", "logp" (for -log10 of the p-value), "cor" (for Pearson correlation) and "abscor" (for the absolute value of Pearson correlation).
#' @param heatmap0 Only applicable if the dcmatrix object contains dependencies between objects that are not one-dimensional. In this case, two heatmaps are plotted, the first using the dependence meausure specified in heatmap0 and the dependencies between the groups in dcmatrix, the second using the dependence masure specified in heatmap and between the corresponding univariate variables.
#' @param size Passed to the options of the scatter plot.
#' @param opacity Passed to the options of the scatter plot.
#' @param discX Numeric vector specifying which columns in X should be interepreted as discrete variables. Factors and characters are always interpreted as discrete.
#' @param discY see discX
#' @param jitt.disc Jitter added to discrete variables.
#' @param jitt.cont Jitter added to continuous variables.
#' @param smooths specifies smooths that fitted to the scatter plots. Options include "none", "xtoy" (x as predictor, y as dependent variable), "ytox" (y as predictor, x as dependent variable) and "both".
#' @param smooth.type Type of the plotted smooths. Options are "loess", "spline" or "regline".
#' @param palette.heatmap Passed to the options of the heatmap.
#' @param palette.heatmap0 Only applicable if the dcmatrix object contains dependencies between objects that are not one-dimensional. In this case, passed to the options of the corresponding heatmap.
#' @param cgtabmax If set to a numeric value, a contingency table is shown for any pair of variables for which each of the variable has not more than cgtabmax unique values.
#' @param subdat optional; Data.frame with the same number of observations as the original data.
#' @param subs Names of some factor(!) variables in subdat. This allows to color the observations corresponding to these variables.
#' @param ... passed to the functions used for the smooths.
#' @return A linked chart object opening in the browser. 
#' 
#' If the dcmatrix object contains only associations between univariate observations, this will consist of one heatmap with corresponding scatterplots.
#' 
#' If the dcmatrix contains associations between groups of variables, there are three windows. The first window displays a fixed heatmap contains the associations between the groups. The second window contains a changing heatmap containing the associations between single observations in selected groups. The third windows display a scatter plot between selected single observations.
#' @examples A linked chart object opening in the browser.
#' X <-  matrix(rnorm(1000), ncol = 10)
#' dcm <- dcmatrix(X, return.data = TRUE)
#' linked(dcm)
#' 
#' X <-  matrix(rnorm(1000), ncol = 10)
#' dcm <- dcmatrix(X, return.data = TRUE, calc.dcor.pw = TRUE, group.X = c(rep("group1",3),rep("group2",3),rep("group3",4)))
#' linked(dcm)
#' 
#' Y <- matrix(rnorm(600), ncol = 6)
#' Y[,5] <- rbinom(100, 3, 0.5)
#' Y[,6] <- rbinom(100, 2, 0.3)
#' dcm <- dcmatrix(X, Y, return.data = TRUE, calc.dcor.pw = TRUE, group.X = c(rep("group1",3),rep("group2",3),rep("group3",4)), group.Y = c(rep("group1",4),rep("group2",1),rep("group3",1)), metr.Y = c("group1" = "euclidean", "group2" = "discrete", "group3" = "discrete"))
#' linked(dcm, discY = c(5,6))
#' @export
linked <- function(dcmatrix, X = NULL, Y =NULL, heatmap = "dcor", heatmap0 = "dcor", size = 2, opacity = 1, discX =NULL, discY =NULL, jitt.disc = 0.1, jitt.cont = 0, smooths ="none", smooth.type = "loess", cgtabmax = NULL, subdat = NULL, subs =NULL, palette.heatmap = RColorBrewer::brewer.pal(9, "Reds"), palette.heatmap0 =RColorBrewer::brewer.pal(9, "YlGnBu"), coldom = NULL, coldom0 = NULL, ...){

  
  beeswarm = FALSE
  docg <- !is.null(cgtabmax)
  dosub <- !is.null(subdat) & !is.null(subs)

  #jittX <- jitt
  #jittY <- jitt
  
  if (dosub) {
    for (i in 1:length(subs))
      subdat[,subs[i]] <- as.factor(subdat[,subs[i]])
  }
 
  
  if (docg) {
    tabnew <- function(X,Y) {
      if (max(length(unique(X)),length(unique(Y))) <= cgtabmax) {
        A <- as.matrix(table(X,Y))
        class(A) <- "matrix"
        return(A)
      }
      else
        as.character("Too many levels.")
    }
  }
  
 

  
  if (is.null(X) & is.null(dcmatrix$X) )
    stop("Linked Chart requires original 
         data. Run dcmatrix again or enter data via argument X.")
  
  if (is.null(X))
    X <- dcmatrix$X
  
  #if (is.vector(X))
  #  X <- as.matrix(X)
  
  X <- as.data.frame(X)

  
  n <- nrow(as.matrix(X))
  colgr <- rep(1,n)
 
  is.facX <- sapply(1:ncol(X), function(t) is.factor(X[,t])|is.character(X[,t])) 
  is.facX[discX] <- TRUE
  for (t in 1:ncol(X))
    if (is.facX[t]) 
      X[,t] <- as.character(X[,t])
      
 
  if (beeswarm)
    lc_points <- lc_beeswarm
  else
    lc_points <- lc_scatter
  
  triplot <- FALSE
  group.X <- group.Y <-  dcmatrix$group.X
 
  
  if (anyDuplicated(group.X)) {
    triplot <- TRUE
  }
     gX <- gY <- dcmatrix$groupslistX


  
  
  if (!is.null(Y) | !is.null(dcmatrix$Y)) {
    group.Y <- dcmatrix$group.Y
    if (is.null(Y))
      Y <- dcmatrix$Y
    if (anyDuplicated(dcmatrix$group.Y)) {
      triplot <- TRUE
    }
      gY <- dcmatrix$groupslistY
      Y <- as.data.frame(Y)
    #if (is.vector(Y))
    #  Y <- as.matrix(Y)
  }
  else
    Y <- X
     

  
  is.facY <- sapply(1:ncol(Y), function(t) is.factor(Y[,t])|is.character(Y[,t]))
  is.facY[discY] <- TRUE
  for (t in 1:ncol(Y))
    if (is.facY[t]) 
      Y[,t] <- as.character(Y[,t])
  
  if (dosub) {
    if (nrow(as.matrix(subdat)) != nrow(as.matrix(X)))
      stop("Subgroup data does not match with size of X and Y")
  }
  
  
  splxtoy <- splytox <- FALSE
  if (smooths == "both") {
    splxtoy <- TRUE
    splytox <- TRUE
  }
  
    if (smooths == "xtoy")
    splxtoy <- TRUE
  
  if (smooths == "ytox")
    splytox <- TRUE
  
  
  if (smooth.type == "loess") {
    smoothfunc <- function (X,Y,...) {
      form <- as.formula(Y~X)
      result <- loess(form,...)
      oX <- order(result$x)
      result$x <- result$x[oX]
      result$y <- result$fitted[oX]
      return(result)
    }
  } else if (smooth.type == "spline") {
    smoothfunc <- function (X,Y,...) {
      smooth.spline(X,Y,...)
    }
  } else if (smooth.type == "regline") {
      smoothfunc <- function (X,Y,...) {
        form <- as.formula(Y~X)
        result <- lm(form,...)
        oX <- order(X)
        result$x <- X[oX]
        result$y <- result$fitted[oX]
        return(result)
      }
  } else {
      smoothfunc <- function (X,Y,...) {
        result <- match.fun(smooth.type)(X,Y,...)
        return(result)
      }
    }

 
   if (triplot) {
     
     
    hm <- switch(heatmap,dcov = dcmatrix$dcov.pw, dcor = dcmatrix$dcor.pw, logp = -log10(dcmatrix$pvalue.pw), abscor = abs(dcmatrix$cor), cor = dcmatrix$cor)
     
    if (is.null(hm))
      stop("Pairwise dcor/dcov/cor/p-values required for linked chart.")
    
  
   hm0 <- switch(heatmap0,dcov = dcmatrix$dcov, dcor = dcmatrix$dcor, logp = -log10(dcmatrix$pvalue))
     
    if (is.null(hm0))
      stop("dcmatrix object does not contain the desired entries specified in heatmap0.")
    
  
   openPage( useViewer=FALSE, layout="table2x3" ) 
     
      if (is.null(coldom0))
     coldom0 <- switch(heatmap0, dcov= c(min(min(hm0),0),max(hm0)), dcor = c(min(min(hm0),0),1), logp= c(0,8))

   
   
   if (heatmap0 == "logp")
     hm0[which(hm0>8,arr.ind=TRUE)] <- 8
   
   if (anyNA(hm0) | any(hm0 == Inf)) {
     warning("NAs and Infs were set to 0 for plotting.")
     hm0[is.na(hm0)] <- 0
     hm0[hm0 == Inf] <- 0
   }
     
    lc_heatmap(
      dat(
        colourDomain = coldom0,
        value = hm0,
        palette = palette.heatmap0,
        on_click = function(k) {      #  \  
            grpX <<- gX[[k[1]]]
            grpY <<- gY[[k[2]]]
            updateCharts("heatmap")       #  /
        }
      ),
      place = "A1",
      chartId ="heatmap0"
    )

    
    if (anyNA(hm) | any(hm == Inf)) {
      warning("NAs and Infs were set to 0 for plotting.")
      hm[is.na(hm)] <- 0
      hm[hm == Inf] <- 0
    }
    
    if (is.null(coldom))
      coldom <- switch(heatmap, dcov= c(min(min(hm),0),max(hm)), dcor = c(min(min(hm),0),1), logp= c(0,8), abscor = c(0,1), cor = c(-1,1))

    

    grpX <- gX[[1]]
    grpY <- gY[[1]]
  
    lc_heatmap(
      dat(
        colourDomain = coldom,
        value = matrix(hm[grpX,grpY],nrow=length(grpX),ncol=length(grpY)),
        palette = palette.heatmap,
        on_click = function(k) {        
          sampleX <<- grpX[k[1]]        
          sampleY <<- grpY[k[2]]
          cc <<- intersect(which(complete.cases(X[,sampleX])), which(complete.cases(Y[,sampleY])))
          jittX <<- is.facX[sampleX] * jitt.disc + (1-is.facX[sampleX]) * jitt.cont
          jittY <<- is.facY[sampleY] * jitt.disc + (1-is.facY[sampleY]) * jitt.cont
          if (docg) {
            tab <<- tabnew(X[cc,sampleX],Y[cc,sampleY])
            updateCharts("cgtable")
          }
          if (splxtoy) {
            spl <- smoothfunc(X[cc,sampleX], Y[cc,sampleY],...)
            splX <<- spl$x
            splY <<- spl$y
          }
          if (splytox) {
            spl2 <- smoothfunc(Y[cc,sampleY], X[cc,sampleX],...)
            spl2X <<- spl2$y
            spl2Y <<- spl2$x
          }
          updateCharts( "A3" )      
        }
      ),
      place = "A2",
      chartId ="heatmap"
    )
    
    sampleX <- 1
    sampleY <- 1
    cc <- intersect(which(complete.cases(X[,1])), which(complete.cases(Y[,1])))
    jittX <- is.facX[1] * jitt.disc + (1-is.facX[1]) * jitt.cont
    jittY <- is.facY[1] * jitt.disc + (1-is.facX[1]) * jitt.cont
    
    if (docg)
      tab <- tabnew(X[cc,1],Y[cc,1])
    
    if (splxtoy) {
      spl <- smoothfunc(X[cc,sampleX],Y[cc,sampleY],...)
      splX <<- spl$x
      splY <<- spl$y
    }
    
    
    if (splytox) {
      spl2 <- smoothfunc(Y[cc,sampleY], X[cc,sampleX],...)
      spl2X <<- spl2$y
      spl2Y <<- spl2$x
    }
    


      lc_points(
        dat(
          x = X[cc,sampleX],
          y = Y[cc,sampleY] ,
          size = size,
          opacity = opacity,
          colourValue = colgr[cc],
          jitterX = jittX,
          jitterY = jittY
        ),
        place = "A3"
      )

    
      if (docg)
        lc_html(dat(content = tab), place = "B1", chartId = "cgtable")  
      
      if (dosub)
        lc_input(type = "radio", labels = c("none",subs), on_click = function(value) {
          if (value == 1) {colgr <<- rep(1,n)}
          else {colgr <<- subdat[,subs[value-1]]
          }
          updateCharts( "A2")
        }, value = 1,place = "B2")
    

    
    if (docg)
      lc_html(dat(content = tab), place = "B1", chartId = "cgtable")  
    
    if (dosub)
      lc_input(type = "radio", labels = c("none",subs), on_click = function(value) {
        if (value == 1) {colgr <<- rep(1,n)}
        else {colgr <<- subdat[,subs[value-1]]
        }
        updateCharts( "A2")
      }, value = 1, place = "B2")
    
    colsplx  <- "blue"
    
    if (splxtoy) {  
      
      lc_line(
        dat(
          x = splX,
          y = splY,
          color = colsplx
          # on_mouseover = function(d) {
          #   colsplx <<- "red"
          #   updateCharts("A3")
          # },
          # on_mouseout = function() {
          #  colsplx <<- "blue"
          #  updateCharts("A3")
          # },
          # on_click = function(d) {
          #   splX <<- min(X[,sampleX])
          #   splY <<- min(Y[,sampleY])
          #   updateCharts("A3")
          # }
        ),
        place = "A3",
        addLayer = T
      )
    } 
    
    colsply <- "goldenrod"
    
    if (splytox) {  
      
      lc_path(
        dat(
          color = colsply,
          #on_mouseover = function(d) {
          #  colsply <<- "red"
          #  updateCharts("A3")
          #},
          #on_mouseout = function() {
          #  colsply <<- "goldenrod"
          #  updateCharts("A3")
          #},
          #on_click = function(d) {
          #  spl2X <<- min(X[,sampleX])
          #  spl2Y <<- min(Y[,sampleY])
          #  updateCharts("A3")
          #},
          x = spl2X,
          y = spl2Y
        ),
        place = "A3",
        addLayer = T
      )
    } 
    
  } else {  
    
    #X <- X[,group.X]
    #Y <- Y[,group.Y]
   
    hm <- switch(heatmap,dcov = dcmatrix$dcov, dcor = dcmatrix$dcor, logp = -log10(dcmatrix$pvalue), abscor = abs(dcmatrix$cor), cor = dcmatrix$cor)
    
    
    if (is.null(hm))
      stop("dcmatrix object does not contain the desired entries specified in heatmap0.")
    
    if (anyNA(hm) | any(hm == Inf)) {
      warning("NAs and Infs were set to 0 for plotting.")
      hm[is.na(hm)] <- 0
      hm[hm == Inf] <- 0
    }
      
  
    if (is.null(coldom))
      coldom <- switch(heatmap, dcov= c(min(min(hm),0),max(hm)), dcor = c(min(min(hm),0),1), logp= c(0,8), abscor = c(0,1), cor = c(-1,1))

    openPage( useViewer=FALSE, layout="table2x2" )
    
    lc_heatmap(
      dat(
        colourDomain = coldom,
        value = hm,
        palette = palette.heatmap,
        on_click = function(k) {     
          sampleX <<- k[1]         
          sampleY <<- k[2]
          cc <<- intersect(which(complete.cases(X[,sampleX])), which(complete.cases(Y[,sampleY])))
          jittX <<- is.facX[sampleX] * jitt.disc + (1-is.facX[sampleX]) * jitt.cont
          jittY <<- is.facY[sampleY] * jitt.disc + (1-is.facY[sampleY]) * jitt.cont
          if (docg) {
            tab <<- tabnew(X[cc,sampleX],Y[cc,sampleY])
            updateCharts("cgtable")
          }
          if (splxtoy) {
            spl <- smoothfunc(X[cc,sampleX], Y[cc,sampleY],...)
            splX <<- spl$x
            splY <<- spl$y
          }#  |  charts
          if (splytox) {
            spl2 <- smoothfunc(Y[cc,sampleY], X[cc,sampleX],...)
            spl2X <<- spl2$y
            spl2Y <<- spl2$x
           }
          
         # oX <- order(splX)
          #splX <- splX[oX]
          #splY <- splY[oX]
          updateCharts( "A2" )       #  /
        }
      ),
      place = "A1",
      chartId ="heatmap"
    )
    
  
    #lc_colourSlider(chart = "heatmap")
    
    sampleX <- 1
    sampleY <- 1
    cc <- intersect(which(complete.cases(X[,1])), which(complete.cases(Y[,1])))
    
    jittX <- is.facX[1] * jitt.disc + (1-is.facX[1]) * jitt.cont
    jittY <- is.facY[1] * jitt.disc + (1-is.facX[1]) * jitt.cont
 
    if (docg)
      tab <- tabnew(X[cc,1],Y[cc,1])
   
      
    
    if (splxtoy) {
      spl <- smoothfunc(X[cc,sampleX],Y[cc,sampleY],...)
      splX <<- spl$x
      splY <<- spl$y
    }
    
    
    if (splytox) {
      spl2 <- smoothfunc(Y[cc,sampleY], X[cc,sampleX],...)
      spl2X <<- spl2$y
      spl2Y <<- spl2$x
    }
    
  # #  if (jitt > 0) {
  #   mX <- min(X[cc,sampleX])
  #   MX <- max(X[cc,sampleX]) 
  #   mY <- min(Y[cc,sampleY])
  #   MY <- max(Y[cc,sampleY])
  #   if (jittrel) {
  #     jittX <- (MX-mX)*jitt
  #     jittY <- (MY-mY)*jitt
  #   }
  #   }

    # if (jittdom) {
    #   lc_points(
    #     dat(
    #       x = X[cc,sampleX],
    #       y = Y[cc,sampleY] ,
    #       size = size,
    #       opacity = opacity,
    #       colourValue = colgr[cc],
    #       jitterX = jittX,
    #       jitterY = jittY,
    #     layerDomainX = c(mX-jittX,MX+jittX),
    #       layerDomainY = c(mY-jittY,MY+jittY)
    #     ),
    #     place = "A2"
    #   )
    # } else {
    #   lc_points(
    #     dat(
    #       x = X[cc,sampleX],
    #       y = Y[cc,sampleY] ,
    #       size = size,
    #       opacity = opacity,
    #       colourValue = colgr[cc],
    #       jitterX = jittX,
    #       jitterY = jittY
    #     ),
    #     place = "A2"
    #   )
   #}
    
    
    
      lc_points(
        dat(
          x = X[cc,sampleX],
          y = Y[cc,sampleY] ,
          size = size,
          opacity = opacity,
          colourValue = colgr[cc],
          jitterX = jittX,
          jitterY = jittY
        ),
        place = "A2"
      )

      if (docg)
        lc_html(dat(content = tab), place = "B1", chartId = "cgtable")  
      
      if (dosub)
        lc_input(type = "radio", labels = c("none",subs), on_click = function(value) {
          if (value == 1) {colgr <<- rep(1,n)}
          else {colgr <<- subdat[,subs[value-1]]
          }
          updateCharts( "A2")
        }, value = 1,place = "B2")
      
    
  colsplx <- "black"
    
  if (splxtoy) {  
    
    lc_line(
      dat(
        x = splX,
        y = splY,
        color = colsplx,
        on_click = function(i) {
          colsplx <<- "white"
          updateCharts("A2")
        }
      ),
      place = "A2",
      addLayer = T
    )
    
  } 
  
  
  if (splytox) {  
    
    lc_path(
      dat(
        x = spl2X,
        y = spl2Y
      ),
      place = "A2",
      addLayer = T
    )
  }
}
    
}