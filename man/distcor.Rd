% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/distcov_function.R
\name{distcor}
\alias{distcor}
\title{Calculates the distance correlation \insertCite{szekely2007,szekely2009brownian}{dcortools}.}
\usage{
distcor(
  X,
  Y,
  affine = FALSE,
  standardize = FALSE,
  bias.corr = FALSE,
  type.X = "sample",
  type.Y = "sample",
  metr.X = "euclidean",
  metr.Y = "euclidean",
  use = "all",
  algorithm = "auto"
)
}
\arguments{
\item{X}{contains either the first  sample or its corresponding distance matrix.

In the first case, X can be provided either as a vector (if one-dimensional), a matrix or a data.frame (if two-dimensional or higher). 

In the second case, the input must be a distance matrix corresponding to the sample of interest.

If X is a sample, type.X must be specified as "sample". If X is a distance matrix, type.X must be specified as "distance".}

\item{Y}{see X.}

\item{affine}{logical; specifies if the affinely invariant distance correlation \insertCite{dueck2014affinely}{dcortools} should be calculated or not.}

\item{standardize}{logical; specifies if X and Y should be standardized dividing each component by its standard deviations. No effect when affine = TRUE.}

\item{bias.corr}{logical; specifies if the bias corrected version of the sample distance correlation \insertCite{huo2016fast}{dcortools} should be calculated.}

\item{type.X}{For "distance", X is interpreted as a distance matrix. For "sample", X is intepreted as a sample.}

\item{type.Y}{see type.X.}

\item{metr.X}{specifies the metric which should be used to compute the distance matrix for X (ignored when type.X = "distance").

 Options are "euclidean", "discrete", "alpha", "minkowski", "gauss", "gaussauto", "boundsq" or user-specified metrics (see examples).
 
 For "alpha", "minkowski", "gauss", "gaussauto" and "boundsq", the corresponding parameters are specified via "c(metric,parameter)", c("gaussian",3) for example uses a Gaussian metric with bandwith parameter 3; the default parameter is 2 for "minkowski" and "1" for all other metrics.
 
 See \insertCite{lyons2013distance,sejdinovic2013equivalence,bottcher2017detecting;textual}{dcortools} for details.}

\item{metr.Y}{see metr.X.}

\item{use}{specifies how to treat missing values. "complete.obs" excludes NA's, "all" uses all observations.}

\item{algorithm}{specifies the algorithm used for calculating the distance correlation. 

"fast" uses an O(n log n) algorithm if the observations are one-dimensional and metr.X and metr.Y are either "euclidean" or "discrete", see also \insertCite{huo2016fast;textual}{dcortools}. 

"memsave" uses a memory saving version of the standard algorithm with computational complexity O(n^2) but requiring only O(n) memory. 

"standard" uses the classical algorithm. User-specified metrics always use the classical algorithm.

"auto" chooses the best algorithm for the specific setting using a rule of thumb.}
}
\value{
numeric; the distance correlation between samples X and Y.
}
\description{
Calculates the distance correlation \insertCite{szekely2007,szekely2009brownian}{dcortools}.
}
\examples{
X <- rnorm(200)
Y <- rnorm(200)
Z <- X + rnorm(200)
dim(X) <- dim(Y) <- dim(Z) <- c(20,10)

# Demonstration that biased-corrected distance correlation is often more meaningful than without using bias-correction
distcor(X,Y) 
distcor(X,Z) 
cor(X,Y,bias.corr=T)
distcor(X,Z,bias.corr=T)

# For more examples of the different option, see the documentation of distcov.
}
\references{
\insertRef{bottcher2017detecting}{dcortools}

\insertRef{dueck2014affinely}{dcortools}

\insertRef{huo2016fast}{dcortools}

\insertRef{lyons2013distance}{dcortools}

\insertRef{sejdinovic2013equivalence}{dcortools}

\insertRef{szekely2007}{dcortools}

\insertRef{szekely2009brownian}{dcortools}
}
