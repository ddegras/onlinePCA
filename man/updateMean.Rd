\name{updateMean}
\alias{updateMean}

\title{
Update the sample mean vector
}

\description{
This functions updates the sample mean vector with respect to new data.
}

\usage{
updateMean(xbar, x, n, ff, byrow = TRUE)
}

\arguments{
  \item{xbar}{sample mean vector.}
  \item{x}{new data.}
  \item{n}{sample size prior to the observation of new data.}
  \item{ff}{
forgetting factor: a number between 0 and 1, set to \code{1/n} by default.}
\item{byrow}{
logical; should be set to TRUE if the new data vectors are in rows (i.e., variables in columns) and FALSE otherwise.}
}

\details{
The forgetting factor \code{ff} determines the balance between past and present observations in the PCA update: the closer it is to 1 (resp. to 0), the more weight is placed on current (resp. past) observations. At least one of the arguments \code{n} and \code{ff} must be specified. If \code{ff} is specified, its value overrides the argument \code{n}; otherwise, \code{ff} is set to \code{1/n} which corresponds to the assumption of a stationnary observation process.\cr
If unspecified, the argument \code{byrow} defaults to TRUE if \code{x} is a matrix and to FALSE if \code{x} is a vector.
}

\value{
The updated sample mean vector. 
}

\seealso{
\code{\link{updateCovariance}}
}

\examples{
n <- 1e4
n0 <- 5e3
p <- 10
mat <- matrix(runif(n*p), n, p)

## Direct computation of the mean
xbar <- colMeans(mat)

## Recursive computation of the covariance
xbar0 <- colMeans(mat[1:n0,])
xbarrec <- updateMean(xbar0, mat[(n0+1):n,], n0)

## Check equality
all.equal(xbar, xbarrec)
}

