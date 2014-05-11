\name{perturbationRpca}
\alias{perturbationRpca}

\title{
 Update the PCA using a rank 1 perturbation method}

\description{
This function recursively updates the PCA with respect to a single new data vector, using the (fast) perturbation method of Hegde et al. (2006).}

\usage{perturbationRpca(d, Q, x, n, ff, center)}

\arguments{
  \item{d}{vector of eigenvalues.}
  \item{Q}{matrix of principal components (eigenvectors of the covariance matrix) stored in columns.}
  \item{x}{new data vector.}
  \item{n}{current sample size (prior to observation of new data)}
  \item{ff}{forgetting factor: a number between 0 and 1, set to \code{1/n} by default.}
\item{center}{centering vector (optional).}
}

\details{
The forgetting factor \code{ff} determines the balance between past and present observations in the PCA update: the closer it is to 1 (resp. to 0), the more weight is placed on current (resp. past) observations. At least one of the arguments \code{n} and \code{ff} must be specified. If \code{ff} is specified, its value overrides the argument \code{n}; otherwise, \code{ff} is set to \code{1/n} which corresponds to the assumption of a stationnary observation process. If \code{center} is not specified, then the data \code{x} is not centered prior to PCA. 
}

\value{
A list with components 
  \item{values}{the updated eigenvalues.}
  \item{vectors}{the updated eigenvectors.}
}

\references{
Hegde et al. (2006) Perturbation-Based Eigenvector Updates for On-Line Principal Components Analysis and Canonical Correlation Analysis. 
\emph{Journal of VLSI Signal Processing}.}

\note{
This perturbation method is based on large sample approximations. It can be highly inaccurate for small/medium sized samples and should not be used in this case.}

\seealso{
\code{\link{secularRpca}}
}

\examples{
n <- 1e4
n0 <- 5e3
p <- 10

mat <- matrix(runif(n*p), n, p)
mat <- mat * rep(sqrt(12 * (1:p)), each = n)
# The eigenvalues of mat should be close to 1, 2, ..., p
# and the corresponding eigenvectors should be close to 
# the canonical basis of R^p

## Perturbation-based PCA
xbar <- colMeans(mat[1:n0,])
pca <- batchpca(cov(mat[1:n0,]))
for (i in (n0+1):n) {
  xbar <- updateMean(xbar, mat[i,],  i - 1)
  pca <- perturbationRpca(pca$values, pca$vectors, mat[i,], 
    n = i - 1, center = xbar)
}
pca

# Compare to batch PCA 
eigen(cov(mat), TRUE)
}

