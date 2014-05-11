\name{batchpca}
\alias{batchpca}

\title{
Compute the PCA of a covariance matrix}

\description{
This function performs the PCA of a covariance matrix, returning an arbitrary number of principal components (eigenvectors) and corresponding eigenvalues.}

\usage{
batchpca(C, k)
}

\arguments{
  \item{C}{covariance matrix.}
  \item{k}{number of PC and eigenvalues to compute.}
}

\details{
If the number of PC to be computed is small (\eqn{k <= ncol(C)/10}), the eigendecomposition is efficiently carried out by the function \code{\link[rARPACK]{eigs_sym}} of the package \code{rARPACK}, which implements an Implicitly Restarted Lanczos Method (IRLM). Otherwise the eigendecomposition is performed with the R base function \code{\link[base]{eigen}}.}

\value{
A list with components
\item{values}{the first \code{k} eigenvaluesof \code{C}.}
\item{vectors}{the first \code{k} PC/eigenvectors of \code{C}.}
}

\references{
\url{http://www.caam.rice.edu/software/ARPACK/}
}

\examples{
n <- 1e4
p <- 500
k <- 10

mat <- matrix(runif(n*p), n, p)
mat <- mat * rep(sqrt(12 * (1:p)), each = n)
C <- cov(mat)
# The eigenvalues of mat should be close to 1, 2, ..., p
# and the corresponding eigenvectors should be close to 
# the canonical basis of R^p

## k first eigenpairs of C using base function eigen
system.time({
  eigenC <- eigen(C, TRUE)
  eigenC$values <- eigenC$values[1:k]
  eigenC$vectors <- eigenC$vectors[,1:k]
  })

## k first eigenpairs of C using function batchpca
system.time(eigenCBatch <- batchpca(C, k))

## Check equality
all.equal(eigenC$values, eigenCBatch$values)
# Reorient eigenvectors if needed
for (i in 1:k) {
  if (sum(eigenCBatch$vectors[,i] * eigenC$vectors[,i]) < 0) 
    eigenCBatch$vectors[,i] <- - eigenCBatch$vectors[,i]
}
all.equal(eigenC$vectors, eigenCBatch$vectors)
}

