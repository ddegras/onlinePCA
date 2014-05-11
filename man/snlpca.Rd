\name{snlpca}
\alias{snlpca}

\title{Subspace Network Learning PCA}

\description{This function updates the PCA with respect to a new data vector using the Subspace Network Learning algorithm of Oja (1992).}

\usage{
snlpca(Q, x, gamma, center, type = c("exact", "nn"))
}

\arguments{
 \item{Q}{matrix of principal components stored in columns.}
  \item{x}{new data vector.}
   \item{gamma}{gain parameter.}
 \item{center}{centering vector for \code{x} (optional).}
  \item{type}{string specifying the type of implementation: "exact" or "nn" (neural network).}
}

\details{
The gain parameter \code{gamma} is a positive number that determines the weight to be placed on the new data vector \code{x} in the PCA update. 
For larger values of \code{gamma}, more weight is placed on \code{x} and less on \code{Q}. A common choice for \code{gamma} is of the form \code{c/n}, with \code{n} the sample size and \code{c} a suitable positive constant. \cr 
The Subspace Network Learning PCA can be implemented exactly or through a neural network. The latter is less accurate but faster.}

\value{A matrix of principal components/eigenvectors stored in columns.}

\references{
Oja (1992). Principal components, Minor components, and linear neural networks. \emph{Neural Networks.}
}


\note{
Unlike the Stochastic Gradient Ascent algorithm which consistently estimates the principal components (i.e., the eigenvectors of the theoretical covariance matrix), the Subspace Network Learning algorithm only provides the linear space spanned by the PCs and not the PCs themselves. 
}

\seealso{
\code{\link{sgapca}}
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

## SNL-PCA
xbar <- colMeans(mat[1:n0,])
pca <- batchpca(cov(mat[1:n0,]))$vectors
for (i in (n0+1):n) {
  xbar <- updateMean(xbar, mat[i,],  i - 1)
  pca <- snlpca(pca, mat[i,], 1 / i, center = xbar)
}
pca
# Corresponding projection matrix
snlP <- tcrossprod(pca)

# Compare to batch PCA 
(batch <- eigen(cov(mat), TRUE)$vectors)
batchP <- tcrossprod(batch)

# Relative  distance between the two projectors (Frobenius metric)
norm(snlP - batchP, "2") / norm(batchP, "2")
}

