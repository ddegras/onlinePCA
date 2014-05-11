\name{sgapca}
\alias{sgapca}

\title{Stochastic Gradient Ascent PCA}

\description{This function updates the PCA with respect to a new data vector, based on the Stochastic Gradient Ascent algorithm of Oja (1992).}

\usage{
sgapca(Q, x, gamma, center, type = c("exact", "nn"))
}

\arguments{
 \item{Q}{matrix of principal components stored in columns.}
  \item{x}{new data vector.}
   \item{gamma}{gain parameter(s).}
 \item{center}{centering vector for \code{x} (optional).}
  \item{type}{string specifying the type of implementation: "exact" or "nn" (neural network).}
}

\details{
The gain scalar or vector \code{gamma} determines the weight placed on the new data in updating each principal component. It is specified as a single positive number or as a vector of the same length as \code{x}. It  For larger values of \code{gamma}, more weight is placed on \code{x} and less on \code{Q}. A common choice for (the components of) \code{gamma} is of the form \code{c/n}, with \code{n} the sample size and \code{c} a suitable positive constant. \cr 
The Stochastic Gradient Ascent PCA can be implemented exactly or through a neural network. The latter is less accurate but faster.
}

\value{
A matrix of updated principal components of the same dimension as \code{Q}. 
}

\references{
Oja (1992). Principal components, Minor components, and linear neural networks. \emph{Neural Networks.}
}

%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{snlpca}}
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

## SGA-PCA
xbar <- colMeans(mat[1:n0,])
pca <- batchpca(cov(mat[1:n0,]))$vectors
for (i in (n0+1):n) {
  xbar <- updateMean(xbar, mat[i,],  i - 1)
  pca <- sgapca(pca, mat[i,], 2 / i, center = xbar)
}
pca

# Compare to batch PCA 
eigen(cov(mat), TRUE)$vectors
}

