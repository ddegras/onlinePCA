\name{secularRpca}
\alias{secularRpca}

\title{
Update the PCA using secular equations}

\description{
This function recursively updates a PCA with respect to a new data vector. The PCA eigenvalues are computed as the roots of a secular equations. The principal components are then  either deduced by direct calculation or computed with the approximation method of Gu and Eisenstat (1994).  
}

\usage{
secularRpca(d, Q, x, n, ff, center, tol = 1e-10, reortho = FALSE) 
}

\arguments{
  \item{d}{vector of eigenvalues.} 
  \item{Q}{matrix of corresponding principal components (eigenvectors of the covariance matrix) stored in columns.}
 \item{x}{new data vector.} 
 \item{n}{sample size before new observation.}  
  \item{ff}{forgetting factor; a real number between 0 and 1 set to \code{1/n} by default.}
  \item{center}{centering vector for the new data (optional).}
  \item{tol}{error tolerance for the computation of eigenvalues.}  
\item{reortho}{logical; FALSE by default. If TRUE, the approximation method of Gu and Eisenstat (1994) is used to compute the eigenvectors. 
}
}

\details{
The method of secular equations provides accurate eigenvalues in all but pathological cases. On the other hand, the perturbation method implemented by \code{\link{perturbationRpca}} typically runs much faster but is only accurate for large sample sizes \code{n}. \cr
The default method for computing PC's/eigenvectors (\code{reortho = FALSE}) is that of Li et al. (2000). It is very accurate for the first few PC's but loss of accuracy and orthogonality may occur for the next PC's. In contrast the method of Gu and Eisenstat (1994) is robust against small errors in the computation of eigenvalues. It provides PC's that may be less accurate than the default method  but for which strict orthogonality is guaranteed.\cr
Lower values of the forgetting factor \code{ff} place more weight on the current PCA decomposition while higher values place more weight on the new data. The default value \code{ff = 1/n} corresponds to the stationarity assumption. 
}

\value{
A list with components  
\item{values}{updated eigenvalues in decreasing order}
\item{vectors}{updated principal components}
}

\references{
Gu, M. and Eisenstat, S. C. (1994). A stable and efficient algorithm for the rank-one modification of the symmetric eigenproblem. \emph{SIAM Journal of Matrix Analysis and Applications.}\cr
Li, W., Yue, H. H., Valles-Cervantes, S. and Qin, S. J. (2000). Recursive PCA for adaptive process monitoring. \emph{Journal of Process Control}.
}

\seealso{
\code{\link{perturbationRpca}}
}

%\examples{
%}

