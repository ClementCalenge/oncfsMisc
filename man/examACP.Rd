\name{examACP}
\alias{examACP}
\title{
  Monte Carlo test to select the number of principal components in Principal
  component analysis.
}
\description{
  \code{examACP} carries out a randomization test to select the number
  of principal component in Principal Component analysis.
}
\usage{
examACP(pc, ncolmax = 3, niter = 100)
}
\arguments{
  \item{pc}{An object of class "pca" returned by \code{dudi.pca}}
  \item{ncolmax}{The number of principal components for which the test
    is to be carried out.}
  \item{niter}{The number of iterations of the randomization.}
}
\details{
  This test is based on the reconstruction of the centered and scaled
  table X analyzed with the PCA (e.g. Gabriel 1971).  Let V be the
  matrix containing the normed eigenvectors of the matrix (X^t)DXQ
  (where D is the diagonal matrix containing the row weights, and Q is
  the diagonal matrix containing the column weights).  Similarly, let U
  be the normed eigenvectors XQ(X^t)D.  It is well known that X = V
  (sqrt(Lambda)) U'
  
  In other words, X = Sum_j (v_j) sqrt(lambda_j) (u_j'). We can
  therefore consider that the table X is a sum of "partial tables" A_j =
  (v_j) sqrt(lambda_j) (u_j'). Thus, each "partial table" A_j
  corresponds to the contribution of the principal component j to the
  reconstruction of the table X.  The function \code{examACP} is based
  on this observation.
  
  Now, let X_(r) = Sum_(j=r)^p (a_j) sqrt(lambda_j) (u_j').
  Therefore, the first eigenvalue of the principal component analysis of
  X_(2) will be identical to the rth eigenvalue of the principal
  component analysis of X.

  To test whether the pattern identified by the rth axis of the PCA is
  significantly different from what would have been obtained if the
  columns of X_(r) were uncorrelated, we can randomize the order of the
  values in each column of X_(r), and perform again a PCA on the
  randomized X_(r).  Repeating this process \code{niter} times, we
  obtain \code{niter} values of the rth eigenvalue under the hypothesis
  of uncorrelated columns of X (conditional on the values of the sth
  eigenvalue with s<r).  If the observed value of lambda_j is outside
  the randomized distribution, this indicates that the j_th principal
  component expresses a nonrandom pattern.\cr

  This function returns a plot showing the actual eigenvalues (black
  dots), as well as the following limits of the distribution of
  randomized values: (i) range (min - max) = light grey area, (ii) 2.5\% 
  and 97.5\% quantiles = grey area, (iii) 5\% and 95\% quantiles = dark
  grey area.
}
\value{
  A list with two elements:
  \item{obs}{a vector containing the \code{ncolmax}  first observed eigenvalues of the PCA}
  \item{sim}{a matrix containing \code{niter} rows and \code{ncolmax}
    columns, containing the randomized eigenvalues of the PCA}
}
\references{
  Gabriel, K.R. (1971) The biplot graphic display of matrices with
  application to principal component analysis. \emph{Ecology},
  \bold{70}, 164-168.
}
\author{
  Clement Calenge \email{clement.calenge@oncfs.gouv.fr}
}
\seealso{
  \code{\link{brokenstick}} for another test.
}
\examples{
library(ade4)

data(deug)

dudi.deug <- dudi.pca(deug$tab, scannf = FALSE, nf = 4)

examACP(dudi.deug)

}
\keyword{multivariate}
\keyword{nonparametric}
