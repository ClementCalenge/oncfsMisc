\name{brokenstick}
\alias{brokenstick}
\title{
  A randomization test based on the broken stick distribution to select
  the number of principal axes in factor analyses
}
\description{
  This function performs a randomization test based on the broken stick
  distribution to select the number of principal axes in factor analyses.
}
\usage{
brokenstick(pc, niter = 1000)
}
\arguments{
  \item{pc}{any object of class \code{"dudi"} returned by the functions
    of the package ade4.}
  \item{niter}{The number of iterations of the randomization test.}
}
\details{
  All factor analyses rely consist in the search of a small number of
  "meaningful directions" in a multidimensional space.  How to measure
  the information brought by a given direction depends on the analysis
  (e.g. the principal component analysis measures the information with
  the variance, the linear discriminant analysis measures the
  information with the ratio between-classes variance/total variance,
  etc.).  The first axes of a factor analysis maximize this
  information (the amount of information expressed by a given axis is
  measured by the corresponding eigenvalue). The broken stick approach
  allows to select the best number of axes from a factor analysis.
  
  Consider a factor analysis returning p eigenvalues. These eigenvalues
  are compared with a broken stick distribution to select the number of
  axes (Legendre and Legendre 1998, p. 410). This distribution is
  simulated \code{niter} times, and is represented for each axis.
}
\value{
  A list with two elements:
  \item{obs}{a vector containing the \code{ncolmax}  first observed eigenvalues of the PCA}
  \item{sim}{a matrix containing \code{niter} rows and \code{ncolmax}
    columns, containing the randomized eigenvalues of the PCA}
}
\references{
  Legendre, P. and Legendre, L. (1998)
  \emph{Numerical Ecology}. Elsevier Science: Amsterdam.
}
\author{
  Clement Calenge \email{clement.calenge@oncfs.gouv.fr}
}
\seealso{
  \code{\link{examACP}} for another selection approach of the number of
  PCA axes.
}
\examples{
library(ade4)

data(deug)

dudi.deug <- dudi.pca(deug$tab, scannf = FALSE, nf = 4)

brokenstick(dudi.deug)

}
\keyword{multivariate}
