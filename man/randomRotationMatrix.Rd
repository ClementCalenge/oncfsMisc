\name{randomRotationMatrix}
\alias{randomRotationMatrix}
\title{
  Random rotation matrix
}
\description{
  This function samples a random matrix from the group of orthogonal
  matrices, using the methods described by Diaconis and Shahshahani (1987).
}
\usage{
randomRotationMatrix(dim = 3, method = c("A","B"))
}
\arguments{
  \item{dim}{The dimension of the rotation matrix.}
  \item{method}{The method to be used (see details).} 
}
\details{
  Two methods are implemented to sample a random rotation matrix.  The
  method A of Diaconis and Shahshahani (1987) is the classical method: a
  \code{dim x dim} matrix M is filled with random numbers drawn from a
  normal distribution.  Then, a QR decomposition of this matrix is
  performed, factorized such that R has positive diagonal entries (as
  required by theorem 3.2 in Stewart 1980). The matrix Q has the
  distribution mu on On, where On is the group of orthogonal matrices of
  order n and mu is the Haar measure (equivalent of the uniform
  distribution for groups). Ozols (2009) gives a detailed account (and
  accessible to non-mathematicians) of this approach (both why the use
  of a normal distribution leads to a random orthogonal matrix
  characterized by a Haar measure, and how to standardize the QR
  decomposition so that R has positive diagonal entries).

  The method B of Diaconis and Shahshahani (1987) is based on the
  subgroup algorithm.  We refer the reader to this paper for further
  details.
}
\value{
  a matrix with \code{n} columns and \code{n} rows, containing the
  rotation matrix.
}
\references{
  Stewart, G. 1980. The efficient generation of random orthogonal
  matrices with an application to condition estimators. \emph{SIAM
  Journal on Numerical Analysis}, 17, 403--409.
  
  Diaconis P. and Shahshahani, M. 1987. The Subgroup Algorithm for Generating
  Uniform Random Variables. \emph{Probability in the Engineering and
    Informational Sciences}, 1, 15--32.

  Ozols, M. 2009. How to generate a random unitary matrix. Report
  available at the following URL:
  \url{http://home.lu.lv/~sd20008/papers/essays/Random\%20unitary\%20\%5Bpaper\%5D.pdf}
}
\author{
  Clement Calenge \email{clement.calenge@oncfs.gouv.fr}
}
\examples{

## Random matrix in 12 dimensions
mat <- randomRotationMatrix(12)

}
\keyword{multivariate}
