\name{pp2.criteria}
\alias{pp2.criteria}
\alias{pp2.holes}
\alias{pp2.centralmass}
\alias{pp2.FriedmanTukey}
\title{
  Several interesting criteria for two-dimensional projection pursuit
}
\description{
  These functions implement several interesting criteria for the two
  dimensional projection pursuit (see Details).
}
\usage{

pp2.holes(x, p)

pp2.centralmass(x, p)

pp2.FriedmanTukey(x)

}
\arguments{
  \item{x}{a list containing two components, each component being a
    vector containing the scores of individuals on a given direction.}
  \item{p}{The dimension of the space in which the criterion is to be
    maximized (i.e. the number of columns of the data.frame passed to
    the function \code{posse1995}).}
}
\details{
  \code{pp.holes} has been used by Cook and Swayne (2007, p. 30) to
  identify the projections with few points in the centre of the plot
  (this criterion is quickly optimized, but will miss a structure if
  groups are not separated at the origin -- this is the case of the
  dataset \code{deug} actually, see examples).\cr
  
  \code{pp.centralmass} has been used by Cook and Swayne (2007, p. 30)
  to identify projections with a lot of point at the origin of the
  plot (this criterion is quickly optimized, and is useful to identify
  outliers).\cr

  \code{pp.FriedmanTukey} is the original criterion of Friedman and
  Tukey (1973).
}
\value{
  A numeric value containing the value of the criterion.
}
\references{
  Friedman, J. and Tukey, J. 1973. A projection pursuit for exploratory
  data analysis. \emph{IEEE Transactions on Computers}, c23: 881--890.

  Friedman, J. 1987. Exploratory projection pursuit. \emph{Journal of
  the American Statistical Association}, 82: 249--266.  
}
\author{
  Clement Calenge \email{clement.calenge@ofb.gouv.fr}
}
\seealso{
  \code{\link{projpurs}} for additional details on the one-dimensional
  projection pursuit algorithm.
}
\examples{
\dontrun{
## Example with the deug dataset
data(deug)
tab <- deug$tab


## Example with pp2.centralmass (p = 9) (par is unused here)
foo <- function(x, par)
     pp2.centralmass(x, 9)
set.seed(12) ## to make the calculations reproducible
prp <- posse1995(tab, foo, sphere=TRUE)
plot(prp)


## Example with pp2.holes
foo <- function(x, par)
     pp2.holes(x, 9)
set.seed(12) ## to make the calculations reproducible
prp <- posse1995(tab, foo, sphere=TRUE)
plot(prp)


## Example with pp2.FriedmanTukey
set.seed(12) ## to make the calculations reproducible
foo <- function(x, par)
     pp2.FriedmanTukey(x)
prp <- posse1995(tab, foo, sphere=TRUE)


}
}
\keyword{multivariate}
