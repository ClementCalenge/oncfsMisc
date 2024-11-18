\name{posse1995}
\alias{posse1995}
\title{
  Algorithm of Posse (1995) for two-dimensional projection pursuit
}
\description{
  \code{projpurs} implements two-dimensional projection pursuit for
  exploratory data analysis.  
}
\usage{

posse1995(df, foo, par = NULL, centsca = TRUE, sphere = FALSE, c = 1,
          tol = 1e-8, maxhalf = 50, maxit = 100000,
          nrep = 500, violent = FALSE,
          triesIfViolent = 10)

}
\arguments{
  \item{df}{The data.frame to be analysed.}
  \item{foo}{A function that takes two arguments named \code{x} and
    \code{par}, implementing the criterion to be maximized. Note that
    when the function foo is programmed, it should consider that
    \code{x} is a unnamed *list* with two components, each component
    being the scores of the rows of \code{df} on a given direction.}
  \item{par}{The object \code{par} to be passed to the function \code{foo}.}
  \item{centsca}{logical: whether the data.frame \code{df} should be
    centered and scaled prior to the search.}
  \item{sphere}{logical: whether the cloud of points stored in the
    data.frame \code{df} should be sphericized prior to the search (see
    details).}
  \item{c, tol, maxhalf}{The parameters of the algorithm (see details).}
  \item{maxit}{an integer value. The maximum number of iterations of the algorithm.}
  \item{nrep}{an integer value. To evaluate the information in the planes found
    by the algorithm, the function samples \code{nrep} random planes
    in the multidimensional space and calculates the criterion
    \code{foo} on these directions.  The actual value of the optimal plane
    is then compared to these \code{nrep} simulations.}
  \item{violent}{logical: whether the "violent" version should be used
    (see details)}
  \item{triesIfViolent}{if \code{violent} is TRUE, the number of random
    vectors used as starting values for the approaches.}
}
\details{
  This function implements the algorithm of Posse (1995) for
  two-dimensional projection pursuits.  Whereas the one-dimensional
  projection pursuit consists in the search of *directions* on which a
  given function is maximized, the two-dimensional projection pursuit
  consists in the search of *two-dimensional planes* on which a given
  criterion is maximized.

  The algorithm consists in the random search of a plane in the
  C-dimensional space (where C is the number of columns of the
  data.frame \code{df}).  This algorithm starts with a global search, and
  refines the search progressively.  It starts with a random plane
  defined by two orthogonal vectors alpha and beta.  At each step, a
  random vector v is sampled on the C-dimensional unit sphere (see
  \code{?mvrspherunif}), and two candidates planes (a1,b1) and (a2,b2)
  are computed, based on this random vector: first the "composite"
  vectors a1r = (alpha + c*v), a2r = (alpha - c*v), b1r = (beta - (a1^t
  beta)a1) and b2r = (beta - (a2^t beta)a2) are calculated. Then the
  vectors a1, a2, b1 and b2 are derived from these vectors after
  standardization to unit length.  Note that the value of c determines
  the weight of the random vector in the calculation of the candidate
  planes (it measures the size of the neighborhood).  At each step, the
  algorithm selects the best plane among the planes (alpha,beta),
  (a1,b1) and (a2,b2). After \code{maxhalf} iterations without any
  improvement of the value of the criterion, the value of \code{c} is
  halved and the search for best planes is performed more locally. This
  process is continued until \code{c < tol}.

  Note that the function implements a randomization approach to
  evaluate the importance of the criterion \code{foo} on the plane found
  by the algorithm.  Thus \code{nrep} planes are randomly sampled in the
  C-dimensional space, and the function compares the observed value of
  the criterion for the "best" plane with this distribution of
  randomized values.
}
\note{
  Note that this function does not implement the structure removal
  approach advocated by Posse (1995).
}
\value{
  A list of class \code{"oprpu"} (see \code{?projpurs}).
}
\references{
  Posse, C. 1995. Tools for Two-Dimensional Exploratory Projection
  Pursuit. \emph{Journal of Computational and Graphical Statistics}, 4:
  83--100.
}
\author{
  Clement Calenge \email{clement.calenge@ofb.gouv.fr}
}
\seealso{
  \code{\link{projpurs}} for one dimensional projection pursuit, and for
  plotting functions. \code{\link{pp2.criteria}} for a list of
  interesting criteria
}
\examples{

\dontrun{
################################
##
## First example: we simulate a random dataset:

## To make the calculations reproducible
set.seed(12)

## Consider the following plane:
xy <- cbind(runif(5000, -1,1), runif(5000,-1,1))
xy <- xy[abs(xy[,1])>0.5|abs(xy[,2])>0.5,]
plot(xy)


## We can see that there is a hole in the dataset. Now, add 6 other
## random dimensions
dat <- cbind(xy,do.call(cbind,lapply(1:6, function(x) runif(nrow(xy)))))


## Perform a random rotation of these data
a <- randomRotationMatrix(ncol(dat))
dat2 <- as.data.frame(dat\%*\%a)
dat <- as.data.frame(dat)

## To identify the hole: we define a criterion
## giving a very large weight to scores very different
## from 0
foo <- function(x,par) sum(-1/(x[[1]]^2+x[[2]]^2))

## The projection pursuit
aa <- posse1995(dat2, foo)
plot(aa, label=FALSE)


################################
##
## Second example: we work on the deug dataset (package ade4)
## and we try to reproduce a principal component analysis

data(deug)
tab <- deug$tab

## The criterion: sum of the squared scores
foo <- function(x, par)
    sum(x[[1]]^2+x[[2]]^2)


## We then use the function posse1995 to find the directions that
## maximize this criterion:
set.seed(12) ## to make the calculations reproducible
prp <- posse1995(tab, foo)

## Value of the criterion
prp$val.crit
## 408.93

## Compare with the criterion for the PCA:
pc <- dudi.pca(tab, scannf = FALSE)
foo(pc$li)

## We can have better results. Use the "violent" version
set.seed(12) ## to make the calculations reproducible
prp <- posse1995(tab, foo, violent=TRUE)
prp$val.crit
## 451

## The two approaches give pretty similar results:
par(mfrow=c(1,2))
s.label(prp$li)
s.label(pc$li)


}

}
\keyword{multivariate}
