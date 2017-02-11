\name{projpurs}
\alias{projpurs}
\alias{plot.oprpu}
\alias{print.oprpu}
\title{
  One dimensional projection pursuit
}
\description{
  \code{projpurs} implements one-dimensional projection pursuit for
  exploratory data analysis.  
}
\usage{

projpurs(df, foo, par = NULL, centsca = TRUE, sphere = FALSE,
         nrep = 100, maxit = 10000,
         violent = FALSE, triesIfViolent=10)

\method{plot}{oprpu}(x, xax = 1, yax = 2, nam="sprpu", summarize = 3, \dots)

\method{print}{oprpu}(x, \dots)

}
\arguments{
  \item{df}{The data.frame to be analysed.}
  \item{foo}{A function that takes two arguments named \code{x} and
    \code{par}, implementing the criterion to be maximized (see details
    and examples).}
  \item{par}{The object \code{par} to be passed to the function \code{foo}.}
  \item{centsca}{logical: whether the data.frame \code{df} should be
    centered and scaled prior to the search.}
  \item{sphere}{logical: whether the cloud of points stored in the
    data.frame \code{df} should be sphericized prior to the search (see
    details).}
  \item{nrep}{an integer value. To evaluate the information in the directions/planes found
    by the algorithm, the function samples \code{nrep} random directions
    in the multidimensional space and calculates the criterion
    \code{foo} on these directions.  The actual value of the direction
    is then compared to these \code{nrep} simulations.}
  \item{maxit}{an integer value. The maximum number of iterations of the algorithm.}
  \item{violent}{logical: whether the "violent" version should be used
    (see details)}
  \item{triesIfViolent}{if \code{violent} is TRUE, the number of random
    vectors used as starting values for the approaches.}
  \item{x}{an object of class \code{"oprpu"} returned by the function
    \code{projpurs}.}
  \item{xax,yax}{The indices of the directions to be plotted as the x
    and y axes.}
  \item{nam}{The name of the grid viewport on which the results will be
    plotted}
  \item{summarize}{By default, the variable names are truncated to
    \code{summarize} characters before plotting.  Alternatively, if
    \code{summarize} is negative, the original variable names will be
    plotted.}
  \item{\dots}{Additionnal arguments to be passed to or from other
    functions.}
}
\details{
  The projection pursuit approach is a dimension reduction approach:
  a given data.frame with N rows and C columns defines a cloud of N
  points in a C-dimensional space.  The projection pursuit algorithm
  seeks to find a linear projections of multivariate data that express
  most "information" in the data.  Of course, the results of the
  algorithm will depend on the criteria used to measure "information"
  (many interesting criteria are described on the help page of
  \code{?pp.criteria}).

  The function \code{projpurs} implements a one-dimensional projection
  pursuit.  In other words, this algorithm seeks to find a
  one-dimensional projection of the data which maximizes the
  information measured by the criterion \code{foo} (see examples for an
  illustration of how these criteria can be defined).  Once a direction
  has been found, the algorithm searches another direction, orthogonal
  to the first one that also maximises the criterion \code{foo}, and a
  third, etc.  In total C-1 dimensions will be found by the algorithm
  (the last one can be deduced from the C-1 first dimensions.

  To find the direction that maximize a criterion, the function
  \code{projpurs} implement the solid angle transform (SAT) described by
  Friedman and Tukey (1973, appendix).  A given direction in the
  C-dimensional space is defined by a normed vector, i.e. a point on a
  sphere.  To find the direction that maximizes the criterion is a
  constrained problem (the solution should be of length one).  Friedman
  and Tukey (1973) indicate that the SAT approach reversibly maps the
  surface of a unit sphere in a C-dimensional space to a
  (C-1)-dimensional unit hypercube (this approach is implemented in the
  functions \code{SAT} and \code{invSAT}).  We use a logit transform to
  reversibly map point in a (C-1)-dimensional infinite Euclidean space
  to a point in a (C-1)-dimensional unit hypercube.  Therefore, the
  constrained problem is transformed into an unconstrained problem,
  which can be solved using "classical" techniques.  Actually, the
  function \code{projpurs} implements the Nelder-Mead algorithm (see
  \code{?optim}).

  Note that the maximum identified by the algorithm may be a local
  maximum (there is no guarantee that the solution found is the actual
  global maximum); however this may still be useful in an exploratory
  context (the fact that a criterion measuring the "information" is
  maximized, even only locally, may bring interesting information on the
  process). To maximize the information on the solution, the
  pre-optimization approach described by Friedman (1987) is carried out
  by the function.

  The starting values for the coefficients a1, a2, ..., aP  that compose
  this vector may therefore have an influence on the solution. By
  default, these starting values are randomly sampled on the unit sphere
  (using the same approach implemented in the function
  \code{mvrspherunif}).  Note that, to maximize the information on the
  resulting axis, it may be interesting to randomly sample several
  directions in the C-dimensional space, then to apply the algorithm,
  and finally to keep the best direction among the tested
  directions.  This is what the function \code{projpurs} does when the
  argument \code{violent} is set to TRUE: it randomly samples
  \code{triesIfViolent} directions in the C-dimensional space, applies
  the PP algorithm on it, and returns the best of all tested
  directions.

  Note that for several criteria, it might be interesting to first
  sphere the data (i.e. to "destroy" the correlation structure of the
  cloud of points). This is done by performing a preliminary principal
  component analysis of the data.frame, and using the criterion
  \code{foo} on the cloud of points defined by the table of standardized
  row coordinates (e.g. see Friedman 1987).
  
  Note that the function implements a randomization approach to evaluate
  the importance of the criterion \code{foo} on the found directions
  with what would have been obtained by sampling randomly \code{nrep}
  directions in the C-dimensional space, and by comparing the observed
  values with the distribution of randomized values.  This test is
  performed by considering the current orthogonality constraints
  (i.e. to evaluate the value of the criterion on the first axis, the
  first value is compared to the distribution of values calculated on
  directions sampled in the C-dimensional space; to evaluate the value
  of the criterion on the second axis, the second value is compared to
  the distribution of values corresponding to the vectors sampled in the
  (C-1)-dimensional space orthogonal to the first axis, etc.).

  The function \code{plot.oprpu} uses functions from the grid package to
  plot a summary of the results.  This function presents, from top to
  bottom and from left to right: (i) violin plots of the randomization
  tests (the violin plots show the distributions of randomized values,
  and the red points correspond to the observed value); (ii) The
  correlations between the variables and the first axis; (iii) the
  elements of the vectors found by the analysis; (iv) the quantile plots
  (see \code{?quplot}) of the two axes; (v) the main window, showing the
  cloud of points.  Note that it is possible to use other
  grid-based functions to customize these results (see the help page of
  \code{?miscgrid} for several possibilities -- see also the examples).
  
}
\value{
  The function \code{projpurs} returns a list of class \code{"oprpu"},
  containing the following elements:
  \item{tab}{The data.frame analyzed by the function (if \code{sphere =
      TRUE}, the data.frame \code{tab} is the data.frame containing the
    standardized row coordinates on the principal components)}
  \item{c1}{The orthonormal vectors found by the analysis}
  \item{li}{The coordinates of the rows of the data.frame \code{tab} on
    the directions defined by \code{c1}}
  \item{cor}{The correlation coefficients between the original variables
    and the row scores on the axes of the analysis}
  \item{criterion}{The function \code{foo}}
  \item{par}{The parameters of the function \code{foo}}
  \item{call}{The call to the function}
  \item{sphere}{Whether the data have been sphered prior to the
    analysis.}
  \item{maxit}{The maximum number of iterations of the Nelder-Mead
    algorithm}
  \item{randtest}{A list with as many elements as there are axes
    (C-1).  Each element is a list with two elements: one element named
    obs containing the value of the criterion on the corresponding
    direction, and one element named sim containing the values of the
    random directions sampled in the subspace explored by the
    algorithm.}
  \item{val.crit}{a vector containing the value of the criterion on the
    axes of the analysis.}
}
\references{
  Friedman, J. and Tukey, J. 1973. A projection pursuit for exploratory
  data analysis. \emph{IEEE Transactions on Computers}, c23: 881--890.

  Friedman, J. 1987. Exploratory projection pursuit. \emph{Journal of
  the American Statistical Association}, 82: 249--266.  
}
\author{
  Clement Calenge \email{clement.calenge@oncfs.gouv.fr}
}
\seealso{
  \code{\link{pp.criteria}} for a list of interesting criteria for the
  projection pursuit.  \code{\link{SAT}} for a further description of
  the solid angle transform.  \code{\link{optim}} for further details on
  the Nelder-Mead approach, \code{\link{miscgrid}} for grid-based
  functions that can be used with the function \code{plot}.
}
\examples{

\dontrun{
################################
##
## First example: we work on the deug dataset (package ade4)

data(deug)
tab <- deug$tab

## We want to find the directions where the Median absolute deviation
## is maximal, i.e. the vector x for which median(abs(x-median(x))) is
## maximal.
## We first write a function foo, depending on two parameters named x
## and par, that calculates this criterion. Note that we need no other
## parameter than x to calculate the MAD, but we still use an argument
## par, since it is a requirement of projpurs
foo <- function(x, par)
     median(abs(x-median(x)))

## We then use the function projpurs to find the directions that
## maximize this criterion:
set.seed(12) ## to make the calculations reproducible
prp <- projpurs(tab, foo)

## Have a look at the results
plot(prp)

## Another interesting criterion is the criterion pp.Friedman1987
## (described in Friedman, 1987). Very slow, but may give
## interesting results (only has a meaning when "sphere=TRUE")
## We use this criterion with, say, J=5
crit <- function(x, par)
    pp.Friedman1987(x, 5)

## projection pursuit again:
prp2 <- projpurs(tab, crit, sphere=TRUE)
plot(prp2)


################################
##
## Second example: we work on the iris dataset

## We want to find the direction where the ratio
## between-class variance/total variance is maximal. We program this
## criterion. This time, we use the Species as parameter par to calculate
## this criterion:

foo2 <- function(x,par) {
    sum(tapply(x,par,function(z) sum(rep((mean(z)-mean(x))^2, length(z)))))/
    sum(tapply(x,par, function(z) sum((z-mean(z))^2)))
}


## And then, we use the projection pursuit, with the argument
## sphere=TRUE
prp3 <- projpurs(iris[,1:4], foo2, par=iris$Species, sphere=TRUE)

## Show the results
plot(prp3)

## This plot is not very informative. It would be more interesting
## to show the classes as the main graph (use of s.classg: this
## plot will be added on the main window
s.classg(prp3$li, prp3$par)

## We can see a clear discrimination between the three species

## By the way, note that this projection pursuit is exactly identical
## to a discriminant analysis:

dis <- discrimin(dudi.pca(iris[,1:4], scannf = FALSE), iris$Species,
                 scannf = FALSE) 

## The row coordinates are identical:
## first axis:
plot(dis$li[,1], prp3$li[,1])

## second axis:
plot(dis$li[,2], prp3$li[,2])

## And the loadings are the same
plot(dis$fa[,1], prp3$c1[,1])
}

}
\keyword{multivariate}
