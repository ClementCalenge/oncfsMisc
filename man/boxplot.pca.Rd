\name{boxplot.pca}
\alias{boxplot.pca}
\title{Quick boxplot of variables prior to a PCA}
\description{
  This function allows to draw a boxplot of the variables in a
  data.frame before a PCA (allows to see the symmetry of the variable).
}
\usage{
\method{boxplot}{pca}(x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{an object of class "pca" returned by the function \code{dudi.pca}}
  \item{\dots}{additional arguments to be passed to and from other functions}
}
\author{
  Clement Calenge \email{clement.calenge@oncfs.gouv.fr}
}
\seealso{
  \code{\link[ade4]{dudi.pca}}, ~~~
}
\examples{
library(ade4)

data(deug)

dudi.deug <- dudi.pca(deug$tab, scannf = FALSE, nf = 4)

boxplot(dudi.deug)
}
\keyword{ hplot }
\keyword{ multivariate }
