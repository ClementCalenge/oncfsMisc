\name{spm}
\alias{spm}
\title{
  Scatterplot Matrix
}
\description{
  This function shows a scatterplot matrix of a data.frame (quantile
  plot on the diagonal).
}
\usage{
spm(x, density = FALSE, ajoutQ = FALSE, add0 = TRUE, avaxes = TRUE)
}
\arguments{
  \item{x}{
    A data.frame
  }
  \item{density}{
    logical value. Whether a 2D kernel smoothing of the scatterplot
    should be added.
  }
  \item{ajoutQ}{
    logical value. Whether the quartiles should be added on the plot.
  }
  \item{add0}{
    logical value. Whether the axes x=0 and y=0 should be plotted.
  }
  \item{avaxes}{
    logical value. Whether axes labels should be plotted.
  }
}
\details{
  This function shows a scatterplot matrix.  The quantile plots (with
  values of quantiles 25, 50 and 75 percent shown) are displayed on the
  diagonal. 
}
\value{
  a plot.
}
\author{
  Clement Calenge \email{clement.calenge@ofb.gouv.fr}
}
\examples{

spm(iris[,1:4])

spm(iris[,1:4], density=TRUE, avaxes=FALSE)

}
\keyword{hplot }
