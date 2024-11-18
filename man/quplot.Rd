\name{quplot}
\alias{quplot}
\title{
  Draw a quantile plot
}
\description{
  This function draws a quantile plot of a variable.
}
\usage{
quplot(x, ylab = "", main = "", \dots)
}
\arguments{
  \item{x}{a numeric vector for which the quantile plot is desired.}
  \item{ylab}{character string.  The y label.}
  \item{main}{character string.  The title of the plot}
  \item{\dots}{additional arguments to be passed to the function \code{plot}}
}
\details{
  Although too rarely used, quantile plots are a simple yet effective
  way to display the distribution of a variable (see Cleveland 1993, section 2.1).
}
\references{
  Cleveland, W. (1993) \emph{Visualizing data}. Hobart Press.
}
\author{
  Clement Calenge \email{clement.calenge@ofb.gouv.fr}
}
\examples{

## Quantile plot characterizing different distribution
par(mfrow=c(2,2))

## Normal distribution (symmetric distribution)
quplot(rnorm(1000), ylab="Value of the variable",
       main="Normal distribution")

## Student distribution: notice the difference with
## the normal distribution (leptokurtic = more "outliers"
quplot(rt(1000, 3), ylab="Value of the variable",
       main="Student distribution (3 df; leptokurtic)")

## Right-skewed distribution: chi-squared with 1 df
quplot(rchisq(1000, 1), ylab="Value of the variable",
       main="Chi-squared distribution (1 df)")

## Uniform distribution (platykurtic)
quplot(runif(1000), ylab="Value of the variable",
       main="Uniform distribution")


}
\keyword{distribution}
\keyword{aplot}
