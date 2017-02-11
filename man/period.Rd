\name{period}
\alias{period}
\alias{fitmod}
\title{Calculates a periodogram}
\description{
  \code{period} calculates the periodogram associated with a time series.
}
\usage{
period(yt, ncons = 1L)
fitmod(per, increaseres = 1L)
}
\arguments{
  \item{yt}{the time series for which the periodogram is desired.}
  \item{ncons}{the number of frequencies to be considered (see details).}
  \item{per}{an object of class "period".}
  \item{increaseres}{an integer greater than 0, defining the resolution
    of the model (see details).}
}
\details{
  The function \code{period} implements the periodogram as described by Diggle
  (1990).  The periodogram can be used to identify the frequencies
  explaining most of the variability of the variable \code{yt} (this
  variable is centred and scaled prior to the calculations).

  It is possible to model the time series \code{yt} as the sum of the
  sine and cosine of the \code{ncons} frequencies explaining most
  variability with the function \code{fitmod}.  More precisely, the
  function fits the model:
  y_t = Sum_i (a_i *  cos(omega_i * t) + b_i*sin(omega_i * t)) +
  epsilon_t
  Where the sum is over the \code{ncons} most important Fourier
  frequencies.

  Note that \code{yt} is supposed to be a regular
  time-series. Successive values are separated by 1 unit of time.  If
  \code{increaseres = 1}, the prediction of \code{fitmod} will be
  carried out for every moment corresponding to a measure in \code{yt},
  i.e. for each unit of time. It is possible to have a prediction of
  \code{yt} at a smaller resolution increasing the parameter
  \code{increaseres} (the prediction will be calculated every
  1/increaseres unit of time).  Note that the function \code{fitmod} is
  provided for exploratory plotting purposes only.
  
}
\value{
  The function \code{period} returns an object of class \code{period},
  i.e. a data.frame with two columns:
  \item{fre}{The frequency}
  \item{spe}{the associated spectrum}
  The function fitmod returns a data.frame with
}
\references{
  Diggle, P. (1990) \emph{Time series. A biostatistical
    introduction}. Oxford University Press.
}
\author{
  Clement Calenge \email{clement.calenge@oncfs.gouv.fr}
}
\examples{
## Use the dataset lh and identify the most important frequencies
pe <- period(lh)

## focus on the two most important frequencies
pe <- period(lh, 2)

## plot the corresponding model
plot(fitmod(pe), ty="l")

## Increase the resolution for a smoother prediction
plot(fitmod(pe, increaseres = 100), ty="l", lwd=2)

## add the centred and scaled series on this plot:
lines(1:length(lh), (lh-mean(lh))/sd(lh), col="red")

}
\keyword{ts}
