\name{clock}
\alias{clock}
\alias{startclock}
\alias{stopclock}
\title{
  Displays a clock
}
\description{
  \code{startclock} indicates that we want to start a
  clock.\code{clock} displays the clock.  \code{stopclock} stops the
  clock. The clock is working on the background.
}
\usage{
startclock()
clock()
stopclock()
}
\author{
  Clement Calenge \email{clement.calenge@oncfs.gouv.fr}
}
\examples{
\dontrun{
## Starts a clock
startclock()
clock()

## Stops the clock
stopclock()
}
}
\keyword{ misc }
