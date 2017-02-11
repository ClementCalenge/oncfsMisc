\name{oncoptions}
\alias{oncoptions}
\title{Setting options for the oncfsMisc package }
\description{
  The function \code{oncoptions} defines the options for the package
  \code{oncfs}.
}
\usage{
oncoptions(\dots)
}
\arguments{
  \item{\dots}{a list of named elements giving the value of options for
    the package \code{oncfs}}
}
\details{
  This functions sets the value of components of the internal list
  \code{.oncoptions}, which contains the values of the options.   One
  option is currently implemented:
  
  \code{runha}{  logical value indicating whether a clock is running
    (see the help page of the function \code{clock}.}
}
\value{
  a list of options for the package \code{oncfs}
}
\author{ Clement Calenge \email{clement.calenge@oncfs.gouv.fr} }
\examples{
## No examples


}
\keyword{environment}
