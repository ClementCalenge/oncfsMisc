\name{exploreClust}
\alias{exploreClust}
\alias{calculeR2}
\title{
  Function to Explore Hierarchical Cluster Analysis
}
\description{
  \code{exploreClust} carries out different hierarchical clustering
  algorithms and explores different related statistics.  The function
  \code{calculeR2} calculates the R2 of a given partition of a
  \code{hclust} object.
}
\usage{

calculeR2(P, hc)

exploreClust(ds, P=NULL, meth=c("clust","R2","break","cut"),
             n=NULL,  hang=-1, \dots)

}
\arguments{
  \item{ds}{The distance matrix to be passed to \code{hclust} }
  \item{P}{Optionally, the data.frame/matrix used to calculate \code{ds}}
  \item{meth}{A character string indicating what is to be displayed (see
    details).}
  \item{hc}{An object of class \code{hclust}}
  \item{n}{The number of groups for the method \code{"cut"}}
  \item{hang}{The parameter \code{hang} to be passed to the function
    \code{plot.hclust}}
  \item{\dots}{Other parameters to be passed to the function
    \code{plot.hclust}} 
}
\details{
  The function \code{exploreClust} uses different clustering algorithms
  (single, complete, average, centroid, mcquitty, Ward, see the help
  page of \code{hclust}) to carry out a hierarchical cluster analysis of
  the dataset P, and display for each algorithm: either the clustering
  (\code{meth="clust"}), how the R-squared changes with the number of
  clusters (\code{meth="R2"}), how the height varies with fusion of
  clusters, to identify a "break" in the decrease (\code{meth="break"}),
  or a plot of the first two columns of \code{P} with the corresponding
  partition of \code{n} groups (\code{meth="cut"}).

  The function \code{calculeR2} calculates the R-squared associated with
  the number of groups in a clustering.
}
\author{
  Clement Calenge \email{clement.calenge@ofb.gouv.fr}
}
\seealso{
  \code{\link{hclust}} for a description of the different clustering
  algorithm.
}
\examples{

## Data
serp <- structure(c(0.047, 0.051, 0.08, 0.068, 0.109, 0.103, 0.146, 0.133, 
0.164, 0.158, 0.146, 0.189, 0.199, 0.174, 0.203, 0.174, 0.182, 
0.199, 0.209, 0.244, 0.278, 0.225, 0.234, 0.254, 0.287, 0.35, 
0.313, 0.34, 0.313, 0.342, 0.373, 0.362, 0.348, 0.366, 0.375, 
0.391, 0.391, 0.42, 0.397, 0.389, 0.381, 0.411, 0.407, 0.434, 
0.409, 0.424, 0.416, 0.276, 0.281, 0.299, 0.328, 0.305, 0.317, 
0.34, 0.332, 0.295, 0.375, 0.35, 0.35, 0.371, 0.405, 0.448, 0.403, 
0.395, 0.428, 0.469, 0.489, 0.448, 0.512, 0.522, 0.487, 0.538, 
0.548, 0.557, 0.563, 0.581, 0.585, 0.596, 0.587, 0.573, 0.581, 
0.557, 0.532, 0.624, 0.655, 0.718, 0.673, 0.696, 0.688, 0.626, 
0.643, 0.622, 0.716, 0.767, 0.741, 0.776, 0.73, 0.751, 0.761, 
0.771, 0.639, 0.674, 0.661, 0.637, 0.622, 0.65, 0.62, 0.596, 
0.587, 0.55, 0.576, 0.506, 0.561, 0.539, 0.519, 0.582, 0.484, 
0.458, 0.491, 0.452, 0.467, 0.498, 0.471, 0.484, 0.478, 0.498, 
0.502, 0.498, 0.482, 0.519, 0.591, 0.517, 0.559, 0.545, 0.569, 
0.526, 0.569, 0.569, 0.596, 0.63, 0.62, 0.622, 0.654, 0.665, 
0.676, 0.63, 0.696, 0.661, 0.68, 0.678, 0.715, 0.709, 0.683, 
0.713, 0.755, 0.741, 0.794, 0.741, 0.778, 0.757, 0.798, 0.765, 
0.787, 0.759, 0.787, 0.77, 0.75, 0.802, 0.748, 0.713, 0.785, 
0.704, 0.672, 0.639, 0.606, 0.565, 0.58, 0.617, 0.72, 0.728, 
0.665, 0.691, 0.752, 0.545, 0.513, 0.493, 0.535, 0.495, 0.508, 
0.545, 0.574, 0.585, 0.532, 0.508, 0.783, 0.761, 0.746, 0.759, 
0.715, 0.739), .Dim = c(100L, 2L), .Dimnames = list(NULL, c("x", 
"y")))

## Distance
ds <- dist(serp)

## hclust
exploreClust(ds, serp)

## Cut at three groups
exploreClust(ds, serp, "cut", 3)

## Show how the height decreases
exploreClust(ds, serp, "break", xlim=c(0,30))

## Shows R-squared
exploreClust(ds, serp, "R2", xlim=c(0,30))


}
\keyword{multivariate}
\keyword{nonparametric}
