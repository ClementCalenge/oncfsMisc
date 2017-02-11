\name{maxent}
\alias{maxent}
\alias{maxentcore}
\alias{prepmaxent}
\alias{print.maxent}
\alias{coef.maxent}
\alias{predict.maxent}

\title{
  Maximum entropy modelling of species distribution
}
\description{
  These functions can be used to fit MaxEnt models to species
  distribution data.  \code{prepmaxent} converts the environment
  variables into "features" (see details).  \code{maxentcore} fits the
  MaxEnt model per se.  The function \code{predict} can be used to
  predict the relative probability of presence in new units, and the
  function \code{coef} can be used to retrieve the coefficients
  associated to the features.
}
\usage{
prepmaxent(dfo, thr = NULL, lithr = NULL, nthr = 10, rare = 0.1,
           PreviousObject = NULL)

maxentcore(df, y, beta = 0, betaj = NULL,
           verbose = FALSE, critconv = 0.0000000001,
           maxIter = 10000, typealgo = c("BFGS","Sequential"))

\method{print}{maxent}(x, \dots)

\method{predict}{maxent}(object, newdata = NULL,
                         type=c("link","response"), \dots)

\method{coef}{maxent}(object, \dots)
}
\arguments{
  \item{dfo}{A data.frame containing the explanatory variables to be
    used in the model.}
  \item{thr}{Optionally, a vector containing the column numbers of the variables in
    \code{dfo} that should be transformed into threshold
    features (if null, no variable will be transformed).}
  \item{lithr}{Optionnally, a list with the same length as \code{thr},
    where each element contains a vector storing the thresholds to be
    used for the corresponding element in \code{thr}.}
  \item{nthr}{numeric value.  The number of thresholds to define for each variable.} 
  \item{rare}{numeric value.  The function will delete the binary
    features which are equal to 1 for less than \code{100*rare}\% units.} 
  \item{PreviousObject}{an object returned by a previous call to
    \code{prepmaxent}, containing all the required elements to allow the
    transformation of variables into features (see details and examples).} 
  \item{df}{The data.frame containing the P features to be used in the
    fit (only numeric variables should be contained in this data.frame)
    for each one of the N spatial units (presence + context points).}
  \item{y}{A numeric vector with containing the type of units in
    \code{df}. This vector should contain only 1 (presence) and 0
    (context points)}
  \item{betaj}{a numeric vector with the same length as \code{df},
    containing the values of the penalty parameters (see details)}
  \item{beta}{The value of beta to be used as a penalty (see details).
    This parameter will be taken into account only if \code{betaj = NULL}.}
  \item{verbose}{logical.  Whether information on the fitting algorithm
           should be printed.}
  \item{critconv}{The tolerance used to identify convergence
    (when the penalized likelihood at step i differs from the penalized
    likelihood at step i-1 by less than \code{critconv}, we consider
    that there is convergence) }
  \item{maxIter}{Maximum number of iterations of the algorithm}
  \item{typealgo}{a character string indicating the type of algorithm to
    use for the fit (see details)}
  \item{x,object}{an object of class \code{"maxent"}}
  \item{newdata}{optionally, a data frame in which to look for features
    with which to predict.  If omitted, the fitted linear predictors
    are used.}
  \item{type}{the type of prediction required.  The default is on the scale
    of the linear predictors; the alternative "response" is on
    the scale of the response variable.}
  \item{\dots}{additional arguments to be passed to or from other functions}
}
\details{
  The Maxent algorithm has been proposed to estimate the probability
  p(x) that a given resource unit x (e.g. a pixel) is used by a given
  species based on a random sample of presence of the species, as well
  as a sample of available sites on the study area (Phillips et al.,
  2006).  For each site/presence, we know the value of the V
  environmental variables describing the environment (slope, elevation,
  etc.). 
 
  The Maxent approach consist in the search of the maximum entropy
  probability distribution p(x) over all available units x such that:
  
  |e[fj] - m[fj]| <= betaj

  Where m[fj] the expectation of the feature j under the probability
  distribution p(x), e[fj] the sample mean of the variable j calculated
  on the species presence, and betaj a tuning parameter fixed by the
  user (see below).
  
  This distribution is the distribution with the largest uncertainty
  (i.e. placing as few constraints as possible on the distribution)
  satisfying the constraint described above. It is a Gibbs distribution:

  p(x)  = exp(lambda\%*\%fx)/Z
  
  Where \code{p(x)} is the probability that an individual randomly
  sampled in the population is localised in the unit \code{x},
  \code{lambda} is a vector of parameters (with one parameter per
  environmental feature), and \code{fx} is a vector containing the value
  of the environmental features in the unit x.

  In theory, one coefficient betaj is needed for each feature j. If
  \code{betaj} is \code{NULL}, the function \code{maxentcore} will
  calculate one coefficient betaj per feature with a general coefficient
  \code{beta} given by the user, using the simplest formula proposed by
  Phillips and Dudik (2008, p. 163):

  betaj = beta*sd(fj)/sqrt(n)

  where sd(fj) is the standard deviation of the feature j measured on
  the sample of presences and n is the number of presences. The
  coefficient \code{beta} will generally be selected using a
  cross-validation approach (see examples).  Note that Phillips and
  Dudik (2008) propose an alternative approach, consisting to define a
  value of \code{betaj} for each feature type.  Our implementation makes it
  possible to use this approach (the user can pass a vector of
  coefficients \code{betaj} to the function \code{maxentcore}).

  Two possible implementations of the Maxent algorithm are available. On
  one hand, the \code{"Sequential"} algorithm corresponds to the
  approach detailed by Dudik et al. (2004).  The algorithm \code{"BFGS"}
  is an implementation of the BFGS algorithm (see \code{?optim}).  From
  our experience, the BFGS algorithm is generally faster.
  
  Note the difference between "features" and "variables": the function
  \code{prepmaxent} transform environmental variables into features. For
  each continuous variable \code{x} in the data.frame \code{dfo}, the
  function \code{prepmaxent} will include in the data.frame used for
  prediction the following features: (i) the variable itself ("linear"
  feature), (ii) the squared variable (this allows to constrain the
  maximum entropy distribution to have a variance on the variable
  \code{x} corresponding to the variance of \code{x} on observed
  presences, see Phillips et al. 2006), (iii) if specified, threshold
  features, i.e. a set of dummy variables indicating whether the value of \code{x}
  is above a given threshold (the list of thresholds is contained in
  \code{lithr}). Note that, by default, 10 equally spaced thresholds are
  defined for each numeric variable in \code{dfo}. The function
  \code{prepmaxent} will also include among the features "product
  features", i.e. the product between all possible pairs of continuous
  variables.  Note that factor variables are transformed as a set of
  dummy variables (one dummy variable per level of the factor, taking
  the value 1 if the resource unit is characterized by this level of the
  factor, and 0 otherwise) at the beginning of the function
  \code{prepmaxent}, and are then treated as other continuous variables.
  Finally, the function \code{prepmaxent} will check if there are any
  replicated features among the calculated features (e.g., for dummy
  variables, the squared feature and the linear feature will be
  identical, and one of them will be removed).  Note also that, each
  feature will be rescaled between 0 and 1 (i.e., xp =
  (x-min(x))/(max(x)-min(x))) by the function \code{prepmaxent}. The
  function \code{prepmaxent} will include all the required elements to
  reproduce this transformation variable -> features in the attribute
  \code{"ForReuse"} of the resulting object.  Although 
  this attribute is not of direct interest for the user, it will be
  required if the user wishes to perform prediction based on a new
  dataset (see examples).  

}
\value{
  The function \code{prepmaxent} returns a data.frame containing the
  features that will be used by the maxent algorithm, with one
  attribute named \code{ForReuse} containing all the elements required
  to allow the transformation of the variables into features (see
  details). 
  
  The function \code{maxentcore} return an object of class
  \code{"maxent"}, containing the results of the fit.
  
}
\references{
  Phillips, S., Anderson, R. and Schapire, R. 2006. Maximum entropy
  modeling of species geographic distributions. \emph{Ecological
    modelling}, 190: 231--259.
  
  Phillips, S.J. and Dudik, M. 2008. Modeling of species distributions
  with Maxent: new extensions and a comprehensive
  evaluation. \emph{Ecography}, 31: 161--175. 

  Dudik, M., Phillips, S. and Schapire, R. 2004. Performance guarantees
  for regularized maximum entropy density estimation. \emph{Proceedings
  of the 17th Annual Conference on Computational Learning Theory},
  472--486. 
}
\author{
  Clement Calenge \email{clement.calenge@oncfs.gouv.fr}
}
\examples{
\dontrun{
## Load the data
data(lynxjura)

## We store the data in objects
map <- lynxjura$map
locs <- lynxjura$locs

## To illustrate the approach, we transform
## one of the variables as a factor
map@data[,1] <- cut(map@data[,1], 3)
levels(map@data[,1]) <- c("low","medium","high")

## Identify the environmental characteristics where the species is present:
pres <- over(locs, map)
pres <- pres[!is.na(pres[,1]),]

## Sample of 2000 context points
contex <- map@data[sample(1:nrow(map), 2000),]

## rbind the two, and create a vector with 1 (presences) and 0 (context
## points)
df <- rbind(pres,contex)
resp <- c(rep(1,nrow(pres)), rep(0, nrow(contex)))

## Calculate the features from the variables. All the numeric variables
## are transformed into threshold features
pr <- prepmaxent(df, thr=2:4)
head(pr)

## fit maxent model: warning!! This may take several minutes
me <- maxentcore(pr, resp, beta=0.01, verbose=TRUE)

## the model
me

## vector of coefficients associated to the features
coef(me)

## predict and map the results
pr2 <- prepmaxent(map@data, PreviousObject = pr)
prediction <- predict(me, newdata=pr2, type="response")
prediction <- data.frame(prediction)
coordinates(prediction) <- coordinates(map)
gridded(prediction) <- TRUE
image(prediction)

#######################################################################
##
## Here, the coefficient beta has been selected with
## a cross-validation approach.

## We illustrate this approach below
## WARNING: this approach is very slow. We provide the code as it might
## be interesting for the user to be able to reproduce this
## cross-validation approach, but this code took nearly 20 minutes to be
## executed on a PC with an Intel core i5 (note that we also diminished
## the convergence criterion).

## We define 5 random groups of 663 points (the 5th
## group contains 2 additionnal points)
set.seed(9280)
gr <- sample(c(rep(1:5, each=663),5,5))

## A sample of beta values to be tested (here, we consider
## only 5 values, but in practical applications, a larger
## number of values might be tested)
betatest <- c(0,0.01, 0.02,0.03,0.04)

## for each possible value of beta
lires <- lapply(1:5, function(ibeta) {
   ## Information on the value of ibeta
   cat("##### value of beta:", betatest[ibeta], "\n")

   beta <- betatest[ibeta]

   ## and for each group:
   do.call("rbind", lapply(1:5, function (i) {
      cat("group number ", i, "\n")
      ## get the data excluding group i
      prb <- pr[gr!=i,]
      respb <- resp[gr!=i]

      ## fit the model with these data and this value of beta
      meb <- maxentcore(prb, respb, beta=beta, critconv=1e-7)

      ## Calculates the prediction for each unit of the group i, that
      ## has not been used to fit the model
      pred <- predict(meb, newdata=pr[gr==i,])
      data.frame(obs=resp[gr==i], pred=pred)
   }))
})

## summing the log-likelihood for each value of beta:
ll <- sapply(lires, function(x) sum(x[x[,1]>0.5,2]))
plot(betatest, ll, xlab="Value of beta",
     ylab="predicted Log-likelihood", ty="b")

## This indicates that the highest predicted log-likelihood would be
## obtained for beta = 0.03 (which is the value we used above)


## Of course, we could go much further, trying different values for the
## different types of features, as in Phillips and Dudik 2008.

}
}
\keyword{multivariate}
