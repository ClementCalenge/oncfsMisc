\name{kuomallick}
\alias{kuomallick}
\alias{summary.kuomal}
\alias{plot.kuomal}
\alias{print.kuomal}
\alias{removeBurnin}
\alias{evalconv.effects}
\title{
  Kuo-Mallick approach to the variable selection problem in a Bayesian
  generalized linear model. 
}
\description{
  The function \code{kuomallick} implements the Kuo-Mallick approach to
  the problem of variable selection in a Bayesian generalized linear
  model.  The function \code{removeBurnin} can be used to remove the
  first iterations before convergence of the MCMC ("burn-in" sample).
  The functions \code{plot} and \code{evalconv.effects} can be used to
  assess the convergence of the MCMC.
}
\usage{
kuomallick(y, X, nsim = 10000, theta = "auto", sd = 100,
          type = c("binomial", "poisson", "gaussian"),
          gaussSigma2 = var(y), listevar = NULL)

\method{summary}{kuomal}(object, \dots)

\method{print}{kuomal}(x, \dots)

removeBurnin(x, n=500)

\method{plot}{kuomal}(x, \dots)

evalconv.effects(x, ngroups = 5)
}
\arguments{
  \item{y}{a numerical vector containing the response variable.}
  \item{X}{a data.frame containing the explanatory variables.}
  \item{nsim}{the number of MCMC iterations.}
  \item{theta}{either the character string \code{"auto"} or a numeric
    vector.  If \code{"auto"}, the starting values of the coefficients
    ki (see details) are set to 1, and the starting values of the
    coefficients ai are randomly sampled in the prior distribution of
    these parameters.  Alternatively, the starting values of these
    parameters can be passed in a vector containing the starting values in the
    following order: (i) the coefficients ki in the same
    order as in \code{X} (or in \code{listevar}, if specified), (ii)
    the intercept, and (iii) the coefficients ai (see details).}
  \item{sd}{the standard deviation of the prior normal distribution on
    the coefficients ai (see details).}
  \item{type}{the type of fitted model (binomial, poisson and gaussian
    models are for the moment the only families available)}
  \item{gaussSigma2}{when \code{type="gaussian"}, the starting value for
    the residual variance}
  \item{listevar}{How the columns of X define the variables to be
    considered in the modeling approach (see details).}
  \item{x,object}{an object of class \code{"kuomal"} returned by the
    function \code{kuomallick}}
  \item{n}{the number of iterations to be removed (part of the burn-in
    sample)}
  \item{ngroups}{the number of groups to be considered to evaluate the
    convergence of the chain (see details)}
  \item{\dots}{additionnal arguments to be passed to or from other
    functions}
}
\details{
  The approach of Kuo and Mallik allows to account for the uncertainty
  on model selection when selecting the variables. This approach consist
  to model the response variable using the following model: 
  
  link(y) = a0 + a1*k1*X1 + a2*k2*X2 + a3*k3*X3 + ... + ap*kp*Xp

  where link() is a link function (logit for binomial type, log for
  poisson type, identity for gaussian). The coefficients a1, a2, a3,
  ..., ap are the coefficients of the logistic regression associated
  with the variables X1...Xp stored in the table \code{X}. The
  coefficients k1,...kp are parameters that can take the value 1 or 0:
  when the variable Xi is in the model, the coefficient ki takes the
  value 1, and 0 otherwise.
  
  What is interesting with this approach is that it can be used to
  estimate the probability (sensu Bayes) that a given variable affects
  the response variable.  With the MCMC (Markov chain Monte
  Carlo), it is indeed possible to simulate values for the coefficients
  ai and ki from the probability distribution of these coefficients. The
  coefficients ki can take only two values (0/1), then MCMC allows to
  estimate the probability that a given variable Xi has an effect on the
  response (for a given variable Xi, this probability is estimated with
  the proportion of values ki=1 generated with the MCMC).  This approach
  also allows to identify the association between the variables (i.e. if
  a given variable Xi has an effect on the response only if the variable
  Xj is also in the model, see examples).

  This model is a Bayesian model.  A normal prior is defined for the
  coefficients a1, ..., ap, with a standard deviation defined by the
  parameter \code{sd}. Equal prior probabilities of 0.5 are defined for
  the two possible outcomes of the coefficients ki (either ki=0 or
  ki=1). The model is fitted with a Gibbs sampler, using adaptive
  regression sampling to sample in the complete conditional
  distribution.

  To check for the convergence of the chain, the usual methods can be
  used for the coefficients ai (plotting, \code{raftery.diag} in the
  package coda, etc.).  However, this is more difficult for the
  coefficients ki.  We provide the function \code{evalconv.effects} to
  evaluate this convergence. This function splits the chain into
  \code{ngroups} groups of consecutive iterations, and calculates, for
  each variable i and each group, the proportion of iterations for which
  ki=1. If the mixing properties are slow, this will appear as a slow
  change of these proportions with the iteration number (see
  examples).

  The argument \bold{X} of the function \code{kuomallick} is a
  data.frame that contains the variables X1, X2, ..., Xp in
  columns. By default, the function \code{kuomallick} first preprocess
  this table to transform it in a suitable way for the fit. The
  resulting data.frame \code{Xdesign} is then used for the fit.  Thus,
  if a variable Xi of the table \code{x} is a factor, each category j of
  this variables is automatically transformed into a dummy variable Xij
  prior to the fit. The whole set of coefficients aij associated to
  the variables Xij derived from the factor Xi is multiplied by a unique
  coefficient ki in the model. Internally, the function calculates the
  element \code{listevar} associated to this data.frame \code{Xdesign}:
  each element of \code{listevar} correspond to one column of \code{X}
  and is a vector containing the indices of the columns of
  \code{Xdesign} corresponding to this variable. Note that it is
  possible to build the data.frame \code{Xdesign} and the list
  \code{listevar} prior to the use of the functions (see examples).
}
\value{
  The function \code{kuomallick} returns an object of class
  \code{"kuomal"}, i.e. a list with the following components:
  \item{resu}{a matrix with \code{nsim} rows and as many columns as
    there are parameters (for each variable Xi, two coefficients ai and
    ki -- as well as the standard deviation of the response in the case
    of a Gaussian model)}
  \item{y}{The response variable}
  \item{Xdesign}{The matrix of variables used in the model. When the
  variables in X are all numerical, this matrix is identical to
  \code{X}. When X contains factor variables, each category is
  associated to a dummy variable in X}
  \item{thetaori}{The starting values of the parameters}
  \item{sd}{The standard deviation of the prior on the parameters ai}
  \item{burnin}{The number of iterations removed from \code{resu} as
    part of the burnin sample}
  \item{listevar}{a list with each element corresponding to one variable
  in the model. Each element contains the index of the column(s)
  defining the variable in \code{Xdesign}. When the variable is numeric,
  the corresponding element contains only one number. When the variable
  is a factor, the corresponding element contains as many numbers as
  there are dummy variables coding for the factor (in general, the
  number of categories minus 1).}
}
\references{
  Kuo, L. and Mallick, B. 1998. Variable selection for regression
  models. \emph{Sankya: The Indian Journal of Statistics}, 60: 65--81.
}
\author{
  Clement Calenge \email{clement.calenge@oncfs.gouv.fr}
}
\examples{

\dontrun{

#######################################################
##
##  First example.
##
## To make the calculation reproducible
set.seed(123)

## We set up the data
fert <- swiss$Fertility
X <- swiss
X$Fertility <- NULL

## We fit the model
lm1 <- kuomallick(fert, X, type="gaussian")

## Remove 500 observations as burnin sample
lm1 <- removeBurnin(lm1, 500)

## Convergence of the coefficients ai
plot(lm1)

## Infant mortality is weird: a use of the raftery approach:
raftery.diag(lm1$resu)

## Do not pay attention to the coefficients ki (name starting with
## "Effect."), as the raftery diagnostic is not valid in this case. 
## For Intercept, Education, Infant.Mortality, the mixing is slow.

## Examinate the convergence of the coefficients ki
evalconv.effects(lm1)

## It seems ok, except for Agriculture (slow decrease in the
## probability).
## In theory, we should increase the number of iterations
## in the function kuomallick... Here, we just demonstrate this 
## function, so that we do not really care.

lm1

## Based on the results obtained, we can conclude that:
## The Infant mortality has a probability very close to 1 to be in
## the "true" model predicting the fertility (estimated to be equal to
## 1, actually). The Education has 87\% of chances to be in the "true"
## model.
##
## Looking at the "Best models", we can see that the Best model, having
## 63\% of chances to be the true model, includes only these two
## variables. Note that the next "Best" model also contains the
## proportion of catholic, but has only 13\% of chances to be the "true"
## model (the variable "catholic" alone has 16\% to be in the true
## model).

#######################################################
##
##  Same example, but we illustrate the case where one
##  variable is categorical

set.seed(123)

## We set up the data
fert <- swiss$Fertility
X <- swiss
X$Fertility <- NULL

## We transform the Agriculture into a factor with four categories
X$Agriculture <- cut(X$Agriculture, 4)

## We fit the model
lm2 <- kuomallick(fert, X, type="gaussian")

## Show the results
head(lm2$resu)

## Note that there is only one coefficient ki associated to the
## agriculture, even if there are 3 coefficients ai
## The coefficients ai correspond to the variables in X
## Notice the difference between X:
head(X)

## And Xdesign in the results
head(lm2$Xdesign)

## The indices of the columns in Xdesign corresponding to each variable
## are stored in the list listevar
lm2$listevar


#######################################################
##
##  Same example, but we illustrate the case where we 
##  prepare the data prior to the fit

set.seed(123)

## We set up the data
fert <- swiss$Fertility
X <- swiss
X$Fertility <- NULL

## We transform the Agriculture into a factor with four categories
X$Agriculture <- cut(X$Agriculture, 4)

## So we consider, again the table X where the variable agriculture
## is categorical
head(X)

## we use the function model.matrix to build the design matrix
Xdesi <- model.matrix(glm(fert~., data=X))
head(Xdesi)

## We remove the intercept (it will be added by the function)
Xdesi <- Xdesi[,-1]

## We then build a list listevar indicating which column of Xdesi
## corresponds to each variable
li <- list(Agri=c(1:3), Exam=4, Educ=5, Catho=6, Mort=7)

## We fit the model
lm3 <- kuomallick(fert, as.data.frame(Xdesi), type="gaussian", listevar=li)

## You can compare the results with lm2 (they are the same)

## Note that we can define starting values for the parameters, for example:
## 1. The coefficients ki are all set to 1
## 2. The intercept, and all the coefficients ai are set to 0
initval <- c(rep(1, ncol(X)), 0, rep(0, ncol(Xdesi)))

lm4 <- kuomallick(fert, as.data.frame(Xdesi), type="gaussian",
                  listevar=li, theta=initval)

lm4

## By the way, note that changing the initial values changes the
## conclusions in comparison to lm3, lm2, etc.. This confirms that the
## fit of the model used here is problematic (the readeur could try a
## longer run, by increasing nsim in kuomallik)
## reach a better convergence)
## Warning!!!! Very slow approach!!!
lm5 <- kuomallick(fert, as.data.frame(Xdesi), type="gaussian",
                  listevar=li, theta=initval, nsim=100000)

}

}
\keyword{models}
