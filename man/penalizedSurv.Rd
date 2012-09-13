\name{penalizedSurv}
\alias{penalizedSurv}
\alias{penalizedSurv,ExpressionSet,character-method}
\alias{penalizedSurv,ExpressionSet,Surv-method}
\alias{penalizedSurv,data.frame,Surv-method}
\alias{penalizedSurv,matrix,Surv-method}
\title{Penalized regression using the \code{penalized} R package}
\description{Ridge or Lasso penalized regression.
}
\usage{
\S4method{penalizedSurv}{data.frame,Surv}(X, y, learnind, penalty, unpenalizedvarnames=NULL,essential=FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
   \item{X}{Gene expression data. The following formats are available:
  \itemize{
   \item{\code{matrix}} Rows correspond to observations, columns to variables.
   \item{\code{data.frame}}  Rows correspond to observations, columns to variables.
   \item{\code{ExpressionSet}} survHD will extract the expressions using exprs().
  } }
  \item{y}{Survival Response, an object of class:
  \itemize{
   \item{\code{Surv}} if \code{X} is of type data.frame or matrix
   \item{\code{character}} if \code{X} is of type \code{ExpressionSet}. In
   this case \code{y} may be the name of 
   the survival response in the phenoData of \code{X}. If survival time and indicator are stored separately in the phenoData
   one can specify a two-element character vector the first element representing the survival time variable.
          }}
  \item{learnind}{An index vector specifying the observations that
                  belong to the learning set. May be \code{missing};
                  in that case, the learning set consists of all
                  observations and predictions are made on the
                  learning set.}
  \item{penalty}{Type of penalization: "ridge", "lasso", or
    "elasticnet".}
  \item{unpenalizedvarnames}{Character vector which is a subset of
    colnames(X), naming variables which will be unpenalized.}
    \item{essential}{shall the model output (see slot model of class \code{\link{ModelLearned}})
                  include the bare essentials only (especially for tune this option might save a lot of memory)?}
  \item{\dots}{Further arguments that shall be passed to downstream
    functions.  Potential downstream functions include \code{optL1}, \code{optL2}, and
    \code{penalized} from the penalized package, and \code{opt2D} from the pensim
    package.  If lambda is not specified in the ..., then a tuning
    will be performed using \code{penalized:optL1} (lasso),
    \code{penalized:optL2} (ridge), or \code{pensim:opt2D}.  For ridge
    or lasso, lambda (if specified at all) should be a single numeric
    value.  For elasticnet, lambda (if specified at all) should be a
    numeric vector of length two.  Whether or not tuning is performed,
    regression for the final model is performed using the
    \code{penalized} function from the penalized R package.}
    
}

\value{An object of class \code{\link{ModelLearned}}.}

\references{
Goeman, J.J. (2010) L1 Penalized Estimation in the Cox Proportional Hazards Model, Biometrical Journal, 52, 70-84

}

\author{
  Levi Waldron \email{lwaldron@hsph.harvard.edu}
  Christoph Bernau \email{bernau@ibe.med.uni-muenchen.de}
  Martin Slawski \email{ms@cs.uni-sb.de}
  Anne-Laure Boulesteix \email{boulesteix@ibe.med.uni-muenchen.de}
}

\seealso{\code{\link{learnSurvival}}, \code{\link{penalized}}}

\examples{
\dontrun{
##load example data
if(require(pensim)){
require(survHD)
data(beer.exprs,package="pensim")
data(beer.survival,package="pensim")
##create Surv object (y)
beerY <- with(beer.survival,Surv(os,status))
##create mini X matrix (only using first 200 genes)
beerX <- t(as.matrix(beer.exprs))
beerX <- beerX[,1:200]
colnames(beerX) <- make.names(colnames(beerX),unique=TRUE)
##scale everything for testing.
beerX <- scale(beerX)
##define training and test sets:
set.seed(4)
allind <- 1:nrow(beerX)
learnind <- sample(allind,size=round(length(allind)*9/10))
testind <- allind[-learnind]
##fit a model with ridge regression, and with L2 penalty equal to 100:
fit.penalized <- penalizedSurv(X=beerX,y=beerY, learnind=learnind, penalty="ridge", lambda=100)
##linear risk score
predict(fit.penalized@model,newdata=beerX[testind,],type='lp')
##survival probabilities at each uncensored time
mytimegrid <- sort(c(0.0,beerY[beerY[,2]==1,1]))
SurvivalProbs <- predict(fit.penalized,beerX[testind,],type='SurvivalProbs',timegrid=mytimegrid,gbm=TRUE)@SurvivalProbs
plot(SurvivalProbs,xlab="Time", ylab="Probability of Survival")
}
}
}
\keyword{multivariate}
