\name{plusMinusSurv}
\alias{plusMinusSurv}
\alias{plusMinusSurv,ExpressionSet,character-method}
\alias{plusMinusSurv,ExpressionSet,Surv-method}
\alias{plusMinusSurv,data.frame,Surv-method}
\alias{plusMinusSurv,matrix,Surv-method}
\title{Simple average learner}
\description{The plusminus learner assigns coefficients of +/-1 to each
  feature, optionally scaled by the standard deviation of the feature
  for equal contributions to the continuous risk score.  Signs are
  determined from the sign of the univariate Cox regression for each
  feature.
}
\usage{
\S4method{plusMinusSurv}{data.frame,Surv}(X, y, learnind,essential=FALSE, modeltype="plusminus", ...)
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
   \item{\code{character}} if \code{X} is of type \code{ExpressionSet}. In this case y is the name of 
   the survival response in the phenoData of \code{X}. If survival time and indicator are stored separately in the phenoData
   one can specify a two-element character vector the first element representing the survival time variable.
          }}
 
  \item{learnind}{An index vector specifying the observations that
                  belong to the learning set. May be \code{missing};
                  in that case, the learning set consists of all
                  observations and predictions are made on the
                  learning set.}
                  \item{essential}{shall the model output (see slot model of class \code{\link{ModelLearned}})
                  include the bare essentials only (especially for tune
		  this option might save a lot of memory)?}
		\item{modeltype}{Can be "plusminus",
		  "compoundcovariate", "tscore", "voting",
		  "positiveriskvoting", or "negativeriskvoting".  See
		  \code{plusMinus} for details.}
  \item{\dots}{Further arguments that shall be passed to \code{plusMinus} function}

}

\value{An object of class \code{\link{ModelLearned}}.}

\author{
  Levi Waldron \email{lwaldron@hsph.harvard.edu}
  Christoph Bernau \email{bernau@ibe.med.uni-muenchen.de}
}

\seealso{\code{\link{learnSurvival}}, \code{\link{plusMinus}}}

\examples{
\dontrun{
##load example data
if(require(pensim)){
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
##fit a model with ten non-zero coefficients:
fit.plusminus <- plusMinusSurv(X=beerX,y=beerY, learnind=learnind, lambda=10)
##linear risk score
predict(fit.plusminus@model,newdata=beerX[testind,],type='lp')
##survival probabilities at each uncensored time
mytimegrid <- sort(c(0.0,beerY[beerY[,2]==1,1]))
SurvivalProbs <- predict(fit.plusminus,beerX[testind,],type='SurvivalProbs',timegrid=mytimegrid,gbm=TRUE)@SurvivalProbs
plot(SurvivalProbs,xlab="Time", ylab="Probability of Survival")
}
}
}
\keyword{multivariate}
