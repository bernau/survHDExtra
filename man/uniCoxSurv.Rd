\name{uniCoxSurv}
\alias{uniCoxSurv}
\alias{uniCoxSurv,ExpressionSet,character-method}
\alias{uniCoxSurv,ExpressionSet,Surv-method}
\alias{uniCoxSurv,data.frame,Surv-method}
\alias{uniCoxSurv,matrix,Surv-method}
\title{Univarate shrinkage prediction in the Cox model}
\description{See
  http://cran.r-project.org/web/packages/uniCox/index.html for more details
}
\usage{
\S4method{uniCoxSurv}{data.frame,Surv}(X, y, learnind, essential=FALSE,...)
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
                  \item{essential}{shall the model output (see slot model of class \code{\link{ModelLearned}})
                  include the bare essentials only (especially for tune this option might save a lot of memory)?}
  \item{\dots}{Further arguments that shall be passed to \code{uniCox}
    function from the uniCox R package (http://cran.r-project.org/web/packages/uniCox/index.html)
  }
  }

\value{An object of class \code{\link{ModelLearned}}.}

\references{Tibshirani RJ: Univariate shrinkage in the cox model for high
  dimensional data. \emph{Stat Appl Genet Mol Biol 2009, 8:Article21.}}

\author{
  Levi Waldron \email{lwaldron@hsph.harvard.edu}
  Christoph Bernau \email{bernau@ibe.med.uni-muenchen.de}
}

\seealso{\code{\link{learnSurvival}}}

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
##fit a model with penalty equal to 10:
fit.unicox <- uniCoxSurv(X=beerX,y=beerY, learnind=learnind, lambda=10)
##linear risk score
predict(fit.unicox@model,newdata=beerX[testind,],type='lp')
##survival probabilities at each uncensored time
mytimegrid <- sort(c(0.0,beerY[beerY[,2]==1,1]))
SurvivalProbs <- predict(fit.unicox,beerX[testind,],type='SurvivalProbs',timegrid=mytimegrid,gbm=TRUE)@SurvivalProbs
plot(SurvivalProbs,xlab="Time", ylab="Probability of Survival")
}
}
}
\keyword{multivariate}

