\name{customUniCox}
\alias{customUniCox}
\title{Univarate shrinkage prediction in the Cox model}
\description{
Usually, this function is not called directly but embedded in a call to \code{customSurv} from the main package \code{survHD} (see examples).

See
  http://cran.r-project.org/web/packages/uniCox/index.html for more details
}
\usage{
customUniCox(Xlearn, Ylearn, learnind,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{Xlearn}{Gene expression data (a \code{data.frame})}
  
  \item{Ylearn}{Survival Response, an object of class \code{Surv}.}
 
  \item{learnind}{An index vector specifying the observations that
                  belong to the learning set. May be \code{missing};
                  in that case, the learning set consists of all
                  observations and predictions are made on the
                  learning set.}
                  
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
##load example data
if(require(pensim)){
library(survHD)
library(survHDExtra)
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
fit.unicox <- customSurv(X=beerX,y=beerY, learnind=learnind, lambda=10,customSurvModel=customUniCox)
##linear risk score
predict(fit.unicox@model,newdata=beerX[testind,],type='lp')
##survival probabilities at each uncensored time
mytimegrid <- sort(c(0.0,beerY[beerY[,2]==1,1]))
SurvivalProbs <- predict(fit.unicox,beerX[testind,],type='SurvivalProbs',timegrid=mytimegrid,gbm=TRUE)@SurvivalProbs
plot(SurvivalProbs,xlab="Time", ylab="Probability of Survival")

}
}
\keyword{multivariate}

