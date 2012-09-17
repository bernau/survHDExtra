\name{customRSF}
\alias{customRSF}
\title{Random Survival Forests using the \code{randomSurvivalForest} R package}
\description{Constructs Random Survival Forests using the function \code{rsf}.
Usually, this function is not called directly but embedded in a call to \code{customSurv} from the main package \code{survHD} (see examples).
}
\usage{
customRSF(Xlearn, Ylearn, learnind, ...)
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
 
  \item{\dots}{Further arguments that shall be passed to function \code{rsf}
  from package{randomSurvivalForest}.}
    
}

\value{An object of class \code{\link{ModelLearned}}.}

\references{
Ishwaran, H., et al. (2008) Random Survival Forests, The Annals of Applied Statistics, 2, 841-860

}

\author{
 
  Christoph Bernau \email{bernau@ibe.med.uni-muenchen.de}

}

\seealso{\code{\link{learnSurvival}}}

\examples{
##load example data
if(require(pensim)){
require(survHD)
require(survHDExtra)
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
fit.rsf <- customSurv(X=beerX,y=beerY, learnind=learnind, ntree=50,customSurvModel=customRSF)
##linear risk score
##survival probabilities at each uncensored time
mytimegrid <- sort(c(0.0,beerY[beerY[,2]==1,1]))
SurvivalProbs <- predict(fit.rsf,beerX[testind,],type='SurvivalProbs',timegrid=mytimegrid,gbm=FALSE)@SurvivalProbs
plot(SurvivalProbs,xlab="Time", ylab="Probability of Survival")
}
}

\keyword{multivariate}