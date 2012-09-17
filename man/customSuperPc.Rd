\name{customSuperPc}
\alias{customSuperPc}

\title{Univarate shrinkage prediction in the Cox model}
\description{See
  http://cran.r-project.org/web/packages/superpc/index.html for more
  details

Usually, this function is not called directly but embedded in a call to \code{customSurv} from the main package \code{survHD} (see examples).

  Note that "lambda", if specified, corresponds to the "threshold"
  argument in superpc.  
}
\usage{
customSuperPc(Xlearn, Ylearn, learnind,...)
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
        
  \item{\dots}{Further arguments that shall be passed to \code{superpc.train}
    function from the superpc R package (http://cran.r-project.org/web/packages/superpc/index.html)
  }
  }

\value{An object of class \code{\link{ModelLearned}}.}

\references{Bair E, Tibshirani R: Semi-Supervised Methods to Predict Patient Survival from Gene Expression Data. \emph{PLoS Biology 2004, 2:e108.}}

\author{
  Levi Waldron \email{lwaldron@hsph.harvard.edu}
  Christoph Bernau \email{bernau@ibe.med.uni-muenchen.de}
  Markus Riester <riester@jimmy.harvard.edu>
}

\seealso{\code{\link{learnSurvival}}}

\examples{
require(survHD)
require(survival)
require(survHDExtra)	
	set.seed(20)
	nsamples <- 100
	X <- matrix(rnorm(nsamples*1000),nrow=nsamples)
	colnames(X) <- make.names(1:ncol(X))
	rownames(X) <- make.names(1:nrow(X))
	time <- rexp(nsamples)
	cens <- sample(0:1,size=nsamples,replace=TRUE)
	y <- Surv(time,cens)
ModelLearned<-customSurv(X=X,y=y,learnind=1:75,customSurvModel=customSuperPc)	
##linear risk score
predict(ModelLearned@model,newdata=X[-(1:75),],type='lp')
##survival probabilities at each uncensored time
mytimegrid <- sort(c(0.0,y[y[,2]==1,1]))
SurvivalProbs <- predict(ModelLearned,X[-(1:75),],type='SurvivalProbs',timegrid=mytimegrid,gbm=TRUE)@SurvivalProbs
plot(SurvivalProbs,xlab="Time", ylab="Probability of Survival")

}
\keyword{multivariate}

