\name{coxBoostSurv}
\alias{coxBoostSurv}
\alias{coxBoostSurv,matrix,Surv-method}
\alias{coxBoostSurv,data.frame,Surv-method}
\alias{coxBoostSurv,ExpressionSet,character-method}
\alias{coxBoostSurv,ExpressionSet,Surv-method}
\title{Componentwise Boosting}
\description{Roughly speaking, Boosting combines 'weak learners'
            in a weighted manner in a stronger ensemble.
}
\usage{
\S4method{coxBoostSurv}{data.frame,Surv}(X, y, learnind, unpenalizedvars=NULL, essential=FALSE,...)
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
  \item{unpenalizedvars}{matrix of unpenalized variables (same structure as \code{X})}
  \item{essential}{shall the model output (see slot model of class \code{\link{ModelLearned}})
                  include the bare essentials only (especially for tune this option might save a lot of memory)?}
  \item{\dots}{Further arguments that shall be passed to function \code{CoxBoost} (package: \code{CoxBoost}).}
}
\value{An object of class \code{\link{ModelLearned}}.}
\references{   Binder, H., et al. (2009) Boosting for high-dimensional time-to-event data with
competing risks, Bioinformatics, 25 (7), 890-896
}
\author{Christoph Bernau \email{bernau@ibe.med.uni-muenchen.de}}
\seealso{\code{\link{learnSurvival}}} 
\examples{
\dontrun{
require(survHD)
require(survival)
	
set.seed(20)
nsamples <- 100
X <- matrix(rnorm(nsamples*1000),nrow=nsamples)
colnames(X) <- make.names(1:ncol(X))
rownames(X) <- make.names(1:nrow(X))
time <- rexp(nsamples)
cens <- sample(0:1,size=nsamples,replace=TRUE)
y <- Surv(time,cens)

ModelLearned<-coxBoostSurv(X=X,y=y,learnind=1:75)
}
}
\keyword{multivariate}
