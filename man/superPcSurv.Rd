\name{superPcSurv}
\alias{superPcSurv}
\alias{superPcSurv,ExpressionSet,character-method}
\alias{superPcSurv,ExpressionSet,Surv-method}
\alias{superPcSurv,data.frame,Surv-method}
\alias{superPcSurv,matrix,Surv-method}
\title{Univarate shrinkage prediction in the Cox model}
\description{See
  http://cran.r-project.org/web/packages/superpc/index.html for more
  details

  Note that "lambda", if specified, corresponds to the "threshold"
  argument in superpc.  
}
\usage{
\S4method{superPcSurv}{data.frame,Surv}(X, y, learnind, essential=FALSE ,...)
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

ModelLearned<-superPcSurv(X=X,y=y,learnind=1:75)	

}
}
\keyword{multivariate}

