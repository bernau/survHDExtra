\name{Kaplan-Meier plot (ModelLearned)}
\alias{plot,ModelLearned-method}
\alias{plot,ModelLearned,missing-method}
\alias{plot,ModelBase,missing-method}
\title{Kaplan-Meier plot}
\description{ The plot visualizes a Kaplan-Meier plot from a ModelLearned
or ModelLinear object}

\arguments{
  \item{x}{An object of class \code{\link{ModelLearned}} or
  \code{\link{ModelLinear}}.}
  \item{newdata}{gene expressions for the test data}
  \item{newy}{survival response for the test data}
  \item{method}{A method for patient stratification supported by
  \code{\link{plotKMStratifyBy}}. Default "median". Ignored if cutpoints is not
  NULL.}
  \item{cutpoints}{Specifies cutpoints. Default NULL, which means that the
  specified method is used to calculate cutpoints.}
  \item{labels}{
  The labels of the groups specified by the cutpoints. Default is "High Risk"
  and "Low Risk" for one cutpoint, and "High Risk", "Mid Risk", and "Low Risk"
  for two cutpoints (when method="tertiles").}
  \item{\dots}{Further graphical options passed to \code{\link{plotKM}}.}
}
\value{Return value of \code{\link{plotKMStratifyBy}}.}

\examples{
\dontrun{
library(survHD)
	data(beer.exprs)
	set.seed(123)
	X<-t(as.matrix(beer.exprs))
	data(beer.survival)
	y<-Surv(beer.survival[,2],beer.survival[,1])
	
	ls<-generateLearningsets(y=y[,1],method='CV',fold=3,niter=1)
	
	svaggr<-learnSurvival(X=X,y=y,survmethod='coxBoostSurv',LearningSets=ls)
	plot(svaggr@ModelLearnedlist[[2]],newdata=svaggr@X[svaggr@ModelLearnedlist[[2]]@learnind,],newy=svaggr@ModelLearnedlist[[2]]@y)

}
}
\author{Markus Riester  \email{markus@jimmy.harvard.edu}}
\seealso{\code{\link{plotKM}}, \code{\link{plotKMStratifyBy}}}
