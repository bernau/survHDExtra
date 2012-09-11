\name{Kaplan-Meier plot (survaggr)}
\alias{plot,survaggr,missing-method}
\title{Kaplan-Meier plot for survaggr-objects}
\description{The plot visualizes a Kaplan-Meier plot from a survaggr
object}
\arguments{
  \item{x}{An object of class \code{\link{survaggr}}}
  \item{method}{A method for patient stratification supported by
  \code{\link{plotKMStratifyBy}}. Default "median". Ignored if
  \code{\link{plotKMStratifyBy}} argument \code{cutpoints} is not
  NULL.}
  \item{labels}{
  The labels of the groups specified by the cutpoints. Default is "High Risk"
  and "Low Risk" for one cutpoint, and "High Risk", "Mid Risk", and "Low Risk"
  for two cutpoints (when \code{method="tertiles"}).}
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

svaggr<-learnSurvival(X=X,y=y,survmethod='coxBoostSurv',learningsets=ls)
plot(svaggr)
}
}
\author{Markus Riester  \email{markus@jimmy.harvard.edu}}
\seealso{\code{\link{plotKM}}, \code{\link{plotKMStratifyBy}}}
