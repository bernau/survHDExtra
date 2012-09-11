\name{plotKMStratifyBy}
\alias{plotKMStratifyBy}
\title{
Kaplan-Meier plot
}
\description{
 A simple but feature-complete Kaplan-Meier plot with sane defaults.
}
\usage{
plotKMStratifyBy(method = "median", y, linearriskscore, cutpoints=NULL, labels
= NULL, plot = TRUE, ...)
}
\arguments{
  \item{method}{
  A method for stratifying patients based on their linear risk score. Default
  is stratify based on the cohort median. Valid values are "median" and
  "tertile" as well as the methods provided by the maxstat package for
  calculating optimal cutpoints. Ignored if cutpoints is not NULL.
}
  \item{y}{
  Surv object with the survival times for each patient.
}
  \item{linearriskscore}{
  Vector of linear risk scores for each patient.
}
  \item{cutpoints}{
  Vector of cutpoints. Used to stratify patients based on their linear risk
  scroes. Use "method" to calculate cutpoints.
}
  \item{labels}{
  The labels of the groups specified by the cutpoints. Default is "High Risk"
  and "Low Risk" for one cutpoint, and "High Risk", "Mid Risk", and "Low Risk"
  for two cutpoints.
}
  \item{plot}{
  Plot the Kaplan-Meier. If FALSE, only cutpoints are returned.
}
  \item{\dots}{
  Additional arguments passed to plotKM.
}
}
\details{
This is a convenient wrapper around plotKM for the common task of stratifying
patients based on a linear risk score.
}
\value{
\item{cutpoints }{The optimal cutpoint of the linear risk score.}
\item{p.value }{The probability of observing this maximal log-rank statistic.}
\item{hr }{The hazard ratio of this cutpoint} 
}
\author{Markus Riester \email{markus@jimmy.harvard.edu}}
\seealso{\code{\link{plotKM}}}
\examples{
\dontrun{
require(survHD)
# stratify patients by median of a linear risk score
plotKMStratifyBy("median",y=Surv(bladder$stop, bladder$event),linearriskscore=1/bladder$enum)
# find optimal cutoffs with the maxstat package, using the Lau92 method
plotKMStratifyBy("Lau92",y=Surv(bladder$stop, bladder$event),linearriskscore=1/bladder$enum)
}
}

