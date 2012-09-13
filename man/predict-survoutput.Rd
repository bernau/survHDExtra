\name{Prediction from ModelLearned objects}
\docType{methods}
\alias{predict,ModelLearned-method}
\alias{predictsurvhd,ModelLearned,data.frame-method}

\title{Predictions from a survival model}
\description{General predictor function for survival models in survHD}
\usage{
\S4method{predictsurvhd}{ModelLearned,data.frame}(object, newdata,type,timegrid=NULL,gbm=TRUE,...)
\S4method{predict}{ModelLearned}(object,...)
}
\arguments{

  \item{object}{object of class \code{ModelLearned}}
\item{newdata}{Gene expression data for which predictions shall be returned. The following formats are available:
  \itemize{
   \item{\code{matrix}} Rows correspond to observations, columns to variables.
   \item{\code{data.frame}}  Rows correspond to observations, columns to variables.
   \item{\code{ExpressionSet}} survHD will extract the expressions using exprs().
  } }
  \item{type}{class \code{character}, can be one of the following:
  \itemize{
  \item "lp", if the linear predictor is to be computed
  \item "SurvivalProbs", if the survival probabilities are to be computed
  }}
  \item{timegrid}{the timepoints for which the survival probabilities are to be computed}
  \item{gbm}{shall the default estimation procedure from package \code{gbm} be used for estimating the basehazard?}
  \item{\dots}{currently unused argument}
}

\seealso{\code{\link{gbmbasehaz}}}
\keyword{multivariate}