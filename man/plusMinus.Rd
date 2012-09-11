\name{plusMinus}
\alias{plusMinus}
\title{Simple average learner}
\description{Creates a linear risk score to predict y from X using one
  of several simple methods:

  plusminus: a very simple method: 1) Fit a univariate Cox model
  for each feature.
  2) Set the beta (coefficients of the multivariate linear risk score)
  for a specified number of features, starting at the top of a list
  ranked by magnitude of the univariate Cox regression coefficient, to
  +1 if the univariate regression coefficient is positive and -1 if it
  is negative.  3) Remaining betas are set to zero.  The number of
  non-zero coefficients is specified as a tuning parameter.
  
  compoundcovariate:  Coefficients are equal to the coefficient
  of the univariate Cox regression

  tscore: features are defined as "good prognosis" or "bad
  prognosis", and the risk score is equal to the t-statistic for a
  pooled variance t-test for good vs. bad prognosis genes.  This method
  was used by The Cancer Genome Atlas (Integrated genomic analyses of
  ovarian carcinoma.  Nature 2011, 474:609-615.)

  voting: a threshold Ti and a coefficient Bi are defined for each
  feature, and the contribution of each feature i is sign(Xi - Ti) * Bi.
  In other words: in the typical case of +-1 coefficients, if a
  bad-prognosis gene has expression above its threshold, its vote is +1,
  and below its threshold, its vote is -1.  Conversely if a
  good-prognosis gene has expression above its threshold, its vote is
  -1, and above its threshold, its vote is +1.  The risk score is the
  sum of the votes of each feature.

  positiveriskvoting: Similar to the voting model, but features can
  only make votes for positive contribution to the risk.

  negativeriskvoting: Similar to the voting model, but features can only
  make votes for negative contribution to the risk, ie positive
  contributions to good prognosis.  An example of this voting scheme in
  use is found in Kang et al., A DNA repair pathway-focused score for
  prediction of outcomes in ovarian cancer treated with platinum-based
  chemotherapy. \emph{Journal of the National Cancer Institute}. 2012
  104(9):670-81.
  
  This function provides the core functionality, plusMinusSurv
  integrates it with the survHD package.
}
\usage{
plusMinus(X, y, lambda=NULL, tuningpar="nfeatures", standardize=FALSE,
  directionality="posneg", ties.method="average", votingthresholdquantile=0.5, modeltype="plusminus")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{Gene expression data. A \code{matrix}. Rows correspond to observations, columns to variables.}           
  \item{y}{Survival Response (an object of class \code{Surv})}
  \item{lambda}{tuning parameter, by default all available variables
  will be included.}
  \item{tuningpar}{if tuningpar="nfeatures", lambda corresponds to the
  number of non-zero coefficients in the model.  If tuningpar="pval",
  lambda corresponds to the p-value from Wald Test of the univariate Cox
  model, above which coefficients will be set to zero.}
  \item{standardize}{Whether to standardize the coefficients by variable
  variance, so that each variable carries an equal weight.  Default is
  FALSE.}
  \item{directionality}{posneg: allow positive and negative
  coefficients.  pos: positive coefficients only.  neg: negative
  coefficients only.}
  \item{ties.method}{How to settle ties when rankingcoefficients, see ?rank.}
  \item{votingthresholdquantile}{Feature thresholds are set at this
  quantile for each feature.  ie, votingthresholdquantile=0.5, each
  feature will make a positive risk vote if its value is above the
  median in the training set and it is a poor-prognosis feature
  (positive coefficient), and make a negative risk vote if its value
  is below the median.  These votes are reversed if it is a
  good-prognosis feature (negative coefficient).}
  \item{modeltype}{Currently supported model types are:  "plusminus",
  "compoundcovariate", "tscore", "voting", "positiveriskvoting", and
  "negativeriskvoting".  Note that positiveriskvoting is not the same as
  setting directionality="pos", for example: with
  modeltype="positiveriskvoting", no feature will ever make a negative
  contribution to risk.  With modeltype="voting" and
  directionality="pos", poor prognosis features only will be selected
  (positive coefficients), but these can make a negative contribution to
  the risk score if their value is below the threshold.}
%  \item{...}{Arguments passed on to \code{\link{rowCoxTests}}.}
}

\value{An object of class \code{\link{coxph}}.}

\author{
  Levi Waldron \email{lwaldron@hsph.harvard.edu}
}

\seealso{\code{\link{learnSurvival}}, \code{\link{plusMinusSurv}}}
\keyword{multivariate}