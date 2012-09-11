\name{Barplot}
\alias{plot,genesel-method}
\alias{plot,genesel,missing-method}
\title{Barplot of variable importance}
\description{ The plot visualizes
             variable importance by a barplot. The height
             of the barplots correspond to variable importance.
             What variable importance exactly means depends
             on the method chosen when calling \code{\link{geneSelection}},
             s. \code{\link{genesel}}.}

\arguments{
  \item{x}{An object of class \code{\link{genesel}}}
  \item{top}{Number of top genes whose variable importance should be displayed.
             Defaults to 10.}
  \item{iter}{Iteration number (\code{learningset}) for which variable
              importance should be displayed.}
  \item{\dots}{Further graphical options passed to \code{barplot}.}
}
\note{Note the following
      \itemize{
      \item If \code{scheme = "multiclass"}, only one plot will be made.
            Otherwise, one plot will be made for each binary scenario
            (depending on whether \code{"scheme"} is \code{"one-vs-all"}
            or \code{"pairwise"}).
      \item Variable importance do not make sense for variable selection
            (ranking) methods that are essentially discrete, such as
            the Wilcoxon-Rank sum statistic or the Kruskal-Wallis statistic.
      \item For the methods \code{"lasso", "elasticnet", "boosting"}
            the number of nonzero coefficients can be very small, resulting
            in bars of height zero if \code{top} has been chosen too
	    large.
	    }
      }
\value{No return.}

\author{
Christoph Bernau \email{bernau@ibe.med.uni-muenchen.de}

Martin Slawski \email{ms@cs.uni-sb.de}

        Anne-Laure Boulesteix \email{boulesteix@ibe.med.uni-muenchen.de}}
\seealso{\code{\link{genesel}}, \code{\link{geneSelection}}}
\keyword{multivariate}
