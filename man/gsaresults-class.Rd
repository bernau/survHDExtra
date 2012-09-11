\name{gsaresults-class}
\docType{class}
\alias{gsaresults-class}
\alias{gsaresults}

\title{"gsaresults"}
\description{Object returned by \code{\link{gsaWilcoxSurv}}}
\section{Slots}{
	 \describe{
    \item{\code{geneset.names}}{\code{character} vector of gene sets names}
    \item{\code{p.values}}{The p-values of the gene sets, as calculated
    \code{wilcoxGST} from the limma package.}
    \item{\code{statistics}}{numeric vector giving the values of statistics
    that were used to rank genes in the GSA.}
    \item{\code{genesets.used}}{The actual gene sets used in the GSA. This
    includes only symbols found in \code{X}.}
    \item{\code{geneset.cluster}}{If the \code{cluster} option was set, this
    includes a factor assigning each gene set to a cluster. This is useful for
    combining similar gene sets.}
  }
}


\author{
        Markus Riester \email{markus@jimmy.harvard.edu}, Levi Waldron \email{lwaldron@hsph.harvard.edu},Christoph Bernau \email{bernau@ibe.med.uni-muenchen.de}
}

\seealso{
        \code{\link{gsaWilcoxSurv}} }
\keyword{GSA}
