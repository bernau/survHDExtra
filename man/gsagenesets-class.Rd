\name{gsagenesets-class}
\docType{class}
\alias{gsagenesets-class}
\alias{gsagenesets}

\title{"gsagenesets"}
\description{Container for gene sets, as used by \code{\link{gsaWilcoxSurv}}}
\section{Slots}{
	 \describe{
    \item{\code{genesets}}{List of gene sets, which are lists of
    \code{character} vectors of gene identifiers.}
    \item{\code{geneset.names}}{\code{character} vector of gene sets names}
    \item{\code{geneset.descriptions}}{An optional vector of \code{character}
    descriptions of the gene
    sets.}
    \item{\code{geneset.direction}}{Can be used to flag gene sets as up- or
    down-regulated genes.} 
  }
}


\author{
        Markus Riester \email{markus@jimmy.harvard.edu}, Levi Waldron \email{lwaldron@hsph.harvard.edu},Christoph Bernau \email{bernau@ibe.med.uni-muenchen.de}
}

\seealso{
        \code{\link{gsaWilcoxSurv}} }
\keyword{GSA}
