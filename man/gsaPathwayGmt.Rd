\name{gsaPathwayGmt}
\alias{gsaPathwayGmt}
\title{
    Retrieve pathway genesets.
}
\description{
    Provides the Biocarta, Reactome, KEGG and NCI pathway databases as geneset
    objects. Utilizes the \code{graphite} package.
}
\usage{
    gsaPathwayGmt(pathway,type="symbol")
}
\arguments{
    \item{pathway}{
        Either a string with the pathway name (reactome, kegg, nci, biocarta)
        or a list of pathway objects as defined in the \code{graphite}
        package.
    }
    \item{type}{
        Genesets can bei either official gene \code{symbols} (default) or
        \code{entrez} ids.
    }
}
\value{
    An object of class \code{\link{gsagenesets}}.
}
\details{
    Provides quick access to the main pathway databases.
}
\author{
    \email{markus@jimmy.harvard.edu}, Levi Waldron
    \email{lwaldron@hsph.harvard.edu}, Christoph Bernau
    \email{bernau@ibe.med.uni-muenchen.de}
}
\seealso{
    \code{\link{gsaWilcoxSurv}}, \code{\link{gsaReadGmt}}
}
\examples{
    require(survHD)
    gmt.reactome <- gsaPathwayGmt("reactome")
}
\keyword{GSA}
