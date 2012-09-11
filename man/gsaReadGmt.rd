\name{gsaReadGmt}
\alias{gsaReadGmt}
\title{
    Read in a gene set collection from a .gmt file
}
\description{
    Read in a gene set collection from a .gmt file. Code was adapted from the
    GSA package by Tibshirani and Efron.
}
\usage{
    gsaReadGmt(filename, direction)
}
\arguments{
    \item{filename}{
        The name of a file to read data values from. Should be a tab-separated text
        file, with one row per gene set. Column 1 has gene set names (identifiers),
        column 2 has gene set descriptions, remaining columns are gene ids for genes
        in that geneset.
    }
    \item{direction}{
        Optionally, if the gene set consists of only up- or down-regulated genes, this
        direction ("up" or "down").
    }
}
\value{
    An object of class \code{\link{gsagenesets}}.
}
\details{
    Reads gene sets in the GMT format, for example from MSigDB.
}
\author{
    Robert Tibshirani (adapted from the GSA package), Markus Riester
    \email{markus@jimmy.harvard.edu}, Levi Waldron
    \email{lwaldron@hsph.harvard.edu}, Christoph Bernau
    \email{bernau@ibe.med.uni-muenchen.de}
}
\seealso{
    \code{\link{gsaWilcoxSurv}}
}
\examples{
    require(survHD)
    gmt <- gsaReadGmt(system.file("extdata/ovarian_gene_signatures.gmt",
    package = "survHD"))
}
\keyword{GSA}
