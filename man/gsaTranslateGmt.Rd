\name{gsaTranslateGmt}
\alias{gsaTranslateGmt}
\alias{gsaTranslateGmt,gsagenesets,matrix,character-method}
\alias{gsaTranslateGmt,gsagenesets,ExpressionSet,character-method}
\alias{gsaTranslateGmt,gsagenesets,data.frame,character-method}
\title{GSA Translate GMT }
\description{ 
    Converts a GMT file so that symbols match rownames of X.
}
\usage{
\S4method{gsaTranslateGmt}{gsagenesets,matrix,character}(gmt, X, genes)
\S4method{gsaTranslateGmt}{gsagenesets,ExpressionSet,character}(gmt, X, genes)
\S4method{gsaTranslateGmt}{gsagenesets,data.frame,character}(gmt, X, genes)
}
\arguments{
    \item{gmt}{An object of class \code{\link{gsagenesets}} or a list of gene
    sets.}
    \item{X}{An object of class \code{\link{ExpressionSet}}, matrix or
    data.frame.}
    \item{genes}{An \code{character} vector with the gene names of X.}
}
\value{An object of class \code{\link{gsagenesets}}.}
\details{
    Can be used to translate gene sets, for example from official gene symbols
    to manufacturer probe ids.
}
\examples{
    library(survHD)
    library(survHDExtra)
    library(hu6800.db)
    library(annotate)

    data(beer.exprs)
    data(beer.survival)
    gmt <- gsaReadGmt(system.file("extdata/ovarian_gene_signatures.gmt", package = "survHD"))
    # the Gmt file contains gene symbols, so we translate it to Affymetrix 
    genes <- getSYMBOL(rownames(beer.exprs), "hu6800")
    gmt.affy <- gsaTranslateGmt(gmt, beer.exprs, genes)
}
\author{
        Markus Riester \email{markus@jimmy.harvard.edu}, Levi Waldron
        \email{lwaldron@hsph.harvard.edu}, Christoph Bernau
        \email{bernau@ibe.med.uni-muenchen.de} 
}
\seealso{\code{\link{gsaReadGmt}}, \code{\link{gsaTranslateGmt}} }
