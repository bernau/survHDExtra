\name{GSA plots}
\alias{plot,gsaresults-method}
\alias{plot,gsaresults,missing-method}
\alias{plot,gsagenesets,missing-method}
\title{GSA plots}
\description{
    Plots for the visualization of GSA results.
}
\usage{
\S4method{plot}{gsaresults,missing}(x, type="overlap", geneset.id1,
geneset.id2, show.es=TRUE, title=NULL, ...)
\S4method{plot}{gsagenesets,missing}(x, type="overlap", ... )
}
\arguments{
    \item{x}{An object of class \code{\link{gsaresults}} or
    \code{\link{gsagenesets}}.}
    \item{type}{The type of plot. Can be "overlap" for an overlap plot, and
        "barcode" for an GSEA-like barcode plot (the latter is only available
        for \code{\link{gsaresults}}). Defaults to "overlap", which shows
    the similarity of (significant) gene sets.}
    \item{geneset.id1}{For "barcode", the id of the geneset in \code{x} to be visualized.}
    \item{geneset.id2}{For "barcode", optional second id, useful when \code{geneset.id1} are up-regulated
        genes and \code{geneset.id2} down-regulated genes.}
    \item{show.es}{For "barcode", show the running Kolmogorov-Smirnov statistic.}
    \item{title}{For "barcode", print this title.} 
    \item{\dots}{Further graphical options passed to \code{levelplot} for
        "overlap" and \code{barcodeplot} for "barcode".}
}
\value{
    No return.
}
\details{
    Provides two kinds of plots for GSA results. The default is a heatmap that
    visualizes the overlap of significant gene sets. The second "barcode" plot is
    very similar to the GSEA (Subramanian 2005) plots.  The first panel of the plot
    shows an enrichment score, which is the running deviation from the expected
    uniformly distributed ranking of genes. More precisely, if m is the gene set
    size and n the total number of genes in the ranking, the enrichment score at
    rank i is increased by 1-m/n if gene i is in the gene set, otherwise it is
    decreased by m/n. The enrichment score is scaled, so that +1 and -1 correspond
    to 95 percent confidence intervals. The second plot shows the barcode plot from
    the limma package, in which the bars visualize the position of the genes.  It
    is possible to visualize two gene sets in one plot, which makes sense for gene
    sets that are divided in up- and down-regulated genes. This makes it possible
    to not only check for non-randomness, but also for consistency of the
    expression direction.
}
\references{
    Subramanian A, Tamayo P, Mootha VK, Mukherjee S, Ebert BL,
    Gillette MA, Paulovich A, Pomeroy SL, Golub TR, Lander ES, Mesirov
    JP (2005). Gene set enrichment analysis: a knowledge-based
    approach for interpreting genome-wide expression profiles. Proc
    Natl Acad Sci USA, 102, 15545-15550.
}    
\author{
    Markus Riester \email{markus@jimmy.harvard.edu}, Levi Waldron
    \email{lwaldron@hsph.harvard.edu}, Christoph Bernau
    \email{bernau@ibe.med.uni-muenchen.de}
}       
\seealso{\code{\link{gsaresults}}, \code{\link{gsagenesets}}}
\keyword{GSA}
