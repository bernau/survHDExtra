\name{gsaWilcoxSurv}
\alias{gsaWilcoxSurv}
\alias{gsaWilcoxSurv,gsagenesets,missing,missing,character,numeric-method}
\alias{gsaWilcoxSurv,gsagenesets,matrix,Surv,missing,missing-method}
\alias{gsaWilcoxSurv,gsagenesets,ExpressionSet,Surv,missing,missing-method}
\alias{gsaWilcoxSurv,gsagenesets,data.frame,Surv,missing,missing-method}
\alias{gsaWilcoxSurv,list,missing,missing,character,numeric-method}
\alias{gsaWilcoxSurv,list,ANY,Surv,missing,missing-method}
\title{Gene Set Analysis GSA, experimental implementation}
\description{ 
    Test whether a set of genes is highly ranked relative to other
    genes in terms of a given statistic, for example a ranking based on univariate
    Cox coefficients. Genes are assumed to be independent.
}
\usage{
\S4method{gsaWilcoxSurv}{gsagenesets,missing,missing,character,numeric}(
gmt, genenames, statistics,  p.value=0.1, cluster=FALSE,
cluster.threshold=0.3,...)
\S4method{gsaWilcoxSurv}{gsagenesets,ExpressionSet,Surv,missing,missing}(gmt,
X, y, genenames, statistics, ...)
}
\arguments{
    \item{gmt}{An object of class \code{\link{gsagenesets}} or a list of gene
    sets.}
    \item{X}{An object of class \code{\link{ExpressionSet}} or matrix.}
    \item{y}{An object of class \code{\link{Surv}}.}
    \item{genenames}{Alternatively to \code{X} and \code{y}, a
    \code{character} vector of gene names 
    rcorresponding to genes in \code{gmt}.}
    \item{statistics}{\code{numeric} vector for genes specified in
    \code{genenames}, any genewise statistic by which genes can be
              ranked.} 
    \item{p.value}{Only report gene sets with p-value lower than this cutoff.}
    \item{cluster}{If true, then significant gene sets are clustered and ranked
    together.}
    \item{cluster.threshold}{Minimum overlap of gene set clustering.}
    \item{...}{Additional arguments passed from alternative S4 signatures or
    to the wilcoxGST function.}
}
\value{An object of class \code{\link{gsaresults}}.}
\details{
    A convenient wrapper around the \code{wilcoxGST} function from the limma
    package. See the limma package documentation for details. 

    This method makes it a little bit easier to test gene sets from a GMT file
    with the \code{wilcoxGST} function. It further adds a clustering of
    gene sets that passed the \code{p.value} threshold. The cluster parameter
    \code{cluster.threshold} specifies the minimum overlap of genes (default
    is 0.3 or 30 percent) in the clustering.
}
\examples{
    library(survHD)
     library(survHDExtra)
    set.seed(100)
    # create some random data
    x<-matrix(rnorm(1000*20),ncol=20)
    dd<-sample(1:1000,size=100)
    u<-matrix(2*rnorm(100),ncol=10,nrow=100)
    x[dd,11:20]<-x[dd,11:20]+u
    y<-Surv(c(rnorm(10)+1,rnorm(10)+2), rep(TRUE, 20))
    genenames=paste("g",1:1000,sep="")
    rownames(x) = genenames

    # create some random gene sets
    genesets=vector("list",50)
    for(i in 1:50){
        genesets[[i]]=paste("g",sample(1:1000,size=30),sep="")
    }
    geneset.names=paste("set",as.character(1:50),sep="")
    gmt <- new("gsagenesets", genesets=genesets, geneset.names=geneset.names)
    gsa.res <- gsaWilcoxSurv(gmt, X=x,y=y,cluster=FALSE,p.value=0.3 )

    # show the similarity of significant gene sets
    plot(gsa.res)

    # display a barcode of up to two gene sets.
    plot(gsa.res, type="barcode",geneset.id1="set22", geneset.id2="set33")
    
    library(genefilter)
    # use a pre-ranking of genes
    genes.ttest = rowttests(x, as.factor(c(rep(1,10),rep(2,10))))
    gsa.res.tt  <- gsaWilcoxSurv(gmt, genenames=rownames(x),
        statistics=genes.ttest[,1])

    # now test some gene signatures on our Affymetrix example data (Beer et. al
    # 2002)
    data(beer.exprs)
    data(beer.survival)
    library(hu6800.db)
    library(annotate)

    gmt <- gsaReadGmt(system.file("extdata/ovarian_gene_signatures.gmt", package = "survHD"))
    # the Gmt file contains gene symbols, so we translate it to Affymetrix 
    genes <- getSYMBOL(rownames(beer.exprs), "hu6800")
    gmt.affy <- gsaTranslateGmt(gmt, beer.exprs, genes)
    gsa.res <- gsaWilcoxSurv(gmt.affy, beer.exprs, Surv(beer.survival[,2], beer.survival[,1]))
}
\author{
    Markus Riester \email{markus@jimmy.harvard.edu}, Levi Waldron
    \email{lwaldron@hsph.harvard.edu}, Christoph Bernau
    \email{bernau@ibe.med.uni-muenchen.de} 
}
\seealso{\code{\link{gsaReadGmt}}, \code{\link{gsaTranslateGmt}} }
