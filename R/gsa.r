# Filename: gsa.r
# Title: Various GSA utils
# 
# Author: Markus Riester
# Email: <markus@jimmy.harvard.edu>
# Date of creation: 08/07/2012
#
# Brief description:
#   Several function that make the use of the limma GSA functions a little
#   bit easier.
#
###############################################################################




gsaReadGmt <- function (filename, direction=NULL) 
{
    a = scan(filename, what = list("", ""), sep = "\t", quote = NULL, 
        fill = TRUE, flush = TRUE, multi.line = FALSE)
    geneset.names = a[1][[1]]
    geneset.descriptions = a[2][[1]]
    dd = scan(filename, what = "", sep = "\t", quote = NULL)
    nn = length(geneset.names)
    n = length(dd)
    ox = rep(NA, nn)
    ii = 1
    for (i in 1:nn) {
        while ((dd[ii] != geneset.names[i]) | (dd[ii + 1] != 
            geneset.descriptions[i])) {
            ii = ii + 1
        }
        ox[i] = ii
        ii = ii + 1
    }
    genesets = vector("list", nn)
    for (i in 1:(nn - 1)) {
        i1 = ox[i] + 2
        i2 = ox[i + 1] - 1
        geneset.descriptions[i] = dd[ox[i] + 1]
        genesets[[i]] = dd[i1:i2]
    }
    geneset.descriptions[nn] = dd[ox[nn] + 1]
    genesets[[nn]] = dd[(ox[nn] + 2):n]
    out = new("gsagenesets",genesets = genesets, geneset.names = geneset.names, 
        geneset.descriptions = geneset.descriptions)
    if (!is.null(direction)) {
        if (!length(grep("^(up|down)$",direction))) 
            stop("Unknown direction, must be 'up' or 'down'")
        out@geneset.direction <- direction
    }
    return(out)
}

gsaClusterGenesets <- function(gmt, cluster.threshold) {
    universe = unique(unlist(gmt@genesets))
    imat = t(mapply(cbind,lapply(gmt@genesets, function(x) sapply(universe, function(y) ifelse(y %in% x,1,0)))))
    colnames(imat) = universe
    Amx = imat %*% t(imat)
    minS=outer(diag(Amx),diag(Amx), pmin)
    x=Amx/minS
    dd.row = as.dendrogram(hclust(dist(x)))
    row.ord = order.dendrogram(dd.row)
    label <- rep(0, nrow(x))
    idx = sapply(2:nrow(x), function(i) { label[i] <<- label[i-1];
    if(x[row.ord[i-1], row.ord[i]]< cluster.threshold) label[i] <<- label[i-1]+1 })
    label <- label[order(row.ord)]
    as.factor(label)
}

gsaPathwayGmt <- function(pathway, type="symbol") {
    require(graphite)
    if (class(pathway) == "character") {
        pathway <- switch(pathway, 
            "biocarta" = biocarta,
            "reactome" = reactome,
            "kegg"     = kegg,
            "nci"      = nci,
            stop("Unknown pathway (not biocarta, reactome, kegg, nci)")
        )
    } else {
        if (class(pathway) != "list" || class(pathway[[1]]) != "pathway")
            stop("Unknown pathway (not biocarta, reactome, kegg, nci)")
    }
    tmp <- lapply(pathway,function(x) convertIdentifiers(x,type)@nodes)
    new("gsagenesets", genesets = tmp, geneset.names=names(tmp))
}

setMethod("gsaTranslateGmt", signature(gmt="gsagenesets", X="ExpressionSet", genes="character"), function(gmt, 
X, genes) {
 gsaTranslateGmt(gmt,exprs(X), genes)  
})

setMethod("gsaTranslateGmt", signature(gmt="gsagenesets", X="data.frame", genes="character"), function(gmt, 
X, genes) {
 gsaTranslateGmt(gmt,as.matrix(X), genes)  
})

setMethod("gsaTranslateGmt", signature(gmt="gsagenesets", X="matrix", genes="character"), function(gmt, 
X, genes) {
    if (length(rownames(X)) != length(genes)) 
        stop("length of rownames of X and genes not equal")
    if (sum(duplicated(genes)) > 0)
        warning("genes not unique, will pick first matching probeset.")
    gmt@genesets = lapply(gmt@genesets, function(x)
    rownames(X)[na.omit(match(x, genes))])
    gmt    
})


setMethod("gsaWilcoxSurv", signature(gmt="gsagenesets", X="missing",
y="missing",genenames="character",statistics="numeric"), function(gmt, 
genenames, statistics,  p.value=0.1, cluster=FALSE, cluster.threshold=0.3, ...) {
    require(limma)
    
    # now use the wilcoxGST function to do the GSA
    res <- lapply(gmt@genesets, function(gs) {
        genes <- gs[gs %in% genenames]
        idx <- match( genes, genenames )
        list(p=wilcoxGST(idx,statistics,...), genes=genes)
    })
    
    # get the names of the genesets 
    if (is.null(names(res))) names(res) <- 1:length(res)
    if (length(gmt@geneset.names)) names(res) <- gmt@geneset.names
    
    # filter results by p-value
    res <- res[sapply(res, function(xx) xx$p < p.value)]
    
    if (!length(res)) return(new("gsaresults"))
    df <- 
    data.frame(geneset.names=names(res), p.values=sapply(res, function(x) x$p),
        stringsAsFactors=FALSE)

    if (cluster) {
        # Combine similar gene sets
        used.genesets <- new("gsagenesets", genesets=lapply(res, function(xx) xx$genes))
        cat("Clustering", length(used.genesets@genesets), "gene sets with p < ",
        p.value)
        df$geneset.cluster <- gsaClusterGenesets(used.genesets, cluster.threshold)
        df.o <- df[order(df$p.values),]
        cluster.order <- unique(df.o$geneset.cluster)
        df.o <- df.o[order(match(df.o$geneset.cluster,cluster.order)),]
    } else {
        df.o <- df[order(df$p.values),]
    }

    res <- new("gsaresults",
        geneset.names   = df.o$geneset.names,
        p.values        = df.o$p.values,
        genesets.used   = lapply(res[match(df.o$geneset.names,names(res))],
        function(xx) xx$genes),
        statistics      = statistics)
    names(res@statistics) = genenames    
    if (!is.null(df.o$geneset.cluster)) res@geneset.cluster <-
        df.o$geneset.cluster
    res
})

setMethod("gsaWilcoxSurv",
signature(gmt="list",X="ANY",y="Surv",genenames="missing",
statistics="missing"), function(gmt,X,y,...) {
    gsaWilcoxSurv(new("gsagenesets", genesets = gmt), X,y,...)
})

setMethod("gsaWilcoxSurv",
signature(gmt="list",X="missing",y="missing",genenames="character",statistics="numeric"),
function(gmt,genenames,statistics,...) {
    gsaWilcoxSurv(new("gsagenesets", genesets = gmt), genenames,statistics,...)
})

setMethod("gsaWilcoxSurv", signature(gmt="gsagenesets", X="matrix", y="Surv",
genenames="missing", statistics="missing"), function(gmt, X, y,
... ) {

        if (is.null(rownames(X))) 
            stop("X must have valid rownames.")
   
        idx <- complete.cases(y)
        if ( sum(!idx) > 0) {
            warning("Samples with unknown y removed.")
            X <- X[,idx]
            y <- y[idx]
        }
        genes.rowtest <- rowCoxTests(X,y)
        gsaWilcoxSurv(gmt,
        genenames=rownames(genes.rowtest),statistics=(1-genes.rowtest[,3])*sign(genes.rowtest[,1]),...)
    })

setMethod("gsaWilcoxSurv", signature(gmt="gsagenesets",
X="ExpressionSet",y="Surv",genenames="missing",statistics="missing"), function(gmt, X, y,
... ) {  gsaWilcoxSurv(gmt, X=exprs(X),y=y, ...) })

setMethod("gsaWilcoxSurv", signature(gmt="gsagenesets",
X="data.frame",y="Surv",genenames="missing", statistics="missing"), function(gmt, X, y,
... ) {  gsaWilcoxSurv(gmt, X=as.matrix(X),y=y, ...) })

