setMethod("plot", signature(x = "gsaresults", y = "missing"),
function(x, type="overlap", geneset.id1, geneset.id2=NULL, show.es=TRUE,
title=NULL,...) {
   switch(type,
              overlap = .overlapPlot(x@genesets.used,...),
              barcode = .barcodePlot(x, geneset.id1, geneset.id2, show.es,
                title, ...),
              stop("Unknown plot type for object gsaresults"))
})

setMethod("plot", signature(x = "gsagenesets", y = "missing"),
function(x, type="overlap", ...) {
   switch(type,
              overlap = .overlapPlot(x@genesets,...),
              stop("Unknown plot type for object gsaresults"))
})


.characterPlot <- function(char, x, y, deltax, deltay, cex=1){
	spltchar <- unlist(strsplit(char, ""))
	for(s in seq(along=spltchar)) points(x+(s-1)*deltax, y+deltay, pch=spltchar[s], cex=cex)
}

.barcodePlot <- function(x, geneset.id1, geneset.id2=NULL,show.es=TRUE, title=NULL,...) {
	index  <- x@genesets.used[[geneset.id1]]
	index2 <- NULL 
	if (!is.null(geneset.id2)) index2 <- x@genesets.used[[geneset.id2]]
	if (is.null(index)) stop("geneset.id1 not found.")
	require(limma)
	if (show.es) {
		par(mfrow=c(2,1))
		par(mar=c(4, 4, 2.5, 2))
		.ksEnrichmentPlot(x@statistics, index=index, index2=index2, main=title)
		barcodeplot(x@statistics, index=index, index2=index2, ... )
	} else {
		barcodeplot(x@statistics, index=index, index2=index2, ... )
	}
}

.overlapPlot <- function(genesets, xlab="", ylab="" , cex=0.6,...) {
	require(lattice)
	universe = unique(unlist(genesets))
	imat = t(mapply(cbind,lapply(genesets, function(x) sapply(universe, function(y) ifelse(y %in% x,1,0)))))
	colnames(imat) = universe
	Amx = imat %*% t(imat)
	minS=outer(diag(Amx),diag(Amx), pmin)
	x=Amx/minS
	dd.row = as.dendrogram(hclust(dist(x)))
	row.ord = order.dendrogram(dd.row)
	dd.col = as.dendrogram(hclust(dist(t(x))))
	col.ord = order.dendrogram(dd.col)
	levelplot(x[row.ord,col.ord], scales=list(x = list(rot = 90),cex=cex),
			xlab=xlab, ylab=ylab,
			col.regions=colorRampPalette(c("white","darkblue"))(256), ...)
}

.ksEnrichmentPlot <- function(statistics, index, index2=NULL, col.bars=NULL,...) {
	es1 <- .ksEnrichmentScore(statistics, index)
	es2 <- .ksEnrichmentScore(statistics, index2)
	if (is.null(col.bars)) 
		if (is.null(index2)) 
			col.bars = c("black", "black","red")
		else col.bars = c("red", "blue", "black")
	plot(.ksEnrichmentScoreSimplify(es1), col=col.bars[1], type="l", lwd=2, 
			ylab="Enrichment Score",xlab="Rank in Ordered Dataset", ylim=c(min(es1,es2),max(es1,es2)), ...)    
	if (!is.null(index2)) {
		points(.ksEnrichmentScoreSimplify(es2), col=col.bars[2], type="l", lwd=2)    
	}
	abline(h=0,col=col.bars[3],lwd=2)
}

# calculate the Kolomogorov-Smirnov enrichment (st statistic, gt genset)
.ksEnrichmentScore <- function(st, gt) {
	if (is.null(gt)) return(0)
	n <- length(st)
	m <- length(gt)
	dc <- m/n
	rnk <- match(st[gt], sort(st, decreasing=TRUE))
	es <- cumsum(1:n %in% rnk) -(m*(1:n))/n
	
	# permutation test to scale the enrichment score, there must be an
	# analytical way 
	r <- lapply(1:500,function(i) quantile(cumsum(1:n %in% sample(n,m)) -
								(m*(1:n))/n,p=c(0,1)))
	sf <- mean(abs(quantile(unlist(r), p=c(0.025,0.975))))
	es/sf
}

# simplify for printing
.ksEnrichmentScoreSimplify <- function(es) {
	idx = c(TRUE,sapply(2:(length(es)-1), function(i) es[i]>es[i-1]), TRUE)
	idx[which(idx)-1] = TRUE
	list(x = (1:length(es))[idx],y=es[idx])
}

