
testGSA <- function(){
    require(survHD)
    require(survival)

    set.seed(2)
    
    # test reading
    gmt <- gsaReadGmt(system.file("extdata/ovarian_gene_signatures.gmt",
    package = "survHD"))
    
    checkEquals(length(gmt@geneset.descriptions), 14)
    checkEquals(length(gmt@geneset.names), 14)
    checkEquals(length(gmt@genesets), 14)
    data(beer.exprs)
    data(beer.survival)
    require(hu6800.db)
    require(annotate)

    # the Gmt file contains gene symbols, so we translate it to Affymetrix 
    genes <- getSYMBOL(rownames(beer.exprs), "hu6800")
    gmt.affy <- gsaTranslateGmt(gmt, beer.exprs, genes)

    checkEquals(length(gmt@genesets), length(gmt.affy@genesets))
    checkEquals(gmt@geneset.names, gmt.affy@geneset.names)
    # test if when we back translate to symbols, we match the genes in the original gene set
    checkTrue(sum(sapply(1:length(gmt@genesets), function(i) sum(!getSYMBOL( gmt.affy@genesets[[i]], "hu6800") %in% gmt@genesets[[i]]))) == 0, "All affy probes genes in original gene set.")

    # test wilcoxGST example
    stat <- rnorm(100)
    sel <- 1:10; stat[sel] <- stat[sel]+1

    gmtx <- new("gsagenesets", genesets=list(paste(make.names(sel))))

    ret.survhd <- gsaWilcoxSurv(gmtx, genenames=make.names(1:100), statistics=stat,p.value=1)
    ret.limma <- wilcoxGST(sel,stat)
    checkEqualsNumeric(ret.limma, ret.survhd@p.values)

    gmt.biocarta <- try(gsaPathwayGmt("biocarta"))
    checkTrue(class(gmt.biocarta) == "gsagenesets", 
        "Return object is of class gsagenesets")

    gmt.biocarta <- try(gsaPathwayGmt(biocarta[1:20],type="entrez"))
    checkTrue(class(gmt.biocarta) == "gsagenesets", 
        "Return object is of class gsagenesets")

    gmt.biocarta.error <- try(gsaPathwayGmt("biocartaXX"),silent=TRUE)
    checkTrue(class(gmt.biocarta.error) == "try-error", 
        "Error with invalid pathway")
}
