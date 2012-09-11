###all classes ######

setClass(Class = "gsagenesets", representation(genesets = "list", geneset.names = "character", 
				geneset.descriptions = "character", geneset.direction = "character"))

setClass(Class = "gsaresults", representation(geneset.names = "character", statistics = "numeric", 
				p.values = "numeric", genesets.used = "list", geneset.cluster = "factor"))