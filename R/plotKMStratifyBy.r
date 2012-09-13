plotKMStratifyBy <-
function(method="median", y, ModelLinear, cutpoints=NULL, labels=NULL,
plot=TRUE,...) {
   if (!is.null(cutpoints)) return(.kmStratify(y, ModelLinear, cutpoints,
       labels,plot,...))
   pval.defined = NULL
   cutpoints  = switch(method,
    "median"  = median(ModelLinear, na.rm=TRUE),
    "tertile" = quantile(ModelLinear,p=c(1/3,2/3), na.rm=TRUE),
    NULL)
    if (is.null(cutpoints)) {
        require(maxstat)
        if (!method %in% eval(formals(maxstat)$pmethod))
            stop(paste("Unknown method. Valid methods are: ",
                paste(c("median","tertile",
                eval(formals(maxstat)$pmethod)),collapse=", "))) 
        mstat = maxstat.test(y~ModelLinear,data=data.frame(ModelLinear=ModelLinear),smethod="LogRank", pmethod=method)
        cutpoints = mstat$estimate
        pval.defined = mstat$p.value
    }
    if (is.null(labels)) labels = .defaultLabels(cutpoints)
    res = .kmStratify(y, ModelLinear, cutpoints, labels, plot,
        pval.defined=pval.defined, ...)
    invisible(list(cutpoints=cutpoints, p.value = pval.defined, hr = res$hr,
    strata=res$strata))
}

.kmStratify <-
function(y, ModelLinear, cutpoints,
labels=NULL, plot=TRUE,
...) {
    if (is.null(labels)) labels = .defaultLabels(cutpoints)

    if (length(cutpoints) + 1 != length(labels))  
        stop("Number of cutpoint labels wrong.")

    strata = rep(labels[1], length(ModelLinear))

    for (i in 1:length(cutpoints)) strata[ModelLinear <= sort(cutpoints,decreasing=TRUE)[i]] = labels[i+1]

    hr = NULL
    strata <- as.factor(strata)
    
    if (plot) hr = plotKM(y=y,strata = strata,...)
    list(hr=hr, strata=strata)
}
