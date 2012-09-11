
testplotKM <- function(){
	if(1==0){
    require(survHD)
    require(maxstat)
    data(DLBCL)

    set.seed(2)

    # test the maxstat wrapper
    ret.maxstat <- maxstat.test(Surv(time, cens) ~ MGE, data=DLBCL,
        smethod="LogRank",pmethod="Lau92")

    ret.survhd <- plotKMStratifyBy("Lau92", y=Surv(DLBCL$time, DLBCL$cens),
        linearriskscore=DLBCL$MGE)

    checkEqualsNumeric(ret.maxstat$p.value, ret.survhd$p.value)
    checkEqualsNumeric(ret.maxstat$estimate, ret.survhd$cutpoint)

    ret.survhd <- plotKMStratifyBy("median", y=Surv(DLBCL$time, DLBCL$cens),
        linearriskscore=DLBCL$MGE)
    checkEqualsNumeric(median(DLBCL$MGE), ret.survhd$cutpoint)
    ret.survhd <- plotKMStratifyBy(cutpoints=0, y=Surv(DLBCL$time, DLBCL$cens),
        linearriskscore=DLBCL$MGE)
    checkEquals( "HR 0.412; 95% CI, 0.175 to 0.969; P = 0.0362",
        ret.survhd$hr)
    
    # test linearriskscore object    
    set.seed(2)
    nsamples <- 100
    X <- matrix(rnorm(nsamples*200),nrow=nsamples)
    colnames(X) <- make.names(1:ncol(X))
    rownames(X) <- make.names(1:nrow(X))
    time <- rexp(nsamples)
    cens <- sample(0:1,size=nsamples,replace=TRUE)
    y <- Surv(time,cens)
    lrs <- rnorm(25)
    names(lrs) <- rownames(X)[1:25]
    lrs.obj <- new("linearriskscore", coefficients=lrs, modeltype="compoundcovariate")
    ret.survhd <- plot(lrs.obj, newdata=X, newy=y, show.n.risk=FALSE)
     
    checkEqualsNumeric( 0.803998, ret.survhd$cutpoint, tolerance=0.00001)
    
    # modeltype optional
    lrs.obj <- new("linearriskscore", coefficients=lrs)
    ret.survhd <- plot(lrs.obj, newdata=X, newy=y, show.n.risk=FALSE)
     
    checkEqualsNumeric( 0.803998, ret.survhd$cutpoint, tolerance=0.00001)
}
2+2
}