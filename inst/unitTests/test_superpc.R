testSuperpc <- function(breakme=TRUE){
if(1==0){
    require(superpc)
    require(survHD)

    x<-matrix(rnorm(1000*40),ncol=40)
    ##I added these next two lines to demonstrate that some dimnames
    ##are OK, but other legal and illegal dimnames are not (see
    ##breakme)
    rownames(x) <- make.names(1:nrow(x))
    colnames(x) <- make.names(1:ncol(x))
    ##Illegal characters in column and row names should not break
    ##things, but currently does for superpc:
    if(breakme){
        bad.characters <- c("3", "`", "@", "/", "-", " ", "	", "\n", "\\")
        colnames(x) <- paste(sample(bad.characters, ncol(x), replace=TRUE), 1:ncol(x), sep="")
        rownames(x) <- paste(sample(bad.characters, nrow(x), replace=TRUE), 1:nrow(x), sep="")
        ##The next two lines make the dimnames valid and unique, but
        ##the predictions are _still_ wrong.
        colnames(x) <- make.names(colnames(x), unique=TRUE)
        rownames(x) <- make.names(rownames(x), unique=TRUE)
    }    
    y<-10+svd(x[1:60,])$v[,1]+ .1*rnorm(40)
    censoring.status      <- sample(c(rep(1,30),rep(0,10)))
    censoring.status.test <- sample(c(rep(1,30),rep(0,10)))

    featurenames <- paste("feature",as.character(1:1000),sep="")
    data<-list(x=x,y=y, censoring.status=censoring.status,
               featurenames=featurenames)

    idx = 1:20
    data.test <-list(x=x[,idx],y=y[idx],
                     censoring.status=censoring.status.test[idx],
                     featurenames=featurenames)
##    rownames(data$x) <- make.names(rownames(data$x), unique=TRUE)
##    rownames(data.test$x) <- make.names(rownames(data.test$x), unique=TRUE)
    set.seed(332)
    a<- superpc.train(data, type="survival")
    aa<-superpc.cv(a,data)
    n.components = which.max(apply(aa$scor,1,max))
    threshold = aa$thresholds[which.max(aa$scor[n.components,])]
    fit<- superpc.predict(a, data, data.test, threshold=threshold,
                          n.components=n.components,
                          prediction.type="continuous")

    set.seed(332)
    fit.superpc <- superPcSurv(X=t(x),y=Surv(y,censoring.status))
    risk.superpc <- predict(fit.superpc@model,newdata=t(data.test$x), type="lp")

    ##redictions from package directly should be the same as from survHD:
    checkEquals(fit$v.pred.1df, risk.superpc@lp)
	
}

}
