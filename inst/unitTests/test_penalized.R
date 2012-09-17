testPenalized <- function(){

    require(penalized)
    require(survHD)
    require(survHDExtra)
    set.seed(332)
    x<-matrix(rnorm(1000*40),ncol=40)
    ##Illegal characters in column and row names should not break things:
    bad.characters <- c("3", "`", "@", "/", "-", " ", "	", "\n", "\\")
    rownames(x) <- paste(sample(bad.characters, nrow(x), replace=TRUE), 1:nrow(x), sep="")
	###but duplicated names are not reasonable and may produce errors
    rownames(x)[6]<-'8765646'
	
    y<-10+svd(x[1:60,])$v[,1]+ .1*rnorm(40)
    censoring.status      <- sample(0:1, ncol(x), replace=TRUE)
    y <- Surv(y, censoring.status)

    idx = 1:20

##    rownames(data$x) <- make.names(rownames(data$x), unique=TRUE)
##    rownames(data.test$x) <- make.names(rownames(data.test$x), unique=TRUE)
    set.seed(1)
    fit.lasso <- penalized(response=y[idx], penalized=t(x[ ,idx]), lambda1=1)
    fit.ridge <- penalized(response=y[idx], penalized=t(x[ ,idx]), lambda2=1000)
    cc.lasso <- structure(rep(0, nrow(x)), .Names=rownames(x))
    cc.lasso[ match(names(coefficients(fit.lasso)), names(cc.lasso))] <- coefficients(fit.lasso)
    cc.ridge <- coefficients(fit.ridge)
    pred.penalized.lasso <- (t(x[ ,idx]) %*% cc.lasso)[,1]
    pred.penalized.ridge <- (t(x[ ,idx]) %*% cc.ridge)[,1]
    
    set.seed(1)
    fit.ridge.survhd <- customSurv(X=t(x), y=y, learnind=idx, penalty="ridge", lambda=1000,customSurvModel=customPenalized)
    fit.lasso.survhd <- customSurv(X=t(x), y=y, learnind=idx, penalty="lasso", lambda=1,customSurvModel=customPenalized)

    
    ##get predict using predict function
    pred.lasso.survhd <- predict(fit.lasso.survhd, type="lp")@lp
    pred.ridge.survhd <- predict(fit.ridge.survhd, type="lp")@lp
    ##get predictions manually
	cc.ridge.survhd <- coefficients(fit.ridge.survhd@model@mod)
	cc.lasso.survhd <- structure(rep(0, nrow(x)), .Names=names(cc.ridge.survhd))
	cc.lasso.survhd[match(names(coefficients(fit.lasso.survhd@model@mod)), names(cc.lasso.survhd))] <- coefficients(fit.lasso.survhd@model@mod)
   
    cc.ridge.survhd <- coefficients(fit.ridge.survhd@model@mod)
    pred.lasso.survhdmanual <- (t(x[ ,idx]) %*% cc.lasso.survhd)[,1]
    pred.ridge.survhdmanual <- (t(x[ ,idx]) %*% cc.ridge.survhd)[,1]
	names(pred.lasso.survhdmanual) <-1:20
	names(pred.ridge.survhdmanual)<-1:20

    ##check predict function against doing predictions in this script:
	checkEquals(pred.lasso.survhd,pred.lasso.survhdmanual)
    checkEquals(pred.ridge.survhd, pred.ridge.survhdmanual)

    ##coefficients and predictions from package directly should be the
    ##same as from survHD:
	names(cc.lasso.survhd)<-names(cc.lasso)<-NULL
    checkEquals(cc.lasso.survhd, cc.lasso)
	names(cc.ridge.survhd)<-names(cc.ridge)<-NULL
    checkEquals(cc.ridge.survhd, cc.ridge)
	names(pred.penalized.lasso)<-1:20
	names(pred.penalized.ridge)<-1:20
    checkEquals(pred.lasso.survhd, pred.penalized.lasso)
    checkEquals(pred.ridge.survhd, pred.penalized.ridge)
   

}
