testPlusMinusVoting<-function(){
	if(1==0){
library(survHD)
require(pensim)
require(survHD)
require(survival)
set.seed(1)
##One good variable
simdat <- create.data(nvars=c(9,1), cors=c(0,0), associations=c(3, -3), firstonly=c(FALSE, FALSE), nsamples=200, censoring=c(2,10), response="timetoevent")

X <- as.matrix(simdat$data[,-match(c("time","cens"),colnames(simdat$data))])
y <- Surv(simdat$data$time,simdat$data$cens)
X <- sweep(X, 2, apply(X, 2, median))

checkEquals <- function(x, y) print(all (x == y))

model <- list()
for (modeltype in c("voting", "negativeriskvoting", "positiveriskvoting")){
    model[[modeltype]] <- list()
    for (directionality in c("posneg", "pos", "neg")){
        print(c(modeltype, directionality))
        model[[modeltype]][[directionality]] <- plusMinusSurv(X=X, y=y, modeltype=modeltype, lambda=10, directionality=directionality)
        if(directionality == "posneg")
            checkEquals(sign( model[[modeltype]][[directionality]]@model@coefficients ), sign(simdat$associations) )
        if(directionality == "pos")
            checkEquals(sign( model[[modeltype]][[directionality]]@model@coefficients ), as.integer(simdat$associations > 0) )
        if(directionality == "neg")
            checkEquals(sign( model[[modeltype]][[directionality]]@model@coefficients ), -as.integer(simdat$associations < 0) )        
    }
}

##synthetic data with known correct votes:
newdat <- t(sapply(0:10, function(x) c(rep(1, x), rep(-1, 10-x))))
colnames(newdat) <- colnames(X)[1:ncol(newdat)]

newdat.votes <- list()
for (modeltype in c("voting", "negativeriskvoting", "positiveriskvoting")){
    newdat.votes[[modeltype]] <- list()
    for (directionality in c("posneg", "pos", "neg")){
        print(c(modeltype, directionality))
        newdat.votes[[modeltype]][[directionality]] <- sweep(newdat, 2, sign(model[[modeltype]][[directionality]]@model@coefficients), "*")
        if(modeltype == "positiveriskvoting")
            newdat.votes[[modeltype]][[directionality]][ newdat.votes[[modeltype]][[directionality]] < 0] <- 0
        if(modeltype == "negativeriskvoting")
            newdat.votes[[modeltype]][[directionality]][ newdat.votes[[modeltype]][[directionality]] > 0] <- 0
        checkEquals(predict(model[[modeltype]][[directionality]], newdata=newdat, type="lp")@lp, rowSums(newdat.votes[[modeltype]][[directionality]]) )
    }
}
}
2+2
}
