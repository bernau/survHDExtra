test_plusMinusSurv <- function(){
	if(1==0){
    require(pensim)
    require(survHD)
    require(survival)
    set.seed(1)
    ##One good variable
    simdat <- create.data(nvars=c(5,5,30),
                          cors=c(0,0,0),
                          associations=c(3, -3, 0),
                          firstonly=c(FALSE, FALSE, FALSE),
                          nsamples=200,
                          censoring=c(2,10),
                          response="timetoevent")
    X <- as.matrix(simdat$data[,-match(c("time","cens"),colnames(simdat$data))])
    y <- Surv(simdat$data$time,simdat$data$cens)
    ##loop over several options to plusMinusSurv, and make sure the
    ##coefficients are correct:
    for (i.directionality in c("posneg", "pos", "neg")){
        print(paste("directionality =", i.directionality))
        for (i.modeltype in c("plusminus", "compoundcovariate", "tscore", "voting", "positiveriskvoting", "negativeriskvoting")){
            print(paste("modeltype = ", i.modeltype))            
            cc.truth <- simdat$associations
            if(i.directionality == "pos")
                cc.truth[cc.truth < 0] <- 0
            if(i.directionality == "neg")
                cc.truth[cc.truth > 0] <- 0
            ##signs and names of coefficients should be identical to what we simulated:
            if(i.modeltype == "tscore" & i.directionality != "posneg"){
                cc <- try(plusMinusSurv(X=X, y=y,
                                        modeltype=i.modeltype,
                                        lambda=ifelse(i.directionality == "posneg", 10, 5),
                                        directionality=i.directionality)@model@coefficients, silent=TRUE)
                checkEquals( class(cc) == "try-error", TRUE )
            }else{
                cc <- plusMinusSurv(X=X, y=y,
                                    modeltype=i.modeltype,
                                    lambda=ifelse(i.directionality == "posneg", 10, 5),
                                    directionality=i.directionality)@model@coefficients
                checkEquals( identical(sign(cc), sign(cc.truth)), TRUE )
            }
        }
    }
}
2+2
}
