testPlusminusTuningFull <- function(){
    require(survHD)
    require(pensim)
    require(survival)
    ##values to loop over during test:
    kGoodVar <- c(1, 10)
    kTestMeasures <- c("CvPLogL", "UnoC", "HarrellC")   
    set.seed(1)
    ##test all combinations of # of good variables and tuning/evaluation metric
    for (n.good.var in kGoodVar){
        for (test.measure in kTestMeasures){
            print(c(n.good.var, test.measure))
            ##generate some data
            simdat <- create.data(nvars=c(n.good.var, 100-n.good.var),
                                  cors=c(0,0),
                                  associations=c(2,0),
                                  firstonly=c(FALSE,FALSE),
                                  nsamples=200,
                                  censoring=c(2,10),
                                  response="timetoevent")
            tuneargs <- list()
            tuneargs$X <- as.matrix(simdat$data[,-match(c("time","cens"),colnames(simdat$data))])
            tuneargs$y <- Surv(simdat$data$time,simdat$data$cens)
            tuneargs$survmethod <- "plusMinusSurv"
            tuneargs$fold <- 5
            tuneargs$grids <- list(lambda=c(1, 10, 50))
            tuneargs$measure <- test.measure
            if(test.measure == "UnoC")
                tuneargs$tau <- 10
            getBest.args <- list()
            ##tuning using whole dataset
            getBest.args$tuneres <- do.call(tune, tuneargs)
            ##prepare to evaluate
            getBest.args$res.ind=1
            getBest.args$measure <- test.measure
            if(test.measure == "UnoC")
                getBest.args$tau <- 10
            ##optimum lambda
            lambda.opt <- do.call(getBest, getBest.args)@par$lambda
            ##check that the correct number of variables was selected
            checkEquals(lambda.opt, n.good.var)
        }
    }
}
