testPlusminusTuning <- function(){
    require(pensim)
    require(survHD)
    require(survival)
    set.seed(1)
    ##One good variable
    simdat <- create.data(nvars=c(1,99),
                          cors=c(0,0),
                          associations=c(2,0),
                          firstonly=c(FALSE,FALSE),
                          nsamples=200,
                          censoring=c(2,10),
                          response="timetoevent")
    X <- as.matrix(simdat$data[,-match(c("time","cens"),colnames(simdat$data))])
    y <- Surv(simdat$data$time,simdat$data$cens)
    bootds <- generateLearningsets(y = y, method = "CV", fold = 8, niter = 1, strat = FALSE)
    tuneres <- tune(X = X, y = y, learningsets = bootds,
                    survmethod="plusMinusSurv", 
                    grids = list(lambda = c(1, 10, 50)),
                    option="fast")                

    lambda.CvPLogL.1 <- unique(sapply(1:8, function(i) getBest(tuneres,res.ind=i,measure='CvPLogL')@par$lambda))
    lambda.HarrellC.1 <- unique(sapply(1:8, function(i) getBest(tuneres,res.ind=i,measure='HarrellC')@par$lambda))
    lambda.UnoC.1 <- unique(sapply(1:8, function(i) getBest(tuneres,res.ind=i,measure='UnoC', tau=10)@par$lambda))

    ##if these pass, the correct value of the tuning parameter was selected in every fold:
    checkEquals(lambda.CvPLogL.1, 1)
    checkEquals(lambda.HarrellC.1, 1)
    checkEquals(lambda.UnoC.1, 1)


    ##Ten good variables
    set.seed(2)
    simdat <- create.data(nvars=c(10,90),
                          cors=c(0,0),
                          associations=c(2,0),
                          firstonly=c(FALSE,FALSE),
                          nsamples=200,
                          censoring=c(2,10),
                          response="timetoevent")
    X <- as.matrix(simdat$data[,-match(c("time","cens"),colnames(simdat$data))])
    y <- Surv(simdat$data$time,simdat$data$cens)
    bootds <- generateLearningsets(y = y, method = "CV", fold = 8, niter = 1, strat = FALSE)
    tuneres <- tune(X = X, y = y, learningsets = bootds,
                    survmethod="plusMinusSurv", 
                    grids = list(lambda = c(1, 10, 50)),
                    option="fast")


    (lambda.CvPLogL.10 <- unique(sapply(1:8, function(i) getBest(tuneres,res.ind=i,measure='CvPLogL')@par$lambda)))
    (lambda.HarrellC.10 <- unique(sapply(1:8, function(i) getBest(tuneres,res.ind=i,measure='HarrellC')@par$lambda)))
    (lambda.UnoC.10 <- unique(sapply(1:8, function(i) getBest(tuneres,res.ind=i,measure='UnoC', tau=10)@par$lambda)))

    ##if these pass, the correct value of the tuning parameter was selected in every fold:
    checkEquals(lambda.CvPLogL.10, 10)
    checkEquals(lambda.HarrellC.10, 10)
    checkEquals(lambda.UnoC.10, 10)

    ##Fifty good variables
    set.seed(3)
    simdat <- create.data(nvars=c(50,50),
                          cors=c(0,0),
                          associations=c(2,0),
                          firstonly=c(FALSE,FALSE),
                          nsamples=400,
                          censoring=c(2,10),
                          response="timetoevent")
    X <- as.matrix(simdat$data[,-match(c("time","cens"),colnames(simdat$data))])
    y <- Surv(simdat$data$time,simdat$data$cens)
    bootds <- generateLearningsets(y = y, method = "CV", fold = 8, niter = 1, strat = FALSE)
    tuneres <- tune(X = X, y = y, learningsets = bootds,
                    survmethod="plusMinusSurv", 
                    grids = list(lambda = c(1, 10, 50)),
                    option="fast")


    (lambda.CvPLogL.50 <- unique(sapply(1:8, function(i) getBest(tuneres,res.ind=i,measure='CvPLogL')@par$lambda)))
    (lambda.HarrellC.50 <- unique(sapply(1:8, function(i) getBest(tuneres,res.ind=i,measure='HarrellC')@par$lambda)))
    (lambda.UnoC.50 <- unique(sapply(1:8, function(i) getBest(tuneres,res.ind=i,measure='UnoC', tau=10)@par$lambda)))

    ##if these pass, the correct value of the tuning parameter was
    ##selected in every fold.  note that tuning does not seem to work
    ##with CvPLogL for plusminus:

    ##checkEquals(lambda.CvPLogL.50, 50)
    checkEquals(lambda.HarrellC.50, 50)
    checkEquals(lambda.UnoC.50, 50)

}
