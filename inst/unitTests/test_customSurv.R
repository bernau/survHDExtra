testCustomSurv <- function(){
    library(survHD)    
####random survival forest as a custom survival model function
####Xlearn and Ylearn are obligatory inputs 
    customRSF<-function(Xlearn,Ylearn,learnind,...){
        ##load required packages
        require(randomSurvivalForest)
        ##handle inputs
        ll<-list(...)
        datarsf<-data.frame(Xlearn,time=Ylearn[,1],status=Ylearn[,2])
        ll$data <- datarsf
        ll$formula<-as.formula('Surv(time,status)~.')
        ##call actual model function rsf from randomSurvivalForest
        output.rsf <- do.call("rsf", args = ll)
        ##define prediction function which will be stored in slot predfun
        ##and called by predictsurvhd (signature(survhdcustom))
###
        predfun<-function(object,newdata,type,timegrid=NULL,...){
            require(randomSurvivalForest)
            ##either type lp or type survprobs must be implemented
            ##for typ lp the obligatory return class is linpred
            if (type == "lp") {
                stop("Random Forests don't provide linear predictors, sorry.")
            }
            ##for typ survprobs the obligatory return class is Breslow
            else if (type == "survprobs") {
                modelobj <- object@mod
                if (is.null(timegrid)) {
                    stop("No timegrid specified.")
                }
                ##create data for which predictions are to be performed
                ##function checks for response but does not use it (fake response)
                predsrsf <- predict.rsf(object = modelobj, test=data.frame(newdata,time=rexp(n=nrow(newdata)),status=sample(c(0,1),nrow(newdata),replace=T)))
                ##predict-function provides predictions for training times only
                ##->interpolate for timepoints in timegrid
                curves<-exp(-t(apply(predsrsf$ensemble,1,FUN=function(z) approx(x=predsrsf$timeInterest,y=z,xout=timegrid)$y)))
                ##create breslow-object
                pred <- new("breslow", curves = curves, time = timegrid)
                ##create survprobs-object embedding the breslow-object
                pred <- new("survprobs", survprobs = pred)
            }
            else stop('Invalid "type" argument.')
            return(pred)  
        }
        ##now create customsurvhd-object (which is the obligatory output-object)
        custommod<-new("survhdcustom",mod=output.rsf,predfun=predfun,extraData=list())
        return(custommod)
    }
    
###define function in global envir
## assign(x="customRSF",value=customRSF,envir=.GlobalEnv)
    
###test it
    ##data
    data(beer.exprs)
    set.seed(123)
    X<-t(as.matrix(beer.exprs))[, 1:100]  #just 100 genes for speed
    data(beer.survival)
    y<-Surv(beer.survival[,2],beer.survival[,1])
    ##learningset	
    ls<-generateLearningsets(y=y[,1],method='CV',fold=3,niter=2)
    ##gene selection
    gsel<-geneSelection(X=X,y=y,method='fastCox',learningsets=ls,criterion='coefficient')

    ##learnsurvival	
    set.seed(1)
##    browser()
    svaggr <- learnSurvival(X=X,y=y,genesel=gsel,nbgene=30,survmethod='customSurv',customSurvModel=customRSF,learningsets=ls)

    ###tune
    tuneres<-tune(X=X,y=y,genesel=gsel,nbgene=30,survmethod='customSurv',customSurvModel=customRSF,learningsets=ls,grids = list(ntree = 100*seq(1, 10, 4)))

###use tuneres in learnSurvival
    svaggr<-learnSurvival(X=X,y=y,genesel=gsel,nbgene=30,survmethod='customSurv',customSurvModel=customRSF,learningsets=ls,tuneres=tuneres,measure="PErrC",timegrid=4:10,gbm=FALSE,addtune=list(genesel=gsel,nbgene=30))
###error expected because rsf does not provide lp
    expect.msg1 <- "Error in function (object, newdata, type, timegrid = NULL, ...)  : \n  Random Forests don't provide linear predictors, sorry.\n"
    expect.err1 <- try(evaluate(svaggr,measure='CvPLogL'), silent=TRUE)
    checkEquals(class(expect.err1), "try-error")
    checkEquals(expect.err1[1], expect.msg1)
###error expected because rsf does not provide lp
    expect.msg2 <- "Error in gbmbasehaz(object, timegrid = timegrid) : \n  Model does not provide a linear predictor.\n"
    expect.err2 <- try(evaluate(svaggr,measure='PErrC',timegrid=4:10,gbm=T), silent=TRUE)
    checkEquals(class(expect.err2), "try-error")
    checkEquals(expect.err2[1], expect.msg2)
###works without lp
    perrc.output <- evaluate(svaggr,measure='PErrC',timegrid=4:10,gbm=F)	

###--------------------------
###End of RSF test, start of CoxBoost test
###--------------------------
        
    
####coxboost as a custom survival model function
    ##Xlearn and Ylearn are obligatory inputs 
    customSurvModel<-function(Xlearn,Ylearn,unpenalizedvars=NULL,learnind,...){
        ##load required packages
        require(CoxBoost)
        ##handle inputs
        ll <- list(...)
        ll$time <- Ylearn[, 1]
        ll$status <- Ylearn[, 2]
        ll$trace = FALSE
        ll$standardize = TRUE
###
        ##handle unpenalized var
        if (is.null(unpenalizedvars) == TRUE) Xlearn <- data.frame(Xlearn[, ])
        else Xlearn <- data.frame(unpenalizedvars[learnind, ], Xlearn[, ])
        if (!is.null(unpenalizedvars)) 
            ll$unpen.index <- unpeninds
###        
        ll$x <- as.matrix(Xlearn)
        ll$learnind<-NULL
        ##call actual model function CoxBoost
        output.coxboost <- do.call("CoxBoost", args = ll)
        ##define prediction function which will be stored in slot predfun
        ##and called by predictsurvhd (signature(survhdcustom))
###
        predfun<-function(object,newdata,type,timegrid=NULL,...){
            require(CoxBoost)
                                        #either type lp or type survprobs must be implemented
                                        #for typ lp the obligatory return class is linpred
            if (type == "lp") {
                survobj <- object@mod
                pred <- predict(object = survobj, type = "lp", newdata = newdata)
                pred <- structure(pred[1, ], .Names = colnames(pred))
                pred <- new("linpred", lp = pred)
            }
                                        #for typ survprobs the obligatory return class is Breslow
            else if (type == "survprobs") {
                survobj <- object@mod
                if (is.null(timegrid)) {
                    stop("No timegrid specified.")
                }
                curves <- predict(object = survobj, newdata = newdata, times = timegrid, 
                                  type = "risk")
                pred <- new("breslow", curves = curves, time = timegrid)
                pred <- new("survprobs", survprobs = pred)
            }
            else stop('Invalid "type" argument.')
            return(pred)  
        }
        ##now create customsurvhd-object (which is the obligatory output-object)
        custommod<-new("survhdcustom",mod=output.coxboost,predfun=predfun,extraData=list())
        return(custommod)
    }
###function customSurvModel must be available in all environments!
    assign(x='customSurvModel',value=customSurvModel,envir=.GlobalEnv)
    
    
    ##test it
    set.seed(123)
    ##learningset	
    ls<-generateLearningsets(y=y[,1],method='CV',fold=3,niter=2)
    
    ##gene selection
    gsel<-geneSelection(X=X,y=y,method='fastCox',learningsets=ls,criterion='coefficient')
    ##learnsurvival	
    checkEquals(class(gsel)[1], "genesel")
	
    ##use learnSurvival on customSurv and coxBoostSurv, model should be identical:
    set.seed(1)
    svaggr<-learnSurvival(X=X,y=y,genesel=gsel,nbgene=100,survmethod='customSurv',customSurvModel=customSurvModel,learningsets=ls)
    checkEquals(class(svaggr)[1], "survaggr")
    set.seed(1)
    svaggr.native <- learnSurvival(X=X,y=y,genesel=gsel,nbgene=100,survmethod='coxBoostSurv',learningsets=ls)
    for (i in 1:length(svaggr@survoutputlist)){
        print(i)
        mod <- svaggr@survoutputlist[[i]]@model@mod
        mod.native <- svaggr.native@survoutputlist[[i]]@model@mod
      checkEquals(mod, mod.native) 
    }

    ##tune
	
    set.seed(1)
    tuneres<-tune(X=X,y=y,genesel=gsel,nbgene=100,
                  survmethod='customSurv',customSurvModel=customSurvModel,
                  learningsets=ls,
                  grids = list(stepno = c(10 * seq(1, 10, 4))))
    checkEquals(class(tuneres)[1], "tuningresult")
    set.seed(1)
    tuneres.native <- tune(X=X,y=y,genesel=gsel,nbgene=100,
                  survmethod='coxBoostSurv',
                  learningsets=ls,
                  grids = list(stepno = c(10 * seq(1, 10, 4))))
    for (i in 1:ls@iter){
        checkEquals(getBest(tuneres, res.ind=i), getBest(tuneres.native, res.ind=i))
    }
    
    ##use tuneres in learnSurvival
    rm(svaggr, svaggr.native)
    set.seed(1)
    svaggr <- learnSurvival(X=X,y=y,genesel=gsel,nbgene=100,survmethod='customSurv',customSurvModel=customSurvModel,learningsets=ls,tuneres=tuneres)
    checkEquals(class(svaggr)[1], "survaggr")
    ##natively
    set.seed(1)
    svaggr.native <- learnSurvival(X=X,y=y,genesel=gsel,nbgene=100,survmethod='coxBoostSurv',learningsets=ls,tuneres=tuneres.native)

    for (i in 1:length(svaggr@survoutputlist)){
        print(i)
        mod <- svaggr@survoutputlist[[i]]@model@mod
        mod.native <- svaggr.native@survoutputlist[[i]]@model@mod
        checkEquals(mod, mod.native)
    }

    ##Evaluations match too:
   checkEquals(evaluate(svaggr,measure='CvPLogL'), evaluate(svaggr.native,measure='CvPLogL'))
    checkEquals(evaluate(svaggr,measure='PErrC',timegrid=30:35,gbm=TRUE),
                evaluate(svaggr.native,measure='PErrC',timegrid=30:35,gbm=TRUE))
    checkEquals(evaluate(svaggr,measure='PErrC',timegrid=30:35,gbm=FALSE),
                evaluate(svaggr.native,measure='PErrC',timegrid=30:35,gbm=FALSE))
    
}
