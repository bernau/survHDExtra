####random survival forest as a custom survival model function
##Xlearn and Ylearn are obligatory inputs 
customRSF<-function(Xlearn,Ylearn,learnind,...){
	###load required packages
	require(randomSurvivalForest)
	###handle inputs
	print(class(Xlearn))
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
#either type lp or type survprobs must be implemented
#for typ lp the obligatory return class is linpred
		if (type == "lp") {
			stop("Random Forests don't provide linear predictors, sorry.")
		}
#for typ survprobs the obligatory return class is Breslow
		else if (type == "SurvivalProbs") {
			modelobj <- object@mod
			if (is.null(timegrid)) {
				stop("No timegrid specified.")
			}
			
			###create data for which predictions are to be performed
			###function checks for response but does not use it (fake response)
			predsrsf <- predict.rsf(object = modelobj, test=data.frame(newdata,time=rexp(n=nrow(newdata)),status=sample(c(0,1),nrow(newdata),replace=T)))
			###predict-function provides predictions for training times only
			###->interpolate for timepoints in timegrid
			curves<-exp(-t(apply(predsrsf$ensemble,1,FUN=function(z) approx(x=predsrsf$timeInterest,y=z,xout=timegrid)$y)))
			###create breslow-object
			pred <- new("breslow", curves = curves, time = timegrid)
			###create survprobs-object embedding the breslow-object
			pred <- new("SurvivalProbs", SurvivalProbs = pred)
		}
		
		else stop('Invalid "type" argument.')
		return(pred)  
	}
	
	
	###now create customsurvhd-object (which is the obligatory output-object)
	custommod<-new("ModelCustom",mod=output.rsf,predfun=predfun,extraData=list())
	
	return(custommod)
}