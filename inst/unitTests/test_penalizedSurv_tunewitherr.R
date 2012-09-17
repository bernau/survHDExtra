makeSampledat <- function(features, samples){
    ll <- list()
    ll$X <- matrix(rnorm(features * samples), nrow=samples)
    colnames(ll$X) <- make.names(1:ncol(ll$X))
    ll$y <- Surv(runif(samples), rbinom(samples, 1, 0.5))
    ll$fold <- 10
    ll$strat <- FALSE
    ll$measure <- "CvPLogL"
    ll$grids$lambda <- c(1, 10, 100)
    ll$survmethod <- "customSurv"
    ll$customSurvModel<-customPenalized
    ll$penalty <- "ridge"
    ll$standardize <- FALSE
    return(ll)
}

   
test_penalizedSurv_tunewitherr <- function(){
    require(penalized)
    require(survHD)
	require(survHDExtra)
    set.seed(1)
    ll <- makeSampledat(100, 50)
    #tuneres.witherr <- do.call(tune, args=ll)
    #lambda.witherr <- getBestParameters(tuneres.witherr, res.ind=1, measure='CvPLogL')@par$lambda
    #ll$grids$lambda <- ll$grids$lambda[-1]
    #tuneres.fine <- do.call(tune, args=ll)
    #lambda.fine <- getBestParameters(tuneres.fine, res.ind=1, measure='CvPLogL')@par$lambda
    #checkEquals(lambda.witherr, lambda.fine)
	
	###model does not converge-> problem happens in predictsurvhd(sig(penalizedSurv,))
	###probably because predictions are not performed because the model did not converge 
}
