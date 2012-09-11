# Filename: superpcsurv.r 
# Title: One of many classifiers.
# 
# Author: Markus Riester, adapted from M. Slawski and A-L Boulesteix 
# Email: <markus@jimmy.harvard.edu> 
# Date of creation: May 7, 2012
# 
# Brief description:
#   Returns an object of class survoutput.
# 
###############################################################################


setGeneric("superPcSurv", function(X, y, ...) standardGeneric("superPcSurv"))

.buildData <- function(X, y) {
    list(x = t(X), y = y[, 1], censoring.status = y[, 2], featurenames = rownames(X))
}

setMethod("superPcSurv", signature(X = "data.frame", y = "Surv"), function(X, y, 
    learnind, essential = FALSE, ...) {
    require(superpc)
    nrx <- nrow(X)
    yAll <- y
    if (nrx != nrow(y)) 
        stop("Number of rows of 'X' must agree with length of y \n")
    if (missing(learnind)) 
        learnind <- 1:nrx
    if (length(learnind) > nrx) 
        stop("length of 'learnind' must be smaller than the number of observations. \n")
    
    if (!identical(all.equal(colnames(X), make.names(colnames(X), unique = TRUE)), 
        TRUE)) {
        colnames(X) <- make.names(colnames(X), unique = TRUE)
        warning("converting colnames(X) to valid, unique names")
    }
    
    Ylearn <- y[learnind]
    Xlearn <- data.frame(X[learnind, ])
    ll <- list(...)
    
    ll$data <- .buildData(Xlearn, Ylearn)
    
    good.arg.names.train <- names(formals(superpc.train))
    good.arg.names.cv <- names(formals(superpc.cv))
    good.arg.names.predict <- names(formals(superpc.predict))
    
    if (hasArg(lambda)) {
        names(ll)[match("lambda", names(ll))] <- "threshold"
    } else if (!hasArg(threshold)) {
        print("Tuning the penalty parameters...")
        ll$fit <- do.call(superpc.train, args = ll[names(ll) %in% good.arg.names.train])
        fit.cv <- do.call(superpc.cv, args = ll[names(ll) %in% good.arg.names.cv])
        ll$n.components = which.max(apply(fit.cv$scor, 1, max, na.rm = TRUE))
        stopifnot(ll$n.components > 0 && ll$n.components < 4)
        ll$threshold = fit.cv$thresholds[which.max(fit.cv$scor[ll$n.components, ])]
        ll <- ll[-match("fit", names(ll))]
    }
    ll$object <- do.call("superpc.train", args = ll[names(ll) %in% good.arg.names.train])
    
    ll$newdata = ll$data
    linear.predictor <- do.call("superpc.predict", args = ll[names(ll) %in% good.arg.names.predict])$v.pred.1df
    linear.predictor <- new("linpred", lp = linear.predictor)
    ll$newdata <- NULL
    
    
    model.out <- new("SuperPcSurv", mod = ll[names(ll) %in% good.arg.names.predict])
    
    if (essential == TRUE) {
        2 + 2
        # prune model.out
    }
    
    return( new("survoutput", y = Ylearn, linear.predictor = linear.predictor, 
        learnind = learnind, method = "superPcSurv", model = model.out) )
})

setMethod("superPcSurv", signature(X = "matrix", y = "Surv"), function(X, y, ...) {
    superPcSurv(X = as.data.frame(X, check.names = FALSE), y = y, ...)
})

setMethod("superPcSurv", signature(X = "ExpressionSet", y = "Surv"), function(X, y, ...) {
    superPcSurv(X = as.data.frame(t(exprs(X)), check.names = FALSE), y = y, ...)
})

setMethod("superPcSurv", signature(X = "ExpressionSet", y = "character"), function(X, 
    y, ...) {
    Xdat <- as.data.frame(t(exprs(X)), check.names = FALSE)
    superPcSurv(X = Xdat, y = .fetchyFromEset(X,y), ...)
})

setMethod("predictsurvhd", signature(object = "SuperPcSurv", newdata = "data.frame"), 
    function(object, newdata, type, timegrid = NULL, ...) {
        if (type == "lp") {
            ll <- list(...)
            ll$newdata = list(x = t(newdata))
            ll = c(ll, object@mod)
            pred <- do.call("superpc.predict", args = ll)$v.pred.1df
            pred <- new("linpred", lp = pred)
        } else if (type == "survprobs") {
            stop("Currently no SuperPC-specific method for predicting survival probabilities implemented. Please try \"gbm=TRUE\".")
        }
        
        return(pred)
    }) 
