# Filename: uniCoxSurv.r 
# Title: One of many classifiers.
# 
# Author: Levi Waldron, adapted from M. Slawski and A-L Boulesteix 
# Email: <lwaldron.research@gmail.com> 
# Date of creation: Apr 19, 2012
# 
# Brief description:
#   Returns an object of class survoutput.
# 
# arguments: 
#   -X: matrix of variables (rows are observations,columns are
#       variables) 
#   -y: survival response of class Surv 
#   -learnind: vector indicating which observations are to be used for 
#              training the survival model
# 
# Value: 
#    survoutput
#
###############################################################################



### X=data.frame, y=Surv
setGeneric("uniCoxSurv", function(X, y, ...) standardGeneric("uniCoxSurv"))
setMethod("uniCoxSurv", signature(X = "data.frame", y = "Surv"), function(X, y, learnind, 
    essential = FALSE, ...) {
    require(uniCox)
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
    
    ll$x <- Xlearn
    ll$y <- Ylearn[, 1]
    ll$status <- Ylearn[, 2]
    
    if (hasArg(lambda)) {
        names(ll)[match("lambda", names(ll))] <- "lamlist"
    } else {
        print("Tuning the penalty parameter...")
        good.arg.names <- names(formals(uniCox))
        ll$fit <- do.call(uniCox, args = ll[names(ll) %in% good.arg.names])
        good.arg.names <- names(formals(uniCoxCV))
        fit.cv <- do.call(uniCoxCV, args = ll[names(ll) %in% good.arg.names])
        optimum.iter <- which.min(fit.cv$se[fit.cv$se > 0])
        ll$lamlist <- unlist(ll$fit)[[paste("lamlist", optimum.iter, sep = "")]]
        ll <- ll[-match("fit", names(ll))]
    }
    
    output <- do.call("uniCox", args = ll)
    
    
    if (essential == TRUE) {
        2 + 2
        # prune output here
    }
    
    output2 <- new("UniCoxSurv", mod = output)
    ## Predictions on learning set can be used later for calculating baseline
    ## hazards
    linear.predictor <- predict(output2, newdata = Xlearn, type = "lp")
    ret.obj <- new("survoutput", y = Ylearn, linear.predictor = linear.predictor, 
        learnind = learnind, method = "uniCoxSurv", model = output2)
    
    return(ret.obj)
    
})

setMethod("uniCoxSurv", signature(X = "matrix", y = "Surv"), function(X, y, ...) {
    uniCoxSurv(X = as.data.frame(X, check.names = FALSE), y = y, ...)
})

setMethod("uniCoxSurv", signature(X = "ExpressionSet", y = "Surv"), function(X, y, ...) {
    uniCoxSurv(X = as.data.frame(t(exprs(X)), check.names = FALSE), y = y, ...)
})

setMethod("uniCoxSurv", signature(X = "ExpressionSet", y = "character"), function(X, 
    y, ...) {
    Xdat <- as.data.frame(t(exprs(X)), check.names = FALSE)
    uniCoxSurv(X = Xdat, y = .fetchyFromEset(X,y), ...)
})

setMethod("predictsurvhd", signature(object = "UniCoxSurv", newdata = "data.frame"), 
    function(object, newdata, type, timegrid = NULL, ...) {
        require(uniCox)
        if (type == "lp") {
            unicox.object <- object@mod
            pred <- predict.uniCox(unicox.object, x = newdata, ...)
            pred <- structure(pred[, 1], .Names = rownames(pred))
            pred <- new("linpred", lp = pred)
        }
        
        if (type == "survprobs") {
            stop("Currently no uniCox-specific method for predicting survival probabilities implemented. Please try \"gbm=TRUE\".")
        }
        return(pred)
    }) 
