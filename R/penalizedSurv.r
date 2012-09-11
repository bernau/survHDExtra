# Filename: penalizedSurv.r
# Title: One of many classifiers.
#
# Author: Levi Waldron, adapted from M. Slawski and A-L Boulesteix
# Email: <lwaldron.research@gmail.com>
# Date of creation: Jan 23, 2011
#
# Brief description:
#   Returns an object of class survoutput.
#
# Arguments:
#   -X: matrix of variables (rows are observations,columns are variables)
#   -y: survival response of class Surv 
#   -learnind: vector indicating which observations are to be used for training the survival model
#   -penalty: type of penalization (L1 or L2)
#   -unpenalizedvars: names of variables whose coefficients are not penalized
#   -models: shall the model objects be returned
#
# Value:
#   survoutput
###############################################################################


setGeneric("penalizedSurv", function(X, y, ...) standardGeneric("penalizedSurv"))

setMethod("penalizedSurv", signature(X = "data.frame", y = "Surv"), function(X, y, 
    learnind, penalty = "ridge", unpenalizedvarnames = NULL, essential = FALSE, ...) {
    nrx <- nrow(X)
    yAll <- y
    if (nrx != nrow(y)) 
        stop("Number of rows of 'X' must agree with length of y \n")
    if (missing(learnind)) 
        learnind <- 1:nrx
    if (length(learnind) > nrx) 
        stop("length of 'learnind' must be smaller than the number of observations. \n")
    
    Ylearn <- y[learnind]
    Xlearn <- data.frame(X[learnind, ], check.names = FALSE)
    ll <- list(...)
    
    if (is.null(unpenalizedvarnames)) {
        ll$penalized <- Xlearn
    } else {
        if (!(all(unpenalizedvarnames %in% colnames(Xlearn)) & is.character(unpenalizedvarnames))) 
            stop("unpenalizedvarnames should be a character vector whose elements are found in colnames(X)")
        ll$penalized <- Xlearn[-match(unpenalizedvarnames, colnames(Xlearn))]
        ll$unpenalized <- Xlearn[match(unpenalizedvarnames, colnames(Xlearn))]
    }
    ll$response <- Ylearn
    ll$trace = FALSE
    ## ll$standardize=TRUE ##FEATURE REQUEST: allow user to set this.
    
    if (hasArg(lambda)) {
        if (identical(penalty, "ridge") || identical(penalty, "lasso")) {
            if (length(ll$lambda) != 1 | (class(ll$lambda) != "numeric" && class(ll$lambda) != 
                "integer")) 
                stop("For ridge or lasso penalty, lambda must be a numeric vector of length 1")
        } else if (identical(penalty, "elasticnet")) {
            if (length(ll$lambda) != 2 | (class(ll$lambda) != "numeric" && class(ll$lambda) != 
                "integer")) 
                stop("For elasticnet penalty, lambda must be a numeric vector of length 2")
        }
    } else {
        print("Tuning the penalty parameter...")
        if (identical(penalty, "ridge")) {
            good.arg.names <- names(formals(optL2))
            ll$lambda <- do.call(optL2, args = ll[names(ll) %in% good.arg.names])$lambda
        } else if (identical(penalty, "lasso")) {
            good.arg.names <- names(formals(optL1))
            ll$lambda <- do.call(optL1, args = ll[names(ll) %in% good.arg.names])$lambda
        } else if (identical(penalty, "elasticnet")) {
            good.arg.names <- names(formals(opt2D))
            opt2D.output <- do.call(opt2D, args = ll[names(ll) %in% good.arg.names])
            ll$lambda <- opt2D.output[which.max(optD.output[, "cvl"]), c("L1", "L2")]
        } else {
            stop("Unknown penalty. Must be either ridge, lasso, or elasticnet.")
        }
    }
    
    if (identical(penalty, "ridge")) {
        names(ll)[names(ll) == "lambda"] <- "lambda2"  ##correct naming of tuning parameter
    } else if (identical(penalty, "lasso")) {
        names(ll)[names(ll) == "lambda"] <- "lambda1"  ##correct naming of tuning parameter
    } else if (identical(penalty, "elasticnet")) {
        ll$lambda1 <- ll$lambda[1]
        ll$lambda2 <- ll$lambda[2]
    }
    
    ## filtering out unused arguments allows the user to pass arguments, some of
    ## which will be used for tuning, and others for training the final model:
    good.arg.names <- names(formals(penalized))
    output <- try(do.call("penalized", args = ll[names(ll) %in% good.arg.names]), 
        silent = TRUE)
    ## browser(expr=eval(expression(class(output) == 'try-error'))) Predictions on
    ## learning set can be used later for calculating baseline hazards
    
    if (class(output) == "penfit") {
        output2 <- new("PenalizedSurv", mod = output)
        linear.predictor <- predict(output2, newdata = Xlearn, type = "lp")
        if (essential == TRUE) {
            2 + 2
            ## prune output
        }
    } else {
        ## Produce NULL model and associated linear.predictor if there was a problem
        ## with the penalized call.
        output2 <- new("PenalizedSurv", mod = NULL)
        linear.predictor <- predict(NULL, newdata = Xlearn)
    }
    
    return( new("survoutput", y = Ylearn, linear.predictor = linear.predictor, 
        learnind = learnind, method = "penalizedSurv", model = output2) )
})

setMethod("penalizedSurv", signature(X = "matrix", y = "Surv"), function(X, y, ...) {
    penalizedSurv(X = as.data.frame(X, check.names = FALSE), y = y, ...)
})

setMethod("penalizedSurv", signature(X = "ExpressionSet", y = "Surv"), function(X, y, ...) {
    penalizedSurv(X = as.data.frame(t(exprs(X)), check.names = FALSE), y = y, ...)
})

setMethod("penalizedSurv", signature(X = "ExpressionSet", y = "character"), function(X, 
    y, ...) {
    Xdat <- as.data.frame(t(exprs(X)), check.names = FALSE)
    penalizedSurv(X = Xdat, y = .fetchyFromEset(X,y), ...)
})

setMethod("predictsurvhd", signature(object = "PenalizedSurv", newdata = "data.frame"), 
    function(object, newdata, type, timegrid = NULL, ...) {
        ##
        ll <- list(...)  #dots call is not actually used
        ## browser(expr=eval(expression(is.null(object@mod))))
        cc <- try(coefficients(object@mod, "all"), silent = TRUE)
        if (type == "lp") {
            if (class(cc) == "try-error" | is.null(cc)) {
                warning("This model produced an error, predictions not made.")
                pred <- predict(NULL, newdata = newdata)
            } else {
                if (!object@mod@converged) {
                  ## return NAs if model did not converge
                  warning("This model did not converge, predictions not made.")
                  pred <- predict(NULL, newdata = newdata)
                } else {
                  ## does the model have an intercept?
                  if (grepl("(Intercept)", names(cc)[1], fixed = TRUE)) {
                    intercept <- cc[1]
                    cc <- cc[-1]
                  } else {
                    intercept <- 0
                  }
                  ## do not use columns of newdata for which there are not coefficients
                  newdata <- newdata[, colnames(newdata) %in% names(cc)]
                  cc <- cc[match(colnames(newdata), names(cc))]
                  if (!identical(names(cc), colnames(newdata))) 
                    stop("names of coefficients do not match colnames(newdata)")
                  pred <- (as.matrix(newdata) %*% cc) + intercept
                  pred <- structure(pred[, 1], .Names = rownames(pred))
                  pred <- new("linpred", lp = pred)
                }
            }
        }
        if (type == "survprobs") {
            stop("Currently no Penalized-specific method for predicting survival probabilities implemented. Please try \"gbm=TRUE\".")
        }
        return(pred)
    }) 
