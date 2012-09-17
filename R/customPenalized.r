####random survival forest as a custom survival model function
##Xlearn and Ylearn are obligatory inputs 
customPenalized<-function(Xlearn,Ylearn,learnind,...){
	library(penalized)
	ll<-list(...)
	
	if(is.null(ll$penalty)) ll$penalty<-'ridge'
	
	unpenalizedvarnames<-ll$unpenalizedvarnames
	
	if (is.null(unpenalizedvarnames)) {
		ll$penalized <- Xlearn
	} else {
		if (!(all(unpenalizedvarnames %in% colnames(Xlearn)) & is.character(unpenalizedvarnames))) 
			stop("unpenalizedvarnames should be a character vector whose elements are found in colnames(X)")
		ll$penalized <- Xlearn[-match(unpenalizedvarnames, colnames(Xlearn))]
		ll$unpenalized <- Xlearn[match(unpenalizedvarnames, colnames(Xlearn))]
	}
	
	ll$trace <- FALSE
	ll$response <- Ylearn
	
	
	if (is.null(ll$lambda)==FALSE) {
		if (identical(ll$penalty, "ridge") || identical(ll$penalty, "lasso")) {
			if (length(ll$lambda) != 1 | (class(ll$lambda) != "numeric" && class(ll$lambda) != 
						"integer")) 
				stop("For ridge or lasso penalty, lambda must be a numeric vector of length 1")
		} else if (identical(ll$penalty, "elasticnet")) {
			if (length(ll$lambda) != 2 | (class(ll$lambda) != "numeric" && class(ll$lambda) != 
						"integer")) 
				stop("For elasticnet penalty, lambda must be a numeric vector of length 2")
		}
	} else {
		print("Tuning the penalty parameter...")
		if (identical(ll$penalty, "ridge")) {
			good.arg.names <- names(formals(optL2))
			ll$lambda <- do.call(optL2, args = ll[names(ll) %in% good.arg.names])$lambda
		} else if (identical(ll$penalty, "lasso")) {
			good.arg.names <- names(formals(optL1))
			ll$lambda <- do.call(optL1, args = ll[names(ll) %in% good.arg.names])$lambda
		} else if (identical(ll$penalty, "elasticnet")) {
			good.arg.names <- names(formals(opt2D))
			opt2D.output <- do.call(opt2D, args = ll[names(ll) %in% good.arg.names])
			ll$lambda <- opt2D.output[which.max(opt2D.output[, "cvl"]), c("L1", "L2")]
		} else {
			stop("Unknown penalty. Must be either ridge, lasso, or elasticnet.")
		}
	}
	
	if (identical(ll$penalty, "ridge")) {
		names(ll)[names(ll) == "lambda"] <- "lambda2"  ##correct naming of tuning parameter
	} else if (identical(ll$penalty, "lasso")) {
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

#output <- do.call("penalized", args = ll[names(ll) %in% good.arg.names])
		
	## browser(expr=eval(expression(class(output) == 'try-error'))) Predictions on
	## learning set can be used later for calculating baseline hazards
	
	
	predfun<-function(object, newdata, type, timegrid = NULL, ...) {
		require(penalized)
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
					pred <- new("LinearPrediction", lp = pred)
				}
			}
		}
		if (type == "SurvivalProbs") {
			stop("Currently no Penalized-specific method for predicting survival probabilities implemented. Please try \"gbm=TRUE\".")
		}
		return(pred)
	}
	
	
	if (class(output) == "penfit") {
		###create customsurvhd-object
		custommod<-new("ModelCustom",mod=output,predfun=predfun,extraData=list())
		
	} else {
		
		###now create customsurvhd-object (with null-model)      
		custommod <- new("ModelCustom", mod = NULL,predfun=predfun,extraData=list())  
	}
	
	return(custommod)
}
