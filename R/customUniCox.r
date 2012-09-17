

customUniCox<-function(Xlearn, Ylearn, learnind,...){
			require(uniCox)
			
		
			if (!identical(all.equal(colnames(Xlearn), make.names(colnames(Xlearn), unique = TRUE)), 
					TRUE)) {
				colnames(Xlearn) <- make.names(colnames(Xlearn), unique = TRUE)
				warning("converting colnames(X) to valid, unique names")
			}
			
		
			ll <- list(...)
			
			ll$x <- Xlearn
			ll$y <- Ylearn[, 1]
			ll$status <- Ylearn[, 2]
			
			if (is.null(ll$lambda)==FALSE) {
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
			
			predfun<-function(object, newdata, type, timegrid = NULL, ...) {
				require(uniCox)
				if (type == "lp") {
					unicox.object <- object@mod
					pred <- predict.uniCox(unicox.object, x = newdata, ...)
					pred <- structure(pred[, 1], .Names = rownames(pred))
					pred <- new("LinearPrediction", lp = pred)
				}
				
				if (type == "SurvivalProbs") {
					stop("Currently no uniCox-specific method for predicting survival probabilities implemented. Please try \"gbm=TRUE\".")
				}
				return(pred)
			} 
		
			
			
			custommod <- new("ModelCustom", mod=output,predfun=predfun,extraData=list())
			
			return(custommod)
			
		}



 
		
