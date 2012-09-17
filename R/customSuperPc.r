

customSuperPc<-function(Xlearn, Ylearn, learnind, ...) {
			require(superpc)
			
			.buildData <- function(X, y) {
				list(x = t(X), y = y[, 1], censoring.status = y[, 2], featurenames = rownames(X))
			}
			
		

		
		
			
			if (!identical(all.equal(colnames(Xlearn), make.names(colnames(Xlearn), unique = TRUE)), 
					TRUE)) {
				colnames(Xlearn) <- make.names(colnames(Xlearn), unique = TRUE)
				warning("converting colnames(X) to valid, unique names")
			}
			
			ll <- list(...)
			
			ll$data <- .buildData(Xlearn, Ylearn)
			
			good.arg.names.train <- names(formals(superpc.train))
			good.arg.names.cv <- names(formals(superpc.cv))
			good.arg.names.predict <- names(formals(superpc.predict))
			
			if (is.null(ll$lambda)==FALSE) {
				names(ll)[match("lambda", names(ll))] <- "threshold"
			} else if (is.null(ll$threshold)) {
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
			ll$newdata <- NULL
			
	        output<-ll[names(ll) %in% good.arg.names.predict]
					
					predfun<-function(object, newdata, type, timegrid = NULL, ...) {
						require(superpc)
				if (type == "lp") {
					ll <- list(...)
					ll$newdata = list(x = t(newdata))
					ll = c(ll, object@mod)
					pred <- do.call("superpc.predict", args = ll)$v.pred.1df
					pred <- new("LinearPrediction", lp = pred)
				} else if (type == "SurvivalProbs") {
					stop("Currently no SuperPc-specific method for predicting survival probabilities implemented. Please try \"gbm=TRUE\".")
				}
				
				return(pred)
			}
		
			custommod<-new("ModelCustom",mod=output,predfun=predfun,extraData=list())
			return(custommod)
		}


