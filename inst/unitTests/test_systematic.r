test_all <- function() {
    # generate some example data
    set.seed(2)
    nsamples <- 100
    X <- matrix(rnorm(nsamples * 200), nrow = nsamples)
    colnames(X) <- make.names(1:ncol(X))
    rownames(X) <- make.names(1:nrow(X))
    time <- rexp(nsamples)
    cens <- sample(0:1, size = nsamples, replace = TRUE)
    y <- Surv(time, cens)
    dfy <- data.frame(y = y)
    colnames(dfy) <- "y"
    rownames(dfy) <- rownames(X)
    
    # generate parameter combinations, this should be extended to cover basically everything
    tests <- list(X = list(X, data.frame(X), ExpressionSet(t(X)), ExpressionSet(t(X), 
        phenoData = as(dfy, "AnnotatedDataFrame"))), y = list(y, y, y, "y"))
    
    learningsets <- generateLearningsets(y = y, method = "CV", fold = 5, niter = 1, 
        strat = TRUE)
    
    res.survhd <- lapply(1:length(tests$X), function(i) learnSurvival(X = tests$X[[i]], 
        y = tests$y[[i]], learningsets = learningsets, survmethod = "penalizedSurv", 
        penalty = "ridge", lambda = 1000))
    
    sapply(1:3, function(i) checkEquals(res.survhd[[i]], res.survhd[[i + 1]]))
} 
