\name{gsaWilcoxSurv-methods}
\docType{methods}
\alias{gsaWilcoxSurv-methods}
\title{General method for performing a Gene Set Analysis (GSA), experimental implementation}
\description{ 
    Performs gene selection for the following signatures:
}
\section{Methods}{
    \describe{
        \item{gmt = "gsagenesets",X = "missing", y = "missing", genenames =
            "character", statistics="missing"}{signature 1}
        \item{gmt = "gsagenesets",X = "matrix", y = "Surv", genenames =
            "missing", statistics="missing"}{signature 2}
        \item{gmt = "gsagenesets",X = "ExpressionSet", y = "Surv", 
            genenames = "missing", statistics="missing"}{signature 3}
        \item{gmt = "gsagenesets",X = "data.frame", y = "Surv", genenames =
            "missing", statistics="missing"}{signature 4}
        \item{gmt = "list",X = "missing", y = "missing", genenames =
            "character", statistics="numeric"}{signature 5}
        \item{gmt = "list",X = "ANY", y = "Surv",
            genenames="missing", statistics="missing"}{signature 6}
    }
    For further argument and output information, consult
    \code{\link{gsaWilcoxSurv}}.
}
\keyword{GSA}
