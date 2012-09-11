\name{plotKM}
\alias{plotKM}
\title{Kaplan-Meier plot
}
\description{
 A simple but feature-complete Kaplan-Meier plot with sane defaults.
}
\usage{
plotKM(y, strata, xlab = paste("Time (", .guess.time.unit(y), ")", sep = ""), 
    ylab = "Survival (\%)", main = NULL, censor.at = NULL, show.n.risk = TRUE, 
    n.risk.step = ifelse(show.n.risk, round(ifelse(is.null(censor.at), max(y[, 1], na.rm = TRUE), censor.at)/60) * 10, NULL), 
    pval.defined = NULL, chisq.dist = NULL, 
    show.legend = TRUE, show.n = !show.n.risk, legendpos = "bottomleft", npos
    = "bottomleft", HRpos = ifelse(is.null(main), "margintop", "topright"), 
    show.HR = ifelse(length(levels(strata)) == 2, TRUE, FALSE), show.PV =
    TRUE, inverse.HR = TRUE, col = NULL, cex.base = 1, cex.lab = cex.base,
    cex.legend = cex.base, cex.HR = cex.base, cex.n.risk = cex.base * 0.85,
    cex.n = cex.base, lwd = 2, lty = 1:length(levels(strata)), ...)
}
\arguments{
  \item{y}{
  Surv object with the survival times for each patient.
}
  \item{strata}{
    A factor of group labels for patient stratification.
}
  \item{xlab}{
  The label of the x-axis. Default will be "Time (unit)", where unit is
  guessed from the range of survival times in y (Days, Months, Years).
}
  \item{ylab}{
  The label of the y-axis.
}
  \item{main}{
  The title of the Kaplan-Meier plot.
}
  \item{censor.at}{
  Censor all uncensored patients at the specified time.
}
  \item{show.n.risk}{
  If TRUE, show the numbers of samples at risk for each time step.
}
  \item{n.risk.step}{
    Vector specifying the time to be the steps for displaying the number of individuals at risk.
}
  \item{pval.defined}{
    Print the specified p-value instead of calculating it.
}
  \item{chisq.dist}{
  Provide an empirical chi-square distribution instead of using the log-rank
  test. Useful for leave-one-out cross-validation.
}
  \item{show.legend}{
  Show the legend.
}
  \item{show.n}{
  Show the number of patients. If \code{show.legend == TRUE}, then the number will be
  added to the legend, per default only if \code{show.n.risk == FALSE} (because this
  information is listed there). Without a legend, then the numbers are printed
  at the position specified with npos.
}
  \item{legendpos}{
    The position of the legend with the group labels (topleft, topright,
    bottomleft, bottomright). 
}
  \item{npos}{
    The position of the numbers of patients. See \code{show.n}.
}
  \item{HRpos}{
  Position of the hazard ratio (HR). Default is outside the plot on the top
  if no title ("main" argument) is given ("margintop"), otherwise "topright".
}
  \item{show.HR}{
  Show the HR and its confidence interval. Default TRUE when two groups
  are compared, otherwise FALSE.
}
  \item{show.PV}{
  Show the p-value of observed survival differences between the groups.
}
  \item{inverse.HR}{
  Display the inverse of the hazard ratio (HR).  
}
  \item{col}{
  The color of the Kaplan-Meier curves.
}
  \item{cex.base}{
  Default font size. Can be used to increase all labels at once, for figure panels for example.
}
  \item{cex.lab}{
  Font size of axis labels.
}
  \item{cex.legend}{
  Font size of the legend.
}
  \item{cex.HR}{
  Font size of the HR.
}
  \item{cex.n.risk}{
  Font size of the numbers at risk table
}
  \item{cex.n}{
  Font size of the numbers of patients per strata label (\code{show.n}).
}
  \item{lwd}{
  The line width of the Kaplan-Meier curves.
}
  \item{lty}{
  The line style of the Kaplan-Meier curves.
}
  \item{\dots}{
  Additional arguments passed to the plot() function.
}
}
\details{
This function produces publication quality Kaplan-Meier plots. The numbers at
risk table was stolen from the \code{survcomp} package.
}
\value{
The observed hazard ratio as formatted string.
}
\author{
Markus Riester \email{markus@jimmy.harvard.edu}
}
\seealso{\code{\link{plotKM}}, \code{\link{plotKMStratifyBy}}}
\examples{
\dontrun{
require(survHD)
plotKM(y=Surv(bladder$stop, bladder$event),strata=as.factor(bladder$enum))

plotKM(y=Surv(bladder$stop, bladder$event),strata=as.factor(bladder$enum),
    show.n.risk=FALSE,show.legend=FALSE,cex.base=1.4)
}
}
