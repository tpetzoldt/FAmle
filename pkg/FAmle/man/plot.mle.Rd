\name{plot.mle}
\Rdversion{1.2}
\alias{plot.mle}
\title{Diagnostic Plots for the Fitted Model}
\description{This function returns diagnotic plots for a \code{\link{mle}} object.}
\usage{
\method{plot}{mle}(x, ci = FALSE, alpha = 0.05,...)
}
\arguments{
  \item{x}{\code{\link{mle}} object corresponding to the fitted model.}
  \item{ci}{Whether or not approximate confidence intervals should be added to the return period and quantile plots.}
  \item{alpha}{\code{1-alpha} is the requested coverage probability for the confidence interval.}
\item{\dots}{none...}
}


\seealso{\code{\link{mle}}, \code{\link{Q.conf.int}}}
\examples{
data(yarns)
x <- yarns$x
fit.1 <- mle(x,'weibull',c(.1,.1))
fit.2 <- mle(x,'logis',c(.1,.1))
plot(fit.1,TRUE,.05)
dev.new();plot(fit.2,TRUE,.05)
}

\keyword{aplot}
\keyword{models}
