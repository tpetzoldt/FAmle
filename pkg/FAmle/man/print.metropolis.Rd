\name{print.metropolis}
\Rdversion{1.2}
\alias{print.metropolis}
\title{Bayesian Estimation of Univariate Probability Distributions}
\description{See \code{\link{metropolis}}.}
\usage{\method{print}{metropolis}(x, stats.fun = NULL,...)}
\arguments{
  \item{x}{\code{\link{metropolis}} object corresponding to the fitted model.}
  \item{stats.fun}{An optional function that may be provided by the user in order to obtain a posterior summary (see Example).}
  \item{\dots}{none...}
}


\seealso{
\code{\link{metropolis}}, \code{\link{print}}
}

\examples{

data(yarns)
x <- yarns$x
fit.x <- mle(x,'gamma',c(.1,.1))
bayes.x <- metropolis(fit.x,50,trans.list=
	list(function(x) exp(x), function(x) exp(x)))
print(bayes.x)
print(bayes.x,stats.fun=function(x) c(mean=mean(x),CV=sd(x)/mean(x)))
}

\keyword{print}
\keyword{models}

