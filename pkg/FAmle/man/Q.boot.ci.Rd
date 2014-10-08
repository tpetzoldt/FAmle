\name{Q.boot.ci}
\Rdversion{1.2}
\alias{Q.boot.ci}
\title{Parametric Bootstrap Confidence Intervals for p-th Quantile}
\description{This function can be used to derive parametric bootstrap confidence intervals for the \code{p}-th quantile of the fitted distribution (see \code{\link{mle}}).}
\usage{Q.boot.ci(p,boot,alpha=.1)}
\arguments{
  \item{p}{Vector of probabilities.}
  \item{boot}{An object obtained using \code{\link{boot.mle}}.}
  \item{alpha}{\code{1-alpha} is the interval's coverage probability.}
}

\value{This functions returns two types of bootstrap confidence intervals for the \code{p}-th quantile - one is based on the "percentile" method, while the other corresponds to the basis bootstrap interval or "reflexion" (see References).}

\references{
Davison, A.C., and Hinkley, D.V. (1997) Bootstrap methods and their application. Cambridge University Press.}



\note{See References for other means of deriving bootstrap intervals.}


\seealso{\code{\link{boot.mle}}, \code{\link{mle}}, \code{\link{Q.conf.int}}}
\examples{
data(yarns)
x <- yarns$x
fit.x <- mle(x,'gamma',c(.1,.1))
Q.conf.int(p=c(.5,.9,.95,.99),model=fit.x,alpha=.01,ln=FALSE)
# should be run again with B = 1000, for example...
boot.x <- boot.mle(model=fit.x,B=50)
Q.boot.ci(p=c(.5,.9,.95,.99),boot=boot.x,alpha=.01)
}
\keyword{models}

