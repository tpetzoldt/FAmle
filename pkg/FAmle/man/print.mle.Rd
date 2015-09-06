\name{print.mle}
\Rdversion{1.2}
\alias{print.mle}
\title{Maximum Likelihood Estimation of Univariate Probability Distributions}
\description{See \code{\link{mle}}.}
\usage{\method{print}{mle}(x,...)}
\arguments{
  \item{x}{\code{mle} object corresponding to the fitted model.}
  \item{\dots}{none...}
}

\seealso{
\code{\link{mle}}, \code{\link{print}}
}

\examples{
data(yarns)
x <- yarns$x
fit.x <- mle(x,'gamma',c(.1,.1))
print(fit.x)
print.mle(fit.x)
}

\keyword{print}
\keyword{models}

