\name{Q.conf.int}
\Rdversion{1.2}
\alias{Q.conf.int}
\title{Approximate Confidence Intervals for p-th Quantile}
\description{This function can be used to derive approximate confidence intervals for the \code{p}-th quantile of the fitted distribution (see \code{\link{mle}}).}
\usage{Q.conf.int(p, model, alpha = 0.1, ln = FALSE)}
\arguments{
  \item{p}{Vector of probabilities.}
  \item{model}{\code{mle} object corresponding to the fitted model.}
  \item{alpha}{\code{1-alpha} is the interval's coverage probability.}
  \item{ln}{whether or not the confidence interval of the \code{p}-th quantile should be computed on the natural logarithmic scale (see Details).}
}
\details{
The p-th quantile confidence interval is derived using the observed Fisher's information matrix in conjuction with the well-known delta method. Here, \code{Q.conf.int} allows the user to chose between two types of confidence intervals: one that is computed on the original scale and one that is computed on the quantile's natural logarithmic scale.
}
\value{The function returns a 3-by-\code{length(p)} array containing, for each value of \code{p}, the confidence interval's lower and upper bounds, as well as the quantile point estimate (maximum likelihood).}

\references{Rice, J.A. (2006) Mathematical statistics and data analysis. Duxbury Press, 3rd edition (regarding the Delta method).}


\seealso{\code{\link{plot.mle}}}
\examples{
data(yarns)
x <- yarns$x
fit.x <- mle(x,'gamma',c(.1,.1))
Q.conf.int(p=c(.5,.9,.95,.99),model=fit.x,alpha=.01,ln=FALSE)
Q.conf.int(p=c(.5,.9,.95,.99),model=fit.x,alpha=.01,ln=TRUE)
}
\keyword{models}

