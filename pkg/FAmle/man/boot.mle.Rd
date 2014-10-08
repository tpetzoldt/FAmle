\name{boot.mle}
\Rdversion{1.2}
\alias{boot.mle}
\title{Bootstrap Distribution for Fitted Model}
\description{This function allows the user to obtain draws from the (parametric) bootstrap distribution of the fitted model's parameters.}
\usage{boot.mle(model, B = 200, seed = NULL, start = NULL,
 method = "Nelder-Mead")}
\arguments{
  \item{model}{\code{\link{mle}} object corresponding to the fitted model.}
  \item{B}{Requested number of bootstrap samples.}
  \item{seed}{A seed may be specified (see \code{\link{set.seed}})}
  \item{start}{Starting values for the optimization algorithm (if \code{is.null(start)==TRUE}, the fitted model's parameters are used as starting values).}
  \item{method}{The optimization method to be used (see \code{\link{optim}} and \code{\link{mle}}).}
}
\details{Parametric bootstrap -- see References.}
\value{
  \item{model}{\code{\link{mle}} object corresponding to the fitted model.}
  \item{B}{Requested number of bootstrap samples.}
  \item{seed}{The specified seed (see \code{\link{set.seed}})}
  \item{par.star}{Array containing realized values from the bootstrap distribution of the maximum likelihood parameter estimators.}
  \item{gof}{The bootstrap distributions of two goodness-of-fit statistics: Anderson-Darling statistic and Pearson's correlation coefficient for the pair ("observed quantiles","fitted quantiles").}
  \item{p.value}{Bootstrap p-values for the two goodness-of-fit statistics.}
  \item{failure.rate}{The proportion of bootstrap samples for which optimization failed using the specified starting values.}
  \item{total.time}{The total amount of time required to generate \code{B} bootstrap samples.}
}
\references{Davison, A.C., and Hinkley, D.V. (1997). Bootstrap methods and their application. Cambridge University Press.}



\seealso{\code{\link{mle}}, \code{\link{Q.conf.int}}, \code{\link{Q.boot.ci}}}
\examples{
data(yarns)
x <- yarns$x
fit.x <- mle(x,'weibull',c(.1,.1))
boot.x <- boot.mle(fit.x,B=10)
boot.x$par.star
boot.x$p.value
}
\keyword{models}
\keyword{optimize}
\keyword{robust}

