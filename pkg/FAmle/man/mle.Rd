\name{mle}
\Rdversion{1.2}
\alias{mle}
\title{Maximum Likelihood Estimation of Univariate Probability Distributions}
\description{For a given dataset, this function serves to find maximum likelihood parameter estimates for some specified parametric probability distribution.}
\usage{mle(x, dist, start = NULL, method = "Nelder-Mead")}
\arguments{
  \item{x}{A univariate dataset (a vector).}
  \item{dist}{Distribution to be fitted to \code{x}.}
  \item{start}{Starting parameter values for the optimization algorithm (see \code{\link{optim}}).}
  \item{method}{The optimization method to be used (see \code{\link{optim}}).}
}

\value{
  \item{fit}{\code{optim} output (see \code{\link{optim}}).}
  \item{x.info}{Array that contains the following columns: 
	\cr \code{i}: (\code{1:length(x)}), 
	\cr \code{x}: (original dataset),
	\cr \code{z}: (sorted dataset),
	\cr \code{Fx}: (CDF of \code{x} evaluated at the estimated parameter values),
	\cr \code{Fz}: (sorted values of \code{Fx}),
	\cr \code{Emp}: (\code{i/(length(x)+1)}),
	\cr \code{zF}: (\code{distr(Emp,'dist',par.hat,'q')} evaluated at estimated parameter values (\code{par.hat})),
	\cr \code{fx}: (PDF of \code{x} evaluated at the estimated parameter values),
	\cr \code{fz}: (PDF of \code{z} evaluated at the estimated parameter values)
	}
  \item{dist}{Distribution fitted to \code{x}.}
  \item{par.hat}{Vector of estiamted parameters.}
  \item{cov.hat}{Observed Fisher's information matrix.}
  \item{k}{Number of parameters}
  \item{n}{Number of observations (i.e., \code{length(x)}).}
  \item{log.like}{Log-likelihood value evaluated at the estimated parameter (i.e. \code{par.hat}).}
  \item{aic}{Akaike information criterion computed as \code{2*k - 2*log.like}.}
  \item{ad}{Anderson Darling statistic evaluated at the estimated parameter values.}
  \item{data.name}{Name for \code{x}.}
  \item{rho}{Pearson's correlation coefficient computed as \code{cor(x.info[,'z'],x.info[,'zF'])}.}
}

\seealso{\code{\link{optim}}, \code{\link{distr}}, \code{\link{boot.mle}}, \code{\link{metropolis}}, \code{\link{Q.conf.int}}}

\examples{
data(yarns)
x <- yarns$x
fit.x <- mle(x,'weibull',c(.1,.1))
fit.x
names(fit.x)
#plot(fit.x)
#plot(fit.x,TRUE,alpha=.01)
p <- c(.9,.95,.99)
distr(p,model=fit.x,type='q')
Q.conf.int(p,fit.x,.01)
Q.conf.int(p,fit.x,.01,TRUE)
}
\keyword{optimize}
\keyword{models}

