\name{distr}
\Rdversion{1.2}
\alias{distr}
\title{Distribution functions 4-in-1}
\description{This function can be used to call any of the 4 functions specific to a given probability distribution available in \R.}
\usage{distr(x, dist, param, type = "d", model = NULL, ...)}
\arguments{
  \item{x}{Vector (or array) of quantiles, vector (or array) of probabilities, or number of observations.}
  \item{dist}{Distribution name.}
  \item{param}{Vector (or array) of parameters.}
  \item{type}{Type of function to be called (\code{'d'}, \code{'p'}, \code{'q'}, or \code{'r'}).}
  \item{model}{Object from the class \code{\link{mle}} - may be specified instead of \code{param} and \code{dist}.}
  \item{\dots}{Additional arguments \code{log}, \code{lower.tail}, and \code{log.p}, depending on \code{type}.}
}
\details{For each distribution available in \R, 4 functions can be called. For example, for the normal distribution, the following 4 functions are available: \code{\link{dnorm}}, \code{\link{pnorm}}, \code{\link{qnorm}}, and \code{\link{rnorm}}. For the normal distribution, based on the argument \code{type}, \code{distr} may be used to call any one of the previous four functions.}
\value{Returns the density, the distribution function, the quantile function, or random variates.}

\note{Most functions in \code{\link{FAmle}} rely upon \code{distr}.}

\examples{
## Example 1
dnorm(-4:4,0,1,log=TRUE)
distr(-4:4,'norm',c(0,1),type='d',log=TRUE)

## Example 2
mu.vec <- c(1,100,100)
sigma.vec <- c(1,11,111)
n <- 3
set.seed(123)
rnorm(n,mu.vec,sigma.vec)
set.seed(123)
distr(n,'norm',cbind(mu.vec,sigma.vec),'r')

## Example 3
qnorm(.9,mu.vec,sigma.vec)
distr(.9,'norm',cbind(mu.vec,sigma.vec),'q')

}

\keyword{distribution}

