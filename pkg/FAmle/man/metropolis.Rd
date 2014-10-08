\name{metropolis}
\Rdversion{1.2}
\alias{metropolis}
\title{Bayesian Estimation of Univariate Probability Distributions}
\description{For a given dataset, this function serves to approximate (using a Metropolis algorithm) the posterior distribution of the parameters  for some specified parametric probability distribution.}
\usage{
metropolis(model, iter = 1000, tun = 2, trans.list = NULL,
	start = NULL, variance = NULL, prior = NULL, burn = 0,
	uniroot.interval = c(-100, 100),pass.down.to.C=FALSE)
}
\arguments{
  \item{model}{\code{\link{mle}} object corresponding to the fitted (by maximum likelihood) model. A \code{list(x=dataset,dist=distribution)} object may also be provided, but the user will then have to make sure to specify the arguments \code{start} and \code{variance}. Moreover, the latter two arguments will have to be specified on their transformed scales (see \code{trans.list}).}
  \item{iter}{The requested number of iterations - the Markov Chain's length.}
  \item{tun}{A tuning constant; value by which the covariance matrix of the multivariate normal proposal will be multiplied - see References.}
  \item{trans.list}{A \code{\link{list}} object containing a function for each parameter that is to be estimated. For each parameter, the function must correspond to the inverse transformation that will determine the parametrization for which the simulation will be carried out (see Example and Details).}
  \item{start}{A vector of starting values for the algorithm. If \code{NULL}, the maximum likelihood parameter estimates will be used as starting values for the Markov Chain. If \code{model} is not an object from the class \code{mle}, this argument will have to be specified, along with the argument \code{variance}. Moreover, as already stated above, the user will have to make sure that both \code{start} and \code{variance} are those for the transformed parameters (see \code{trans.list}).}
  \item{variance}{Covariance matrix of the multivariate normal proposal distribution. If \code{NULL}, the observed Fisher's information will be used and multiplied by the specified \code{tun}. As for \code{start}, this argument needs to be specified if \code{model} is not from the class \code{mle}.}
  \item{prior}{A function that corresponds to the joint prior distribution (see Example). Note that the prior distribution will be evaluated on the transformed parameter space(s).}
  \item{burn}{Burn-in period (see References).}
  \item{uniroot.interval}{Default is c(-100,100). This interval is used by \code{R}'s function \code{\link{uniroot}} to search for the inverse of each element in \code{trans.list}.}
  \item{pass.down.to.C}{If \code{TRUE}, the iterative task is passed down to a C program for faster implementation of the MCMC algorithm.}
}
\details{This function uses a single block Metropolis algorithm with multivariate normal proposal. For this function to work properly, all parameters should be defined on the real line - parameter transformation(s) might be required. If \code{trans.list} is not specified, the function will assume that the parameter distributions are all defined on the real line (i.e., \code{function(x) x} will be used for each parameter). If no prior distribution is provided, an improper prior distribution - uniform on the interval )-Inf,+Inf( - will be used for all parameters (i.e., prior distribution proportional to 1 - function(x) 1).

In order to minimize the number of arguments for \code{metropolis}, the function automatically computes the inverse of \code{trans.list}: this suppresses the need for the user to provide both the "inverse transformation" and the "transformation". However, problems may occur, and it is why the user is allowed to alter \code{uniroot.interval}. Depending on the number of errors reported, future versions of this package may end up requesting that a list for both the "inverse transformation" and the "transformation" be provided by the user.

A nice list of references is provided below for more information on topics such as: MCMC algorithms, tuning of Metropolis-Hastings algorithms, MCMC convergence diagnostics, the Bayesian paradigm...

}
\value{
  \item{rate}{MCMC acceptance rate. This value is computed before applying the burn-in; i.e., it is computed for \code{sims.all}.}
  \item{total.time}{Total computation time.}
  \item{sims.all}{Array containing all iterations.}
  \item{sims}{Array containing iterations after burn-in.}
  \item{input}{Inputted \code{\link{mle}} object.}
  \item{iter}{Number of iterations.}
  \item{prior}{Prior distribution.}
  \item{burn}{Integer corresponding to the number of iterations to be discarded - burn-in period.}
  \item{M}{Parameter vector whose elements correspond to the parameter values (on the scales specified by \code{trans.list}) obtained at the last iteration of the Metropolis sampler; i.e. \code{sims[iter,]}.}
  \item{V}{Covariance matrix computed using, after removing the burn-in period, the joint posterior distribution of the parameters (on the scales specified by \code{trans.list}). This matrix might be used to tune the MCMC algorithm.}
}
\references{
Gelman, A., Carlin, J.B., Stern, H.S., and Rubin, D.B. (2004). Bayesian data analysis, 2nd edition, Chapman & Hall/CRC.

Carlin, B.P, and Louis, T.A. (2009). Bayesian methods for data analysis. Chapman & Hall/CRC.

Gamerman, D., and Lopes H.F. (2006). Markov Chain Monte Carlo: Stochastic simulation for Bayesian inference. 2nd edition, Chapman & Hall/CRC.

Gilks, W.R., Richardson, S., and Spiegelhalter, D.J. (1996). Markov Chain Monte Carlo in Practice. Chapman & Hall.
}


\seealso{\code{\link{plot.metropolis}}, \code{\link{mle}}}
\examples{

### These examples should be re-run with, e.g., iter > 2000.
data(yarns)
x <- yarns$x
fit.x <- mle(x,'gamma',c(.1,.1))
bayes.x.no.prior <- metropolis(model=fit.x,iter=150,
	trans.list=list(function(x) x,function(x) exp(x)))
plot(bayes.x.no.prior)

# examples of prior distributions (note that these prior distribution
#	are specified for the transformated parameters;
#	i.e., in this case, 'meanlog' -> 'meanlog' and 'sdlog' -> 'ln.sdlog')
# for the scale parameter only
prior.1 <- function(x) dnorm(x[2],.8,.1)
# for both parameters (joint but independent in this case)
prior.2 <- function(x) dunif(x[1],3.4,3.6)*dnorm(x[2],1,1)

bayes.x.prior.2 <- metropolis(model=fit.x,iter=150,
	trans.list=list(function(x) x,function(x) exp(x)),prior=prior.2)
plot(bayes.x.prior.2)

# Example where 'model' is not from the class 'mle'; i.e. both 'start' and 'variance' need to be specified!
#x <- rweibull(5,2,1)
x <- c(0.9303492,1.0894917,0.9628029,0.6145032,0.4756699)
# Here 'fit.x <- mle(x,'weibull',c(.1,.1))' is not used,
model.x <- list(x=x,dist='weibull')
# and an informative prior distribution is considered to ensure a proper posterior distribution
prior.x <- function(x) dnorm(x[1],log(2),.1)*dnorm(x[2],log(1),.1)
trans.list.x <- list(function(x) exp(x), function(x) exp(x))
bayes.x <- metropolis(model=model.x,iter=150,prior=prior.x,trans.list=trans.list.x,pass.down.to.C=TRUE,start=c(0,0),variance=diag(.1,2,2))

}

\keyword{models}
