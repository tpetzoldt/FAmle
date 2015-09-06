\name{plot.metropolis}
\Rdversion{1.2}
\alias{plot.metropolis}
\title{A Function to Plot metropolis objects}

\description{This function allows to user to call different plots for visual assessment of the posterior distribution(s).}

\usage{\method{plot}{metropolis}(x, plot.type = "carlin", pos = 1:x$iter, ...)}
\arguments{
  \item{x}{\code{\link{mle}} object corresponding to the fitted model.}
  \item{plot.type}{The user may choose betweew:
	\cr \code{carlin} returns the same plot as in Carlin and Louis (2009) (see References);
	\cr \code{ts} returns \code{\link{plot.ts}};
	\cr \code{pairs} returns a \code{\link{pairs}};
	\cr \code{hist} returns an \code{\link{hist}} for each marginal posterior distribution;
	\cr \code{post.pred} returns an histogram of the data's posterior predictive distribution.
}
  \item{pos}{May be used by the user to plot a subset (i.e. a random subset, \code{sample})) of the posterior distribution when \code{pairs} is called. This avoids using too much memory while building the plot.}
  \item{\dots}{Additional arguments pertaining to function \code{\link{plot.default}}.}
}



\references{See list of references for \code{\link{metropolis}}.}



\seealso{\code{\link{metropolis}}}
\examples{
data(yarns)
x <- yarns$x
fit.x <- mle(x,'gamma',c(.1,.1))
bayes.x <- metropolis(model=fit.x,iter=100,
	trans.list=list(function(x) exp(x),function(x) exp(x)))
plot(bayes.x)
plot(bayes.x,'hist',col='cyan')
plot(bayes.x,'pairs',cex=.1,pch=19)
plot(bayes.x,'pairs',pos=sample(1:bayes.x$iter,20),col='red')
plot(bayes.x,'post.pred',col='green')
}

\keyword{aplot}
\keyword{models}

