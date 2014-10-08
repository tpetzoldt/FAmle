\name{FAmle-internal}
\Rdversion{1.2}
\alias{FAmle-internal}
\alias{cdf.plot}
\alias{delta.Q}
\alias{delta.QQ}
\alias{Diff.1}
\alias{Diff.2}
\alias{Diff.3}
\alias{hist.metropolis}
\alias{hist.plot}
\alias{Plot.post.pred}
\alias{post.pred}
\alias{Quantile.plot}
\alias{Return.plot}
\alias{Carlin}
\title{Internal Functions in the FAmle Package}
\description{Internal functions in the FAmle package
.}
\usage{
cdf.plot(z)
delta.Q(p, model, ln = FALSE)
delta.QQ(model, alpha = 0.1, ln = FALSE)
Diff.1(x, f, h = 1e-04)
Diff.2(k, i, model, p, ln = FALSE)
Diff.3(i, model, p, ln = FALSE)
\method{hist}{metropolis}(x, density = TRUE, ...)
\method{hist}{plot}(x,...)
Plot.post.pred(x, ...)
post.pred(z, fun = NULL)
Quantile.plot(z, ci = FALSE, alpha = 0.05)
Return.plot(model, ci = FALSE, alpha = 0.05)
Carlin(x)
}
\arguments{
 \item{z}{A \code{\link{mle}} object.}
 \item{p}{A vector of probabilities.}
 \item{model}{A \code{\link{mle}} object.}
 \item{ln}{Whether or not (TRUE or FALSE) computations should be carried out on the natural logarithmic scale.}
 \item{alpha}{The significance level.}
 \item{x}{Value at which the numerical derivative should be evaluated. For the \code{Carlin} function (see References for \code{\link{metropolis}}), this \code{x} corresponds to an object from the class \code{\link{metropolis}}.}
 \item{f}{A function to be differentiated.}
 \item{h}{Small number representing a small change in \code{x}.}
 \item{k}{Parameter value at which the first derivative should be evaluated.}
 \item{i}{Position of the parameter, within a vector of parameters, with respect to which differentiation should be carried out.}
 \item{density}{Whether or not (TRUE or FALSE) a Kernel density should be added to the histogram - see \code{\link{density}}.}
 \item{\dots}{Additional arguments pertaining to \code{\link{hist}}.}
 \item{fun}{optional argument that may be used to modify the scale on which the histogram will be plotted.}
 \item{ci}{Whether or not (TRUE or FALSE) approximated \emph{100*(1-\code{alpha})} confidence intervals should be added to the plot (either \code{Quantile.plot} or \code{Return.plot}).}
}

\keyword{package}
