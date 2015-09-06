\name{floodsNB}
\docType{data}
\alias{floodsNB}
\title{New Brunswick (Canada) Flood Dataset}
\description{
\code{floodsNB} is a \code{list} object containing the hydrometric stations considered for analysis... Each element from the \code{list} corresponds to an hydrometric station located in the Canadian province of New Brunswick, for which the flow is unregulated. For each station, the following information is available:
	\itemize{
	\item \code{data}: Maximum annual daily mean discharge (in \eqn{m^3/s});
	\item \code{peak}: Maximum annual daily peak discharge (in \eqn{m^3/s});
	\item \code{ln.drain}: Natural logarithm of the drainage area (in \eqn{km^2});
	\item \code{coor}: Coordinates (in latitude and longitude) of the hydrometric station;
	\item \code{status}: Station's status - \code{Active} or \code{Inactive};
	\item \code{Aucoin.2001}: Whether or not (\code{TRUE} or \code{FALSE}) the station is retained for analysis in ....
	}
}
\usage{data(floodsNB)}
\format{A \code{list} object whose elements correspond to distinct hydrometric stations.}
\source{HYDAT database.}
\references{
Environment Canada: \url{http://www.ec.gc.ca/rhc-wsc/default.asp?lang=En&n=894E91BE-1}

}
\keyword{datasets}
