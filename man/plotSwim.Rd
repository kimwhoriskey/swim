\name{plotSwim}
\alias{plotSwim}
\title{Plot a swim object} 
\description{Plots a swim object (the results of fitting the HMMM) as a spatial object with surrounding land boundaries. }

\usage{plotSwim(object)}

\arguments{
\item{object}{An object of class 'swim'.}
}

\author{K. Whoriskey}

\references{
Vector land data were taken from Natural Earth.
}

\examples{
data(blueshark)
fit = fitSwim(blueshark, 6)
plotSwim(fit)
}




