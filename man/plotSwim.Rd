\name{plotSwim}
\alias{plotSwim}
\title{Plot a swim object} 
\description{Plots a swim object (the results of fitting the SHMMM) as a spatial object with surrounding land boundaries.}

\usage{plotSwim(object)}

\arguments{
\item{object}{An object of class 'swim'.}
}

\author{K. Whoriskey}

\examples{
data(blueshark)
fit = fitSwim(blueshark, 6)
plotSwim(fit)
}




