\name{genTrack}
\alias{genTrack}
\title{Simulate animal track with switching behaviors}
\description{Simulates an animal path according to the code of the Switching Hidden Markov Movement Model. Does not incorporate measurement error.}

\usage{genTrack(T=500, theta=c(0,pi), gamma=c(0.7,0.2), 
                Sigma=matrix(c(0.2*0.2,0,0,0.1*0.1),nrow=2), 
                alphas=c(0.8,0.9))}

\arguments{
\item{T}{The length of the animal track.}
\item{theta}{Vector of two turning angles.}
\item{gamma}{Vector of two autocorrelation parameters. Must be between 0 and 1.}
\item{Sigma}{Covariance matrix for the Gaussian process error distribution.}
\item{alphas}{Vector of two switching probabilities: the first is the probability of being in state 1 and staying in state 1, and the second is the probability of being in state 2 and staying in state 2, i.e. alpha_(1,1) and alpha_(2,2) respectively. }
}

\details{...}

\value{A data.frame with three columns: "lon" is the simulated longitude, "lat" is the simulated latitude, and "b" is the simulated behavioral state.}

\references{...}

\author{K. Whoriskey}

\examples{

track = genTrack(500)
plot(track$lat~track$lon, col=ifelse(track$b==1,"red", "black"))

}

