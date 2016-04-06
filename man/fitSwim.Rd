\name{fitSwim}
\alias{fitSwim}
\alias{SHMMM}
\title{Fit a switching movement model (swim)}
\description{Fits a switching movement model to an animal track. Currently the only model available is the Switching Hidden Markov Movement Model, the SHMMM, which is inappropriate for data with large measurement error.}

\usage{fitSwim(data, ts, model_type="SHMMM", regularize=TRUE)}

\arguments{
\item{data}{A data.frame with either two columns if regularize=FALSE, or three columns if regularize=TRUE. The columns must be specifically named. The column "date" should include the date-time step of the animal track as a POSIX class. The columns "lon" and "lat" denote the longitude and latitude, respectively.}
\item{ts}{The time step for interpolating the data, in hours. Fractions of hours are also possible.}
\item{model_type}{The swim model to be fit. Currently the only available option is the SHMMM, the Switching Hidden Markov Movement Model.}
\item{regularize}{Logical. Does your data need to be regularized in time? If so, then use regularize=TRUE. If your data already occurs in discrete-time, then use regularize=FALSE.}
}
\details{...}

\value{A list with four objects. The first, "regData" is the interpolated data used to fit the model. The second, "obj", is the TMB object which has been optimized, i.e. the negative log-likelihood of the fitted switching Markov movement model. The third, "parms", includes the maximum likelihood estimates of the parameters, and their standard errors. The fourth contains the latent behavioral states, estimated via the Viterbi algorithm.}

\references{...}

\author{K. Whoriskey}

\examples{
data(blueshark)
test = fitSwim(blueshark, 6)
test$parms
summary(test)
plot(test)
}

