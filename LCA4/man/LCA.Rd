\name{LCA}
\alias{LCA}
\title{Simulate the Leaky, Competing Accumulator model of perceptual decision-making}
\usage{
LCA(I, kappa, beta, Z, NDT, nTrials=10000, s=0.1, dt=0.001, maxT=5, nonLinear=TRUE, x0=NULL)
}
\description{
Simulates the Leaky, Competing Accumulator model of perceptual decision-making (Usher & McClelland, 2001). Simulation code is written in C for performance.

}
\arguments{
	\item{I}{Vector of inputs of each accumulator}
	\item{kappa}{Leakge/decay parameter}
	\item{beta}{Inhibition parameter (lateral inhibition between accumulators)}
	\item{Z}{Threshold}
	\item{NDT}{Non-decision time in seconds. If a value larger than 1 is provided, the code assumes it to be miliseconds and converts to seconds.}
  \item{nTrials}{Number of trials to simulate. Defaults to 10,000.}
	\item{s}{Variance of Gaussian noise. Defaults to 0.1 (as was used in Miletić et al., 2017)}
	\item{dt}{Temporal resolution of simulation in seconds. Defaults to 0.001, corresponding to a milisecond resolution.}
	\item{maxT}{Maximum decision time to simulate in seconds. Defaults to 5, corresponding to a 5 second decision time. Especially when simulating a large number of trials, set maxT sufficiently high to capture the right tail of the distribution well.}
	\item{nonLinear}{If TRUE, simulates the LCA including the non-linearity that was included in the original formulation (Usher & McClelland, 2001). If FALSE, the non-linearity is ommitted, and accumulators are allowed to obtain values below 0. Defaults to TRUE.}
  \item{x0}{Start points. In the original paper, all start points were assumed to be 0, but it's easy to simulate with different start points. Defaults to NULL, which sets all start points to 0.}
}

\references{
Miletić, S., Turner, B., Forstmann, B. U., and Van Maanen, L. (2017) Parameter recovery for the leaky, competing accumulator model. \emph{Journal of Mathematical Psychology}, \bold{76}, 25--50.

Usher, M., and McClelland, J. L. (2001) The Time Course of Perceptual Choice. \emph{Psychological Review}, \bold{108(3)}, 550--592.
}

\examples{
# Simulate data
dat <- LCA(nTrials=1000, I=c(1.2, 1, 1), kappa=3, beta=3, Z=.12, NDT=.450, s=.1, dt=.001, maxT=5, nonLinear=TRUE, x0=c(.01, .02, .03))
print(head(dat))

# Plot densities
d1 <- density(dat$rt[dat$corr == TRUE])
d2 <- density(dat$rt[dat$corr == FALSE])
with(d1, plot(x=x, y=y, type='l', main='Density estimates of simulated data'))
with(d2, lines(x=x, y=y, col='red'))
}
