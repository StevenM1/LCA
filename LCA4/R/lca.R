# This is a wrapper function (and tester) to call the simulateLCA function built in C.

LCA <- function(I, kappa, beta, Z, NDT, nTrials=1000, s=.1, dt=.001, maxT=5, nonLinear=TRUE, x0=NULL) {
  # Function to simulate the Leaky, Competing Accumulator (LCA, Usher & McClelland, 2001) model
  # Parameters #
  # I (double vector): Vector of input for every accumulator. E.g., c(1.2, 1, 1) simulates 3 accumulators, with inputs 1.2, 1, 1, respectively.
  # kappa (double): leakage
  # beta (double): inhibition
  # Z (double): threshold of accumulation
  # nTrials (integer): Number of trials to simulate.
  # s (double): Variance of Gaussian noise. Defaults to 0.1 (as in Miletic et al. 2017).
  # dt (double): temporal resolution of simulation. Defaults to 0.001, corresponding to a milisecond resolution.
  # maxT (double): maximum time to simulate. Defaults to 5, corresponding to 5 seconds of decision time.
  # nonLinear (bool): if TRUE, simulates the LCA with the non-linearity included in the original paper. If FALSE, the non-linearity is ommitted and negative accumulator values are allowed. Defaults to TRUE.
  # x0 (double vector): Vector of start point values for every accumulator. Defaults to NULL, which sets all start points to 0.

  # Check whether the NDT is larger than 1. If so, assume that the NDT is given in miliseconds, and divide by 1000 to convert to s.
  if(NDT > 1) {
  	warning(paste0("The non-decision time provided is larger than 1. I'll assume you meant ", NDT, " miliseconds"))
  	NDT <- NDT/1000
  }

  # Check whether x0 is NULL. If so, sets x0 to 0. If x0 is provided, check whether the length of x0 and I are the same.
  if(is.null(x0)) {
  	x0 <- rep(0, length(I))
  } else if(length(x0) != length(I)) {
  	stop('The number of accumulators in I is not the same as the number of accumulators in x0')
  }

  # Determine maximum number of time steps to simulate
  maxiter <- length(seq(0, maxT, by=dt))

  # Create empty numeric vectors for the responses and rts
  rts <- resps <- numeric(nTrials)

  # Call C-code
  out = .C('simulateLCA',
           nAcc = as.integer(length(I)), # Pass the number of accumulators to C. All R objects are passed as *pointers*, which makes it surprisingly difficult to obtain the length of the array in C itself.
           I = as.double(I),
           kappa = as.double(kappa),
           beta = as.double(beta),
           Z = as.double(Z),
           s = as.double(s),
           dt = as.double(dt),
           maxiter = as.integer(maxiter),
           resp = as.integer(resps),
           rt = as.double(rts),
           nTrials = as.integer(nTrials),
           nonLinear = as.integer(as.logical(nonLinear)),
           x0 = as.double(x0))

  # Structure output dataframe
  dat <- list(out$rt, out$resp)
  dat <- as.data.frame(dat)
  colnames(dat) <- c('rt','response')

  # Check which answers are correct
  dat$corr <- dat$response == 1

  # Add non-decision time to rt
  dat$rt = dat$rt+NDT

  # Round reaction times to 3 decimals
  dat$rt <- round(dat$rt, 3)
  return(dat)
}

testLCAsim <- function() {
	# Small function to test whether the LCA simulation works properly. Simulates two datasets and plots the data.

	# Simulate some data
	dat1 <- LCA(nTrials=1000, I=c(1.2, 1, 1), kappa=3, beta=3, Z=.2, NDT=.450, s=.1, dt=.001, maxT=5, nonLinear=TRUE, x0=c(.01, .02, .03))
	dat2 <- LCA(nTrials=1000, I=c(1.2, 1, 1), kappa=3, beta=3, Z=.2, NDT=.450, s=.1, dt=.001, maxT=5, nonLinear=FALSE, x0=c(0.01, .02, .03))

	print('Succesfully simulated two sets of data. Example data:')
	print(head(dat1))

	# Plot densities
	par(mfrow=c(1,2))
	plot(density(dat1$rt[dat1$corr == TRUE]), main='', xlab='nTrials = 1000')
	title('Simulated 3 accumulators with non-linearity', sub='I=c(1.2, 1, 1), k=3, b=3, Z=.2, NDT=.450')
	lines(density(dat1$rt[dat1$corr == FALSE]), col='red')

	plot(density(dat2$rt[dat2$corr == TRUE]), main='', xlab='nTrials = 1000')
	title('Simulated 3 accumulators without non-linearity', sub='I=c(1.2, 1, 1), k=3, b=3, Z=.2, NDT=.450')
	lines(density(dat2$rt[dat2$corr == FALSE]), col='red')
	print('Everything seems OK!')
}
