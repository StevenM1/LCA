#include <math.h>
#include <Rmath.h>
#include <stdlib.h>
#include <stdio.h>
#include <R.h>
#include <stdbool.h>

void simulateLCA(int *nAcc,
                 double *I,
                 double *kappa,
                 double *beta,
                 double *Z,
                 double *s,
                 double *dt,
                 int *maxiter,
                 int *resp,
                 double *rt,
                 int *nTrials,
                 int *nonLinear,
                 double *x0)
    /* This function simulates the Leaky, Competing Accumulator (LCA) model (Usher & McClelland, 2001)
    The following arguments are required:
    1. nAcc = Number of accumulators (integer scalar);
    2. I = Input of accumulators (double vector)
    3. kappa = leakage (double scalar)
    4. beta = inhibition (double scalar)
    5. Z = threshold (double scalar)
    6. s = noise variance (double scalar)
    7. dt = step size (0.001 = 1 ms) (double scalar)
    8. maxiter = maximum iterations (double scalar). You can determine a good number for maxiter using length(seq(0, maxT.seconds, by=dt)).
    9. resp = responses vector (integer vector). This vector of length nTrials should be initialized empty.
    10. rt = reaction times vector (double vector). This vector of length nTrials should be initialized empty.
    11. nTrials = number of trials to simulate (integer scalar).
    12. nonLinear = Should the LCA non-linearity be included?
    13. x0 = Vector of start points (double vector).
    */
{
  // initialize some variables
  double rhs, sum;
  int iter, i;
  bool winner;

  // Get random number generator state from R
  GetRNGstate();

  // Get factor to multiply random noise (drawn from N(0, 1)) with. *s = VARIANCE of white noise.
  rhs = sqrt((*dt)*(*s)*(*s));

  // iterate i over nTrials
  for(i = 0; i < *nTrials; i++) {
    winner = false;            // start out without a winner
    iter = 0;                  // time step counter
    resp[i] = (double) -1.0;   // Start trial with resp = -1 (no response given)
    double x[*nAcc];           // initialize accumulator array
    double dummy[*nAcc];       // initialize dummy array (later used to determine inhibition)

    for(int z = 0; z < *nAcc; z++) {
      // Set both accumulator array x and dummy array to 0
      x[z] = x0[z];
      dummy[z] = x0[z];

      // The following statements allow for some debugging (check whether the start point and input are correct)
      // if(i == 0) {
      //  printf("Start point for accumulator %d: %.3f\n", z+1, x[z]);
      //  printf("Input for accumulator %d: %.3f\n", z+1, I[z]);
      // }
    }

    // Within a trial, do a while loop...
    do
    {
      // First, copy current accumulator values into dummy
      for(unsigned int z = 0; z < *nAcc; z++) {
        dummy[z] = x[z];
        dummy[z] = dummy[z]*(*dt)*(*beta);  // Scale by beta to get inhibition!
      }

      // Sum over dummy to get full inhibition of system. Sum = summed inhibition.
      sum = 0;
      for(int a = 0; a < *nAcc; a++)  //
      {
        sum = sum + dummy[a];
      }

      // LCA equations: To each accumulator, add the influences of input, leak, inhibition, and noise.
      for(unsigned int z = 0; z < *nAcc; z++) {
        x[z] = x[z] + (*dt)*(I[z]) - (*kappa)*(x[z])*(*dt) - sum + dummy[z] + rhs*norm_rand(); // A little trick to code easier: subtract all inhibition in the system from each accumulator, and then add it's own contribution to that inhibition back.
      }

      // Add one to the time step counter
      iter = iter+1;

      // Check if there is a winner; ie an accumulator with a value >= threshold
      for(int z = 0; z < *nAcc; z++) {
        if(x[z] >= *Z) {
          resp[i] = z+1;
          winner=true;
        }

        // If simulating a non-linear LCA, check for accumulators with values lower than 0. If these exist, reset to 0.
        if(*nonLinear == 1) {
          if(x[z] < 0 ) {
            x[z] = 0;
          }
        }
      }
    } while(iter < *maxiter && winner==false);  // Stop the while loop if there's a winner, or if time ran out.

    // After while loop is done, set reaction time of current trial. Non-responses are coded as the maximum RT with a response -1.
    rt[i] = ((double) iter)*(*dt) - (*dt)/((double) 2.0);

  }
  PutRNGstate();
}
