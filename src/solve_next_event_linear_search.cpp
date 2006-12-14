/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   contact info, copyright info, etc
------------------------------------------------------------------------- */

#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "solve_next_event_linear_search.h"
#include "spk.h"
#include "random_park.h"
#include "error.h"

using namespace SPPARKS;

/* ---------------------------------------------------------------------- */

SolveNextEventLinearSearch::SolveNextEventLinearSearch(SPK *spk, int narg, char **arg) : Solve(spk, narg, arg)
{
  if (narg != 2) error->all("Illegal solve command");

  int seed = atoi(arg[1]);
  random = new RandomPark(seed);
}

/* ---------------------------------------------------------------------- */

SolveNextEventLinearSearch::~SolveNextEventLinearSearch()
{
  delete random;
}

/* ---------------------------------------------------------------------- */

void SolveNextEventLinearSearch::init(int n, double *propensity)
{
  delete [] prob;
  nevents = n;
  prob = new double[n];
  sum = 0;
  fprintf(screen,"propensity:   ");
  for (int i = 0; i < n; i++) {
    fprintf(screen,"%f ",propensity[i]);
    prob[i] = propensity[i];
    sum += propensity[i];
  }
  fprintf(screen,"\n");
}

/* ---------------------------------------------------------------------- */

void SolveNextEventLinearSearch::update(int n, int *indices, double *propensity)
{
  for (int i = 0; i < n; i++) {
    sum -= prob[indices[i]];
    prob[indices[i]] = propensity[indices[i]];
    sum +=  propensity[indices[i]];
  }
}

/* ---------------------------------------------------------------------- */

int SolveNextEventLinearSearch::event(double *pdt)
{
  int m;

  if (sum == 0.0) return -1;
  double fraction = sum * random->uniform();
  
  double partial = 0.0;
  for (m = 0; m < nevents; m++) {
    partial += prob[m];
    if (partial > fraction) break;
  }

  *pdt = -1.0/sum * log(random->uniform());
  return m;
}
