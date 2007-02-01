/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   contact info, copyright info, etc
------------------------------------------------------------------------- */

#ifndef SOLVE_NEXT_EVENT_GROUP_SEARCH_H
#define SOLVE_NEXT_EVENT_GROUP_SEARCH_H

#include "solve.h"

namespace SPPARKS {

class SolveNextEventGroupSearch : public Solve {
 public:
  SolveNextEventGroupSearch(class SPK *, int, char **);
  ~SolveNextEventGroupSearch();
  void input(int, char **) {}
  void init(int, double *);
  void update(int, int *, double *);
  void update(int, double *);
  void resize(int, double *);
  int event(double *);

 private:
  class RandomPark *random;
  class Groups *groups;
  int nevents;
  double *p;
  double last_size;
  double sum;
  int seed;
  double lo, hi;
};

}

#endif
