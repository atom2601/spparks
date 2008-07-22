/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   http://www.cs.sandia.gov/~sjplimp/spparks.html
   Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories

   Copyright (2008) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level SPPARKS directory.
------------------------------------------------------------------------- */

#ifndef APP_TEST_GROUP_H
#define APP_TEST_GROUP_H

#include "app.h"

namespace SPPARKS_NS {

class AppTestGroup : public App {
 public:
  AppTestGroup(class SPPARKS *, int, char **);
  ~AppTestGroup();
  void init();
  void input(char *, int, char **);
  void run(int, char **);

 private:
  class RandomPark *random;

  int ntimestep;
  double time;
  int nlimit;

  int nevents;               // # of events
  double *propensity;        // propensity of each event
  double pmax,pmin;          // maximum/minimun propensity value
  double tweak;              // percentage propensity tweak
  int seed;                  // random number seed

  double psum;
  int *count;

  int dep_graph;             // 1 if build/store dependency graph, else 0
  int ndep;                  // max number of dependencies from user
  int *ndepends;             // # of events that depend on each event
  int **depends;             // i,j = jth event that depends on ith event
  int *ran_dep;              // random deps for on-the-fly generation

  void iterate();
  void build_dependency_graph();
  double compute_propensity(int);

  void stats(char *);
  void stats_header(char *);
};

}

#endif
