/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   contact info, copyright info, etc
------------------------------------------------------------------------- */

#ifndef DIAG_CLUSTER_H
#define DIAG_CLUSTER_H

#include "stdio.h"
#include <stack>
#include "cluster.h"
#include "diag.h"

namespace SPPARKS {

class DiagCluster : public Diag {
  friend class SweepLattice;

 public:
  DiagCluster(class SPK *, int, char **);
  virtual ~DiagCluster();

  void init(double);
  void compute(double, int);
  void stats(char *);
  void stats_header(char *);

 protected:

  // Functions and Data for Cluster Analysis
  void analyze_clusters(double);
  void write_header();
  void dump_clusters(double);
  void dump_clusters_detailed();
  void generate_clusters();
  void add_cluster(int, int, int, int*);
  void free_clustlist();

  int* cluster_ids;
  int ncluster;
  Cluster* clustlist;
  std::stack<int> cluststack;      // stack for performing cluster analysis

  FILE *fp, *fpdump;
  class AppLattice *applattice;

  class CommLattice *comm;

  double boxxlo,boxxhi,boxylo,boxyhi,boxzlo,boxzhi;    // simulation box bounds

  int *id;                     // global ID (1-N) of site
  double **xyz;                // coords of site

  int nglobal;                 // global # of sites
  int nlocal;                  // # of sites I own
  int nghost;                  // # of ghost sites I store

  enum DumpStyles {STANDARD};
  int dump_style;
  int idump;
};

}

#endif
