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

#include "mpi.h"
#include "math.h"
#include "stdlib.h"
#include "string.h"
#include "output.h"
#include "app.h"
#include "app_lattice.h"
#include "timer.h"
#include "memory.h"
#include "error.h"

using namespace SPPARKS_NS;

enum{GENERAL,LATTICE};     // same as in app.h
enum{INT,DOUBLE};

#define MAXLINE 1024
#define DEFAULT "id lattice x y z"

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

/* ---------------------------------------------------------------------- */

Output::Output(SPPARKS *spk) : Pointers(spk)
{
  MPI_Comm_rank(world,&me);
  MPI_Comm_size(world,&nprocs);

  ndiags = 0;
  diaglist = NULL;
  dump_delta = 0.0;
  idump = 0;
  dump_ilogfreq = 0;
  dump_eps = 1.0e-6;
  stats_delta = 0.0;
  stats_ilogfreq = 0;
  stats_eps = 1.0e-6;

  size_one = 0;
  fp = NULL;
  vtype = NULL;
  vindex = NULL;
  vformat = NULL;
  pack_choice = NULL;
  buf = NULL;
}

/* ---------------------------------------------------------------------- */

Output::~Output()
{
  for (int i = 0; i < ndiags; i++) delete diaglist[i];
  memory->sfree(diaglist);

  if (me == 0 && fp) fclose(fp);

  delete [] vtype;
  delete [] vindex;
  for (int i = 0; i < size_one; i++) delete [] vformat[i];
  delete [] vformat;
  delete [] pack_choice;

  memory->sfree(buf);
}

/* ---------------------------------------------------------------------- */

void Output::init(double time)
{
  // test if dump is defined and propensity is output and doesn't exist

  if (dump_delta > 0.0) {
    int flag = 0;
    for (int i = 0; i < size_one; i++)
      if (pack_choice[i] == &Output::pack_propensity) flag = 1;
    if (flag && !solve)
      error->all("Dumping propensity but no KMC solve performed");
  }

  // setup future stats and dump calls

  if (dump_ilogfreq == 0) {
    dump_time = time + dump_delta;
  } else if (dump_ilogfreq == 1) {
    dump_time = time + dump_delta;
    dump_t0 = time;
    dump_irepeat = 0;
  }

  if (stats_ilogfreq == 0) {
    stats_time = time + stats_delta;
  } else if (stats_ilogfreq == 1) {
    stats_time = time + stats_delta;
    stats_t0 = time;
    stats_irepeat = 0;
  }

  // initialize all diagnostics

  for (int i = 0; i < ndiags; i++)
    diaglist[i]->init(time);

  // print dump file 1st snapshot

  if (dump_delta > 0.0) dump();

  // print stats header and initial stats
  
  stats_header();
  stats(0);
}

/* ---------------------------------------------------------------------- */

void Output::set_stats(int narg, char **arg)
{
  if (narg < 1) error->all("Illegal stats command");
  stats_delta = atof(arg[0]);

  int iarg = 1;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"logfreq") == 0) {
      stats_ilogfreq = 1;
      iarg++;
      if (iarg+1 < narg) {
	stats_nrepeat = atoi(arg[iarg]);
	if (stats_nrepeat < 1) error->all("Illegal stats command");
	iarg++;
	stats_scale = atof(arg[iarg]);
	if (stats_scale <= 1.0) error->all("Illegal stats command");
      } else {
	error->all("Illegal stats command");
      }
    }
    iarg++;
  }
}

/* ---------------------------------------------------------------------- */

void Output::set_dump(int narg, char **arg)
{
  // determine correct kind of app pointer

  if (app->appclass == LATTICE)
    applattice = (AppLattice *) app;
  else error->all("Cannot use dump with off-lattice app");

  // parse dump args

  if (narg < 2) error->all("Illegal dump command");
  dump_delta = atof(arg[0]);
  if (dump_delta <= 0.0) error->all("Illegal dump command");

  if (me == 0) {
    if (fp) fclose(fp);
    fp = fopen(arg[1],"w");
    if (!fp) error->one("Cannot open dump file");
  }

  int iarg = 2;
  if (iarg < narg) {
    if (strcmp(arg[iarg],"logfreq") == 0) {
      dump_ilogfreq = 1;
      iarg++;
      if (iarg+1 < narg) {
	dump_nrepeat = atoi(arg[iarg]);
	if (dump_nrepeat < 1) error->all("Illegal dump command");
	iarg++;
	dump_scale = atof(arg[iarg]);
	if (dump_scale <= 1.0) error->all("Illegal dump command");
	iarg++;
      } else {
	error->all("Illegal dump command");
      }
    }
  }

  // line = one string of concatenated keywords
  // size_one = # of keywords

  char *line = new char[MAXLINE];
  if (iarg == narg) {
    size_one = 5;
    strcpy(line,DEFAULT);
  } else {
    size_one = narg - iarg;
    line[0] = '\0';
    for (int jarg = iarg; jarg < narg; jarg++) {
      strcat(line,arg[jarg]);
      strcat(line," ");
    }
    line[strlen(line)-1] = '\0';
  }

  // allocate per-keyword memory

  vtype = new int[size_one];
  vindex = new int[size_one];
  vformat = new char*[size_one];
  pack_choice = new FnPtrPack[size_one];

  int i = 0;
  char *word = strtok(line," \0");

  while (word) {
    if (strcmp(word,"id") == 0) {
      pack_choice[i] = &Output::pack_id;
      vtype[i] = INT;
    } else if (strcmp(word,"lattice") == 0) {
      pack_choice[i] = &Output::pack_lattice;
      vtype[i] = INT;
      if (app->appclass == LATTICE && applattice->lattice == NULL)
	error->all("Dumping lattice but application does not support it");
    } else if (strcmp(word,"x") == 0) {
      pack_choice[i] = &Output::pack_x;
      vtype[i] = DOUBLE;
    } else if (strcmp(word,"y") == 0) {
      pack_choice[i] = &Output::pack_y;
      vtype[i] = DOUBLE;
    } else if (strcmp(word,"z") == 0) {
      pack_choice[i] = &Output::pack_z;
      vtype[i] = DOUBLE;
    } else if (strcmp(word,"energy") == 0) {
      pack_choice[i] = &Output::pack_energy;
      vtype[i] = DOUBLE;
    } else if (strcmp(word,"propensity") == 0) {
      pack_choice[i] = &Output::pack_propensity;
      vtype[i] = DOUBLE;

    // integer value = iN
    // double value = dN

    } else if (word[0] == 'i') {
      pack_choice[i] = &Output::pack_integer;
      vtype[i] = INT;
      vindex[i] = atoi(&word[1]);
      if (app->appclass != LATTICE)
	error->all("Invalid keyword in dump command");
      if (vindex[i] < 1 || vindex[i] > applattice->ninteger)
	error->all("Invalid keyword in dump command");
      vindex[i]--;
    } else if (word[0] == 'd') {
      pack_choice[i] = &Output::pack_double;
      vtype[i] = DOUBLE;
      vindex[i] = atoi(&word[1]) - 1;
      if (app->appclass != LATTICE)
	error->all("Invalid keyword in dump command");
      if (vindex[i] < 1 || vindex[i] > applattice->ndouble)
	error->all("Invalid keyword in dump command");
      vindex[i]--;

    } else error->all("Invalid keyword in dump command");

    word = strtok(NULL," \0");
    i++;
  }

  delete [] line;

  // setup vformat strings, one per field

  for (i = 0; i < size_one; i++) {
    char *format;
    if (vtype[i] == INT) format = "%d ";
    else if (vtype[i] == DOUBLE) format = "%g ";
    int n = strlen(format) + 1;
    vformat[i] = new char[n];
    strcpy(vformat[i],format);
  }

  // setup dump params and buf

  nglobal = applattice->nglobal;
  nlocal = applattice->nlocal;
  boxxlo = applattice->boxxlo;
  boxxhi = applattice->boxxhi;
  boxylo = applattice->boxylo;
  boxyhi = applattice->boxyhi;
  boxzlo = applattice->boxzlo;
  boxzhi = applattice->boxzhi;

  int nbuf = nlocal*size_one;
  MPI_Allreduce(&nbuf,&maxbuf,1,MPI_INT,MPI_MAX,world);
  buf = (double *) memory->smalloc(maxbuf*sizeof(double),"output:buf");
}

/* ---------------------------------------------------------------------- */

void Output::compute(double time, int done)
{
  int iflag;
  int ntmp;
  double tgoal;

  // check if dump output required

  if (dump_delta > 0.0 && time > dump_time-dump_eps ) {
    dump();
  // calculate new dump time
  // ensure new dump_time exceeds time
    if (dump_ilogfreq == 0) {
      dump_time += dump_delta;
      if (time > dump_time-dump_eps)
	dump_time = ceil(time/dump_delta)*dump_delta;
    } else if (dump_ilogfreq == 1) {
      dump_time += dump_delta;
      dump_irepeat++;
      if (dump_irepeat == dump_nrepeat || time > dump_time-dump_eps) {


	// Calculate next smallest delta that will 
	// reach tgoal within nrepeat steps
	tgoal = time-dump_t0+dump_delta;
	ntmp = MAX(1,static_cast<int>
		   (ceil(log(tgoal/(dump_delta*dump_nrepeat))
			 /log(dump_scale))));
	dump_delta *= pow(dump_scale,ntmp);
	dump_time = ceil(tgoal/dump_delta)*dump_delta;
	dump_irepeat = 0;
      }
    }
  }

  // check if stats output required

  iflag = 0;

  if (stats_delta > 0.0 && time > stats_time-stats_eps) iflag = 1;

  // perform diagnostics

  for (int i = 0; i < ndiags; i++)
    diaglist[i]->compute(time,iflag,done);

  // perform stats (after diagnostics)

  if (iflag || done) stats(1);

  // calculate new stats time
  // ensure new stats_time exceeds time

  if (iflag) {
    if (stats_ilogfreq == 0) {
      stats_time += stats_delta;
      if (time > stats_time-stats_eps)
	stats_time = ceil(time/stats_delta)*stats_delta;
    } else if (stats_ilogfreq == 1) {
      stats_time += stats_delta;
      stats_irepeat++;
      if (stats_irepeat == stats_nrepeat || time > stats_time-stats_eps) {
	// Calculate next smallest delta that will 
	// reach tgoal within nrepeat steps
	tgoal = time-stats_t0+stats_delta;
	ntmp = MAX(1,static_cast<int>
		   (ceil(log(tgoal/(stats_delta*stats_nrepeat))
			 /log(stats_scale))));
	stats_delta *= pow(stats_scale,ntmp);
	stats_time = ceil(tgoal/stats_delta)*stats_delta;
	stats_irepeat = 0;
      }
    }
  }
}

/* ----------------------------------------------------------------------
   print stats
------------------------------------------------------------------------- */

void Output::stats(int init_flag)
{
  char str[2048] = {'\0'};
  char *strpnt = str;

  app->stats(strpnt);
  strpnt += strlen(strpnt);

  if (init_flag) {
    sprintf(strpnt,"%10.3g ",timer->elapsed(TIME_LOOP));
    strpnt += strlen(strpnt);
  } else {
    sprintf(strpnt,"%10.3g ",0.0);
    strpnt += strlen(strpnt);
  }

  for (int i = 0; i < ndiags; i++) {
    diaglist[i]->stats(strpnt);
    strpnt += strlen(strpnt);
  }

  if (me == 0) {
    if (screen)
      fprintf(screen,"%s\n",str);
    if (logfile)
      fprintf(logfile,"%s\n",str);
  }
}

/* ----------------------------------------------------------------------
   print stats header
------------------------------------------------------------------------- */

void Output::stats_header()
{
  char str[2048] = {'\0'};
  char *strpnt = str;

  app->stats_header(strpnt);
  strpnt += strlen(strpnt);

  sprintf(strpnt,"%10s","CPU");
  strpnt += strlen(strpnt);

  for (int i = 0; i < ndiags; i++) {
    diaglist[i]->stats_header(strpnt);
    strpnt += strlen(strpnt);
  }

  if (me == 0) {
    if (screen) fprintf(screen,"%s\n",str);
    if (logfile) fprintf(logfile,"%s\n",str);
  }
}

/* ---------------------------------------------------------------------- */

void Output::add_diag(Diag *diag)
{
  ndiags++;
  diaglist = (Diag **) memory->srealloc(diaglist,ndiags*sizeof(Diag *),
					"output:diaglist");
  diaglist[ndiags-1] = diag;
}

/* ----------------------------------------------------------------------
   dump a snapshot of lattice values as atom coords
------------------------------------------------------------------------- */

void Output::dump()
{
  // proc 0 writes timestep header

  if (me == 0) {
    fprintf(fp,"ITEM: TIMESTEP\n");
    fprintf(fp,"%d\n",idump);
    fprintf(fp,"ITEM: NUMBER OF ATOMS\n");
    fprintf(fp,"%d\n",nglobal);
    fprintf(fp,"ITEM: BOX BOUNDS\n");
    fprintf(fp,"%g %g\n",boxxlo,boxxhi);
    fprintf(fp,"%g %g\n",boxylo,boxyhi);
    fprintf(fp,"%g %g\n",boxzlo,boxzhi);
    fprintf(fp,"ITEM: ATOMS\n");
  }

  idump++;

  // pack my info into buffer

  for (int n = 0; n < size_one; n++) (this->*pack_choice[n])(n);
  int me_size = nlocal*size_one;

  // proc 0 pings each proc, receives it's data, writes to file
  // all other procs wait for ping, send their data to proc 0

  int tmp,nlines;
  MPI_Status status;
  MPI_Request request;

  if (me == 0) {
    for (int iproc = 0; iproc < nprocs; iproc++) {
      if (iproc) {
	MPI_Irecv(buf,maxbuf,MPI_DOUBLE,iproc,0,world,&request);
	MPI_Send(&tmp,0,MPI_INT,iproc,0,world);
	MPI_Wait(&request,&status);
	MPI_Get_count(&status,MPI_DOUBLE,&nlines);
	nlines /= size_one;
      } else nlines = me_size/size_one;
      
      write_data(nlines,buf);
    }
    fflush(fp);

  } else {
    MPI_Recv(&tmp,0,MPI_INT,0,0,world,&status);
    MPI_Rsend(buf,me_size,MPI_DOUBLE,0,0,world);
  }
}

/* ---------------------------------------------------------------------- */

void Output::write_data(int n, double *buf)
{
  int i,j;

  int m = 0;
  for (i = 0; i < n; i++) {
    for (j = 0; j < size_one; j++) {
      if (vtype[j] == INT) fprintf(fp,vformat[j],static_cast<int> (buf[m]));
      else fprintf(fp,vformat[j],buf[m]);
      m++;
    }
    fprintf(fp,"\n");
  }
}

/* ----------------------------------------------------------------------
   one method for every keyword dump can output
   the site quantity is packed into buf starting at n with stride size_one
------------------------------------------------------------------------- */

/* ---------------------------------------------------------------------- */

void Output::pack_id(int n)
{
  int *id;
  id = applattice->id;

  for (int i = 0; i < nlocal; i++) {
    buf[n] = id[i];
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void Output::pack_lattice(int n)
{
  int i,j,k;

  int *lattice = applattice->lattice;
  for (i = 0; i < nlocal; i++) {
    buf[n] = lattice[i];
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void Output::pack_x(int n)
{
  double **xyz = applattice->xyz;

  for (int i = 0; i < nlocal; i++) {
    buf[n] = xyz[i][0];
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void Output::pack_y(int n)
{
  double **xyz = applattice->xyz;

  for (int i = 0; i < nlocal; i++) {
    buf[n] = xyz[i][1];
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void Output::pack_z(int n)
{
  double **xyz = applattice->xyz;

  for (int i = 0; i < nlocal; i++) {
    buf[n] = xyz[i][2];
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void Output::pack_energy(int n)
{
  int i,j,k;

  for (i = 0; i < nlocal; i++) {
    buf[n] = applattice->site_energy(i);
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void Output::pack_propensity(int n)
{
  int i,j,k;

  for (i = 0; i < nlocal; i++) {
    buf[n] = applattice->site_propensity(i);
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void Output::pack_integer(int n)
{
  int *ivec = applattice->iarray[vindex[n]];

  for (int i = 0; i < nlocal; i++) {
    buf[n] = ivec[i];
    n += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void Output::pack_double(int n)
{
  double *dvec = applattice->darray[n];

  for (int i = 0; i < nlocal; i++) {
    buf[n] = dvec[i];
    n += size_one;
  }
}
