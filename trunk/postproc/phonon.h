#ifndef PHONON_H
#define PHONON_H

#include <complex>
#include "dynmat.h"
#include "memory.h"

using namespace std;

class Phonon{
public:
  Phonon(DynMat *);
  ~Phonon();

  DynMat *dynmat;

private:
  int job;
  char *outfile;

  int nq, ndim;
  double **qpts, *wt;
  double **eigs;

  Memory *memory;

  void QMesh();
  void ComputeAll();

  void dos();
  void disp();
  void therm();

  void ldos_egv();
  void ldos_rsgf();

  void dmanyq();
  void vfanyq();
  void DMdisp();

  void smooth(double *, int, double, double);
};

#endif
