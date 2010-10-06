#ifndef LMP_PHONON_H
#define LMP_PHONON_H

#include "stdio.h"
#include "stdlib.h"
#include "memory.h"
#include "interpolate.h"
extern "C"{
#include "f2c.h"
#include "clapack.h"
}

class LMP_PHONON
{
public:

  LMP_PHONON(int, char**);
  ~LMP_PHONON();

  int nx, ny, nz, fftdim;
  char *funit;

  void getDMq(double *);
  void writeDMq(double *);
  void writeDMq(double *, const double, FILE *fp);
  int geteigen(double *);
  void getIntMeth();

private:

  Interpolate *interpolate;
  
  Memory *memory;
  int nucell, sysdim;
  int npt, fftdim2;

  char *binfile, *dmfile, *dmdfile;
  double boltz, eml2f, q[3];

  doublecomplex **DM_all, **DM_q;

};
#endif
