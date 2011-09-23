#ifndef INTERPOLATION_H
#define INTERPOLATION_H

#include "stdio.h"
#include "stdlib.h"
#include "memory.h"
#include <tricubic.h>
extern "C"{
#include "f2c.h"
#include "clapack.h"
}

using namespace std;

class Interpolate{
public:
  Interpolate(int, int, int, int, doublecomplex **);
  ~Interpolate();

  int which;

  void execute(double *, doublecomplex *);

private:
  void tricubic(double *, doublecomplex *);
  void trilinear(double *, doublecomplex *);
  Memory *memory;

  int Nx, Ny, Nz, Npt, ndim;

  doublecomplex **data;
  doublecomplex **Dfdx, **Dfdy, **Dfdz, **D2fdxdy, **D2fdxdz, **D2fdydz, **D3fdxdydz;
  double a[64], f[8], dfdx[8], dfdy[8], dfdz[8], d2fdxdy[8], d2fdxdz[8], d2fdydz[8], d3fdxdydz[8];
  int vidx[8];
};

#endif