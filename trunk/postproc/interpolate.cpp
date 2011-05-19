#include "interpolate.h"
#include "math.h"

/* ----------------------------------------------------------------------------
 * Constructor used to get input from caller, and prepare for tricubic
 * interpolations.
 * ---------------------------------------------------------------------------- */
Interpolate::Interpolate(int nx, int ny, int nz, int ndm, doublecomplex **DM)
{
  Nx = nx;
  Ny = ny;
  Nz = nz;
  Npt = Nx*Ny*Nz;
  ndim = ndm;
  memory = new Memory;

  data = memory->create(data, Npt, ndim, "Interpolate_Interpolate:data");
  Dfdx = memory->create(Dfdx, Npt, ndim, "Interpolate_Interpolate:Dfdx");
  Dfdy = memory->create(Dfdy, Npt, ndim, "Interpolate_Interpolate:Dfdy");
  Dfdz = memory->create(Dfdz, Npt, ndim, "Interpolate_Interpolate:Dfdz");
  D2fdxdy = memory->create(D2fdxdy, Npt, ndim, "Interpolate_Interpolate:D2fdxdy");
  D2fdxdz = memory->create(D2fdxdz, Npt, ndim, "Interpolate_Interpolate:D2fdxdz");
  D2fdydz = memory->create(D2fdydz, Npt, ndim, "Interpolate_Interpolate:D2fdydz");
  D3fdxdydz = memory->create(D3fdxdydz, Npt, ndim, "Interpolate_Interpolate:D2fdxdydz");

  // copy data
  for (int n=0; n<Npt; n++){
    for (int idim=0; idim<ndim; idim++) data[n][idim] = DM[n][idim];
  }

  // get the derivatives
  int n=0;
  const double half = 0.5, one4 = 0.25, one8 = 0.125;
  for (int i=0; i<Nx; i++)
  for (int j=0; j<Ny; j++)
  for (int k=0; k<Nz; k++){
    //facx = facy = facz = 0.5;
    int ip = (i+1)%Nx, jp = (j+1)%Ny, kp = (k+1)%Nz;
    int im = (i-1+Nx)%Nx, jm = (j-1+Ny)%Ny, km = (k-1+Nz)%Nz;

    int p100 = (ip*Ny+j)*Nz+k;
    int p010 = (i*Ny+jp)*Nz+k;
    int p001 = (i*Ny+j)*Nz+kp;
    int p110 = (ip*Ny+jp)*Nz+k;
    int p101 = (ip*Ny+j)*Nz+kp;
    int p011 = (i*Ny+jp)*Nz+kp;
    int pm00 = (im*Ny+j)*Nz+k;
    int p0m0 = (i*Ny+jm)*Nz+k;
    int p00m = (i*Ny+j)*Nz+km;
    int pmm0 = (im*Ny+jm)*Nz+k;
    int pm0m = (im*Ny+j)*Nz+km;
    int p0mm = (i*Ny+jm)*Nz+km;
    int p1m0 = (ip*Ny+jm)*Nz+k;
    int p10m = (ip*Ny+j)*Nz+km;
    int p01m = (i*Ny+jp)*Nz+km;
    int pm10 = (im*Ny+jp)*Nz+k;
    int pm01 = (im*Ny+j)*Nz+kp;
    int p0m1 = (i*Ny+jm)*Nz+kp;
    int p111 = (ip*Ny+jp)*Nz+kp;
    int pm11 = (im*Ny+jp)*Nz+kp;
    int p1m1 = (ip*Ny+jm)*Nz+kp;
    int p11m = (ip*Ny+jp)*Nz+km;
    int pm1m = (im*Ny+jp)*Nz+km;
    int p1mm = (ip*Ny+jm)*Nz+km;
    int pmm1 = (im*Ny+jm)*Nz+kp;
    int pmmm = (im*Ny+jm)*Nz+km;

    for (int idim=0; idim<ndim; idim++){
      Dfdx[n][idim].r = (data[p100][idim].r - data[pm00][idim].r) * half;
      Dfdx[n][idim].i = (data[p100][idim].i - data[pm00][idim].i) * half;
      Dfdy[n][idim].r = (data[p010][idim].r - data[p0m0][idim].r) * half;
      Dfdy[n][idim].i = (data[p010][idim].i - data[p0m0][idim].i) * half;
      Dfdz[n][idim].r = (data[p001][idim].r - data[p00m][idim].r) * half;
      Dfdz[n][idim].i = (data[p001][idim].i - data[p00m][idim].i) * half;
      D2fdxdy[n][idim].r = (data[p110][idim].r - data[p1m0][idim].r - data[pm10][idim].r + data[pmm0][idim].r) * one4;
      D2fdxdy[n][idim].i = (data[p110][idim].i - data[p1m0][idim].i - data[pm10][idim].i + data[pmm0][idim].i) * one4;
      D2fdxdz[n][idim].r = (data[p101][idim].r - data[p10m][idim].r - data[pm01][idim].r + data[pm0m][idim].r) * one4;
      D2fdxdz[n][idim].i = (data[p101][idim].i - data[p10m][idim].i - data[pm01][idim].i + data[pm0m][idim].i) * one4;
      D2fdydz[n][idim].r = (data[p011][idim].r - data[p01m][idim].r - data[p0m1][idim].r + data[p0mm][idim].r) * one4;
      D2fdydz[n][idim].i = (data[p011][idim].i - data[p01m][idim].i - data[p0m1][idim].i + data[p0mm][idim].i) * one4;
      D3fdxdydz[n][idim].r = (data[p111][idim].r-data[pm11][idim].r - data[p1m1][idim].r - data[p11m][idim].r +
                              data[p1mm][idim].r+data[pm1m][idim].r + data[pmm1][idim].r - data[pmmm][idim].r) * one8;
      D3fdxdydz[n][idim].i = (data[p111][idim].i-data[pm11][idim].i - data[p1m1][idim].i - data[p11m][idim].i +
                              data[p1mm][idim].i+data[pm1m][idim].i + data[pmm1][idim].i - data[pmmm][idim].i) * one8;
    }
    n++;
  }
return;
}

/* ----------------------------------------------------------------------------
 * Deconstructor used to free memory
 * ---------------------------------------------------------------------------- */
Interpolate::~Interpolate()
{
  memory->destroy(data);
  memory->destroy(Dfdx);
  memory->destroy(Dfdy);
  memory->destroy(Dfdz);
  memory->destroy(D2fdxdy);
  memory->destroy(D2fdxdz);
  memory->destroy(D2fdydz);
  memory->destroy(D3fdxdydz);
  delete memory;
}

/* ----------------------------------------------------------------------------
 * Tricubic interpolation, by calling the tricubic library
 * ---------------------------------------------------------------------------- */
void Interpolate::tricubic(double *qin, doublecomplex *DMq)
{
  // qin should be in unit of 2*pi/L
  double q[3];
  for (int i=0; i<3; i++) q[i] = qin[i];
  for (int i=0; i<3; i++){
    while (q[i] < 0.)  q[i] += 1.;
    while (q[i] >= 1.) q[i] -= 1.;
  }
  
  int nx = Nx, ny = Ny, nz = Nz;
  int ix = int(q[0]*double(nx));
  int iy = int(q[1]*double(ny));
  int iz = int(q[2]*double(nz));
  double x = q[0]*double(nx)-double(ix);
  double y = q[1]*double(ny)-double(iy);
  double z = q[2]*double(nz)-double(iz);
  int ixp = (ix+1)%nx, iyp = (iy+1)%ny, izp = (iz+1)%nz;
  vidx[0] = (ix*Ny+iy)*Nz+iz;
  vidx[1] = (ixp*Ny+iy)*Nz+iz;
  vidx[2] = (ix*Ny+iyp)*Nz+iz;
  vidx[3] = (ixp*Ny+iyp)*Nz+iz;
  vidx[4] = (ix*Ny+iy)*Nz+izp;
  vidx[5] = (ixp*Ny+iy)*Nz+izp;
  vidx[6] = (ix*Ny+iyp)*Nz+izp;
  vidx[7] = (ixp*Ny+iyp)*Nz+izp;
  for (int idim=0; idim<ndim; idim++){
    for (int i=0; i<8; i++){
      f[i] = data[vidx[i]][idim].r;
      dfdx[i] = Dfdx[vidx[i]][idim].r;
      dfdy[i] = Dfdy[vidx[i]][idim].r;
      dfdz[i] = Dfdz[vidx[i]][idim].r;
      d2fdxdy[i] = D2fdxdy[vidx[i]][idim].r;
      d2fdxdz[i] = D2fdxdz[vidx[i]][idim].r;
      d2fdydz[i] = D2fdydz[vidx[i]][idim].r;
      d3fdxdydz[i] = D3fdxdydz[vidx[i]][idim].r;
    }
    tricubic_get_coeff(&a[0],&f[0],&dfdx[0],&dfdy[0],&dfdz[0],&d2fdxdy[0],&d2fdxdz[0],&d2fdydz[0],&d3fdxdydz[0]); 
    DMq[idim].r = tricubic_eval(&a[0],x,y,z);
    
    for (int i=0; i<8; i++){
      f[i] = data[vidx[i]][idim].i;
      dfdx[i] = Dfdx[vidx[i]][idim].i;
      dfdy[i] = Dfdy[vidx[i]][idim].i;
      dfdz[i] = Dfdz[vidx[i]][idim].i;
      d2fdxdy[i] = D2fdxdy[vidx[i]][idim].i;
      d2fdxdz[i] = D2fdxdz[vidx[i]][idim].i;
      d2fdydz[i] = D2fdydz[vidx[i]][idim].i;
      d3fdxdydz[i] = D3fdxdydz[vidx[i]][idim].i;
    }
    tricubic_get_coeff(&a[0],&f[0],&dfdx[0],&dfdy[0],&dfdz[0],&d2fdxdy[0],&d2fdxdz[0],&d2fdydz[0],&d3fdxdydz[0]); 
    DMq[idim].i = tricubic_eval(&a[0],x,y,z);
  }
}

/* ----------------------------------------------------------------------------
 * method to interpolate the DM at an arbitrary q point;
 * the input q should be a vector in unit of (2pi/a 2pi/b 2pi/c).
 * All q components will be rescaled into [0 1).
 * ---------------------------------------------------------------------------- */
void Interpolate::trilinear(double *qin, doublecomplex *DMq)
{
  // rescale q[i] into [0 1)
  double q[3];
  for (int i=0; i<3; i++) q[i] = qin[i];
  for (int i=0; i<3; i++){
    while (q[i] < 0.)  q[i] += 1.;
    while (q[i] >= 1.) q[i] -= 1.;
  }

  // find the index of the eight vertice
  int ix, iy, iz, ixp, iyp, izp;
  double x, y, z;
  q[0] *= double(Nx);
  q[1] *= double(Ny);
  q[2] *= double(Nz);

  ix = int(q[0])%Nx;
  iy = int(q[1])%Ny;
  iz = int(q[2])%Nz;
  ixp = (ix+1)%Nx;
  iyp = (iy+1)%Ny;
  izp = (iz+1)%Nz;
  x = q[0] - double(ix);
  y = q[1] - double(iy);
  z = q[2] - double(iz);

//--------------------------------------
  int vindex[8];
  vindex[0] = ((ix*Ny)+iy)*Nz + iz;
  vindex[1] = ((ixp*Ny)+iy)*Nz + iz;
  vindex[2] = ((ix*Ny)+iyp)*Nz + iz;
  vindex[3] = ((ix*Ny)+iy)*Nz + izp;
  vindex[4] = ((ixp*Ny)+iy)*Nz + izp;
  vindex[5] = ((ix*Ny)+iyp)*Nz + izp;
  vindex[6] = ((ixp*Ny)+iyp)*Nz + iz;
  vindex[7] = ((ixp*Ny)+iyp)*Nz + izp;
  double fac[8];
  fac[0] = (1.-x)*(1.-y)*(1.-z);
  fac[1] = x*(1.-y)*(1.-z);
  fac[2] = (1.-x)*y*(1.-z);
  fac[3] = (1.-x)*(1.-y)*z;
  fac[4] = x*(1.-y)*z;
  fac[5] = (1.-x)*y*z;
  fac[6] = x*y*(1.-z);
  fac[7] = x*y*z;
  
  // now to do the interpolation
  for (int idim=0; idim<ndim; idim++){
    DMq[idim].r = 0.;
    DMq[idim].i = 0.;
    for (int i=0; i<8; i++){
      DMq[idim].r += data[vindex[i]][idim].r*fac[i];
      DMq[idim].i += data[vindex[i]][idim].i*fac[i];
    }
  }

return;
}

/* ----------------------------------------------------------------------------
 * To invoke the interpolation
 * ---------------------------------------------------------------------------- */
void Interpolate::execute(double *qin, doublecomplex *DMq)
{
  if (which) // 1: tricubic
    tricubic(qin, DMq);
  else       // otherwise: trilinear
    trilinear(qin, DMq);
return;
}
