#include "dynmat.h"
#include "string.h"
#include "math.h"

#define MAXLINE 256

// to intialize the class
DynMat::DynMat(int narg, char **arg)
{
  // get the binary file name from command line option or user input
  char str[MAXLINE];
  if (narg < 2) {
    do printf("\nPlease input the binary file name from fix_phonon: ");
    while (strlen(gets(str)) < 1);

    int n = strlen(str) + 1;
    binfile = new char[n];
    strcpy(binfile, str);
  } else {
    int n = strlen(arg[1]) + 1;
    binfile = new char[n];
    strcpy(binfile, arg[1]);
  }

  // open the binary file
  FILE *fp = fopen(binfile, "rb");
  if (fp == NULL) {
    printf("\nFile %s not found! Programe terminated.\n", binfile);
    exit(1);
  }

  // read data from the binary file
  if ( fread(&sysdim, sizeof(int),    1, fp) != 1) {printf("\nError while reading sysdim from file: %s\n", binfile); fclose(fp); exit(2);}
  if ( fread(&nx,     sizeof(int),    1, fp) != 1) {printf("\nError while reading nx from file: %s\n", binfile); fclose(fp); exit(2);}
  if ( fread(&ny,     sizeof(int),    1, fp) != 1) {printf("\nError while reading ny from file: %s\n", binfile); fclose(fp); exit(2);}
  if ( fread(&nz,     sizeof(int),    1, fp) != 1) {printf("\nError while reading nz from file: %s\n", binfile); fclose(fp); exit(2);}
  if ( fread(&nucell, sizeof(int),    1, fp) != 1) {printf("\nError while reading nucell from file: %s\n", binfile); fclose(fp); exit(2);}
  if ( fread(&boltz,  sizeof(double), 1, fp) != 1) {printf("\nError while reading boltz from file: %s\n", binfile); fclose(fp); exit(2);}

  fftdim = sysdim*nucell; fftdim2 = fftdim*fftdim;
  npt = nx*ny*nz;

  // display info related to the read file
  printf("\n"); for (int i=0; i<60; i++) printf("#"); printf("\n");
  printf("Dynamical matrix is read from file: %s\n", binfile);
  printf("The system size in three dimension: %d x %d x %d\n", nx, ny, nz);
  printf("Number of atoms per unit cell     : %d\n", nucell);
  printf("System dimension                  : %d\n", sysdim);
  printf("Boltzmann constant in used units  : %g\n", boltz);
  for (int i=0; i<60; i++) printf("#"); printf("\n\n");
  if (sysdim<1||sysdim>3||nx<1||ny<1||nz<1||nucell<1){
    printf("Wrong values read from header of file: %s, please check the binary file!\n", binfile);
    fclose(fp); exit(3);
  }
  funit = new char[4];
  strcpy(funit, "THz");
  if (boltz == 1.){eml2f = 1.; delete funit; funit=new char[22]; strcpy(funit,"sqrt(epsilon/(m.sigma^2))");}
  else if (boltz == 0.0019872067) eml2f = 3.256576161;
  else if (boltz == 8.617343e-5) eml2f = 15.63312493;
  else if (boltz == 1.3806504e-23) eml2f = 1.;
  else if (boltz == 1.3806504e-16) eml2f = 1.591549431e-14;
  else {
    printf("WARNING: Because of float precision, I cannot get the factor to convert sqrt(E/ML^2)\n");
    printf("into THz, instead, I set it to be 1; you should check the unit used by LAMMPS.\n");
    eml2f = 1.;
  }

  // now to allocate memory for DM
  memory = new Memory;
  DM_all = memory->create_2d_complex_array(npt, fftdim2, "DynMat:DM_all");
  DM_q   = memory->create_2d_complex_array(fftdim,fftdim,"DynMat:DM_q");

  // read all dynamical matrix info into DM_all
  if ( fread(DM_all[0], sizeof(doublecomplex), npt*fftdim2, fp) != size_t(npt*fftdim2)){
    printf("\nError while reading the DM from file: %s\n", binfile);
    fclose(fp);
    exit(1);
  }

  // now try to read unit cell info from the binary file
  flag_latinfo = 0;
  basis = memory->create_2d_double_array(nucell,sysdim,"DynMat:basis");
  attyp = new int[nucell];
  
  if ( fread(&Tmeasure,   sizeof(double), 1, fp) == 1) flag_latinfo |= 1;
  if ( fread(&basevec[0], sizeof(double), 9, fp) == 9) flag_latinfo |= 2;
  if ( fread(basis[0],    sizeof(double), fftdim, fp) == fftdim) flag_latinfo |= 4;
  if ( fread(&attyp[0],   sizeof(int),    nucell, fp) == nucell) flag_latinfo |= 8;
  fclose(fp);

  if ((flag_latinfo&15) == 15){
    flag_latinfo = 1;
    car2dir();
  } else {
    flag_latinfo = 0;
    Tmeasure = 0.;
  }

  // ask for the interpolation method
  interpolate =  new Interpolate(nx, ny, nz, fftdim2, DM_all);
  getIntMeth();

  return;
}

// to destroy the class
DynMat::~DynMat()
{
 // destroy all memory allocated
 delete []binfile;
 delete []funit;
 if (dmfile) delete []dmfile;
 if (interpolate) delete interpolate;
 if (attyp) delete []attyp;

 memory->destroy_2d_complex_array(DM_all);
 memory->destroy_2d_complex_array(DM_q);
 memory->destroy_2d_double_array(basis);
 delete memory;
}

/* ----------------------------------------------------------------------------
 * method to write DM_q to file, single point
 * ---------------------------------------------------------------------------- */
void DynMat::writeDMq(double *q)
{
  FILE *fp;
  // only ask for file name for the first time
  // other calls will append the result to the file.
  if (dmfile == NULL){
    char str[MAXLINE];
    do  printf("\nPlease input the filename to output the DM at selected q: ");
    while (strlen(gets(str)) < 1);
    int n = strlen(str) + 1;
    dmfile = new char[n];
    strcpy(dmfile, str);
    fp = fopen(dmfile,"w");
  } else {
    fp = fopen(dmfile,"a");
  }
  fprintf(fp,"# q = [%lg %lg %lg]\n", q[0], q[1], q[2]);

  for (int i=0; i<fftdim; i++){
    for (int j=0; j<fftdim; j++) fprintf(fp,"%lg %lg\t", DM_q[i][j].r, DM_q[i][j].i);
    fprintf(fp,"\n");
  }
  fprintf(fp,"\n");
  fclose(fp);
return;
}

/* ----------------------------------------------------------------------------
 * method to write DM_q to file, dispersion-like
 * ---------------------------------------------------------------------------- */
void DynMat::writeDMq(double *q, const double qr, FILE *fp)
{

  fprintf(fp, "%lg %lg %lg %lg ", q[0], q[1], q[2], qr);

  for (int i=0; i<fftdim; i++){
    for (int j=0; j<fftdim; j++) fprintf(fp,"%lg %lg\t", DM_q[i][j].r, DM_q[i][j].i);
  }
  fprintf(fp,"\n");
return;
}

/* ----------------------------------------------------------------------------
 * method to evaluate the eigenvalues of current q-point;
 * return the eigenvalues in egv.
 * cLapack subroutine zheevd is employed.
 * ---------------------------------------------------------------------------- */
int DynMat::geteigen(double *egv, int flag)
{
  char jobz, uplo;
  integer n, lda, lwork, lrwork, *iwork, liwork, info;
  doublecomplex *work;
  doublereal *w = &egv[0], *rwork;

  n     = fftdim;
  if (flag) jobz = 'V';
  else jobz = 'N';

  uplo = 'U';
  lwork = (n+2)*n;
  lrwork = 1 + (5+n+n)*n;
  liwork = 3 + 5*n;
  lda    = n;
  work  = new doublecomplex [lwork];
  rwork = new double [lrwork];
  iwork = new long [liwork];

  zheevd_(&jobz, &uplo, &n, DM_q[0], &lda, w, work, &lwork, rwork, &lrwork, iwork, &liwork, &info);
 
  // to get w instead of w^2; and convert w into v (THz hopefully)
  for (int i=0; i<n; i++){
    if (w[i]>= 0.) w[i] = sqrt(w[i]);
    else w[i] = -sqrt(-w[i]);

    w[i] *= eml2f;
  }

  delete []work;
  delete []rwork;
  delete []iwork;

return info;
}

/* ----------------------------------------------------------------------------
 * method to get the Dynamical Matrix at q
 * ---------------------------------------------------------------------------- */
void DynMat::getDMq(double *q)
{
  interpolate->execute(q, DM_q[0]);

return;
}

/* ----------------------------------------------------------------------------
 * method to select the interpolation method.
 * ---------------------------------------------------------------------------- */
void DynMat::getIntMeth()
{
  char str[MAXLINE];
  int im = 1;
  printf("\n");for(int i=0; i<60; i++) printf("=");
  printf("\nWhich interpolation method would you like to use?\n");
  printf("  1. Tricubic;\n  2. Trilinear;\n");
  printf("Your choice [1]: ");
  if (strlen(gets(str)) >0) im = atoi(strtok(str," \t\n\r\f"));

  im =2-im%2;
  interpolate->which = im;
  printf("Your chose: %d\n", im);
  for(int i=0; i<60; i++) printf("="); printf("\n\n");

return;
}


/* ----------------------------------------------------------------------------
 * private method to convert the cartisan coordinate of basis into fractional
 * ---------------------------------------------------------------------------- */
void DynMat::car2dir()
{
  if (sysdim == 1){
    double h_inv = 1./basevec[0];
    for (int i=0; i<nucell; i++) basis[i][0] *= h_inv;
  } else if (sysdim == 2){
    double h[3], h_inv[3];
    h[0] = basevec[0]; h[1] = basevec[4]; h[2] = basevec[3];
    h_inv[0] = 1./h[0]; h_inv[1] = 1./h[1];
    h_inv[2] = -h[2]/(h[0]*h[1]);
    for (int i=0; i<nucell; i++){
      double x[2];
      x[0] = basis[i][0]; x[1] = basis[i][1];
      basis[i][0] = h_inv[0]*x[0] + h_inv[2]*x[1];
      basis[i][1] = h_inv[1]*x[1];
    }
  } else {
    double h[6], h_inv[6];
    h[0] = basevec[0]; h[1] = basevec[4]; h[2] = basevec[8];
    h[3] = basevec[7]; h[4] = basevec[6]; h[5] = basevec[3];
    for (int i=0; i<3; i++) h_inv[i] = 1./h[i];
    h_inv[3] = -h[3]/(h[1]*h[2]);
    h_inv[4] = (h[3]*h[5]-h[1]*h[4])/(h[0]*h[1]*h[2]);
    h_inv[5] = -h[5]/(h[0]*h[1]);
   
    for (int i=0; i<nucell; i++){
      double x[3];
      x[0] = basis[i][0]; x[1] = basis[i][1]; x[2] = basis[i][2];
      basis[i][0] = h_inv[0]*x[0] + h_inv[5]*x[1] + h_inv[4]*x[2];
      basis[i][1] = h_inv[1]*x[1] + h_inv[3]*x[2];
      basis[i][2] = h_inv[2]*x[2];
    }
  }

return;
}
