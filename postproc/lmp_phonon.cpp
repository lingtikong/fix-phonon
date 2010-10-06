#include "lmp_phonon.h"
#include "string.h"
#include "math.h"

using namespace std;

#define MAX_STRING_LEN 256

// to intialize the class
LMP_PHONON::LMP_PHONON(int argc, char **argv)
{
  // get the binary file name from command line option or user input
  char str[MAX_STRING_LEN];
  if (argc<2){
    int nr = 0;
    while (nr != 1){
      printf("\nPlease input the binary file name from fix_phonon: ");
      nr = scanf("%s", str); while (getchar() != '\n');
    }
    int n = strlen(str) + 1;
    binfile = new char[n];
    strcpy(binfile, str);
  } else {
    int n = strlen(argv[1]) + 1;
    binfile = new char[n];
    strcpy(binfile, argv[1]);
  }

  // open the binary file
  FILE *fp = fopen(binfile, "rb");
  if (fp==NULL ){
    printf("\nFile %s not found! Programe terminated.\n", binfile);
    exit(1);
  }

  // read data from the binary file
  size_t nr;
  nr = fread(&sysdim, sizeof(int), 1, fp); if (nr != 1){printf("\nError while reading from file: %s\n", binfile); fclose(fp); exit(2);}
  nr = fread(&nx,     sizeof(int), 1, fp); if (nr != 1){printf("\nError while reading from file: %s\n", binfile); fclose(fp); exit(2);}
  nr = fread(&ny,     sizeof(int), 1, fp); if (nr != 1){printf("\nError while reading from file: %s\n", binfile); fclose(fp); exit(2);}
  nr = fread(&nz,     sizeof(int), 1, fp); if (nr != 1){printf("\nError while reading from file: %s\n", binfile); fclose(fp); exit(2);}
  nr = fread(&nucell, sizeof(int), 1, fp); if (nr != 1){printf("\nError while reading from file: %s\n", binfile); fclose(fp); exit(2);}
  nr = fread(&boltz,  sizeof(double), 1, fp); if (nr != 1){printf("\nError while reading from file: %s\n", binfile); fclose(fp); exit(2);}

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
  DM_all = memory->create_2d_complex_array(npt, fftdim2, "LMP_PHONON:DM_all");
  DM_q   = memory->create_2d_complex_array(fftdim,fftdim,"LMP_PHONON:DM_q");

  // read all dynamical matrix info into DM_all
  nr = fread(DM_all[0], sizeof(doublecomplex), npt*fftdim2, fp);
  fclose(fp);
  if (nr != size_t(npt*fftdim2)) {printf("\nError while reading from file: %s\n", binfile); exit(4);}

  // ask for the interpolation method
  interpolate =  new Interpolate(nx, ny, nz, fftdim2, DM_all);
  getIntMeth();

  return;
}

// to destroy the class
LMP_PHONON::~LMP_PHONON()
{
 // destroy all memory allocated
 delete []binfile;
 delete []funit;
 if (dmfile) delete []dmfile;
 if (dmdfile) delete []dmdfile;
 if (interpolate) delete interpolate;

 memory->destroy_2d_complex_array(DM_all);
 memory->destroy_2d_complex_array(DM_q);
 delete memory;
}

/* ----------------------------------------------------------------------------
 * method to write DM_q to file, single point
 * ---------------------------------------------------------------------------- */
void LMP_PHONON::writeDMq(double *q)
{
  FILE *fp;
  // only ask for file name for the first time
  // other calls will append the result to the file.
  if (dmfile == NULL){
    char str[MAX_STRING_LEN];
    int nr = 0;
    while (nr != 1){
      printf("\nPlease input the filename to output the DM at selected q: ");
      nr = scanf("%s", str); while ( getchar() != '\n' );
    }
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
void LMP_PHONON::writeDMq(double *q, const double qr, FILE *fp)
{

  fprintf(fp, "%lg %lg %lg %lg ", q[0], q[1], q[2], qr);

  for (int i=0; i<fftdim; i++){
    for (int j=0; j<fftdim; j++) fprintf(fp,"%lg %lg\t", DM_q[i][j].r, DM_q[i][j].i);
  }
  fprintf(fp,"\n");
return;
}

int LMP_PHONON::geteigen(double *egv)
{
  char jobz, uplo;
  integer n, lda, lwork, lrwork, *iwork, liwork, info;
  doublecomplex *work;
  doublereal *w = &egv[0], *rwork;

  jobz = 'N'; uplo = 'U';
  n     = fftdim;
  lwork = n + 1;
  lda    = n;
  lrwork = n;
  liwork = n;
  work = new doublecomplex [lwork];
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

void LMP_PHONON::getDMq(double *q)
{
  interpolate->execute(q, DM_q[0]);

return;
}

void LMP_PHONON::getIntMeth()
{
  char str[MAX_STRING_LEN];
  int im=1;
  printf("\n");for(int i=0; i<60; i++) printf("=");
  printf("\nWhich interpolation method would you like to use?\n");
  printf("  1. Tricubic;\n  2. Trilinear;\n");
  printf("Your choice [1]: ");
  if (strlen(gets(str)) >0) sscanf(str,"%d", &im);
  im =2-im%2;
  interpolate->which = im;
  printf("Your chose: %d\n", im);
  for(int i=0; i<60; i++) printf("="); printf("\n\n");

return;
}
