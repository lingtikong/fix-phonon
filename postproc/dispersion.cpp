#include "dispersion.h"
#include "lmp_phonon.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"

#define MAX_STRING_LEN 256
#define MIN(a,b) ((a)>(b)?(b):(a))
#define MAX(a,b) ((a)>(b)?(a):(b))

/* ----------------------------------------------------------------------------
 * Class Dispersion is the main driver to calculate phonon DOS, phonon
 * dispersion curve and some other things.
 * ---------------------------------------------------------------------------- */
Dispersion::Dispersion(LMP_PHONON *dm)
{
  // pass the class from main
  dynmat = dm;

  // display the menu
  while ( 1 ){
    printf("\n"); for (int i=0; i<27;i++) printf("="); printf(" Menu "); for (int i=0; i<27;i++) printf("="); printf("\n");
    printf("  1. To calculate the phonon DOS;\n");
    printf("  2. To calculate the dispersion curve;\n");
    printf("  3. To get the dynamical matrix at arbitrary q;\n");
    printf("  4. To get the vibration frequency at arbitrary q;\n");
    printf("  5. To get dispersion-like curve for dynamical matrix;\n");
    printf("  6. Reset the interpolation method;\n");
    printf("  0. Exit.\n");
    // read user choice
    char str[MAX_STRING_LEN];
    job = 0;
    printf("Your choice[0]: ");
    if (strlen(gets(str)) > 0) sscanf(str,"%d",&job);
    printf("\nYou chose %d\n", job);
    for (int i=0; i<60;i++) printf("=");printf("\n\n");

    // now to do the job according to user's choice
    if      (job == 1) dos();
    else if (job == 2) disp();
    else if (job == 3) dmanyq(); 
    else if (job == 4) vfanyq(); 
    else if (job == 5) DMdisp(); 
    else if (job == 6) dynmat->getIntMeth();
    else break;
  }
}

/* ----------------------------------------------------------------------------
 * Deconstructor to free memory
 * ---------------------------------------------------------------------------- */
Dispersion::~Dispersion()
{
  dynmat = NULL;
}

/* ----------------------------------------------------------------------------
 * Private method to calculate the phonon DOS
 * ---------------------------------------------------------------------------- */
void Dispersion::dos()
{
  // ask for mesh info
  char str[MAX_STRING_LEN];
  int nx = dynmat->nx, ny = dynmat->ny, nz = dynmat->nz, nr=0;
  printf("\nThe q-mesh size from the read dynamical matrix is: %d x %d x %d\n", nx, ny, nz);
  printf("Please input your desired size to measure the DOS [%d %d %d]: ", nx, ny, nz);
  if (strlen(gets(str)) > 0) sscanf(str,"%d %d %d", &nx, &ny, &nz);

  if (nx<1||ny<1||nz<1) return;
  if (dynmat->nx == 1) nx = 1;
  if (dynmat->ny == 1) ny = 1;
  if (dynmat->nz == 1) nz = 1;
  int nq = nx*ny*nz;
  printf("Your new q-mesh size would be: %d x %d x %d\n", nx,ny,nz);

  // now to calculate the frequencies at all q-points
  int ndim = dynmat->fftdim;
  double egvs[nq][ndim], q[3];
  int idx=0;
  for (int i=0; i<nx; i++){
    q[0] = double(i)/double(nx);
    for (int j=0; j<ny; j++){
      q[1] = double(j)/double(ny);
      for (int k=0; k<nz; k++){
        q[2] = double(k)/double(nz);
        dynmat->getDMq(q);
        dynmat->geteigen(egvs[idx++]);
      }
    }
  }
  // now to get the frequency range
  double fmin, fmax;
  fmin = fmax = egvs[0][0];
  for (int i=0; i<nq; i++){
    for (int j=0; j<ndim; j++){
      fmin = MIN(fmin, egvs[i][j]);
      fmax = MAX(fmax, egvs[i][j]);
    }
  }

  // Now to ask for the output frequency range
  printf("\nThe frequency range of all q-points are: [%g %g]\n", fmin, fmax);
  printf("Please input the desired range to get DOS [%g %g]: ", fmin, fmax);
  if (strlen(gets(str)) > 0) sscanf(str,"%lg %lg", &fmin, &fmax);
  if (fmin > fmax){double swap = fmin; fmin = fmax; fmax = swap;}

  int intv = 100;
  printf("Please input the number of intervals [100]: ");
  if (strlen(gets(str)) > 0) sscanf(str,"%d", &intv);
  intv = MAX(2,intv);

  double finc = (fmax-fmin)/double(intv), finc_inv = 1./finc;;
  double dos[intv];
  for (int i=0; i<intv; i++) dos[i] = 0.;

  // now to calculate the DOS
  int total=0;
  for (int i=0; i<nq; i++){
    for (int j=0; j<ndim; j++){
      idx = int(egvs[i][j]*finc_inv);
      if (idx>=0 && idx<intv){dos[idx] += 1.; total++;}
    }
  }

  // smooth dos ?
  printf("Would you like to smooth the phonon dos? (y/n)[n]: ");
  if (strlen(gets(str)) > 0){
    if (strcmp(str,"y") == 0 || strcmp(str,"Y") == 0){
      smooth(dos, intv, fmin, finc);
    }
  }

  // normalize dos to 1
  double sum = 1./double(total)*finc_inv;
  for (int i=0; i<intv; i++) dos[i] *= sum;

  // now to output the phonon DOS
  nr = 0;
  while (nr !=1){
    printf("Please input the filename to output DOS: ");
    nr = scanf("%s", str); while (getchar() != '\n');
  }

  int n = strlen(str)+1;
  outfile = new char[n];
  strcpy(outfile, str);
  FILE *fp = fopen(outfile, "w");
  fprintf(fp,"# frequency  DOS\n");
  fprintf(fp,"#%s  number\n", dynmat->funit);
  double freq = 0.5*finc;
  for (int i=0; i<intv; i++){
    fprintf(fp,"%lg %lg\n", freq, dos[i]);
    freq += finc;
  }
  fclose(fp);
  delete []outfile;
}

void Dispersion::disp()
{
  // ask the output file name and write the header.
  char str[MAX_STRING_LEN];
  int nr = 0;
  while (nr != 1){
    printf("Please input the filename to output the dispersion data:");
    nr = scanf("%s", str); while (getchar() != '\n');
  }
  int n = strlen(str) + 1;
  outfile = new char[n];
  strcpy(outfile, str);
  FILE *fp = fopen(outfile, "w"); delete []outfile;
  fprintf(fp,"# q     qr    freq\n");
  fprintf(fp,"# 2pi/L  2pi/L %s\n", dynmat->funit);

  // now the calculate the dispersion curve
  double qstr[3], qend[3], q[3], qinc[3], qr=0., dq;
  int nq = MAX(MAX(dynmat->nx,dynmat->ny),dynmat->nz)/2;
  qend[0] = qend[1] = qend[2] = 0.;

  while (1){
    int ndim = dynmat->fftdim;
    double egvs[ndim];
    int quit=0;

    for (int i=0; i<3; i++) qstr[i] = qend[i];
    nr = 0;
    while (nr != 3){
      printf("\nPlease input the start q-point in unit of 2pi/L, q to exit [%g %g %g]: ", qstr[0], qstr[1], qstr[2]);
      if (strlen(gets(str)) > 0){
        if (strcmp(str,"q") ==0 || strcmp(str,"exit") ==0){quit = 1; break;}
        nr = sscanf(str,"%lg %lg %lg", &qstr[0], &qstr[1], &qstr[2]);
      } else  break;
    }
    if (quit) break;

    nr = 0;
    while (nr != 3){
      printf("Please input the end q-point in unit of 2pi/L: ");
      nr = scanf("%lg %lg %lg", &qend[0], &qend[1], &qend[2]); while ( getchar() != '\n' );
    }
    nr = 0;
    while (nr != 1){
      printf("Please input the # of points along the line [%d]: ", nq);
      if (strlen(gets(str)) > 0) nr = sscanf(str,"%d", &nq);
      else break;
    }
    nq = MAX(nq,2);

    for (int i=0; i<3; i++) qinc[i] = (qend[i]-qstr[i])/double(nq-1);
    dq = sqrt(qinc[0]*qinc[0]+qinc[1]*qinc[1]+qinc[2]*qinc[2]);

    for (int i=0; i<3; i++) q[i] = qstr[i];
    for (int ii=0; ii<nq; ii++){
      dynmat->getDMq(q);
      dynmat->geteigen(egvs);
      fprintf(fp,"%lg %lg %lg %lg ", q[0], q[1], q[2], qr);
      for (int i=0; i<ndim; i++) fprintf(fp," %lg", egvs[i]);
      fprintf(fp,"\n");

      for (int i=0; i<3; i++) q[i] += qinc[i];
      qr += dq;
    }
    qr -= dq;
  }
  fclose(fp);
}

/* ----------------------------------------------------------------------------
 * Private method to calculate the phonon DOS
 * ---------------------------------------------------------------------------- */
void Dispersion::dmanyq()
{
  double q[3];
  int nr = 0;
  while (nr != 3){
    printf("Please input the q-point to output the dynamical matrix:");
    nr = scanf("%lg %lg %lg", &q[0], &q[1], &q[2]); while ( getchar() != '\n' );
  }

  dynmat->getDMq(q);
  dynmat->writeDMq(q);
return;
}

void Dispersion::vfanyq()
{
  int ndim = dynmat->fftdim;
  double q[3], egvs[ndim];
  
  int nr = 0;
  while (nr != 3){
    printf("Please input the q-point to output the dynamical matrix:");
    nr = scanf("%lg %lg %lg", &q[0], &q[1], &q[2]); while (getchar() != '\n');
  }
  dynmat->getDMq(q);
  dynmat->geteigen(egvs);
  printf("q-point: [%lg %lg %lg]\n", q[0], q[1], q[2]);
  printf("Vibrational frequencies at this q-point:\n");
  for (int i=0; i<ndim; i++) printf("%lg ", egvs[i]); printf("\n\n");
}

/* ----------------------------------------------------------------------------
 * Private method to get the dispersion-like data for dynamical matrix
 * ---------------------------------------------------------------------------- */
void Dispersion::DMdisp()
{
  // ask the output file name and write the header.
  char str[MAX_STRING_LEN];
  int nr = 0;
  while (nr != 1){
    printf("Please input the filename to output the dispersion-like dynamical matrix data:");
    nr = scanf("%s", str); while (getchar() != '\n');
  }
  int n = strlen(str) + 1;
  outfile = new char[n];
  strcpy(outfile, str);
  FILE *fp = fopen(outfile, "w"); delete []outfile;
  fprintf(fp,"# q     qr    D\n");

  // now the calculate the dispersion-like curve
  double qstr[3], qend[3], q[3], qinc[3], qr=0., dq;
  int nq = MAX(MAX(dynmat->nx,dynmat->ny),dynmat->nz)/2;
  qend[0] = qend[1] = qend[2] = 0.;

  while (1){
    int quit=0;

    for (int i=0; i<3; i++) qstr[i] = qend[i];
    int nr = 0;
    while (nr != 3){
      printf("\nPlease input the start q-point in unit of 2pi/L, q to exit [%g %g %g]: ", qstr[0], qstr[1], qstr[2]);
      if (strlen(gets(str)) > 0){
        if (strcmp(str,"q") ==0 || strcmp(str,"exit") ==0){quit = 1; break;}
        nr = sscanf(str,"%lg %lg %lg", &qstr[0], &qstr[1], &qstr[2]);
      } else  break;
    }
    if (quit) break;

    nr = 0;
    while (nr != 3){
      printf("Please input the end q-point in unit of 2pi/L: ");
      nr = scanf("%lg %lg %lg", &qend[0], &qend[1], &qend[2]); while ( getchar() != '\n' );
    }
    nr = 0;
    while (nr != 1){
      printf("Please input the # of points along the line [%d]: ", nq);
      if (strlen(gets(str)) > 0) nr = sscanf(str,"%d", &nq);
      else break;
    }
    nq = MAX(nq,2);

    for (int i=0; i<3; i++) qinc[i] = (qend[i]-qstr[i])/double(nq-1);
    dq = sqrt(qinc[0]*qinc[0]+qinc[1]*qinc[1]+qinc[2]*qinc[2]);

    for (int i=0; i<3; i++) q[i] = qstr[i];
    for (int ii=0; ii<nq; ii++){
      dynmat->getDMq(q);
      dynmat->writeDMq(q, qr, fp);
      for (int i=0; i<3; i++) q[i] += qinc[i];
      qr += dq;
    }
    qr -= dq;
  }
  fclose(fp);
return;
}

/* ----------------------------------------------------------------------------
 * Private method to smooth the dos
 * ---------------------------------------------------------------------------- */
void Dispersion::smooth(double *array, int npt, double xmin, double delta)
{
  if (npt < 1) return;

  double *tmp, *table;
  tmp   = new double[npt];
  table = new double[npt];
  
  double sigma = 4.*delta, fac = 1./(sigma*sigma);
  double em = xmin + delta*0.5;
  for (int i=0; i<npt; i++){
    tmp[i] = 0.;
    table[i] = exp(-em*em*fac);
    em += delta;
  }
  em = 0.;
  for (int i=0; i<npt; i++){
    double peso = 0., emp = 0.;
    for (int j=0; j<npt; j++){
      double tij = table[abs(i-j)];

      tmp [i] += array[j]*tij;
      peso += tij;
      emp += delta;
    }
    array[i] = tmp[i]/peso;
    em += delta;
  }
  delete []tmp;
  delete []table;

return;
}
