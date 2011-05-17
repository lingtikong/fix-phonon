#include "string.h"
#include "phonon.h"
#include "green.h"
#include "timer.h"

#ifdef UseSPG
extern "C"{
#include "spglib.h"
}
#endif

#define MAXLINE 256
#define MIN(a,b) ((a)>(b)?(b):(a))
#define MAX(a,b) ((a)>(b)?(a):(b))

/* ----------------------------------------------------------------------------
 * Class Phonon is the main driver to calculate phonon DOS, phonon
 * dispersion curve and some other things.
 * ---------------------------------------------------------------------------- */
Phonon::Phonon(DynMat *dm)
{
  // create memory 
  memory = new Memory();

  // pass the class from main
  dynmat = dm;
  sysdim = dynmat->sysdim;
  ndim   = dynmat->fftdim;
  dos  = NULL;
  ldos = NULL;
  qpts = NULL;
  wt   = NULL;
  eigs = NULL;
  locals = NULL;

#ifdef UseSPG
  attyp = NULL;
  atpos = NULL;
#endif

  // display the menu
  char str[MAXLINE];
  while ( 1 ){
    printf("\n"); for (int i=0; i<27;i++) printf("="); printf(" Menu "); for (int i=0; i<27;i++) printf("="); printf("\n");
    printf("  1. Phonon DOS evaluation;\n");
    printf("  2. Phonon dispersion curves;\n");
    printf("  3. Dynamical matrix at arbitrary q;\n");
    printf("  4. Vibrational frequencies at arbitrary q;\n");
    printf("  5. Dispersion-like curve for dynamical matrix;\n");
    printf("  6. Vibrational thermodynamical properties;\n");
    printf("  7. Local phonon DOS from eigenvectors;\n");
    printf("  8. Local phonon DOS by RSGF method;\n");
    printf("  9. Reset the interpolation method;\n");
    printf("  0. Exit.\n");
    // read user choice
    int job = 0;
    printf("Your choice [0]: ");
    if (count_words(gets(str)) > 0) job = atoi(strtok(str," \t\n\r\f"));
    printf("\nYou chose: %d\n", job);
    for (int i=0; i<60;i++) printf("=");printf("\n\n");

    // now to do the job according to user's choice
    if      (job == 1) pdos();
    else if (job == 2) pdisp();
    else if (job == 3) dmanyq(); 
    else if (job == 4) vfanyq(); 
    else if (job == 5) DMdisp(); 
    else if (job == 6) therm(); 
    else if (job == 7) ldos_egv(); 
    else if (job == 8) ldos_rsgf(); 
    else if (job == 9) dynmat->getIntMeth();
    else break;
  }
return;
}

/* ----------------------------------------------------------------------------
 * Deconstructor to free memory
 * ---------------------------------------------------------------------------- */
Phonon::~Phonon()
{
  dynmat = NULL;
  if (qpts) memory->destroy_2d_double_array(qpts);
  if (eigs) memory->destroy_2d_double_array(eigs);
  if (wt)  delete []wt;
  
  if (locals) delete []locals;

  if (dos) delete []dos;
  if (ldos) memory->destroy_3d_double_array(ldos);

#ifdef UseSPG
  if (attyp) delete []attyp;
  if (atpos) memory->destroy_2d_double_array(atpos);
#endif
  if (memory) delete memory;
}

/* ----------------------------------------------------------------------------
 * Private method to calculate the phonon DOS
 * ---------------------------------------------------------------------------- */
void Phonon::pdos()
{
  // get frequencies on a q-mesh
  QMesh();       // generate q-points, hopefully irreducible
  ComputeAll();  // get all eigen values ==> frequencies

  // now to get the frequency range
  char str[MAXLINE];

  fmin = fmax = eigs[0][0];
  for (int iq=0; iq<nq; iq++){
    for (int j=0; j<ndim; j++){
      fmin = MIN(fmin, eigs[iq][j]);
      fmax = MAX(fmax, eigs[iq][j]);
    }
  }

  // Now to ask for the output frequency range
  printf("\nThe frequency range of all q-points are: [%g %g]\n", fmin, fmax);
  printf("Please input the desired range to get DOS [%g %g]: ", fmin, fmax);
  if (count_words(gets(str)) >= 2){
    fmin = atof(strtok(str," \t\n\r\f"));
    fmax = atof(strtok(NULL," \t\n\r\f"));
  }
  if (fmin > fmax){double swap = fmin; fmin = fmax; fmax = swap;}
  printf("The fequency range for your phonon DOS is [%g %g].\n", fmin, fmax);

  ndos = 201;
  printf("Please input the number of intervals [%d]: ", ndos);
  if (count_words(gets(str)) > 0) ndos = atoi(str);
  ndos += (ndos+1)%2; ndos = MAX(2,ndos);

  
  df  = (fmax-fmin)/double(ndos-1);
  rdf = 1./df;
  if (dos) delete []dos;
  dos = new double[ndos];
  for (int i=0; i<ndos; i++) dos[i] = 0.;

  // now to calculate the DOS
  double offset = fmin-0.5*df;
  for (int iq=0; iq<nq; iq++)
  for (int j=0; j<ndim; j++){
    int idx = int((eigs[iq][j]-offset)*rdf);
    if (idx>=0 && idx<ndos) dos[idx] += wt[iq];
  }

  // smooth dos ?
  printf("Would you like to smooth the phonon dos? (y/n)[n]: ");
  if (count_words(gets(str)) > 0){
    if (strcmp(str,"y") == 0 || strcmp(str,"Y") == 0){
      smooth(dos, ndos, fmin, df);
    }
  }

  // normalize dos to 1
  Normalize();

  // output DOS
  writeDOS();

return;
}

/* ----------------------------------------------------------------------------
 * Private method to write the phonon DOS to file
 * ---------------------------------------------------------------------------- */
void Phonon::writeDOS()
{
  if (dos == NULL) return;

  char str[MAXLINE];
  // now to output the phonon DOS
  printf("\nPlease input the filename to output DOS [pdos.dat]: ");
  if (count_words(gets(str)) < 1) strcpy(str, "pdos.dat");
  char *fname = strtok(str," \t\n\r\f");

  FILE *fp = fopen(fname, "w"); fname = NULL;
  fprintf(fp,"# frequency  DOS\n");
  fprintf(fp,"#%s  number\n", dynmat->funit);
  double freq = fmin;
  for (int i=0; i<ndos; i++){
    fprintf(fp,"%lg %lg\n", freq, dos[i]);
    freq += df;
  }
  fclose(fp);

return;
}

/* ----------------------------------------------------------------------------
 * Private method to write the local DOS to files.
 * ---------------------------------------------------------------------------- */
void Phonon::writeLDOS()
{
  if (ldos == NULL) return;

  const double one3 = 1./double(sysdim);
  char str[MAXLINE];
  for (int ilocal=0; ilocal<nlocal; ilocal++){
    sprintf(str,"pldos_%d.dat", locals[ilocal]);
    char *fname = strtok(str," \t\n\r\f");

    FILE *fp = fopen(fname, "w"); fname = NULL;
    fprintf(fp,"#freq xDOS yDOS zDOS total\n");
    double freq = fmin;
    for (int i=0; i<ndos; i++){
      fprintf(fp,"%lg", freq);
      double total = 0.;
      for (int idim=0; idim<sysdim; idim++) {fprintf(fp," %lg",ldos[i][ilocal][idim]); total += ldos[i][ilocal][idim];}
      fprintf(fp," %lg\n", total*one3);
      freq += df;
    }
    fclose(fp);
  }

return;
}

/* ----------------------------------------------------------------------------
 * Private method to calculate the local phonon DOS via the real space Green's
 * function method
 * ---------------------------------------------------------------------------- */
void Phonon::ldos_rsgf()
{
  char str[MAXLINE];
  const double tpi = 8.*atan(1.);
  double **Hessian, scale;
  scale = dynmat->eml2f*tpi; scale *= scale;
  Hessian = memory->create_2d_double_array(ndim, ndim, "dispersion_ldos_Hessian");

  double q0[3];
  q0[0] = q0[1] = q0[2] = 0.;

  dynmat->getDMq(q0);
  for (int i=0; i<ndim; i++)
  for (int j=0; j<ndim; j++) Hessian[i][j] = dynmat->DM_q[i][j].r*scale;

  double *egvs = new double [ndim];
  dynmat->geteigen(egvs, 0);

  fmin = fmax = egvs[0];
  for (int i=1; i<ndim; i++){fmin = MIN(fmin, egvs[i]); fmax = MAX(fmax, egvs[i]);}
  delete []egvs;

  ndos = 201;
  int ik = 0, nit = MAX(ndim*0.1, MIN(ndim,50));
  double eps = 12.; // for Cu with 1000+ atoms, 12 is enough; for small system, eps should be large.

  while (1) {
    // ask for relevant info
    printf("\nThere are %d atoms in each unit cell of your lattice.\n", dynmat->nucell);
    printf("Please input the index of atom in the unit cell to evaluate LDOS [%d]: ", ik);
    if (strlen(gets(str)) > 0) ik = atoi(strtok(str," \t\n\r\f"));
    if (ik < 0 || ik >= dynmat->nucell) break;

    printf("Please input the frequency range to evaluate LDOS [%g %g]: ", fmin, fmax);
    if (count_words(gets(str)) >= 2){
      fmin = atof(strtok(str," \t\n\r\f"));
      fmax = atof(strtok(NULL," \t\n\r\f"));
    }
    if (fmax < fmin) break;
    printf("The frequency range for your LDOS is [%g %g].\n", fmin, fmax);

    printf("Please input the desired number of points in LDOS [%d]: ", ndos);
    if (count_words(gets(str)) > 0) ndos = atoi(strtok(str," \t\n\r\f"));
    if (ndos < 2) break;
    ndos += (ndos+1)%2;

    printf("Please input the maximum # of Lanczos iterations  [%d]: ", nit);
    if (count_words(gets(str)) > 0) nit = atoi(strtok(str," \t\n\r\f"));
    if (nit < 1) break;

    printf("Please input the value of epsilon for delta-function [%g]: ", eps);
    if (count_words(gets(str)) > 0) eps = atof(strtok(str," \t\n\r\f"));
    if (eps <= 0.) break;
    
    // time info
    Timer *time = new Timer();
    printf("\nNow to compute the LDOS for atom %d by Real Space Greens function method ...", ik);

    // run real space green's function calculation
    Green *green = new Green(dynmat->nucell, dynmat->sysdim, nit, fmin, fmax, ndos, eps, Hessian, ik);
    delete green;

    printf("Done!\n");
    time->stop(); time->print(); delete time;
  }
  memory->destroy_2d_double_array(Hessian);

return;
}

/*------------------------------------------------------------------------------
 * Private method to evaluate the phonon dispersion curves
 *----------------------------------------------------------------------------*/
void Phonon::pdisp()
{
  // ask the output file name and write the header.
  char str[MAXLINE];
  printf("Please input the filename to output the dispersion data [pdisp.dat]:");
  if (count_words(gets(str)) < 1) strcpy(str, "pdisp.dat");
  char *fname = strtok(str," \t\n\r\f");

  FILE *fp = fopen(fname, "w"); fname = NULL;
  fprintf(fp,"# q     qr    freq\n");
  fprintf(fp,"# 2pi/L  2pi/L %s\n", dynmat->funit);

  // now the calculate the dispersion curve
  double qstr[3], qend[3], q[3], qinc[3], qr=0., dq;
  int nq = MAX(MAX(dynmat->nx,dynmat->ny),dynmat->nz)/2+1;
  qend[0] = qend[1] = qend[2] = 0.;

  double *egvs = new double [ndim];
  while (1){
    for (int i=0; i<3; i++) qstr[i] = qend[i];

    int quit = 0;
    printf("\nPlease input the start q-point in unit of B1->B3, q to exit [%g %g %g]: ", qstr[0], qstr[1], qstr[2]);
    int nr = count_words(gets(str));
    if ((nr == 1) && (strcmp(str,"q") == 0)) break;
    else if (nr >= 3){
      qstr[0] = atof(strtok(str," \t\n\r\f"));
      qstr[1] = atof(strtok(NULL," \t\n\r\f"));
      qstr[2] = atof(strtok(NULL," \t\n\r\f"));
    }

    do printf("Please input the end q-point in unit of B1->B3: ");
    while (count_words(gets(str)) < 3);
    qend[0] = atof(strtok(str," \t\n\r\f"));
    qend[1] = atof(strtok(NULL," \t\n\r\f"));
    qend[2] = atof(strtok(NULL," \t\n\r\f"));

    printf("Please input the # of points along the line [%d]: ", nq);
    if (count_words(gets(str)) > 0) nq = atoi(strtok(str," \t\n\r\f"));
    nq = MAX(nq,2);

    for (int i=0; i<3; i++) qinc[i] = (qend[i]-qstr[i])/double(nq-1);
    dq = sqrt(qinc[0]*qinc[0]+qinc[1]*qinc[1]+qinc[2]*qinc[2]);

    for (int i=0; i<3; i++) q[i] = qstr[i];
    for (int ii=0; ii<nq; ii++){
      dynmat->getDMq(q);
      dynmat->geteigen(egvs, 0);
      fprintf(fp,"%lg %lg %lg %lg ", q[0], q[1], q[2], qr);
      for (int i=0; i<ndim; i++) fprintf(fp," %lg", egvs[i]);
      fprintf(fp,"\n");

      for (int i=0; i<3; i++) q[i] += qinc[i];
      qr += dq;
    }
    qr -= dq;
  }
  fclose(fp);
  delete []egvs;

return;
}

/* ----------------------------------------------------------------------------
 * Private method to write out the dynamical matrix at selected q
 * ---------------------------------------------------------------------------- */
void Phonon::dmanyq()
{
  char str[MAXLINE];
  double q[3];
  do printf("Please input the q-point to output the dynamical matrix:");
  while (count_words(gets(str)) < 3);
  q[0] = atof(strtok(str," \t\n\r\f"));
  q[1] = atof(strtok(NULL," \t\n\r\f"));
  q[2] = atof(strtok(NULL," \t\n\r\f"));

  dynmat->getDMq(q);
  dynmat->writeDMq(q);

return;
}

/* ----------------------------------------------------------------------------
 * Private method to get the vibrational frequencies at selected q
 * ---------------------------------------------------------------------------- */
void Phonon::vfanyq()
{
  char str[MAXLINE];
  double q[3], egvs[ndim];
  
  do printf("Please input the q-point to output the dynamical matrix:");
  while (count_words(gets(str)) < 3);
  q[0] = atof(strtok(str," \t\n\r\f"));
  q[1] = atof(strtok(NULL," \t\n\r\f"));
  q[2] = atof(strtok(NULL," \t\n\r\f"));

  dynmat->getDMq(q);
  dynmat->geteigen(egvs, 0);
  printf("q-point: [%lg %lg %lg]\n", q[0], q[1], q[2]);
  printf("Vibrational frequencies at this q-point:\n");
  for (int i=0; i<ndim; i++) printf("%lg ", egvs[i]); printf("\n\n");
}

/* ----------------------------------------------------------------------------
 * Private method to get the dispersion-like data for dynamical matrix
 * ---------------------------------------------------------------------------- */
void Phonon::DMdisp()
{
  // ask the output file name and write the header.
  char str[MAXLINE];

  printf("Please input the filename to output the DM data [DMDisp.dat]: ");
  if (count_words(gets(str)) < 1) strcpy(str, "DMDisp.dat");
  char *fname = strtok(str," \t\n\r\f");

  FILE *fp = fopen(fname, "w"); fname = NULL;
  fprintf(fp,"# q     qr    D\n");

  // now the calculate the dispersion-like curve
  double qstr[3], qend[3], q[3], qinc[3], qr=0., dq;
  int nq = MAX(MAX(dynmat->nx,dynmat->ny),dynmat->nz)/2;
  qend[0] = qend[1] = qend[2] = 0.;

  while (1){

    for (int i=0; i<3; i++) qstr[i] = qend[i];

    printf("\nPlease input the start q-point in unit of B1->B3, q to exit [%g %g %g]: ", qstr[0], qstr[1], qstr[2]);
    int n = count_words(gets(str));
    if ((n == 1) && (strcmp(str,"q") == 0)) break;
    else if (n >= 3){
      qstr[0] = atof(strtok(str," \t\n\r\f"));
      qstr[1] = atof(strtok(NULL," \t\n\r\f"));
      qstr[2] = atof(strtok(NULL," \t\n\r\f"));
    }

    do printf("Please input the end q-point in unit of B1->B3: ");
    while (count_words(gets(str)) < 3);
    qend[0] = atof(strtok(str," \t\n\r\f"));
    qend[1] = atof(strtok(NULL," \t\n\r\f"));
    qend[2] = atof(strtok(NULL," \t\n\r\f"));

    printf("Please input the # of points along the line [%d]: ", nq);
    if (count_words(gets(str)) > 0) nq = atoi(strtok(str," \t\n\r\f"));
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
void Phonon::smooth(double *array, int npt, double xmin, double delta)
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

/* ----------------------------------------------------------------------------
 * Private method to calculate the thermal properties
 * ---------------------------------------------------------------------------- */
void Phonon::therm()
{
  // get frequencies on a q-mesh
  QMesh();
  ComputeAll();
 
  // get the filename to output thermal properties
  char str[MAXLINE];

  printf("\nPlease input the filename to output thermal properties [therm.dat]:");
  if (strlen(gets(str)) < 1) strcpy(str, "therm.dat");
  char *fname = strtok(str," \t\n\r\f");
  FILE *fp = fopen(fname, "a"); fname = NULL;
  // header line 
  fprintf(fp,"#Temp   Uvib    Svib     Fvib    ZPE      Cvib\n");
  fprintf(fp,"# K      eV      Kb       eV      eV       Kb\n");

  // constants          J.s             J/K                J
  const double h = 6.62606896e-34, Kb = 1.380658e-23, eV = 1.60217733e-19;

  // first temperature
  double T = dynmat->Tmeasure;
  do {
    // constants under the same temperature; assuming angular frequency in THz
    double h_o_KbT = h/(Kb*T)*1.e12, KbT_in_eV = Kb*T/eV;

    double Uvib = 0., Svib = 0., Fvib = 0., Cvib = 0., ZPE = 0.;
    for (int iq=0; iq<nq; iq++){
      double Utmp = 0., Stmp = 0., Ftmp = 0., Ztmp = 0., Ctmp = 0.;
      for (int i=0; i<ndim; i++){
        if (eigs[iq][i] <= 0.) continue;
        double x = eigs[iq][i] * h_o_KbT;
        double expterm = 1./(exp(x)-1.);
        Stmp += x*expterm - log(1.-exp(-x));
        Utmp += (0.5+expterm)*x;
        Ftmp += log(2.*sinh(0.5*x));
        Ctmp += x*x*exp(x)*expterm*expterm;
        Ztmp += 0.5*h*eigs[iq][i];
      }
      Svib += wt[iq]*Stmp;
      Uvib += wt[iq]*Utmp;
      Fvib += wt[iq]*Ftmp;
      Cvib += wt[iq]*Ctmp;
      ZPE  += wt[iq]*Ztmp;
    }
    Uvib *= KbT_in_eV;
    Fvib *= KbT_in_eV;
    ZPE  /= eV*1.e-12;
    // output result under current temperature
    fprintf(fp,"%lg %lg %lg %lg %lg %lg\n", T, Uvib, Svib, Fvib, ZPE, Cvib);

    printf("Please input the desired temperature (K), enter to exit: ");
    if (count_words(gets(str)) < 1) break;
    T = atof(strtok(str," \t\n\r\f"));
  } while (T > 0.);
  fclose(fp);

return;
}

/* ----------------------------------------------------------------------------
 * Private method to generate the q-points from a uniform q-mesh
 * ---------------------------------------------------------------------------- */
void Phonon::QMesh()
{
  // ask for mesh info
  char str[MAXLINE];
  int nx = dynmat->nx, ny = dynmat->ny, nz = dynmat->nz;
  printf("\nThe q-mesh size from the read dynamical matrix is: %d x %d x %d\n", nx, ny, nz);
  printf("Please input your desired q-mesh size [%d %d %d]: ", nx, ny, nz);
  if (count_words(gets(str)) >= 3){
    nx = atoi(strtok(str," \t\n\r\f"));
    ny = atoi(strtok(NULL," \t\n\r\f"));
    nz = atoi(strtok(NULL," \t\n\r\f"));
  }
  if (nx<1||ny<1||nz<1) return;
  if (dynmat->nx == 1) nx = 1;
  if (dynmat->ny == 1) ny = 1;
  if (dynmat->nz == 1) nz = 1;

#ifdef UseSPG
  // ask method to generate q-points
  int method = 2;
  printf("Please select your method to generate the q-points:\n");
  printf("  1. uniform;\n  2. Monkhost-Pack mesh;\nYour choice [2]: ");
  if (count_words(gets(str)) > 0) method = atoi(strtok(str," \t\n\r\f"));
  method = 2-method%2;
#endif
 
  if (qpts) memory->destroy_2d_double_array(qpts);
  if (wt)  delete []wt;

#ifdef UseSPG
  if (method == 1){
#endif
    nq = nx*ny*nz;
    double w = 1./double(nq);
    qpts = memory->create_2d_double_array(nq, 3, "QMesh_qpts");
    wt = new double [nq];
    int iq = 0;
    for (int i=0; i<nx; i++)
    for (int j=0; j<ny; j++)
    for (int k=0; k<nz; k++){
      qpts[iq][0] = double(i)/double(nx);
      qpts[iq][1] = double(j)/double(ny);
      qpts[iq][2] = double(k)/double(nz);
      wt[iq++] = w;
    }
  }
#ifdef UseSPG
  if ((method == 2) && (atpos == NULL)){
    atpos = memory->create_2d_double_array(dynmat->nucell,3,"phonon_ldos_egv:atpos");
    attyp = new int [dynmat->nucell];
    for (int i=0; i<dynmat->nucell; i++)
    for (int idim=0; idim<3; idim++) atpos[i][idim] = 0.;
    for (int i=0; i<3; i++) latvec[i][i] = 1.;

    int flag_lat_info_read = dynmat->flag_latinfo;

    if ( flag_lat_info_read ){ // get unit cell info from binary file; done by dynmat
      num_atom = dynmat->nucell;
      // set default, in case system dimension under study is not 3.
      for (int i=0; i<dynmat->nucell; i++)
      for (int idim=0; idim<3; idim++) atpos[i][idim] = 0.;
      //for (int i=0; i<3; i++) latvec[i][i] = 1.;

      // get atomic type info
      for (int i=0; i<num_atom; i++) attyp[i] = dynmat->attyp[i];
      // get unit cell vector info
      int ndim = 0;
      for (int idim=0; idim<3; idim++)
      for (int jdim=0; jdim<3; jdim++) latvec[jdim][idim] = dynmat->basevec[ndim++];
      // get atom position in unit cell; fractional
      for (int i=0; i<num_atom; i++)
      for (int idim=0; idim<sysdim; idim++) atpos[i][idim] = dynmat->basis[i][idim];

      // display the unit cell info read
      printf("\n");for (int ii=0; ii<60; ii++) printf("="); printf("\n");
      printf("The basis vectors of the unit cell:\n");
      for (int idim=0; idim<3; idim++) printf("  A%d = %lg %lg %lg\n", idim+1, latvec[0][idim], latvec[1][idim], latvec[2][idim]);
      printf("Atom(s) in the unit cell:\n");
      printf("  No.  type  sx  sy sz\n");
      for (int i=0; i<num_atom; i++) printf("  %d %d %lg %lg %lg\n", i+1, attyp[i], atpos[i][0], atpos[i][1], atpos[i][2]);
      printf("\nIs the above info correct? (y/n)[y]: ");
      if ( (strlen(gets(str)) > 0) && ( (strcmp(str,"y") != 0) && (strcmp(str,"Y") != 0)) ) flag_lat_info_read = 0;
    }

    if (flag_lat_info_read == 0) { // get unit cell info from file or user input
      int latsrc = 1;
      printf("\nPlease select the way to provide the unit cell info:\n");
      printf("  1. By file;\n  2. Read in.\nYour choice [1]: ");
      if (count_words(gets(str)) > 0) latsrc = atoi(strtok(str," \t\n\r\f"));
      latsrc = 2-latsrc%2;
      /*----------------------------------------------------------------
       * Ask for lattice info from the user; the format of the file is:
       * A1_x A1_y A1_z
       * A2_x A2_y A2_z
       * A3_x A3_y A3_z
       * natom
       * Type_1 sx_1 sy_1 sz_1
       * ...
       * Type_n sx_n sy_n sz_n
       *----------------------------------------------------------------*/
      if (latsrc == 1){ // to read unit cell info from file; get file name first
        do printf("Please input the file name containing the unit cell info: ");
        while (count_words(gets(str)) < 1);
        char *fname = strtok(str," \t\n\r\f");
        FILE *fp = fopen(fname,"r"); fname = NULL;
  
        if (fp == NULL) latsrc = 2;
        else {
          for (int i=0; i<3; i++){ // read unit cell vector info; # of atoms per unit cell
            if (fgets(str,MAXLINE,fp) == NULL) {latsrc = 2; break;}
            sscanf(str,"%lg %lg %lg", &latvec[0][i], &latvec[1][i], &latvec[2][i]);
          }
          if (fgets(str,MAXLINE,fp) == NULL) latsrc = 2;
          else {
            num_atom = atoi(strtok(str," \t\n\r\f"));
            if (num_atom > dynmat->nucell){
              printf("\nError: # of atoms read from file (%d) is bigger than that given by the dynamical matrix (%d)!\n", num_atom, dynmat->nucell);
              return;
            }
      
            for (int i=0; i<num_atom; i++){ // read atomic type and fractional positions
              if (fgets(str,MAXLINE,fp) == NULL) latsrc = 2;
              else sscanf(str,"%d %lg %lg %lg",&attyp[i],&atpos[i][0],&atpos[i][1],&atpos[i][2]);
            }
          }
        }
        fclose(fp);
      }
      if (latsrc == 2){
        for (int i=0; i<3; i++){
          do printf("Please input the vector A%d: ", i+1);
          while (count_words(gets(str)) < 3);
          latvec[0][i] = atof(strtok(str," \t\n\r\f"));
          latvec[1][i] = atof(strtok(NULL," \t\n\r\f"));
          latvec[2][i] = atof(strtok(NULL," \t\n\r\f"));
        }
  
        do printf("please input the number of atoms per unit cell: ");
        while (count_words(gets(str)) < 1);
        num_atom = atoi(strtok(str," \t\n\r\f"));
        if (num_atom > dynmat->nucell){
          printf("\nError: # of atoms input (%d) is bigger than that given by the dynamical matrix (%d)!\n", num_atom, dynmat->nucell);
          return;
        }
  
        for (int i=0; i<num_atom; i++){
          do printf("Please input the type, and fractional coordinate of atom No.%d: ", i+1);
          while (count_words(gets(str)) < 4);
          attyp[i] = atoi(strtok(str," \t\n\r\f"));
  
          atpos[i][0] = atof(strtok(NULL," \t\n\r\f"));
          atpos[i][1] = atof(strtok(NULL," \t\n\r\f"));
          atpos[i][2] = atof(strtok(NULL," \t\n\r\f"));
        }
      }
    } // end of read from file or input
  } // end of if (method == 2 && ...

  if (method == 2){
    int mesh[3], shift[3], is_time_reversal = 0;
    mesh[0] = nx; mesh[1] = ny; mesh[2] = nz;
    shift[0] = shift[1] = shift[2] = 0;
    int num_grid = mesh[0]*mesh[1]*mesh[2];
    int grid_point[num_grid][3], map[num_grid];
    double symprec = 1.e-4, pos[num_atom][3];

    for (int i=0; i<num_atom; i++)
    for (int j=0; j<3; j++) pos[i][j] = atpos[i][j];

    nq = spg_get_ir_reciprocal_mesh(grid_point, map, num_grid, mesh, shift,
                               is_time_reversal, latvec, pos, attyp, num_atom, symprec);
    qpts = memory->create_2d_double_array(nq,3,"qpts");
    wt = new double[nq];

    int *iq2idx = new int[num_grid];
    int numq = 0;
    for (int i=0; i<num_grid; i++){
      int iq = map[i];
      if (iq == i) iq2idx[iq] = numq++;
    }
    for (int iq=0; iq<nq; iq++) wt[iq] = 0.;
    numq = 0;
    for (int i=0; i<num_grid; i++){
      int iq = map[i];
      if (iq == i){
        qpts[numq][0] = double(grid_point[i][0])/double(mesh[0]);
        qpts[numq][1] = double(grid_point[i][1])/double(mesh[1]);
        qpts[numq][2] = double(grid_point[i][2])/double(mesh[2]);
        numq++;
      }
      wt[iq2idx[iq]] += 1.;
    }
    delete []iq2idx;

    double wsum = 0.;
    for (int iq=0; iq<nq; iq++) wsum += wt[iq];
    for (int iq=0; iq<nq; iq++) wt[iq] /= wsum;
    
  }
#endif
  printf("Your new q-mesh size would be: %d x %d x %d => %d points\n", nx,ny,nz,nq);

return;
}

/* ----------------------------------------------------------------------------
 * Private method to calculate the local phonon DOS and total phonon DOS based
 * on the eigenvectors
 * ---------------------------------------------------------------------------- */
void Phonon::ldos_egv()
{
  // get local position info
  char str[MAXLINE], *ptr;
  printf("\nThe # of atoms per cell is: %d, please input the atom IDs to compute\n", dynmat->nucell);
  printf("local PDOS, IDs begin with 0: ");
  int nmax = count_words(gets(str));
  if (nmax < 1) return;
  if (locals) delete []locals;
  locals = new int[nmax];

  nlocal = 0;
  ptr = strtok(str," \t\n\r\f");
  while (ptr != NULL){
    int id = atoi(ptr);
    if (id >= 0 && id < dynmat->nucell) locals[nlocal++] = id;

    ptr = strtok(NULL," \t\n\r\f");
  }
  if (nlocal < 1) return;

  printf("Local PDOS for atom(s):");
  for (int i=0; i<nlocal; i++) printf(" %d", locals[i]);
  printf("  will be computed.\n");

  fmin = 0.; fmax = 10.;
  printf("Please input the freqency (nv, THz) range to compute PDOS [%g %g]: ", fmin, fmax);
  if (count_words(gets(str)) >= 2) {
    fmin = atof(strtok(str," \t\n\r\f"));
    fmax = atof(strtok(NULL," \t\n\r\f"));
  }
  if (fmax < 0. || fmax < fmin) return;

  ndos = 201;
  printf("Please input your desired # of points in PDOS [%d]: ", ndos);
  if (count_words(gets(str)) > 0) ndos = atoi(strtok(str," \t\n\r\f"));
  if (ndos < 2) return;
  ndos += (ndos+1)%2;

  df = (fmax-fmin)/double(ndos-1);
  rdf = 1./df;

  // get the q-points
  QMesh();

  // allocate memory for DOS and LDOS
  if (dos) delete []dos;
  if (ldos) memory->destroy_3d_double_array(ldos);

  dos = new double[ndos];
  ldos = memory->create_3d_double_array(ndos,nlocal,sysdim,"phonon_ldos_egv:ldos");

  for (int i=0; i<ndos; i++){
    dos[i] = 0.;
    for (int j=0; j<nlocal; j++)
    for (int k=0; k<sysdim; k++) ldos[i][j][k] = 0.;
  }

  int nprint;
  if (nq > 10) nprint = nq/10;
  else nprint = 1;
  Timer *time = new Timer();

  // memory and pointer for eigenvalues and eigenvectors
  double egval[ndim], offset=fmin-0.5*df;
  doublecomplex **egvec = dynmat->DM_q;

  printf("\nNow to compute the phonons and DOSs "); fflush(stdout);
  for (int iq=0; iq<nq; iq++){
    if ((iq+1)%nprint == 0) {printf("."); fflush(stdout);}

    dynmat->getDMq(qpts[iq]);
    dynmat->geteigen(&egval[0], 1);

    for (int idim=0; idim<ndim; idim++){
      int hit = int((egval[idim] - offset)*rdf);
      if (hit >= 0 && hit <ndos){
        dos[hit] += wt[iq];

        for (int ilocal=0; ilocal<nlocal; ilocal++){
          int ipos = locals[ilocal]*sysdim;
          for (int k=0; k<sysdim; k++){
            double dr = egvec[idim][ipos+k].r, di = egvec[idim][ipos+k].i;
            double norm = dr * dr + di * di;
            ldos[hit][ilocal][k] += wt[iq] * norm;
          }
        }
      }
    }
  }
  egvec = NULL;
  printf("Done!\nNow to normalize the DOSs ...");

  // normalize the measure DOS and LDOS
  Normalize();
  printf("Done! ");
  time->stop(); time->print(); delete time;

  // to write the DOSes
  writeDOS();
  writeLDOS();

return;
}

/* ----------------------------------------------------------------------------
 * Private method to normalize the DOS and/or Local DOS.
 * Simpson's rule is used for the integration.
 * ---------------------------------------------------------------------------- */
void Phonon::Normalize()
{
  double odd, even, sum;
  if (dos){
    odd = even = 0.;
    for (int i=1; i<ndos-1; i +=2) odd  += dos[i];
    for (int i=2; i<ndos-1; i +=2) even += dos[i];
    sum = dos[0] + dos[ndos-1];
    sum += 4.*odd + 2.*even;
    sum = 3.*rdf/sum;
    for (int i=0; i<ndos; i++) dos[i] *= sum;
  }

  if (ldos){
    for (int ilocal=0; ilocal<nlocal; ilocal++)
    for (int idim=0; idim<sysdim; idim++){
      odd = even = 0.;
      for (int i=1; i<ndos-1; i +=2) odd  += ldos[i][ilocal][idim];
      for (int i=2; i<ndos-1; i +=2) even += ldos[i][ilocal][idim];
      sum = ldos[0][ilocal][idim] + ldos[ndos-1][ilocal][idim];
      sum += 4.*odd + 2.*even;
      sum = 3.*rdf/sum;
      for (int i=0; i<ndos; i++) ldos[i][ilocal][idim] *= sum;
    }
  }

return;
}

/* ----------------------------------------------------------------------------
 * Private method to calculate vibrational frequencies for all q-points
 * ---------------------------------------------------------------------------- */
void Phonon::ComputeAll()
{
  int nprint;
  if (nq > 10) nprint = nq/10;
  else nprint = 1;
  Timer *time = new Timer();

  printf("\nNow to compute the phonons "); fflush(stdout);
  // now to calculate the frequencies at all q-points
  if (eigs) memory->destroy_2d_double_array(eigs);
  eigs = memory->create_2d_double_array(nq,ndim,"QMesh_eigs");
  
  for (int iq=0; iq<nq; iq++){
    if ((iq+1)%nprint == 0) {printf("."); fflush(stdout);}

    dynmat->getDMq(qpts[iq]);
    dynmat->geteigen(eigs[iq], 0);
  }
  printf("Done!\n");
  time->stop(); time->print(); delete time;

return;
}

/*------------------------------------------------------------------------------
 * Method to count # of words in a string, without destroying the string
 *----------------------------------------------------------------------------*/
int Phonon::count_words(const char *line)
{
  int n = strlen(line) + 1;
  char *copy = (char *) memory->smalloc(n*sizeof(char),"copy");
  strcpy(copy,line);

  char *ptr;
  if (ptr = strchr(copy,'#')) *ptr = '\0';

  if (strtok(copy," \t\n\r\f") == NULL) {
    memory->sfree(copy);
    return 0;
  }
  n = 1;
  while (strtok(NULL," \t\n\r\f")) n++;

  memory->sfree(copy);
  return n;
}
