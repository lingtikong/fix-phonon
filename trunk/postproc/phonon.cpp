#include "phonon.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "green.h"

#ifdef UseSPG
extern "C"{
#include "spglib.h"
}
#endif

#define MAX_STRING_LEN 256
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
  ndim = dynmat->fftdim;

  // display the menu
  while ( 1 ){
    printf("\n"); for (int i=0; i<27;i++) printf("="); printf(" Menu "); for (int i=0; i<27;i++) printf("="); printf("\n");
    printf("  1. To compute the phonon DOS;\n");
    printf("  2. To compute the phonon dispersion;\n");
    printf("  3. To get the DM at arbitrary q;\n");
    printf("  4. To get the frequencies at arbitrary q;\n");
    printf("  5. To get dispersion-like curve for dynamical matrix;\n");
    printf("  6. To calculate the vibrational thermal properties;\n");
    printf("  7. To compute LDOS from eigenvectors;\n");
    printf("  8. To compute LDOS by RSGF method;\n");
    printf("  9. Reset the interpolation method;\n");
    printf("  0. Exit.\n");
    // read user choice
    char str[MAX_STRING_LEN];
    job = 0;
    printf("Your choice [0]: ");
    if (strlen(gets(str)) > 0) job = atoi(str);
    printf("\nYou chose: %d\n", job);
    for (int i=0; i<60;i++) printf("=");printf("\n\n");

    // now to do the job according to user's choice
    if      (job == 1) dos();
    else if (job == 2) disp();
    else if (job == 3) dmanyq(); 
    else if (job == 4) vfanyq(); 
    else if (job == 5) DMdisp(); 
    else if (job == 6) therm(); 
    else if (job == 7) ldos_egv(); 
    else if (job == 8) ldos_rsgf(); 
    else if (job == 9) dynmat->getIntMeth();
    else break;
  }
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
  if (memory) delete memory;
}

/* ----------------------------------------------------------------------------
 * Private method to calculate the phonon DOS
 * ---------------------------------------------------------------------------- */
void Phonon::dos()
{
  // get frequencies on a q-mesh
  QMesh();
  ComputeAll();

  // now to get the frequency range
  char str[MAX_STRING_LEN];
  double fmin, fmax;

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
  if (strlen(gets(str)) > 0) sscanf(str,"%lg %lg", &fmin, &fmax);
  if (fmin > fmax){double swap = fmin; fmin = fmax; fmax = swap;}

  int ndos = 201;
  printf("Please input the number of intervals [%d]: ", ndos);
  if (strlen(gets(str)) > 0) ndos = atoi(str);
  ndos += (ndos+1)%2; ndos = MAX(2,ndos);

  double finc = (fmax-fmin)/double(ndos-1), finc_inv = 1./finc;;
  double *dos;
  dos = new double[ndos];
  for (int i=0; i<ndos; i++) dos[i] = 0.;

  // now to calculate the DOS
  for (int iq=0; iq<nq; iq++)
  for (int j=0; j<ndim; j++){
    int idx = int((eigs[iq][j]-fmin)*finc_inv+0.5);
    if (idx>=0 && idx<ndos) dos[idx] += wt[iq];
  }

  // smooth dos ?
  printf("Would you like to smooth the phonon dos? (y/n)[n]: ");
  if (strlen(gets(str)) > 0){
    if (strcmp(str,"y") == 0 || strcmp(str,"Y") == 0){
      smooth(dos, ndos, fmin, finc);
    }
  }

  // normalize dos to 1
  double odd = 0., even = 0.;
  for (int i=1; i<ndos-1; i +=2) odd  += dos[i];
  for (int i=2; i<ndos-1; i +=2) even += dos[i];
  double sum = dos[0] + dos[ndos-1];
  sum += 4.*odd + 2.*even;
  sum = 3.* finc_inv / sum;
  for (int i=0; i<ndos; i++) dos[i] *= sum;

  // now to output the phonon DOS
  if (outfile) delete []outfile;
  printf("Please input the filename to output DOS [pdos.dat]: ");
  if (strlen(gets(str)) > 0) {
    int n = strlen(str)+1;
    outfile = new char[n];
    strcpy(outfile, str);
  } else {
    outfile = new char[9];
    strcpy(outfile, "pdos.dat");
  }
  FILE *fp = fopen(outfile, "w");
  fprintf(fp,"# frequency  DOS\n");
  fprintf(fp,"#%s  number\n", dynmat->funit);
  double freq = fmin;
  for (int i=0; i<ndos; i++){
    fprintf(fp,"%lg %lg\n", freq, dos[i]);
    freq += finc;
  }
  fclose(fp);
  delete []outfile;
  delete []dos;

return;
}

/* ----------------------------------------------------------------------------
 * Private method to calculate the local phonon DOS via the real space Green's
 * function method
 * ---------------------------------------------------------------------------- */
void Phonon::ldos_rsgf()
{
  char str[MAX_STRING_LEN];
  double q0[3], egvs[ndim];

  q0[0] = q0[1] = q0[2] = 0.;
  dynmat->getDMq(q0);
  dynmat->geteigen(egvs);

  double vmin, vmax;
  vmin = vmax = egvs[0];
  for (int i=1; i<ndim; i++){vmin = MIN(vmin, egvs[i]); vmax = MAX(vmax, egvs[i]);}

  q0[0] = q0[1] = q0[2] = 0.;
  dynmat->getDMq(q0);
  const double tpi = 8.*atan(1.);
  double **Hessian, scale;
  scale = dynmat->eml2f*tpi; scale *= scale;

  Hessian = memory->create_2d_double_array(ndim, ndim, "dispersion_ldos_Hessian");
  for (int i=0; i<ndim; i++)
  for (int j=0; j<ndim; j++) Hessian[i][j] = dynmat->DM_q[i][j].r*scale;

  int ik = 0, ndos = 201, nit = MAX(ndim*0.1, 50);
  double eps = vmax*vmax*0.1;
  while (1) {
    printf("\nThere are %d atoms in each unit cell of your lattice.\n", dynmat->nucell);
    printf("Please input the index of atom in the unit cell to evaluate LDOS [%d]: ", ik);
    if (strlen(gets(str)) > 0) ik = atoi(str); 
    if (ik < 0 || ik >= dynmat->nucell) break;

    printf("Please input the frequency range to evaluate LDOS [%g %g]: ", vmin, vmax);
    if (strlen(gets(str)) > 0) sscanf(str,"%lg %lg", &vmin, &vmax); vmin = MAX(vmin,0.);
    if (vmax < vmin) break;

    printf("Please input the desired number of points in LDOS [%d]: ", ndos);
    if (strlen(gets(str)) > 0) ndos = atoi(str); ndos += (ndos+1)%2;
    if (ndos < 2) break;

    printf("Please input the maximum # of Lanczos iterations  [%d]: ", nit);
    if (strlen(gets(str)) > 0) nit = atoi(str);
    if (nit < 1) break;

    printf("Please input the value of epsilon for delta-function [%g]: ", eps);
    if (strlen(gets(str)) > 0) eps = atof(str);
    if (eps <= 0.) break;
    
    Green *green = new Green(dynmat->nucell, dynmat->sysdim, nit, vmin, vmax, ndos, eps, Hessian, ik);
    delete green;
  }

return;
}

void Phonon::disp()
{
  // ask the output file name and write the header.
  char str[MAX_STRING_LEN];
  printf("Please input the filename to output the dispersion data [pdisp.dat]:");
  if (outfile) delete []outfile;
  if (strlen(gets(str)) > 0){
    int n = strlen(str) + 1;
    outfile = new char[n];
    strcpy(outfile, str);
  } else {
    outfile = new char [10];
    strcpy(outfile, "pdisp.dat");
  }
  FILE *fp = fopen(outfile, "w"); delete []outfile;
  fprintf(fp,"# q     qr    freq\n");
  fprintf(fp,"# 2pi/L  2pi/L %s\n", dynmat->funit);

  // now the calculate the dispersion curve
  double qstr[3], qend[3], q[3], qinc[3], qr=0., dq;
  int nq = MAX(MAX(dynmat->nx,dynmat->ny),dynmat->nz)/2+1;
  qend[0] = qend[1] = qend[2] = 0.;

  while (1){
    double egvs[ndim];
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

return;
}

/* ----------------------------------------------------------------------------
 * Private method to calculate the phonon DOS
 * ---------------------------------------------------------------------------- */
void Phonon::dmanyq()
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

void Phonon::vfanyq()
{
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
void Phonon::DMdisp()
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
  char str[MAX_STRING_LEN];
  if (outfile) delete []outfile;

  printf("\nPlease input the filename to output thermal properties [therm.dat]:");
  if (strlen(gets(str)) > 0){
    int n = strlen(str) + 1;
    outfile = new char[n];
    strcpy(outfile,str);
  } else {
    outfile = new char[10];
    strcpy(outfile,"therm.dat");
  }

  FILE *fp = fopen(outfile, "w");
  // header line 
  fprintf(fp,"#Temp  Uvib Svib Fvib ZPE Cvib\n");
  fprintf(fp,"# K    eV   Kb   eV   eV  Kb\n");

  // constants          J.s             J/K                J
  const double h = 6.62606896e-34, Kb = 1.380658e-23, eV = 1.60217733e-19;

  // first temperature
  double T;
  while (1){
    printf("Please input the desired temperature (K), enter to exit: ");
    if (strlen(gets(str)) < 1) break;
    sscanf(str,"%lg",&T);
    if (T <= 0.) break;

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
  }
  fclose(fp);

return;
}

/* ----------------------------------------------------------------------------
 * Private method to calculate the phonon frequencies on a q-mesh
 * ---------------------------------------------------------------------------- */
void Phonon::QMesh()
{
  // ask for mesh info
  char str[MAX_STRING_LEN];
  int nx = dynmat->nx, ny = dynmat->ny, nz = dynmat->nz, nr=0;
  printf("\nThe q-mesh size from the read dynamical matrix is: %d x %d x %d\n", nx, ny, nz);
  printf("Please input your desired q-mesh size [%d %d %d]: ", nx, ny, nz);
  if (strlen(gets(str)) > 0) sscanf(str,"%d %d %d", &nx, &ny, &nz);
  if (nx<1||ny<1||nz<1) return;
  if (dynmat->nx == 1) nx = 1;
  if (dynmat->ny == 1) ny = 1;
  if (dynmat->nz == 1) nz = 1;

#ifdef UseSPG
  // ask method to generate q-points
  int method = 2;
  printf("Please select your method to generate the q-points:\n");
  printf("  1. uniform;\n  2. Monkhost-Pack mesh;\n");
  printf("Your choice [2]: ");
  if (strlen(gets(str)) > 0) sscanf(str,"%d", &method);
#endif
 
  if (qpts) memory->destroy_2d_double_array(qpts);
  if (wt)  delete []wt;

#ifdef UseSPG
  if (method == 1){
#endif
    nq = nx*ny*nz;
    qpts = memory->create_2d_double_array(nq, 3, "QMesh_qpts");
    wt = new double [nq];
    int iq = 0;
    for (int i=0; i<nx; i++)
    for (int j=0; j<ny; j++)
    for (int k=0; k<nz; k++){
      qpts[iq][0] = double(i)/double(nx);
      qpts[iq][1] = double(j)/double(ny);
      qpts[iq][2] = double(k)/double(nz);
      wt[iq++] = 1.;
    }
#ifdef UseSPG
  } else {
    double lattice[3][3];
    int num_atom, *types;
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
    printf("Please input the name of file containing the unit cell info: ");
    int n = strlen(gets(str))+1;
    char *fname = new char[n];
    strcpy(fname, str);
    FILE *fp = fopen(fname,"r");
    if (fp == NULL) return;

    for (int i=0; i<3; i++){
      if (fgets(str,MAX_STRING_LEN,fp) == NULL) return;
      sscanf(str,"%lg %lg %lg", &lattice[i][0], &lattice[i][1], &lattice[i][2]);
    }
    if (fgets(str,MAX_STRING_LEN,fp) == NULL) return;
    sscanf(str,"%d",&num_atom);
    if (num_atom < 1) return;
    types = new int[num_atom];
    double position[num_atom][3];
    
    for (int i=0; i<num_atom; i++){
      if (fgets(str,MAX_STRING_LEN,fp) == NULL) return;
      sscanf(str,"%d %lg %lg %lg",&types[i],&position[i][0],&position[i][1],&position[i][2]);
    }
    fclose(fp);
    if (fname) delete []fname;

    int mesh[3], shift[3], is_time_reversal = 0;
    mesh[0] = nx; mesh[1] = ny; mesh[2] = nz;
    shift[0] = shift[1] = shift[2] = 0;
    int num_grid = mesh[0]*mesh[1]*mesh[2];
    int grid_point[num_grid][3], map[num_grid];
    double symprec = 1.e-5;

    nq = spg_get_ir_reciprocal_mesh(grid_point, map, num_grid,
                               mesh, shift, is_time_reversal,
                               lattice, position, types,
                               num_atom, symprec);
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
    
    if (types) delete []types;
  }
#endif
  printf("Your new q-mesh size would be: %d x %d x %d => %d points\n", nx,ny,nz,nq);

return;
}

void Phonon::ldos_egv()
{

}

void Phonon::ComputeAll()
{
  // now to calculate the frequencies at all q-points
  if (eigs) memory->destroy_2d_double_array(eigs);
  eigs = memory->create_2d_double_array(nq,ndim,"QMesh_eigs");
  
  for (int iq=0; iq<nq; iq++){
    dynmat->getDMq(qpts[iq]);
    dynmat->geteigen(eigs[iq]);
  }

return;
}
