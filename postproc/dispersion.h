#include "lmp_phonon.h"
#include <complex>
#include "memory.h"

using namespace std;

class Dispersion{
public:
  Dispersion(LMP_PHONON *);
  ~Dispersion();

  LMP_PHONON *dynmat;

private:
  int job;
  char *outfile;

  int nq;
  double **qpts, *wt;

  Memory *memory;

  void dos();
  void disp();
  void dmanyq();
  void vfanyq();
  
  void DMdisp();

  void smooth(double *, int, double, double);
};
