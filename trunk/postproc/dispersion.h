#include "lmp_phonon.h"
#include <complex>

using namespace std;

class Dispersion{
public:
  Dispersion(LMP_PHONON *);
  ~Dispersion();

  LMP_PHONON *dynmat;

private:
  int job;
  char *outfile;

  void dos();
  void disp();
  void dmanyq();
  void vfanyq();
  
  void DMdisp();

  void smooth(double *, int, double, double);
};
