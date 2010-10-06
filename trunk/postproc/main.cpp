#include "stdio.h"
#include "stdlib.h"
#include "lmp_phonon.h"
#include "dispersion.h"

using namespace std;

int main(int argc, char** argv)
{

  LMP_PHONON *dynmat = new LMP_PHONON(argc, argv);
  Dispersion *driver = new Dispersion(dynmat);

  delete driver;
  delete dynmat;

return 0;
}
