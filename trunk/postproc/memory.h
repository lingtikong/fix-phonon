#ifndef MEMORY_H
#define MEMORY_H

#include "stdio.h"
#include "stdlib.h"
extern "C"{
#include "f2c.h"
#include "clapack.h"
}

class Memory {
public:
  double **create_2d_double_array(int, int, const char *);
  void destroy_2d_double_array(double **);

  double ***create_3d_double_array(int, int, int, const char *);
  void destroy_3d_double_array(double ***);

  doublecomplex **create_2d_complex_array(int, int, const char *);
  void destroy_2d_complex_array(doublecomplex **);
  void *smalloc(int, const char *);
  void sfree(void *);
};

#endif
