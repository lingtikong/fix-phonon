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
  doublecomplex **create_2d_complex_array(int, int, const char *);
  void destroy_2d_complex_array(doublecomplex **);
  void *smalloc(int, const char *);
  void sfree(void *);
};

#endif
