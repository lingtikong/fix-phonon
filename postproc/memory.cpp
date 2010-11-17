#include "memory.h"

/* ----------------------------------------------------------------------------
 * method to allocate/free memory in a safe mode
 * ---------------------------------------------------------------------------- */
void *Memory::smalloc(int n, const char *name)
{
  if (n == 0) return NULL;
  void *ptr = malloc(n);
  if (ptr == NULL) {
    printf("Failed to allocate %d bytes for array %s\n",n,name);
    exit(1);
  }
  return ptr;
}

void Memory::sfree(void *ptr)
{
  if (ptr == NULL) return;
  free(ptr);
}

/* ----------------------------------------------------------------------------
 * method to create 2d double array
 * ---------------------------------------------------------------------------- */
double **Memory::create_2d_double_array(int n1, int n2, const char *name)
{
  double *data   = (double  *) smalloc(n1*n2*sizeof(double),name);
  double **array = (double **) smalloc(n1*sizeof(double *),name);

  int n = 0;
  for (int i = 0; i < n1; i++) {
    array[i] = &data[n];
    n += n2;
  }
return array;
}

/* ----------------------------------------------------------------------------
 * method to destroy 2d double array
 * ---------------------------------------------------------------------------- */
void Memory::destroy_2d_double_array(double **array)
{
  if (array == NULL) return;
  sfree(array[0]);
  sfree(array);
}

/* ----------------------------------------------------------------------
   create a 3d double array 
------------------------------------------------------------------------- */
double ***Memory::create_3d_double_array(int n1, int n2, int n3, const char *name)
{
  int i,j;
  if (n1 < 1 || n2 < 1 || n3 < 1) return NULL;

  double *data = (double *) smalloc(n1*n2*n3*sizeof(double),name);
  double **plane = (double **) smalloc(n1*n2*sizeof(double *),name);
  double ***array = (double ***) smalloc(n1*sizeof(double **),name);

  int n = 0;
  for (i = 0; i < n1; i++) {
    array[i] = &plane[i*n2];
    for (j = 0; j < n2; j++) {
      plane[i*n2+j] = &data[n];
      n += n3;
    }
  }

  return array;
}

/* ----------------------------------------------------------------------
   free a 3d double array 
------------------------------------------------------------------------- */
void Memory::destroy_3d_double_array(double ***array)
{
  if (array == NULL) return;
  sfree(array[0][0]);
  sfree(array[0]);
  sfree(array);
}
/* ----------------------------------------------------------------------------
 * method to create 2d complex array
 * ---------------------------------------------------------------------------- */
doublecomplex **Memory::create_2d_complex_array(int n1, int n2, const char *name)
{
  doublecomplex *data   = (doublecomplex  *) smalloc(n1*n2*sizeof(doublecomplex),name);
  doublecomplex **array = (doublecomplex **) smalloc(n1*sizeof(doublecomplex *),name);

  int n = 0;
  for (int i = 0; i < n1; i++) {
    array[i] = &data[n];
    n += n2;
  }
return array;
}

/* ----------------------------------------------------------------------------
 * method to destroy 2d complex array
 * ---------------------------------------------------------------------------- */
void Memory::destroy_2d_complex_array(doublecomplex **array)
{
  if (array == NULL) return;
  sfree(array[0]);
  sfree(array);
}
