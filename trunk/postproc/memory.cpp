#include "memory.h"

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
