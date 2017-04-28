#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#define ALIGNMENT 2*1024*1024

double * alloc(int *len) {
  assert(*len > 0);
  double * p = (double *)aligned_alloc(ALIGNMENT, sizeof(double)*(*len));
  assert(p != NULL);
  return p;
}

void alloc_free(double ** p)  {
  printf("About to free %p\n", *p);
  free(*p);
}

