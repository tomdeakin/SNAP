#include <stdlib.h>
#include <stdio.h>

#define ALIGNMENT 2*1024*1024

double * alloc(int *len) {
  double * p = (double *)aligned_alloc(ALIGNMENT, sizeof(double)*(*len));
  printf("Allocated %p size %zu. len was %d\n", p, sizeof(double)*(*len), (*len));
  return p;
}

void alloc_free(double ** p)  {
  printf("About to free %p\n", *p);
  free(*p);
}

