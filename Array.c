
#include "Array.h"

void freeArray(Array *array)
{
  free(array->arrayTable);
  free(array);
}

