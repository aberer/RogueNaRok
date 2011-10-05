#ifndef ARRAY_H
#define ARRAY_H

#include <stdlib.h>

typedef struct 
{
  void *arrayTable;
  void *commonAttributes; 
  unsigned int length;
} Array;


#define CLONE_ARRAY_FLAT(FUNCNAME, TYPE, TYPEATTR)                                                \
Array *FUNCNAME (const Array *array)                                                                    \
{                                                                                                 \
    Array *result = CALLOC(1,sizeof(Array));                                                      \
    result->length = array->length;                                                               \
    result->arrayTable = CALLOC(result->length, sizeof(TYPE));                             \
    memcpy(result->arrayTable, array->arrayTable, array->length * sizeof(TYPE) );                 \
												  \
    if( array->commonAttributes )                                                                 \
    {                                                                                             \
	result->commonAttributes = CALLOC(1 , sizeof(TYPEATTR));                                  \
        memcpy(result->commonAttributes, array->commonAttributes, sizeof(TYPEATTR) ) ;            \
    }                                                                                             \
  return result;                                                                                  \
}

void freeArray(Array *array);


#endif
