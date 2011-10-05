#ifndef PROFILEELEM_H
#define PROFILEELEM_H

#include <string.h>
#include <assert.h>

#include "Array.h"
#include "HashTable.h"
#include "common.h"
#include "BitVector.h"
/* #include "Tree.h" */


typedef struct 
{
  BitVector bitVectorLength; 
  BitVector treeVectorLength;  
  BitVector *randForTaxa;	/* random numbers to hash the vectors */
  BitVector lastByte;		/* the padding bits */
} ProfileElemAttr;


typedef struct profile_elem
{
  BitVector *bitVector;
  BitVector *treeVector;
  int treeVectorSupport;
  boolean isInMLTree;
  BitVector id;
  int numberOfBitsSet;
} ProfileElem;

#define GET_PROFILE_ELEM(array,index) (((ProfileElem**)array->arrayTable)[(index)])
#define GET_DROPSET_ELEM(array,index) (((Dropset**)array->arrayTable)[(index)])

int *createNumBitIndex(Array *bipartitionProfile, int mxtips, boolean useInsertionSort);
int sortById(const void *a, const void *b);
int sortBySupport(const void *a, const void *b);
int sortBipProfile(const void *a, const void *b);
Array *cloneProfileArrayFlat(const Array *array);
void addElemToArray(ProfileElem *elem, Array *array);
void freeProfileElem(ProfileElem *elem);
#endif
