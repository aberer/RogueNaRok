/*  RogueNaRok is an algorithm for the identification of rogue taxa in a set of phylogenetic trees. 
 *
 *  Moreover, the program collection comes with efficient implementations of 
 *   * the unrooted leaf stability by Thorley and Wilkinson
 *   * the taxonomic instability index by Maddinson and Maddison
 *   * a maximum agreement subtree implementation (MAST) for unrooted trees 
 *   * a tool for pruning taxa from a tree collection. 
 * 
 *  Copyright October 2011 by Andre J. Aberer
 * 
 *  Tree I/O and parallel framework are derived from RAxML by Alexandros Stamatakis.
 *
 *  This program is free software; you may redistribute it and/or
 *  modify its under the terms of the GNU General Public License as
 *  published by the Free Software Foundation; either version 2 of the
 *  License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  General Public License for more details.
 *
 *  For any other inquiries send an Email to Andre J. Aberer
 *  andre.aberer at googlemail.com
 * 
 *  When publishing work that is based on the results from RogueNaRok, please cite:
 *  Andre J. Aberer, Denis Krompaß, Alexandros Stamatakis. RogueNaRok: an Efficient and Exact Algorithm for Rogue Taxon Identification. (unpublished) 2011. 
 * 
 */

#ifndef PROFILEELEM_H
#define PROFILEELEM_H

#include <string.h>
#include <assert.h>

#include "Array.h"
#include "HashTable.h"
#include "common.h"
#include "BitVector.h"


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

int *createNumBitIndex(Array *bipartitionProfile, int mxtips);
int sortById(const void *a, const void *b);
int sortBySupport(const void *a, const void *b);
int sortBipProfile(const void *a, const void *b);
Array *cloneProfileArrayFlat(const Array *array);
void addElemToArray(ProfileElem *elem, Array *array);
void freeProfileElem(ProfileElem *elem);
#endif
