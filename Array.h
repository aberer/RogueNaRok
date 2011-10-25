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
