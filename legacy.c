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

#include "legacy.h"


void freeHashTable(hashtable *h)
{
  unsigned int
    i,
    entryCount = 0;   

  for(i = 0; i < h->tableSize; i++)
    {
      if(h->table[i] != NULL)
	{
	  entry *e = h->table[i];
	  entry *previous;	 

	  do
	    {
	      previous = e;
	      e = e->next;

	      if(previous->bitVector)
		free(previous->bitVector);

	      if(previous->treeVector)
		free(previous->treeVector);

	      if(previous->supportVector)
		free(previous->supportVector);
	      
	      free(previous);	      
	      entryCount++;
	    }
	  while(e != NULL);	  
	}

    }

  assert(entryCount == h->entryCount);
 
  free(h->table);
}

hashtable *initHashTable(unsigned int n)
{
  static const  unsigned int initTable[] = {64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384, 32768, 65536, 131072, 262144, 524288, 1048576, 2097152, 4194304, 8388608, 16777216, 33554432, 67108864, 134217728, 268435456, 536870912, 1073741824, 2147483648U};
  
  hashtable *h = (hashtable*)CALLOC(1,sizeof(hashtable));
  
  unsigned int
    tableSize,
    i,

#ifndef NDEBUG
    maxSize = (unsigned int)-1,
#endif
    primeTableLength = sizeof(initTable)/sizeof(initTable[0]);

#ifndef NDEBUG
  assert(n <= maxSize);
#endif

  i = 0;

  while(initTable[i] < n && i < primeTableLength)
    i++;
  
  assert(i < primeTableLength);

  tableSize = initTable[i];

  /* printf("Hash table init with size %u\n", tableSize); */

  h->table = (entry**)CALLOC(tableSize, sizeof(entry*));
  h->tableSize = tableSize;  
  h->entryCount = 0;  

  return h;
}


BitVector **initBitVector(All *tr, BitVector *vectorLength)
{
  BitVector **bitVectors = (BitVector **)CALLOC(2 * tr->mxtips, sizeof(BitVector*));
  int i;

  if(tr->mxtips % MASK_LENGTH == 0)
    *vectorLength = tr->mxtips / MASK_LENGTH;
  else
    *vectorLength = 1 + (tr->mxtips / MASK_LENGTH); 
  
  for(i = 1; i <= tr->mxtips; i++)
    {
      bitVectors[i] = (BitVector *)CALLOC(*vectorLength, sizeof(BitVector));
      bitVectors[i][(i - 1) / MASK_LENGTH] |= mask32[(i - 1) % MASK_LENGTH];
    }
  
  for(i = tr->mxtips + 1; i < 2 * tr->mxtips; i++) 
    bitVectors[i] = (BitVector *)CALLOC(*vectorLength, sizeof(BitVector));

  return bitVectors;
}


ProfileElem *addProfileElem(entry *helem, int vectorLength, int treeVectorLength, int numberOfTrees) 
{
  ProfileElem *result = CALLOC(1,sizeof(ProfileElem));
  result->isInMLTree = FALSE; 
  result->bitVector = CALLOC(vectorLength, sizeof(BitVector));
  result->treeVector = CALLOC(treeVectorLength, sizeof(BitVector));
  result->bitVector = memcpy(result->bitVector, helem->bitVector, vectorLength * sizeof(BitVector));
  result->treeVector = memcpy(result->treeVector, helem->treeVector, treeVectorLength * sizeof(BitVector));

  if(NTH_BIT_IS_SET(result->treeVector, numberOfTrees))
    {
      result->isInMLTree = TRUE;
      UNFLIP_NTH_BIT(result->treeVector, numberOfTrees);
    }
  
  result->treeVectorSupport = genericBitCount(result->treeVector, treeVectorLength);

  return result; 
}
