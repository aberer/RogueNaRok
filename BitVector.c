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

#include "BitVector.h"

BitVector *mask32;

void initializeMask()
{
  int i;
  mask32 = CALLOC(MASK_LENGTH, sizeof(BitVector));
  mask32[0] = 1; 
  
  for(i = 1; i < MASK_LENGTH; ++i)
    mask32[i] = mask32[i-1] << 1; 
}


void printBitVector(BitVector *bv, int length)
{
  int i ;  
  for(i = 0; i < length * 32; ++i)
    printf("%d", NTH_BIT_IS_SET(bv, i) ? 1 : 0);
}


void freeBitVectors(unsigned int **v, int n)
{
  int i;

  for(i = 1; i < n; i++)
    free(v[i]);
}


boolean areSameBitVectors(BitVector *a, BitVector *b, int bitVectorLength)
{
  int 
    i; 

  FOR_0_LIMIT(i,bitVectorLength)
    if(a[i] != b[i])
      return FALSE;

  return TRUE;
}


BitVector genericBitCount(BitVector* bitVector, int bitVectorLength)
{
  BitVector 
    i, 
    result = 0;

  for(i = 0; i < bitVectorLength; i++)
    result += BIT_COUNT(bitVector[i]);
  
  return result; 
}


static int iterated_bitcount(BitVector n)
{
    int 
      count=0;    
    
    while(n)
      {
        count += n & 0x1u ;    
        n >>= 1 ;
      }
    
    return count;
}

void compute_bits_in_16bits(void)
{
    BitVector i;    
    
    for (i = 0; i < (0x1u<<16); i++)
        bits_in_16bits[i] = iterated_bitcount(i);
    
    return ;
}

BitVector precomputed16_bitcount (BitVector n)
{
  /* works only for 32-bit int*/
    return bits_in_16bits [n         & 0xffffu]
        +  bits_in_16bits [(n >> 16) & 0xffffu] ;
}



BitVector *copyBitVector(BitVector *bitVector, int bitVectorLength)
{
  BitVector *result = CALLOC(bitVectorLength, sizeof(BitVector));
  memcpy(result, bitVector,bitVectorLength * sizeof(BitVector)); 
  return result;
}



