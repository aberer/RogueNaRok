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


#ifndef BITVECTOR_H
#define BITVECTOR_H

#include <stdlib.h>
#include <assert.h>
#include <string.h>

#include "common.h"

typedef unsigned int BitVector; 

#define BIT_COUNT(x)  precomputed16_bitcount(x)
#define NUMBER_BITS_IN_COMPLEMENT(bipartition) (mxtips - dropRound - bipartition->numberOfBitsSet)
#define GET_BITVECTOR_LENGTH(x) (((x) % MASK_LENGTH) ? ((x) / MASK_LENGTH + 1) : ((x) / MASK_LENGTH))
#define FLIP_NTH_BIT(bitVector,n) (bitVector[(n) / MASK_LENGTH] |= mask32[ (n) % MASK_LENGTH ])
#define UNFLIP_NTH_BIT(bitVector,n) (bitVector[(n) / MASK_LENGTH] &= ~mask32[ (n) % MASK_LENGTH ])
#define NTH_BIT_IS_SET(bitVector,n) (bitVector[(n) / MASK_LENGTH] & mask32[(n) % MASK_LENGTH])
#define NTH_BIT_IS_SET_IN_INT(integer,n) (integer & mask32[n])
#define MASK_LENGTH 32

char bits_in_16bits [0x1u << 16];
extern BitVector *mask32;

void initializeMask();
BitVector genericBitCount(BitVector* bitVector, int bitVectorLength);
BitVector precomputed16_bitcount (BitVector n);
void compute_bits_in_16bits(void);
void printBitVector(BitVector *bv, int length);
void freeBitVectors(BitVector **v, int n);
BitVector *copyBitVector(BitVector *bitVector, int bitVectorLength);
void printBitVector(BitVector *bv, int length);

#endif



