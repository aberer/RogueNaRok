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



