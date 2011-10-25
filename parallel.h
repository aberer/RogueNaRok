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

#ifndef PARALLEL_H
#define PARALLEL_H
#include "HashTable.h"


extern volatile int numberOfThreads; 
extern volatile int numberOfJobs;
extern volatile int jobCycle;
extern volatile int threadJob;
extern volatile char *barrierBuffer;
extern  pthread_mutex_t mutex;

#define THREAD_GET_EVENTS 1 
#define THREAD_COMBINE_EVENTS 2 
#define THREAD_MRE 3 
#define THREAD_EVALUATE_EVENTS 4

typedef struct _parArgs 
{
  HashTable *mergingHash;
  BitVector *candidateBips; 
  Array *bipartitionsById; 
  Array *bipartitionProfile; 
  int *indexByNumberBits;   
  boolean firstMerge; 
  Array *allDropsets;   
  List *consensusBipsCanVanish;
} parallelArguments ; 


parallelArguments *globalPArgs; 

typedef struct
{
  parallelArguments *pArgs;
  int threadNumber;
} threadData;


void masterBarrier(int jobType, parallelArguments *pArgs);
void startThreads();
void *workerThreadWait(void *tData);
void execFunction(parallelArguments *pArgs, int tid, int n);


#endif
