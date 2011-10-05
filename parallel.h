#ifndef PARALLEL_H
#define PARALLEL_H
#include "HashTable.h"


extern volatile int numberOfThreads; 
extern volatile int numberOfJobs;
extern volatile int jobCycle;
extern volatile int threadJob;
extern volatile char *barrierBuffer;

/* TODO shouldn't it be volative too? */
extern  pthread_mutex_t mutex;

#define THREAD_GET_EVENTS 1 
#define THREAD_COMBINE_EVENTS 2 
#define THREAD_MRE 3 

typedef struct _parArgs 
{
  HashTable *mergingHash;
  BitVector *candidateBips; 
  Array *bipartitionsById; 
  Array *bipartitionProfile; 
  int *indexByNumberBits;   
  boolean firstMerge; 
  Array *allDropsets;   
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
