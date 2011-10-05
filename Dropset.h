#ifndef DROPSET_H
#define DROPSET_H

#include "common.h"
#include "List.h"
#include "ProfileElem.h"
#include "HashTable.h"

typedef union _mergeBips
{
  int pair[2];
  IndexList *many; 
} MergingBipartitions;

typedef struct _mergingEvent
{
  MergingBipartitions mergingBipartitions; 
  boolean isComplex;
  int supportLost; 
  int supportGained; 
  boolean computed; 
}  MergingEvent;


typedef struct dropset
{
  IndexList *taxaToDrop;
  int improvement;
  List *primeEvents;
  List *combinedEvents;  
/* #ifdef PARALLEL */
/*   pthread_mutex_t *lock;  */
/* #endif */
} Dropset;


boolean dropsetEqual(HashTable *hashtable, void *entryA, void *entryB);
unsigned int dropsetHashValue(HashTable *hashTable, void *value);
void removeMergingEvent(Dropset *dropset, List *toBeRemoved);  
List *freeMergingEventReturnNext(List *elem); 
void removeDropsetAndRelated(HashTable *mergingHash, Dropset *dropset);
void addEventToDropsetPrime(Dropset *dropset, int a, int b);
List *addEventToDropsetCombining(List *complexEvents, MergingBipartitions primeEvent);
void freeDropsetDeepInHash(void *value);
void addEventToDropsetForCombining(Dropset *dropset, IndexList *mergingBips);
void initializeRandForTaxa(int mxtips);
void freeDropsetDeep(void *values, boolean freeCombinedM);
IndexList *getDropset(ProfileElem *elemA, ProfileElem *elemB, boolean complement, BitVector *neglectThose);
#endif
