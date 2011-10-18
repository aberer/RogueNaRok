#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <unistd.h>
#include <limits.h>

#include "Tree.h"
#include "sharedVariables.h"
#include "Dropset.h"
#include "legacy.h"
#include "newFunctions.h"
#include "Node.h"

#ifdef PARALLEL
#include "parallel.h"
#include <pthread.h> 
#endif

#define PROG_NAME "RogueNaRok"
#define PROG_VERSION "1.0"
#define PROG_RELEASE_DATE "2011-10-12"

/* #define PRINT_VERY_VERBOSE */
/* #define MYDEBUG */

#define PRINT_DROPSETS
/* #define PRINT_TIME */

/* try to produce minimal dropsets */
/* #define MIN_DROPSETS */

#define VANILLA_CONSENSUS_OPT  0
#define ML_TREE_OPT 1
#define MRE_CONSENSUS_OPT 2

#define HASH_TABLE_SIZE_CONST 100

extern unsigned int *randForTaxa;

int bitVectorLength,
  treeVectorLength,
  maxDropsetSize = 1, 
  rogueMode = 0,
  dropRound = 0, 
  taxaDropped = 0, 
  thresh, 
  numberOfTrees,   
  bestLastTime, 
  *cumScores, 
  cumScore = 0, 
  bestCumEver = 0, 

  numBips,
  mxtips; 

Dropset **dropsetPerRound; 

boolean computeSupport = TRUE;

BitVector *droppedTaxa,
  *neglectThose, 
  *paddingBits; 

double labelPenalty = 0., 
  timeInc; 

#ifdef MYDEBUG
void debug_dropsetConsistencyCheck(HashTable *mergingHash)
{
  HashTableIterator *htIter; 
  FOR_HASH(htIter, mergingHash)
    {
      Dropset *ds = getCurrentValueFromHashTableIterator(htIter);
      
      if( NOT ds)
	break;
      
      HashTableIterator *htIter2;
      boolean hasNext2 = TRUE;
      for(htIter2 = createHashTableIterator(mergingHash); htIter2 && hasNext2 ; hasNext2 = hashTableIteratorNext(htIter2))
	{
	  Dropset *ds2 = getCurrentValueFromHashTableIterator(htIter2);
	  if(ds == ds2)
	    continue;

	  if(indexListEqual(ds->taxaToDrop, ds2->taxaToDrop))
	    {
	      PR("duplicate dropset: ");
	      printIndexList(ds->taxaToDrop);
	      PR(" and ");
	      printIndexList(ds2->taxaToDrop);
	      PR("\n");
	      exit(-1);
	    }
	}
      free(htIter2);
    }
  free(htIter);
}
#endif


boolean isCompatible(ProfileElem* elemA, ProfileElem* elemB, BitVector *droppedTaxa)
{
  unsigned int i;
  
  unsigned int 
    *A = elemA->bitVector,
    *C = elemB->bitVector;
  
  FOR_0_LIMIT(i,bitVectorLength)
    if(A[i] & C[i]  & ~ (droppedTaxa[i] | paddingBits[i]) )
      break;
          
  if(i == bitVectorLength)
    return TRUE;
  
  FOR_0_LIMIT(i,bitVectorLength)
    if( ( A[i] & ~C[i]  ) & ~ (droppedTaxa[i] | paddingBits[i]) )
      break;
   
  if(i == bitVectorLength)  
    return TRUE;  
  
  FOR_0_LIMIT(i,bitVectorLength)
    if( ( ~A[i] & C[i] )  & ~ (droppedTaxa[i] | paddingBits[i]) )
      break;
  
  if(i == bitVectorLength)
    return TRUE;  
  else
    return FALSE;
}


#ifdef MYDEBUG
/* ensures, no merging event occurs twice per dropset */
void debug_mergingHashSanityCheck(HashTable *mergingHash, int totalNumberOfBips)
{
  HashTableIterator *htIter; 
  FOR_HASH(htIter, mergingHash)
    {
      Dropset
	*ds = getCurrentValueFromHashTableIterator(htIter);

      if(NOT ds)
	break;
      
      BitVector
	*bv = CALLOC(totalNumberOfBips, sizeof(BitVector));
      
      List
	*meIter = ds->primeEvents;

      FOR_LIST(meIter)
      {
	MergingEvent
	  *me = meIter->value;  
	if(me->isComplex)
	  {
	    IndexList *il =  me->mergingBipartitions.many; 
	    FOR_LIST(il)
	    {
	      assert(NOT NTH_BIT_IS_SET(bv, il->index));
	      FLIP_NTH_BIT(bv, il->index);
	    }
	  }
	else
	  {
	    assert(NOT NTH_BIT_IS_SET(bv, me->mergingBipartitions.pair[0]));
	    assert(NOT NTH_BIT_IS_SET(bv, me->mergingBipartitions.pair[1]));
	    FLIP_NTH_BIT(bv, me->mergingBipartitions.pair[0]);
	    FLIP_NTH_BIT(bv, me->mergingBipartitions.pair[1]);
	  }
      }
      free(bv);
    }
  free(htIter);
}
#endif


boolean myBitVectorEqual(ProfileElem *elemA, ProfileElem *elemB)
{
  boolean normalEqual = TRUE,
    complement = TRUE;

  int i ;
  FOR_0_LIMIT(i,bitVectorLength)
    {
      normalEqual = normalEqual && (  elemA->bitVector[i] == elemB->bitVector[i]);
      complement  =  complement && (elemA->bitVector[i] == ~(elemB->bitVector[i] | droppedTaxa[i] | paddingBits[i]));
    }

  return normalEqual || complement;
}


boolean canMergeWithComplement(ProfileElem *elem)
{
  return mxtips - taxaDropped - 2 * elem->numberOfBitsSet <= 2 * maxDropsetSize;
}


boolean bitVectorEqual(ProfileElem *elemA, ProfileElem *elemB)
{
  boolean normalEqual = TRUE,
    complement = TRUE;

  int i ; 
  FOR_0_LIMIT(i,bitVectorLength)
    {
      normalEqual = normalEqual && (  elemA->bitVector[i] == elemB->bitVector[i]);
      complement  =  complement && (elemA->bitVector[i] == ~(elemB->bitVector[i] | droppedTaxa[i] | paddingBits[i]));
    }

  return normalEqual || complement;
}


boolean mergedBipVanishes(MergingEvent *me, Array *bipartitionsById, IndexList *taxaToDrop)
{
  int vanBits = 0; 

  IndexList
    *iter = taxaToDrop;  
  
  ProfileElem
    *elem = me->isComplex ? GET_PROFILE_ELEM(bipartitionsById, (me->mergingBipartitions).many->index)  : GET_PROFILE_ELEM(bipartitionsById, (me->mergingBipartitions).pair[0]);     
  
  FOR_LIST(iter)
    if(NTH_BIT_IS_SET(elem->bitVector,iter->index))
      vanBits++;
  
  return elem->numberOfBitsSet - vanBits < 2; 
}


/* insert dropset. If it is a multi-taxa dropset, gather all merging
   events of sub-dropsets */
Dropset *insertOrFindDropset(HashTable *hashtable, Dropset *dropset, unsigned int hashValue) 
{
  void
    *result = searchHashTable(hashtable, dropset, hashValue);

  if(result)
    {
      freeDropsetDeep(dropset, TRUE);
      return result;
    }
  else
    {
      insertIntoHashTable(hashtable, dropset, hashValue);      
      return dropset;
    }
}

boolean checkForMergerAndAddEvent(boolean complement, ProfileElem *elemA, ProfileElem *elemB, HashTable *mergingHash)
{
  IndexList
    *dropsetTaxa = getDropset(elemA,elemB,complement, neglectThose);
  
  if(dropsetTaxa)
    {
      Dropset
	*dropset,
	*tmp = CALLOC(1,sizeof(Dropset));
      tmp->taxaToDrop = dropsetTaxa;      
      
      unsigned int hashValue = 0; 
      IndexList *iter =  dropsetTaxa; 
      FOR_LIST(iter)  
      {
	assert(iter->index < mxtips);
	hashValue ^= randForTaxa[ iter->index ];
      }
      
#ifdef PARALLEL
      int position  = hashValue % mergingHash->tableSize;      
      pthread_mutex_lock(mergingHash->lockPerSlot[position]);
#endif
      dropset = insertOrFindDropset(mergingHash, tmp, hashValue);
      addEventToDropsetPrime(dropset, elemA->id, elemB->id);
#ifdef PARALLEL
      pthread_mutex_unlock(mergingHash->lockPerSlot[position]);
#endif
      return TRUE;
    } 
  else 
    return FALSE;
}


/* can i trust you tiny function? (practical relevant -> 0) */
/* boolean bothDropsetsRelevant(ProfileElem *elemA) */
boolean bothDropsetsRelevant(int numBits)
{
  return numBits <= maxDropsetSize && numBits >= mxtips - taxaDropped - maxDropsetSize; 
}


int cleanup_applyOneMergerEvent(MergingEvent *mergingEvent, Array *bipartitionsById, BitVector *mergingBipartitions)
{
  int 
    j;
  
  ProfileElem
    *resultBip, *elem; 

  resultBip = 
    mergingEvent->isComplex
    ? GET_PROFILE_ELEM(bipartitionsById, mergingEvent->mergingBipartitions.many->index) 
    : GET_PROFILE_ELEM(bipartitionsById, mergingEvent->mergingBipartitions.pair[0]);

  if(mergingEvent->isComplex)
    {
	IndexList
	  *iterBip = mergingEvent->mergingBipartitions.many->next; 
	FOR_LIST(iterBip)
	{
	  elem = GET_PROFILE_ELEM(bipartitionsById, iterBip->index);
	  FLIP_NTH_BIT(mergingBipartitions, elem->id);
	  resultBip->isInMLTree |= elem->isInMLTree;
	  FOR_0_LIMIT(j,treeVectorLength)
	    resultBip->treeVector[j] |= elem->treeVector[j]; 
	}
	
	freeIndexList(mergingEvent->mergingBipartitions.many);
	free(mergingEvent);
    }
  else
    {      
      elem = GET_PROFILE_ELEM(bipartitionsById,mergingEvent->mergingBipartitions.pair[1]);
      FLIP_NTH_BIT(mergingBipartitions, elem->id);      
      resultBip->isInMLTree |= elem->isInMLTree;
      FOR_0_LIMIT(j,treeVectorLength)
	resultBip->treeVector[j] |=  elem->treeVector[j];
    }

  resultBip->treeVectorSupport = genericBitCount(resultBip->treeVector, treeVectorLength);
  return resultBip->id;
}


int getSupportOfMRETreeHelper(Array *bipartitionProfile, Dropset *dropset) 
{
  int
    result = 0,
    i,j; 

  BitVector
    *taxaDroppedHere = copyBitVector(droppedTaxa, bitVectorLength); 

  if(dropset)
    {
      IndexList
	*iter = dropset->taxaToDrop;
      FOR_LIST(iter)	
	FLIP_NTH_BIT(taxaDroppedHere, iter->index);  
    }

  qsort(bipartitionProfile->arrayTable, bipartitionProfile->length, sizeof(ProfileElem**), sortBySupport);

  Array *mreBips = CALLOC(1,sizeof(Array)); 
  mreBips->arrayTable = CALLOC((mxtips-3), sizeof(ProfileElem*)); 
  mreBips->length = 0;

#ifdef MYDEBUG
  for(i = 1; i < bipartitionProfile->length; i ++)
    assert(GET_PROFILE_ELEM(bipartitionProfile,i-1)->treeVectorSupport >= GET_PROFILE_ELEM(bipartitionProfile,i)->treeVectorSupport);
#endif

  FOR_0_LIMIT(i, bipartitionProfile->length)
    if(GET_PROFILE_ELEM(bipartitionProfile,i)->treeVectorSupport > thresh)
      addElemToArray(GET_PROFILE_ELEM(bipartitionProfile,i), mreBips);
    else
      break;

  for(; i < bipartitionProfile->length && mreBips->length < mxtips-3; ++i)
    {
      ProfileElem
	*elemA = GET_PROFILE_ELEM(bipartitionProfile, i);
      boolean compatibleP = TRUE; 

      FOR_0_LIMIT(j,mreBips->length)
	{
	  ProfileElem
	    *elemB  = GET_PROFILE_ELEM(mreBips,j);

	  compatibleP &= isCompatible(elemA, elemB, taxaDroppedHere); 

	  if( NOT compatibleP)
	    break;
	}

      if(compatibleP)
	addElemToArray(GET_PROFILE_ELEM(bipartitionProfile,i), mreBips);
    }

  if(computeSupport)
    FOR_0_LIMIT(i,mreBips->length)      
      result +=  GET_PROFILE_ELEM(mreBips,i)->treeVectorSupport;
  else
    result = mreBips->length; 

  free(mreBips->arrayTable);  free(mreBips);
  free(bipartitionProfile->arrayTable);  free(bipartitionProfile);

  return result; 
}


void getSupportGainedThreshold(MergingEvent *me, Array *bipartitionsById)
{
  int
    i; 
  me->supportGained = 0; 
  BitVector
    *tmp;
  boolean isInMLTree = FALSE; 

  if(me->isComplex)
    {
      IndexList
	*iI = me->mergingBipartitions.many;  
      
      int bestPossible = 0; 
      FOR_LIST(iI)
      {	
	ProfileElem
	  *elem = GET_PROFILE_ELEM(bipartitionsById, iI->index);
	bestPossible += elem->treeVectorSupport;
	isInMLTree |= elem->isInMLTree;
      }

      if( rogueMode == VANILLA_CONSENSUS_OPT && bestPossible < thresh)
	return ; 
      if( rogueMode == ML_TREE_OPT && NOT isInMLTree)
	return ;

      tmp = CALLOC(treeVectorLength, sizeof(BitVector));
      
      /* create new bip vector */
      iI = me->mergingBipartitions.many;  
      FOR_LIST(iI)
      {
	ProfileElem
	  *elem = GET_PROFILE_ELEM(bipartitionsById, iI->index);
    
	FOR_0_LIMIT(i, treeVectorLength)
	  tmp[i] |= elem->treeVector[i];
      }
    }
  else
    {
      ProfileElem
	*elemA = GET_PROFILE_ELEM(bipartitionsById, me->mergingBipartitions.pair[0]),
	*elemB = GET_PROFILE_ELEM(bipartitionsById, me->mergingBipartitions.pair[1]);      
   
      if(rogueMode == VANILLA_CONSENSUS_OPT && elemA->treeVectorSupport + elemB->treeVectorSupport < thresh)
      	return;
      
      isInMLTree = elemA->isInMLTree || elemB->isInMLTree;       
      if(rogueMode == ML_TREE_OPT && NOT isInMLTree)
	return; 
      
      tmp = CALLOC(treeVectorLength, sizeof(BitVector));      
      FOR_0_LIMIT(i,treeVectorLength)
	tmp[i] = elemA->treeVector[i] | elemB->treeVector[i];
    }

  int newSup = genericBitCount(tmp, treeVectorLength);
  switch (rogueMode)
    {
    case MRE_CONSENSUS_OPT:
      {
	me->supportGained = computeSupport ? newSup : 1; 
	break; 
      }
    case VANILLA_CONSENSUS_OPT: 
      {
	if(rogueMode == VANILLA_CONSENSUS_OPT  && newSup > thresh)
	  me->supportGained = computeSupport ? newSup : 1 ;
	break; 
      }
    case ML_TREE_OPT: 
      {
	if(isInMLTree )
	  me->supportGained = computeSupport ? newSup : 1 ;
	break; 
      }
    default : 
      assert(0);
    }
  
  free(tmp);
}


int getSupportOfMRETree(Array *bipartitionsById,  Dropset *dropset)
{
  List
    *mergingEvents = NULL; 
  if(dropset)
    {
      if(maxDropsetSize == 1)
	mergingEvents =  dropset->ownPrimeE; 
      else
	{
	  List *iter = dropset->acquiredPrimeE ; 
	  FOR_LIST(iter)
	    APPEND(iter->value, mergingEvents);
	  iter = dropset->complexEvents; 
	  FOR_LIST(iter)
	    APPEND(iter->value, mergingEvents);
	}
    }
  int
    i;

  /* initial case  */
  if( NOT dropset)
    {
      Array *array = cloneProfileArrayFlat(bipartitionsById);
      int tmp = getSupportOfMRETreeHelper( array, dropset);
      return tmp; 
    }

  Array
    *emergedBips = CALLOC(1,sizeof(Array)),
    *finalArray  = CALLOC(1,sizeof(Array)),
    *tmpArray = cloneProfileArrayFlat(bipartitionsById);
  emergedBips->arrayTable = CALLOC(lengthOfList(mergingEvents), sizeof(ProfileElem*));
  emergedBips->length = 0; 

  finalArray->arrayTable = CALLOC(tmpArray->length, sizeof(ProfileElem*));
  finalArray->length = 0;

  /* kill merging bips from array */
  List *meIter = mergingEvents ; 
  FOR_LIST(meIter)
  {
    MergingEvent *me = meIter->value; 
    if(me->isComplex)
      {
	IndexList *iter = me->mergingBipartitions.many; 
	FOR_LIST(iter)	  
	  GET_PROFILE_ELEM(tmpArray, iter->index) = NULL; 
	
	/* create emerged bips in other array */
	ProfileElem *elem = CALLOC(1,sizeof(ProfileElem)); 
	getSupportGainedThreshold(me,bipartitionsById);
	elem->treeVectorSupport = me->supportGained; 
	elem->bitVector = GET_PROFILE_ELEM(bipartitionsById, me->mergingBipartitions.many->index)->bitVector;
	GET_PROFILE_ELEM(emergedBips, emergedBips->length) = elem;
	emergedBips->length++;
      }
    else
      {
	int a = me->mergingBipartitions.pair[0],
	  b = me->mergingBipartitions.pair[1];
	GET_PROFILE_ELEM(tmpArray, a) = NULL; 
	GET_PROFILE_ELEM(tmpArray, b) = NULL; 

	/* create emerged bips in other array */
	ProfileElem *elem = CALLOC(1,sizeof(ProfileElem)); 
	getSupportGainedThreshold(me,bipartitionsById);
	elem->treeVectorSupport = me->supportGained; 
	elem->bitVector = GET_PROFILE_ELEM(bipartitionsById, a)->bitVector;
	GET_PROFILE_ELEM(emergedBips, emergedBips->length) = elem;
	emergedBips->length++;
      }
  }

  /* kill vanishing bips from array */
  FOR_0_LIMIT(i, tmpArray->length)
    {
      if( GET_PROFILE_ELEM(tmpArray,i) ) 
	{
	  ProfileElem *elem = GET_PROFILE_ELEM(tmpArray,i);
	  int remainingBits = elem->numberOfBitsSet; 
	  IndexList *iter = dropset->taxaToDrop; 
	  FOR_LIST(iter)
	    if(NTH_BIT_IS_SET(elem->bitVector,  iter->index))
	      remainingBits--;	  
	  if(remainingBits > 1) 
	    addElemToArray(elem, finalArray); 
	}
    }

  FOR_0_LIMIT(i, emergedBips->length)
    addElemToArray(GET_PROFILE_ELEM(emergedBips, i), finalArray);

  free(tmpArray->arrayTable);  free(tmpArray);

  int result = getSupportOfMRETreeHelper(finalArray, dropset);

  if(maxDropsetSize > 1 )
    freeListFlat(mergingEvents);    

  FOR_0_LIMIT(i,emergedBips->length)
    free(GET_PROFILE_ELEM(emergedBips,i));
  free(emergedBips->arrayTable);  free(emergedBips);
  
  return result; 
}


boolean bipartitionVanishesP(ProfileElem *elem, Dropset *dropset)
{
  IndexList *iter = dropset->taxaToDrop;
  int result = elem->numberOfBitsSet;

  FOR_LIST(iter)
    if(NTH_BIT_IS_SET(elem->bitVector, iter->index))
      result--;  

  return result < 2; 
}


void removeMergedBipartitions(Array *bipartitionsById, Array *bipartitionProfile, BitVector *mergingBipartitions)
{
  int 
    i;

  FOR_0_LIMIT(i,bipartitionProfile->length)
    {
      ProfileElem
	*elem = GET_PROFILE_ELEM(bipartitionProfile,i);

      if( NOT elem )
	continue;

      if(NTH_BIT_IS_SET(mergingBipartitions,elem->id)) 
	{
	  GET_PROFILE_ELEM(bipartitionProfile, i) = NULL;
	  GET_PROFILE_ELEM(bipartitionsById, elem->id) = NULL;
	  freeProfileElem(elem);
#ifdef PRINT_VERY_VERBOSE
	  PR("CLEAN UP: removing %d from bip profile because of merger\n", elem->id);
#endif
	}
    }
}


boolean eventMustBeRecomputed(MergingEvent *meIter, BitVector *mergingBipartitions, BitVector *newCandidates)
{
  boolean
    mustBeRecomputed = FALSE;

  IndexList
    *mergingBipsIter = meIter->mergingBipartitions.many;
  FOR_LIST(mergingBipsIter)
    mustBeRecomputed |= 
    NTH_BIT_IS_SET(mergingBipartitions,mergingBipsIter->index) 
    | NTH_BIT_IS_SET(newCandidates, mergingBipsIter->index);

  return mustBeRecomputed;
}


boolean checkValidityOfEvent(BitVector *obsoleteBips, List *elem)
{
  MergingEvent
    *me = elem->value; 
  boolean
    killP = FALSE;

  if(me->isComplex)
    {      
      IndexList
	*iter  = me->mergingBipartitions.many;      
      FOR_LIST(iter)
	killP |= NTH_BIT_IS_SET(obsoleteBips, iter->index);       
      if(killP)
      	freeIndexList(me->mergingBipartitions.many);
    }
  else
    killP = NTH_BIT_IS_SET(obsoleteBips, me->mergingBipartitions.pair[0]) || NTH_BIT_IS_SET(obsoleteBips, me->mergingBipartitions.pair[1]) ;   

  if(killP)
    {
      free(me);
      return FALSE;
    }
  else
    return TRUE; 
}


#ifdef LATER
void printDropset(Dropset *dropset)
{
  IndexList *iter; 
  iter = dropset->taxaToDrop;
  boolean isFirst = TRUE;
  FOR_LIST(iter)
  {
    PR(isFirst ? "%d" : ",%d", iter->index);
    isFirst = FALSE;
  }
  PR(" [%d] : ", dropset->improvement);
  List *lIter = dropset->primeEvents;
  FOR_LIST(lIter)
  {
    MergingEvent *me = lIter->value;    
    PR("{ %d,%d },", me->mergingBipartitions.pair[0], me->mergingBipartitions.pair[1]);
  }
  PR("\n");  
  if(dropset->combinedEvents)
    {
      PR("\t`->");
      lIter = dropset->combinedEvents;
      FOR_LIST(lIter)
      {
	if(((MergingEvent*)lIter->value)->isComplex)
	  {	
	    isFirst = TRUE;	
	    PR("{ ");
	    iter = ((MergingEvent*)lIter->value)->mergingBipartitions.many;	
	    FOR_LIST(iter)
	    {
	      PR(isFirst ? "%d" : ",%d" ,iter->index);
	      isFirst = FALSE;
	    }
	    PR(" },");
	  }
	else
	  {
	    MergingBipartitions mb = ((MergingEvent*)lIter->value)->mergingBipartitions ; 
	    PR("[ %d,%d ],", mb.pair[0], mb.pair[1]); 
	  }
      }
      PR("\n");
    }
}
#endif


#ifdef LATER
void printMergingHash(HashTable *mergingHash)
{
  if(mergingHash->entryCount < 1)
    {
      PR("Nothing in mergingHash.\n");
      return;
    }

  HashTableIterator *htIter;  
  FOR_HASH(htIter, mergingHash)
    {
      Dropset
	*dropset = getCurrentValueFromHashTableIterator(htIter); 
      printDropset(dropset);
    }
}
#endif


#ifdef MYDEBUG
void debug_assureCleanStructure(HashTable *hashtable, BitVector *mergingBipartitions)
{
  HashTableIterator *htIter;

  FOR_HASH(htIter, hashtable)
    {
      Dropset *ds = (Dropset*)getCurrentValueFromHashTableIterator(htIter); 

      if( NOT ds)
	break;

      List *meIter = ds->primeEvents;
      FOR_LIST(meIter)
      {
	MergingEvent *me = meIter->value; 
	if(me->isComplex)
	  {	    
	    IndexList *iter = me->mergingBipartitions.many;
	    FOR_LIST(iter)
	      if(NTH_BIT_IS_SET(mergingBipartitions, iter->index))
		{
		  PR("%d from merging bipartitions still present. \n", iter->index);
		  printMergingHash(hashtable);
		  exit(-1);
		}
	  }
	else
	  {
	    if(NTH_BIT_IS_SET(mergingBipartitions, me->mergingBipartitions.pair[0]) )
	      {
		PR("%d from merging bipartitions still present. \n", me->mergingBipartitions.pair[0]);
		printMergingHash(hashtable);
		exit(-1);
	      }
	    if(NTH_BIT_IS_SET(mergingBipartitions, me->mergingBipartitions.pair[1]))
	      {
		PR("%d from merging bipartitions still present. \n", me->mergingBipartitions.pair[1]);
		printMergingHash(hashtable);
		exit(-1);
	      }
	  }
      }
    }
  free(htIter);
}
#endif



void cleanup_mergingEvents(HashTable *mergingHash, BitVector *mergingBipartitions, BitVector *candidateBips, int length)
{
  HashTableIterator
    *htIter;

  int i; 
  FOR_0_LIMIT(i,GET_BITVECTOR_LENGTH(length))
    mergingBipartitions[i] |= candidateBips[i];

  if( NOT mergingHash->entryCount)
    return;

  FOR_HASH(htIter, mergingHash)
    {
      Dropset
	*dropset = getCurrentValueFromHashTableIterator(htIter);
    
      /* always remove combined events */
      List *iter = dropset->complexEvents;
      FOR_LIST(iter)
      {
	MergingEvent *me = iter->value; 
	assert(me);
	if(me && me->isComplex)
	  {
	    freeIndexList(me->mergingBipartitions.many);
	    free(me);
	  }	
      }
      freeListFlat(dropset->complexEvents);
      
      /* always remove acquired elems */
      freeListFlat(dropset->acquiredPrimeE);
    }
  free(htIter);

  /* reduce own elements */
  FOR_HASH_2(htIter, mergingHash)
    {
      Dropset
	*dropset = getCurrentValueFromHashTableIterator(htIter);

      assert(dropset);

      /* prime events */
      List 
	*iter = dropset->ownPrimeE, 
	*start = NULL; 
      while(iter)
	{
	  List *next = iter->next;
	  if( checkValidityOfEvent(mergingBipartitions, iter) ) 
	    APPEND(iter->value, start); 
	  free(iter);
	  iter = next; 
	}
      dropset->ownPrimeE = start; 
    }
  free(htIter);  
  
#ifdef MYDEBUG
  debug_assureCleanStructure(mergingHash, mergingBipartitions);
#endif

  free(mergingBipartitions);  
}


/* 
   • inverses the bit vector, if more than half of the remaining taxa
   bits is set. This is better anyway, is the bit vector do not get
   that physically heavy (assuming that a 1 weighs more than a 0)
   • assuming, we already know the numbers of bits set 
   • assuming, bits are unflipped, if the taxon was dropped  
*/
void unifyBipartitionRepresentation(Array *bipartitionArray,  BitVector *droppedTaxa)
{
  int
    i,j,
    bvLen = GET_BITVECTOR_LENGTH(mxtips),
    remainingTaxa = mxtips - genericBitCount(droppedTaxa, bvLen);

#ifdef PRINT_VERY_VERBOSE
  PR("remaining taxa: %d\n", remainingTaxa);
  PR("inverting bit vectors: ");
#endif

  FOR_0_LIMIT(i,bipartitionArray->length)
    {
      ProfileElem
	*elem = GET_PROFILE_ELEM(bipartitionArray,i);

      if( elem
	  && elem->numberOfBitsSet > remainingTaxa / 2)	
	{

#ifdef PRINT_VERY_VERBOSE
	  PR("%d (%d bits set), ", elem->id, elem->numberOfBitsSet);
#endif
	  FOR_0_LIMIT(j,bvLen)
	    elem->bitVector[j] = ~(elem->bitVector[j] | paddingBits[j] |  droppedTaxa[j]);
	  elem->numberOfBitsSet = remainingTaxa - elem->numberOfBitsSet;
	}
    }
#ifdef PRINT_VERY_VERBOSE
  PR("\n");
#endif
}


void printBipartitionProfile(Array *bipartitionProfile)
{
  int i; 
  FOR_0_LIMIT(i,bipartitionProfile->length)
    {
      ProfileElem *elem = GET_PROFILE_ELEM(bipartitionProfile, i);      
      if(elem)
	{
	  PR("%d (%d):\t\t", elem->id, elem->numberOfBitsSet);
	  printBitVector(elem->bitVector, GET_BITVECTOR_LENGTH(mxtips));
	}
      else
	break;
      PR("\n");
    }
}


int getNumberOfBipsPresent(Array *bipartitionArray)
{
  int result = 0; 
  int i; 

  FOR_0_LIMIT(i, bipartitionArray->length)
    {
      if(GET_PROFILE_ELEM(bipartitionArray,i))
	result++;
    }

  return result; 
}


int getInitScore(Array *bipartitionProfile)
{
  int 
    score = 0 , i; 

  if(rogueMode == MRE_CONSENSUS_OPT)
    return getSupportOfMRETree(bipartitionProfile, NULL);

  FOR_0_LIMIT(i,bipartitionProfile->length)
    {
      ProfileElem
	*elem = GET_PROFILE_ELEM(bipartitionProfile,i);

      switch(rogueMode)
	{
	case VANILLA_CONSENSUS_OPT:
	  if(elem->treeVectorSupport > thresh)
	    score += computeSupport ? elem->treeVectorSupport : 1 ; 
	  break;

	case ML_TREE_OPT:
	  if(elem->isInMLTree)
	    score += computeSupport ? elem->treeVectorSupport : 1;
	  break;

	case MRE_CONSENSUS_OPT:
	default:
	  assert(0);
	}
    }

  return score; 
}


void printDropsetImprovement(Dropset *dropset, All *tr, int cumScore)
{  
#ifndef PRINT_DROPSETS
  return ; 
#endif
  
  IndexList
    *iter = dropset->taxaToDrop;
  boolean isFirst = TRUE;

  FOR_LIST(iter)
  {
    PR( isFirst ? ">%d" : ",%d", iter->index);
    isFirst = FALSE;    
  }

  isFirst = TRUE;
  PR("\t");
  iter = dropset->taxaToDrop; 
  FOR_LIST(iter)
  {
    PR(isFirst ? "%s" : ",%s" , tr->nameList[iter->index+1]);
    isFirst = FALSE;
  }

  PR("\t");
  PR("%f\t%f\n",
     (double)dropset->improvement /(computeSupport ?   (double)numberOfTrees : 1) ,  
     (double)cumScore / (double)( (computeSupport ? numberOfTrees : 1)  * (mxtips-3)));
}


void fprintRogueNames(All *tr, FILE *file, IndexList *list)
{
  boolean isFirst = TRUE; 

  FOR_LIST(list)
  {
    if(isFirst)
      {
	fprintf(file, "%s", tr->nameList[list->index+1]);
	isFirst = FALSE;
      }
    else
      fprintf(file, ",%s", tr->nameList[list->index+1]); 
  }
}


void printRogueInformationToFile( All *tr, FILE *rogueOutput, int bestCumEver, int *cumScores, Dropset **dropsetInRound)
{
  int
    i=1,j; 

  boolean reached = bestCumEver == cumScores[0];
  while ( NOT reached)
    {
      fprintf(rogueOutput, "%d\t", i);       
      printIndexListToFile(rogueOutput, dropsetInRound[i]->taxaToDrop); 
      fprintf(rogueOutput, "\t");
      fprintRogueNames(tr, rogueOutput, dropsetInRound[i]->taxaToDrop);
      fprintf(rogueOutput, "\t%f\t%f\n", (double)(cumScores[i]  - cumScores[i-1] )/ (double)tr->numberOfTrees,(double)cumScores[i] / (double)((computeSupport ? numberOfTrees : 1 ) * (mxtips-3)) ); 
      reached = bestCumEver == cumScores[i];
      ++i;
    }

  FOR_0_LIMIT(j,mxtips)
    if(NOT NTH_BIT_IS_SET(neglectThose,j))
      {
	fprintf(rogueOutput, "%d\t%d\t%s\t%s\t%s\n", i, j, tr->nameList[j+1], "NA", "NA");
	i++;
      }
}


void findCandidatesForBip(HashTable *mergingHash, ProfileElem *elemA, boolean firstMerge, Array *bipartitionsById, Array *bipartitionProfile, int* indexByNumberBits)
{
  ProfileElem 
    *elemB;
  int indexInBitSortedArray; 

  boolean
    compMerge = canMergeWithComplement(elemA);

  if(firstMerge)
    {
      if(NOT compMerge && maxDropsetSize == 1)
	indexInBitSortedArray = indexByNumberBits[elemA->numberOfBitsSet +1]; 
      else	
	indexInBitSortedArray = indexByNumberBits[elemA->numberOfBitsSet];
    }
  else
    indexInBitSortedArray = 
      elemA->numberOfBitsSet - maxDropsetSize < 0 ?
      indexByNumberBits[0]
      : indexByNumberBits[elemA->numberOfBitsSet-maxDropsetSize];
	
  for( ;
       indexInBitSortedArray < bipartitionProfile->length
	 && (elemB = GET_PROFILE_ELEM(bipartitionProfile,indexInBitSortedArray))
	 && elemB->numberOfBitsSet - elemA->numberOfBitsSet <= maxDropsetSize ;
       indexInBitSortedArray++)
    { 
      if(
	 maxDropsetSize == 1 && 
	 NOT compMerge && 
	 elemA->numberOfBitsSet == elemB->numberOfBitsSet)
	continue;

      boolean foundOne = FALSE;
      if(compMerge)
	foundOne = checkForMergerAndAddEvent(TRUE,elemA, elemB, mergingHash); 
      
      if(NOT foundOne || bothDropsetsRelevant(elemA->numberOfBitsSet))
	checkForMergerAndAddEvent(FALSE, elemA, elemB, mergingHash);	    
    }
}


/* HashTable * */
void createOrUpdateMergingHash(All *tr, HashTable *mergingHash, Array *bipartitionProfile, Array *bipartitionsById, BitVector *candidateBips, boolean firstMerge, int *indexByNumberBits) 
{
  int  i; 
  FOR_0_LIMIT(i,bipartitionProfile->length)
    if(NTH_BIT_IS_SET(candidateBips, i))
      findCandidatesForBip(mergingHash, GET_PROFILE_ELEM(bipartitionsById, i),  firstMerge, bipartitionsById, bipartitionProfile, indexByNumberBits);  
  
  free(candidateBips);
}


void combineEventsForOneDropset(Array *allDropsets, Dropset *refDropset, Array *bipartitionsById)
{
  List *allEventsUncombined = NULL;
  refDropset->acquiredPrimeE = NULL; 
  refDropset->complexEvents = NULL;
  int eventCntr = 0; 

  if(NOT refDropset->taxaToDrop->next)
    {
      List *iter =  refDropset->ownPrimeE;
      FOR_LIST(iter)
       APPEND(iter->value, refDropset->acquiredPrimeE);
      return; 
    }

  /* gather all events */
  int i; 
  FOR_0_LIMIT(i,allDropsets->length)
    {      
      Dropset *currentDropset = GET_DROPSET_ELEM(allDropsets, i);
      if( isSubsetOf(currentDropset->taxaToDrop, refDropset->taxaToDrop) )
	{
	  List
	    *iter = currentDropset->ownPrimeE; 
	  FOR_LIST(iter)
	  {
	    APPEND(iter->value, allEventsUncombined);
	    eventCntr++;
	  }
	}
    }

  /* transform the edges into nodes */
  HashTable *allNodes = createHashTable(eventCntr * 10, NULL, nodeHashValue, nodeEqual);
  Node *found;
  List *iter = allEventsUncombined; 
  FOR_LIST(iter)
  {
    MergingEvent *me = (MergingEvent*)iter->value;
    int a = me->mergingBipartitions.pair[0]; 
    int b = me->mergingBipartitions.pair[1]; 


    if( ( found = searchHashTableWithInt(allNodes, a) ) )
      APPEND_INT(b,found->edges); 
    else
      {
	Node *node = CALLOC(1,sizeof(Node));
	node->id = a; 
	APPEND_INT(b,node->edges);
	insertIntoHashTable(allNodes, node, a);
      }

    if( ( found = searchHashTableWithInt(allNodes, b) ) )
      APPEND_INT(a,found->edges); 
    else
      {
	Node *node = CALLOC(1,sizeof(Node));
	node->id = b; 
	APPEND_INT(a,node->edges);
	insertIntoHashTable(allNodes, node, b);
      }
  }

  iter = allEventsUncombined;
  FOR_LIST(iter)
  {
    MergingEvent *me = (MergingEvent*)iter->value;
    int a = me->mergingBipartitions.pair[0],
      b = me->mergingBipartitions.pair[1]; 
    
    Node *foundA = searchHashTableWithInt(allNodes,a),
      *foundB = searchHashTableWithInt(allNodes,b);

    if(NOT foundA->edges->next
       && NOT foundB->edges->next) 
      {
	assert(foundA->edges->index == foundB->id ); 
	assert(foundB->edges->index == foundA->id ); 
	APPEND(me, refDropset->acquiredPrimeE); 
      }
    else
      {	
	IndexList
	  *component = findAnIndependentComponent(allNodes,foundA);
	if( component)
	  {
	    MergingEvent *complexMe  = CALLOC(1,sizeof(MergingEvent));
	    complexMe->mergingBipartitions.many = component; 
	    complexMe->isComplex = TRUE;
	    APPEND(complexMe,refDropset->complexEvents);
	  }
      }
  }
  
  destroyHashTable(allNodes, freeNode);
  freeListFlat(allEventsUncombined);  
}


HashTable *combineMergerEvents(HashTable *mergingHash, Array *bipartitionsById) 
{  
  /* hash to array  */
  Array *allDropsets =  CALLOC(1,sizeof(Array));
  allDropsets->arrayTable = CALLOC(mergingHash->entryCount, sizeof(Dropset**));
  
  HashTableIterator *htIter; 
  int cnt = 0; 
  FOR_HASH(htIter, mergingHash)
    {
      GET_DROPSET_ELEM(allDropsets, cnt) = getCurrentValueFromHashTableIterator(htIter);
      cnt++;
    }
  free(htIter);
  assert(cnt == mergingHash->entryCount);
  allDropsets->length = cnt; 

#ifdef PARALLEL
  globalPArgs->allDropsets = allDropsets; 
  globalPArgs->bipartitionsById =bipartitionsById; 
  numberOfJobs = allDropsets->length;
  masterBarrier(THREAD_COMBINE_EVENTS, globalPArgs); 
#else
  int i; 
  FOR_0_LIMIT(i,allDropsets->length)   
    combineEventsForOneDropset(allDropsets, GET_DROPSET_ELEM(allDropsets,i), bipartitionsById);
#endif

  free(allDropsets->arrayTable);
  free(allDropsets);

  return mergingHash; 
}


void getLostSupportThreshold(MergingEvent *me, Array *bipartitionsById)
{
  ProfileElem *elemA, *elemB ; 
  me->supportLost = 0; 
  
  if(me->isComplex)
    {
      IndexList *iI = me->mergingBipartitions.many; 
      
      FOR_LIST(iI)
      {
	elemA = GET_PROFILE_ELEM(bipartitionsById, iI->index);
	switch (rogueMode)
	{
	case VANILLA_CONSENSUS_OPT : 
	  {
	    if(elemA->treeVectorSupport > thresh)
	      me->supportLost += computeSupport ? elemA->treeVectorSupport : 1; 
	    break ;
	  }
	case ML_TREE_OPT: 
	  {
	    if(elemA->isInMLTree)
	      me->supportLost += computeSupport ? elemA->treeVectorSupport : 1  ; 
	    break; 
	  }
	default : 
	  assert(0);
	}
      }
    }
  else
    {       
      elemA = GET_PROFILE_ELEM(bipartitionsById, me->mergingBipartitions.pair[0]);
      elemB = GET_PROFILE_ELEM(bipartitionsById, me->mergingBipartitions.pair[1]);
      
      switch(rogueMode)
	{
	case MRE_CONSENSUS_OPT: 
	case VANILLA_CONSENSUS_OPT: 
	  {
	    if(elemA->treeVectorSupport > thresh)
	      me->supportLost += computeSupport ? elemA->treeVectorSupport : 1 ;
	    if(elemB->treeVectorSupport > thresh)
	      me->supportLost += computeSupport ? elemB->treeVectorSupport : 1;
	    break; 
	  }
	case ML_TREE_OPT:
	  {
	    if(elemA->isInMLTree)
	      me->supportLost += computeSupport ? elemA->treeVectorSupport : 1 ; 
	    if(elemB->isInMLTree)
	      me->supportLost += computeSupport ? elemB->treeVectorSupport : 1 ; 
	  }
	}
    }
}


void evaluateDropset(HashTable *mergingHash, Dropset *dropset,Array *bipartitionsById, List *consensusBipsCanVanish )
{
  int result = 0; 
  List
    *allElems = NULL, 
    *elemsToCheck = NULL;
  
  if(maxDropsetSize == 1)
    elemsToCheck = dropset->ownPrimeE ;
  else
    {
      List *otherIter = dropset->acquiredPrimeE; 
      FOR_LIST(otherIter)
	APPEND(otherIter->value, elemsToCheck);
      otherIter = dropset->complexEvents;
      FOR_LIST(otherIter)
	APPEND(otherIter->value, elemsToCheck);
      allElems = elemsToCheck; 
    }

  BitVector
    *bipsSeen = CALLOC(GET_BITVECTOR_LENGTH(bipartitionsById->length), sizeof(BitVector));

  FOR_LIST(elemsToCheck)    
  {
    MergingEvent *me = (MergingEvent*)elemsToCheck->value;
    
    if(NOT me->computed)
      {
	getLostSupportThreshold(me, bipartitionsById);
	getSupportGainedThreshold(me, bipartitionsById);
	me->computed = TRUE; 
      }
    
    result -= me->supportLost;
    if(  me->supportGained
	 &&  NOT mergedBipVanishes(me, bipartitionsById, dropset->taxaToDrop) )
      result += me->supportGained;   
    
    if(me->isComplex)
      {
	IndexList *iI =  me->mergingBipartitions.many ;	
	FOR_LIST(iI)
	{
	  assert(NOT NTH_BIT_IS_SET(bipsSeen, iI->index));
	  if(NTH_BIT_IS_SET(bipsSeen, iI->index))
	    {
	      PR("problem:");
	      printIndexList(me->mergingBipartitions.many);
	      PR("at ");
	      printIndexList(dropset->taxaToDrop);		
	      PR("\n");
	      exit(0);
	    }
	  FLIP_NTH_BIT(bipsSeen, iI->index);	
	}
      }
    else
      {
	assert( NOT NTH_BIT_IS_SET(bipsSeen, me->mergingBipartitions.pair[0]));
	assert( NOT NTH_BIT_IS_SET(bipsSeen, me->mergingBipartitions.pair[1]));
	FLIP_NTH_BIT(bipsSeen,me->mergingBipartitions.pair[0]);
	FLIP_NTH_BIT(bipsSeen,me->mergingBipartitions.pair[1]);
      }
  }
  freeListFlat(allElems);
  
  
  /* handle vanishing bip */
  List *iter = consensusBipsCanVanish; 
  FOR_LIST(iter)
  {
    ProfileElem *elem = iter->value;
    
    switch(rogueMode)
      {
      case VANILLA_CONSENSUS_OPT :
	{
	  if(elem->treeVectorSupport > thresh
	     && NOT NTH_BIT_IS_SET(bipsSeen, elem->id)
	     && bipartitionVanishesP(elem,dropset))
	    result -= computeSupport ? elem->treeVectorSupport : 1;
	  break; 	
	}
      case ML_TREE_OPT: 
	{
	  if(elem->isInMLTree 
	     && NOT NTH_BIT_IS_SET(bipsSeen, elem->id)
	     && bipartitionVanishesP(elem,dropset))
	    result -= computeSupport ? elem->treeVectorSupport : 1; 
	  break; 
	}
      default: 
	assert(0);
      }
  }  

  free(bipsSeen);
  dropset->improvement = result;
}


List *getConsensusBipsCanVanish(Array *bipartitionProfile)
{
  List
    *consensusBipsCanVanish = NULL; 

  if(rogueMode == VANILLA_CONSENSUS_OPT
     || rogueMode == MRE_CONSENSUS_OPT)
    {
      int i; 
      FOR_0_LIMIT(i,bipartitionProfile->length)
	{
	  ProfileElem
	    *elem = GET_PROFILE_ELEM(bipartitionProfile, i);

	  if(NOT elem)
	    break;

	  if(elem->numberOfBitsSet - maxDropsetSize > 1 )
	    break;
      
	  if(elem->treeVectorSupport > thresh)
	    APPEND(elem,consensusBipsCanVanish);
	}
    }
  else if(ML_TREE_OPT)
    {
      int i; 
      FOR_0_LIMIT(i,bipartitionProfile->length)
	{
	  ProfileElem
	    *elem = GET_PROFILE_ELEM(bipartitionProfile, i);

	  if(NOT elem)
	    break;
      
	  if(elem->isInMLTree)
	    APPEND(elem,consensusBipsCanVanish);
	}
    }

  return consensusBipsCanVanish;
}


Dropset *evaluateEvents(HashTable *mergingHash, Array *bipartitionsById, Array *bipartitionProfile)
{
  Dropset 
    *result = NULL; 

  int i ; 

  List
    *consensusBipsCanVanish = getConsensusBipsCanVanish(bipartitionProfile);

  if( NOT mergingHash->entryCount)
    return NULL ;

  /* gather dropsets in array  */
  Array *allDropsets = CALLOC(1,sizeof(Array)) ;
  allDropsets->length = mergingHash->entryCount; 
  allDropsets->arrayTable = CALLOC(mergingHash->entryCount, sizeof(Dropset*));
       
  int cnt = 0; 
  HashTableIterator *htIter;
  FOR_HASH(htIter, mergingHash)
    {    
      GET_DROPSET_ELEM(allDropsets,cnt) = getCurrentValueFromHashTableIterator(htIter);
      cnt++; 
    }
  free(htIter); 
  assert(cnt == mergingHash->entryCount);
  
  /* compute MRE stuff  */
  if(rogueMode == MRE_CONSENSUS_OPT)
    { 
      
#ifdef PARALLEL
      numberOfJobs = allDropsets->length;
      globalPArgs->bipartitionsById =  bipartitionsById; 
      globalPArgs->allDropsets = allDropsets; 
      masterBarrier(THREAD_MRE, globalPArgs); 
#else
      
      FOR_0_LIMIT(i,allDropsets->length)	  
	{
	  Dropset *dropset =  GET_DROPSET_ELEM(allDropsets, i);
	  dropset->improvement =  getSupportOfMRETree(bipartitionsById, dropset) - cumScore;
	}
#endif     
    }
  

  /* evaluate dropsets */
  if(rogueMode != MRE_CONSENSUS_OPT)
    {
#ifdef PARALLEL
      numberOfJobs = allDropsets->length; 
      globalPArgs->mergingHash = mergingHash; 
      globalPArgs->allDropsets = allDropsets; 
      globalPArgs->bipartitionsById = bipartitionsById; 
      globalPArgs->consensusBipsCanVanish = consensusBipsCanVanish;
      masterBarrier(THREAD_EVALUATE_EVENTS, globalPArgs); 
#else
      FOR_0_LIMIT(i, allDropsets->length)
	{
	  Dropset *dropset =  GET_DROPSET_ELEM(allDropsets, i);   
	  evaluateDropset(mergingHash, dropset, bipartitionsById, consensusBipsCanVanish); 
	}
#endif
    }
  
  FOR_0_LIMIT(i,allDropsets->length)
    {      
      Dropset
	*dropset =  GET_DROPSET_ELEM(allDropsets, i);

      if(NOT result)
	result = dropset;
      else
	{
	  int drSize =  lengthIndexList(dropset->taxaToDrop),
	    resSize = lengthIndexList(result->taxaToDrop);
	  
	  double oldQuality =  labelPenalty == 0.0  
	    ? result->improvement * drSize 
	    :  (double)(result->improvement / (double)(computeSupport ?  numberOfTrees : 1.0)) - (double)resSize;
	  double newQuality = labelPenalty == 0.0 
	    ? dropset->improvement * resSize
	    : (double)(dropset->improvement / (double)(computeSupport ? numberOfTrees : 1.0)) - (double)drSize; 
	  
	  if( (newQuality  >  oldQuality)
	      || ((newQuality == oldQuality)
#ifdef MIN_DROPSETS
		  && drSize < resSize
#else
		  && drSize > resSize
#endif
		  ))
	    result = dropset;
	}
    }
  freeListFlat(consensusBipsCanVanish);

  free(allDropsets->arrayTable);
  free(allDropsets);

  if(labelPenalty == 0.0 && result->improvement > 0)
    return result;
  else if(labelPenalty != 0.0 && (result->improvement / (computeSupport ? numberOfTrees : 1.0) - lengthIndexList(result->taxaToDrop))  > 0.0 )
    return result; 
  else 
    return NULL;
}



void cleanup_updateNumBitsAndCleanArrays(Array *bipartitionProfile, Array *bipartitionsById, BitVector *mergingBipartitions, BitVector *newCandidates, Dropset *dropset)
{
  int profileIndex; 

  FOR_0_LIMIT(profileIndex,bipartitionProfile->length)
    {
      ProfileElem
	*elem = GET_PROFILE_ELEM(bipartitionProfile,profileIndex);
	      
      if( NOT elem )
	continue;
      
      /* check if number of bits has changed  */
      if(NOT NTH_BIT_IS_SET(mergingBipartitions,elem->id)) 
	{	  
	  if( mxtips - taxaDropped - 2 * elem->numberOfBitsSet <= 2 * maxDropsetSize )	  
	    FLIP_NTH_BIT(newCandidates, elem->id);
	  IndexList *iter = dropset->taxaToDrop;
	  boolean taxonDroppedP = FALSE;      
	  FOR_LIST(iter)
	  {
	    if(NTH_BIT_IS_SET(elem->bitVector, iter->index)) 
	      {
		taxonDroppedP = TRUE;
		UNFLIP_NTH_BIT(elem->bitVector, iter->index);
		elem->numberOfBitsSet--;
	      }
	  }

	  if(taxonDroppedP)
	    {
	      if(elem->numberOfBitsSet < 2)
		{ 
		  UNFLIP_NTH_BIT(newCandidates, elem->id);
		  FLIP_NTH_BIT(mergingBipartitions, elem->id);
		}	  
	      else
		FLIP_NTH_BIT(newCandidates, elem->id);
	    }
	}
      
      /* bip has been merged or vanished  */
      if(NTH_BIT_IS_SET(mergingBipartitions,elem->id)) 
	{
	  assert(NOT NTH_BIT_IS_SET(newCandidates, elem->id));
	  GET_PROFILE_ELEM(bipartitionProfile, profileIndex) = NULL;
	  GET_PROFILE_ELEM(bipartitionsById, elem->id) = NULL;
	  freeProfileElem(elem);
	}
    }  
}


BitVector *cleanup_applyAllMergerEvents(Array *bipartitionsById, Dropset *bestDropset, BitVector *mergingBipartitions)
{
  BitVector
    *candidateBips = CALLOC(GET_BITVECTOR_LENGTH(bipartitionsById->length), sizeof(BitVector)) ;

  if( bestDropset) 
    {
#ifdef PRINT_VERY_VERBOSE
      printDropset(bestDropset);
#endif
      
      List *iter = NULL ; 
      if(maxDropsetSize == 1)
	iter = bestDropset->ownPrimeE;
      else 
	iter = bestDropset->acquiredPrimeE; 
      FOR_LIST(iter)
      {
	int newBipId = cleanup_applyOneMergerEvent((MergingEvent*)iter->value, bipartitionsById, mergingBipartitions);
	FLIP_NTH_BIT(candidateBips, newBipId);
      }

      if(maxDropsetSize > 1 )
	{
	  iter = bestDropset->complexEvents;
	  FOR_LIST(iter)
	  {
	    int newBipId = cleanup_applyOneMergerEvent((MergingEvent*)iter->value, bipartitionsById, mergingBipartitions);
	    FLIP_NTH_BIT(candidateBips, newBipId);
	  }
	}
    } 
  
  return candidateBips;
}


void cleanup_rehashDropsets(HashTable *mergingHash, Dropset *bestDropset)
{
  if(maxDropsetSize == 1)
    return; 
  
  IndexList
    *taxaToDrop = bestDropset->taxaToDrop;

  List *allDropsets = NULL; 
  HashTableIterator *htIter; 
  FOR_HASH(htIter, mergingHash)
    {
      Dropset
	*dropset = getCurrentValueFromHashTableIterator(htIter);
      allDropsets = appendToList(dropset, allDropsets);
    }
  free(htIter);
  
  List *iter = allDropsets; 
  FOR_LIST(iter)
  {
    Dropset
      *dropset = (Dropset*)iter->value;

    if( NOT dropset)
      break;

    if(NOT dropset->ownPrimeE || isSubsetOf(dropset->taxaToDrop, taxaToDrop) )
      {
	removeElementFromHash(mergingHash, dropset);
	freeDropsetDeep(dropset, FALSE);
      }
    else if(haveIntersection(dropset->taxaToDrop, taxaToDrop)) /* needs reinsert */
      {
	removeElementFromHash(mergingHash, dropset);

#ifdef MYDEBUG 
	int length = lengthIndexList(dropset->taxaToDrop);
#endif    

	dropset->taxaToDrop = setMinusOf(dropset->taxaToDrop, taxaToDrop);

#ifdef MYDEBUG
	assert(length > lengthIndexList(dropset->taxaToDrop));
#endif
	unsigned int hv = mergingHash->hashFunction(mergingHash, dropset);
	Dropset *found = searchHashTable(mergingHash, dropset, hv);
	if( NOT found)
	  insertIntoHashTable(mergingHash,dropset,hv);
	else			/* reuse the merging events */
	  {
	    List
	      *iter, *next; 
	    for(iter = dropset->ownPrimeE; iter; iter = next)
	      {
		/* TODO potential error: double check, if this stuff did not already occur would be great */
		next = iter->next; 
		iter->next = found->ownPrimeE;
		found->ownPrimeE = iter;
	      }
	    freeIndexList(dropset->taxaToDrop);
	    free(dropset);
	  } 	
      }
  }
  freeListFlat(allDropsets);
}

BitVector *cleanup(All *tr, HashTable *mergingHash, Dropset *bestDropset, BitVector *candidateBips, Array *bipartitionProfile, Array *bipartitionsById)
{
  IndexList
    *ilIter;

  BitVector 
    *bipsToVanish = CALLOC(GET_BITVECTOR_LENGTH(bipartitionsById->length), sizeof(BitVector));

  /* apply merging events for best dropset  */
  candidateBips = cleanup_applyAllMergerEvents(bipartitionsById, bestDropset, bipsToVanish);
  
  if(NOT bestDropset)
    {
      free(bipsToVanish);
      return candidateBips;
    }
	  
  /* add to list of dropped taxa */
  ilIter = bestDropset->taxaToDrop;
  FOR_LIST(ilIter)
    FLIP_NTH_BIT(droppedTaxa,ilIter->index);

  /* remove merging bipartitions from arrays (not candidates) */
  cleanup_updateNumBitsAndCleanArrays(bipartitionProfile, bipartitionsById, bipsToVanish,candidateBips,bestDropset );
  removeElementFromHash(mergingHash, bestDropset);
  cleanup_mergingEvents(mergingHash, bipsToVanish, candidateBips, bipartitionProfile->length);

  cleanup_rehashDropsets(mergingHash, bestDropset);
  
#ifdef PRINT_VERY_VERBOSE
  int i; 
  PR("CLEAN UP: need to recompute bipartitions ");
  FOR_0_LIMIT(i, bipartitionProfile->length)
    if(NTH_BIT_IS_SET(candidateBips, i))
      PR("%d,", i);
  PR("\n");
#endif

#ifdef MYDEBUG	  
  debug_dropsetConsistencyCheck(mergingHash);
#endif

#ifdef PRINT_TIME
  PR("[%f] executed the merging events \n", updateTime(&timeInc));
#endif

#ifdef PRINT_VERY_VERBOSE
  PR("bips present %d (id) %d (profile)\n", getNumberOfBipsPresent(bipartitionsById), getNumberOfBipsPresent(bipartitionProfile));
#endif 

  cumScore += bestDropset->improvement;
  if(cumScore > bestCumEver)
    bestCumEver = cumScore;
  bestLastTime += bestDropset->improvement;
  dropsetPerRound[dropRound+1] = bestDropset; 
  cumScores[dropRound+1] = cumScore;

  printDropsetImprovement(bestDropset, tr, cumScore);

  return candidateBips; 
}


void doomRogues(All *tr, char *bootStrapFileName, char *dontDropFile, char *treeFile, boolean mreOptimisation, int rawThresh)
{
  double startingTime = gettime();
  timeInc = gettime();

  int 
    *indexByNumberBits,
    i;  

  FILE
    *bootstrapTreesFile = getNumberOfTrees(tr, bootStrapFileName),
    *rogueOutput = getOutputFileFromString("droppedRogues");

  BitVector
    *candidateBips;

  HashTable
    *mergingHash = NULL;  

  numberOfTrees = tr->numberOfTrees;

  if(strlen(treeFile))
    {
      rogueMode = ML_TREE_OPT;
      if(mreOptimisation)
	{
	  PR("ERROR: Please choose either support in the MRE consensus tree OR the bipartitions in the ML tree for optimization.\n");
	  exit(-1);
	}
      PR("mode: optimization of support of ML tree bipartitions in the bootstrap tree set.\n");
    }
  else if(mreOptimisation)
    {
      rogueMode = MRE_CONSENSUS_OPT;
      thresh = tr->numberOfTrees  * 0.5;
      PR("mode: optimization on MRE consensus tree. \n");
    }
  else 
    {
      rogueMode = VANILLA_CONSENSUS_OPT;
      thresh = tr->numberOfTrees * rawThresh / 100; 
      if(thresh == tr->numberOfTrees)
	thresh--; 
      PR("mode: optimization on consensus tree. Bipartition is part of consensus, if it occurs in more than %d trees\n", thresh); 
    }

  FILE
    *bestTree = (rogueMode == ML_TREE_OPT) ? myfopen(treeFile,"r") : NULL;

  mxtips = tr->mxtips;
  tr->bitVectorLength = GET_BITVECTOR_LENGTH(mxtips);

  Array 
    *bipartitionProfile = getOriginalBipArray(tr, bestTree, bootstrapTreesFile);

  if(maxDropsetSize >= mxtips - 3)
    {
      PR("\nMaximum dropset size (%d) too large. If we prune %d taxa, then there \n\
 will be no bipartitions left and thus such a pruned tree set can never \n\
 have a higher information content (in terms of RBIC) than the original \n\
 tree.\n", maxDropsetSize, mxtips-3);
      exit(-1);
    }

  dropsetPerRound = CALLOC(mxtips, sizeof(Dropset*)); 
  Dropset    
    *bestDropset = NULL;   

  neglectThose = neglectThoseTaxa(tr, dontDropFile);

  initializeRandForTaxa(mxtips);

  treeVectorLength = GET_BITVECTOR_LENGTH(tr->numberOfTrees);
  bitVectorLength = GET_BITVECTOR_LENGTH(tr->mxtips);
  droppedTaxa = CALLOC(bitVectorLength, sizeof(BitVector));

  paddingBits = CALLOC(GET_BITVECTOR_LENGTH(mxtips), sizeof(BitVector));
  for(i = mxtips; i < GET_BITVECTOR_LENGTH(mxtips) * MASK_LENGTH; ++i)
    FLIP_NTH_BIT(paddingBits,i);

  FOR_0_LIMIT(i,bipartitionProfile->length)
    {
      ProfileElem *elem = ((ProfileElem**)bipartitionProfile->arrayTable)[i];
      elem->numberOfBitsSet = genericBitCount(elem->bitVector, bitVectorLength);
    }

  Array
    *bipartitionsById = CALLOC(1,sizeof(Array)); 
  bipartitionsById->arrayTable = CALLOC(bipartitionProfile->length, sizeof(ProfileElem*));
  bipartitionsById->length = bipartitionProfile->length;
  FOR_0_LIMIT(i,bipartitionsById->length)
    GET_PROFILE_ELEM(bipartitionsById,i) = GET_PROFILE_ELEM(bipartitionProfile, i);
  qsort(bipartitionsById->arrayTable, bipartitionsById->length, sizeof(ProfileElem**), sortById);

  numBips = bipartitionProfile->length;

  cumScore = getInitScore(bipartitionProfile);
  cumScores = CALLOC(mxtips-3, sizeof(int));  
  cumScores[0]  = cumScore;
  bestCumEver = cumScore;

  bestLastTime = cumScore;
  fprintf(rogueOutput, "num\ttaxNum\ttaxon\trawImprovement\tRBIC\n");
  fprintf(rogueOutput, "%d\tNA\tNA\t%d\t%f\n", 0, 0, (double)cumScore /(numberOfTrees * (mxtips-3)) );

  PR("[%f] initialisation done (initScore = %f, numBip=%d)\n", updateTime(&timeInc), (double)cumScore / (double)((tr->mxtips-3) * (computeSupport ? tr->numberOfTrees : 1 ) ), bipartitionsById->length);

  boolean firstMerge= TRUE;
  candidateBips = CALLOC(GET_BITVECTOR_LENGTH(bipartitionProfile->length),sizeof(BitVector));
  FOR_0_LIMIT(i,bipartitionProfile->length)
    FLIP_NTH_BIT(candidateBips,i);

  mergingHash = createHashTable(tr->mxtips * maxDropsetSize * HASH_TABLE_SIZE_CONST,
				NULL,
				dropsetHashValue, 
				dropsetEqual); 

   

#ifdef PARALLEL
  globalPArgs = CALLOC(1,sizeof(parallelArguments));   
  startThreads();
#endif

  /* main loop */
  do 
    {
#ifdef PRINT_VERY_VERBOSE
      PR("ROUND %d ================================================================================================================================================================================================================\n",dropRound);
      PR("dropped vector is: ");
      printBitVector(droppedTaxa, GET_BITVECTOR_LENGTH(mxtips));
      PR("\n");
#endif
      
      /***********/
      /* prepare */
      /***********/
      bestDropset = NULL;
      unifyBipartitionRepresentation(bipartitionProfile,droppedTaxa); 
      indexByNumberBits = createNumBitIndex(bipartitionProfile, mxtips);

#ifdef PRINT_TIME
      PR("[%f] sorting bipartition profile\n", updateTime(&timeInc));
#endif

#ifdef PRINT_VERY_VERBOSE
      printBipartitionProfile(bipartitionProfile);
#endif

      /***********************************/
      /* create / update  merging events */
      /***********************************/
#ifdef PARALLEL
      numberOfJobs = bipartitionProfile->length;
      globalPArgs->mergingHash = mergingHash; 
      globalPArgs->candidateBips = candidateBips; 
      globalPArgs->bipartitionsById = bipartitionsById; 
      globalPArgs->bipartitionProfile = bipartitionProfile ; 
      globalPArgs->indexByNumberBits = indexByNumberBits; 
      globalPArgs->firstMerge = firstMerge; 
      masterBarrier(THREAD_GET_EVENTS, globalPArgs);
      free(candidateBips);      
#else 
      createOrUpdateMergingHash(tr, mergingHash, bipartitionProfile, bipartitionsById, candidateBips, firstMerge, indexByNumberBits );
#endif
      firstMerge = FALSE;      

#ifdef MYDEBUG
      debug_dropsetConsistencyCheck(mergingHash);
      debug_mergingHashSanityCheck(mergingHash, bipartitionProfile->length);
#endif

#ifdef PRINT_TIME
      PR("[%f] computed / updated events\n", updateTime(&timeInc));
#endif

      /* clear */

      /******************/
      /* combine events */
      /******************/
      if(maxDropsetSize > 1)
	mergingHash = combineMergerEvents(mergingHash, bipartitionsById);

#ifdef PRINT_TIME
      PR("[%f] combined events\n", updateTime(&timeInc));
#endif

#ifdef PRINT_VERY_VERBOSE
      if(mergingHash->entryCount > 0)
      	printMergingHash(mergingHash);
#endif

      /**********************/
      /* evaluate dropsets  */
      /**********************/
      bestDropset = evaluateEvents(mergingHash, bipartitionsById, bipartitionProfile);
      free(indexByNumberBits);

#ifdef PRINT_TIME
      PR("[%f] calculated per dropset improvement\n", updateTime(&timeInc));
#endif

      /*****************/
      /*  cleanup      */
      /*****************/
      candidateBips = cleanup(tr, mergingHash, bestDropset, candidateBips, bipartitionProfile, bipartitionsById); 

#ifdef MYDEBUG
      int l,m;
      FOR_0_LIMIT(l,bipartitionProfile->length)
	{
	  ProfileElem
	    *elemA = GET_PROFILE_ELEM(bipartitionProfile,l);

	  if(NOT elemA )
	    continue;

	  for(m = l+1; m < bipartitionProfile->length; ++m)
	    {
	      ProfileElem
		*elemB = GET_PROFILE_ELEM(bipartitionProfile,m);

	      if( NOT elemB)
		continue;

	      if(elemA->numberOfBitsSet == elemB->numberOfBitsSet && myBitVectorEqual(elemA,elemB))
		{
		  PR("%d and %d are equal!\n", elemA->id, elemB->id);
		  printBitVector(elemA->bitVector, bitVectorLength);
		  PR("\n");
		  printBitVector(elemB->bitVector, bitVectorLength);
		  PR("\n");
		  /* assert(0); */
		  exit(-1);
		}
	    }
	}
#endif

#ifdef PRINT_VERY_VERBOSE
      PR("dropped vector is: ");
      printBitVector(droppedTaxa, GET_BITVECTOR_LENGTH(mxtips));
      PR("\n");
#endif
      if(bestDropset)
	taxaDropped += lengthIndexList(bestDropset->taxaToDrop);      

      dropRound++;      
    } while(bestDropset);
  
  /* print out result */  
  printRogueInformationToFile(tr, rogueOutput, bestCumEver,cumScores, dropsetPerRound);

  PR("total time elapsed: %f\n", updateTime(&startingTime));

  /* free everything */   
  FOR_0_LIMIT(i, bipartitionProfile->length)
    {
      ProfileElem *elem = GET_PROFILE_ELEM(bipartitionProfile,i);
      if(elem)
  	freeProfileElem(elem);
    }
  free(((ProfileElemAttr*)bipartitionProfile->commonAttributes));
  freeArray(bipartitionProfile);
  freeArray(bipartitionsById);
  destroyHashTable(mergingHash, freeDropsetDeepInHash);
  /* destroyHashTable(mergingHash, NULL); */

  fclose(rogueOutput);
  for(i= 0 ; i < dropRound + 1; ++i)
    {
      Dropset *theDropset = dropsetPerRound[i];
      if(theDropset)
	freeDropsetDeepInEnd(theDropset);
    }
  free(dropsetPerRound);
  free(neglectThose);
  free(cumScores);
  free(paddingBits);
  free(randForTaxa);
  free(droppedTaxa);
  free(candidateBips);
}


void printHelpFile()
{
  printVersionInfo(FALSE);
  printf("This program implements the RogueNaRok algorithm for rogue taxon identification.\n\nSYNTAX: ./%s -i <bootTrees> -n <runId> [-x <excludeFile>] [-c <threshold>] [-b] [-s <dropsetSize>] [-w <workingDir>] [-h]\n", programName);
  printf("\n\tOBLIGATORY:\n");
  printf("-i <bootTrees>\n\tA collection of bootstrap trees.\n");
  printf("-n <runId>\n\tAn identifier for this run.\n");
  printf("\n\tOPTIONAL:\n");
  printf("-t <bestKnownTree>\n\tIf a single best-known tree (such as an ML or MP\n\t\
tree) is provided, RogueNaRok optimizes the bootstrap support in this\n\t\
best-known tree (still drawn from the bootstrap trees). The threshold\n\t\
parameter is ignored.\n");
  printf("-x <excludeFile>\n\ttaxa in this file (one taxon per line) will not be\n\t\
considered for pruning.\n");
  printf("-c <threshold>\n\t A threshold or mode for the consensus tree that is\n\t\
optimized. Specify a value between 50 (majority rule consensus) and\n\t\
100 (strict consensus) or MR (for the extended majority rule\n\t\
consensus). Note that rogue taxa identified with respect to different\n\t\
thresholds can vary substantially. DEFAULT: MR consensus\n");
  printf("-b\n\tInstead of trying to maximize the support in the consensus tree,\n\t\
the RogueNaRok will try to maximize the number of bipartition in the\n\t\
final tree by pruning taxa. DEFAULT: off\n");
  printf("-L <factor>\n\ta weight factor to penalize for dropset size. \n\t\
Factor=1 is Pattengale's criterion. The higher the value, the more \n\t\
conservative the algorithm is in pruning taxa. DEFAULT: 0.0 (=RBIC)\n");
  printf("-s <dropsetSize>\n\tmaximum size of dropset per iteration. If\n\t\
dropsetSize == n, then RogueNaRok will test in each iteration which\n\t\
tuple of n taxa increases optimality criterion the most and prunes\n\t\
taxa accordingly. This improves the result, but runtimes will\n\t\
increase at least linearly. DEFAULT: 1\n");
  printf("-w <workDir>\n\tA working directory where output files are created.\n");
  printf("-T <num>\n\tExecute RogueNaRok in parallel with <num> threads. You need to compile the program for parallel execution first.\n");
  printf("-h\n\tThis help file.\n");
  printf("\nMINIMAL EXAMPLE:\n./%s -i <bootstrapTreeFile> -n run1\n", programName);
}



int main(int argc, char *argv[])
{
  int
    c,
    threshold = 50;

  char
    *excludeFile = "", 
    *bootTrees = "",
    *treeFile = ""; 

  boolean
    mreOptimisation = FALSE;

  if(sizeof(int) != 4)
    {
      printf("I am sorry, RogueNaRok currently does not support your computer architecture. The code assumes that an integer (type int) consists of 4 bytes.\n");
      assert(sizeof(int) == 4);
    }


  programName = PROG_NAME;
  programVersion = PROG_VERSION;
  programReleaseDate  = PROG_RELEASE_DATE;
  
  while ((c = getopt (argc, argv, "i:t:n:x:whc:s:bT:L:")) != -1)
    switch (c)
      {
      case 'i':
	bootTrees = optarg;
	break;
      case 'T':
	{
#ifndef PARALLEL
	  printf("\n\nFor running RogueNaRok in parallel, please compile with \n\n"); 
	  exit(-1);	  
#else
	  numberOfThreads = wrapStrToL(optarg); 
#endif
	  break; 
	}
      case 'b':
	computeSupport = FALSE;
	break;
      case 'n':
	strcpy(run_id, optarg);
	break;
      case 't':
	treeFile = optarg;
	break;
      case 's':
	maxDropsetSize = wrapStrToL(optarg);
	break;
      case 'x': 
	excludeFile = optarg;
	break;
      case 'w':
	strcpy(workdir, optarg) ; 
	break;
      case 'L':
	labelPenalty = wrapStrToDouble(optarg); 
	break; 
      case 'c':
	{
	  if( NOT strcmp(optarg, "MRE"))
	    {
	      mreOptimisation = TRUE;
	      threshold = 50; 
	    }
	  else
	    threshold = wrapStrToL(optarg);
	  break;
	}
      case 'h':
      default:	
	{
	  printHelpFile();
	  abort ();
	}
      }
  
  /* initialize fast bit counting */
  compute_bits_in_16bits();
  initializeMask();

#ifdef PARALLEL
  if(NOT numberOfThreads)
    {
      printf("\n\nPlease specify the number of threads for parallel execution with -T\n\n");
      exit(-1);
    }
  if(numberOfThreads == 1)
    {
      printf("\n\nCalling parallel version of RogueNaRok with 1 thread is deprecated.\n\
Please compile a sequential version of RogueNaRok instead.\n\n");
      exit(-1);
    }
#endif

  if( NOT strcmp(treeFile, ""))
    rogueMode = ML_TREE_OPT;

  if( NOT strcmp(bootTrees, ""))
    {
      printf("ERROR: Please specify a file containing bootstrap trees via -i.\n");
      printHelpFile();
      exit(-1);
    }  

  if( NOT strcmp(run_id, ""))
    {
      printf("ERROR: Please specify a run-id via -n\n");
      printHelpFile();
      exit(-1);
    }

  if(threshold < 50)
    {
      printf("ERROR: Only accepting threshold values between 50 (MR) and 100 (strict).\n");
      exit(-1);
    }

  if(threshold != 50 &&  strcmp(treeFile, "") )    
    {
      printf("ERROR: threshold option -c not available in combination with best-known tree.\n");
      exit(-1);
    }

  All 
    *tr = CALLOC(1,sizeof(All));  
  setupInfoFile();
  if  (NOT setupTree(tr, bootTrees))
    {
      PR("Something went wrong during tree initialisation. Sorry.\n");
      exit(-1);
    }   

  doomRogues(tr,
  	     bootTrees,
  	     excludeFile,
  	     treeFile,
  	     mreOptimisation,
	     threshold);

  freeTree(tr);
  free(mask32);
  free(infoFileName);

  return 0; 
}
