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

#include "HashTable.h"
#ifdef PARALLEL
#include <pthread.h>
#endif

HashTable *createHashTable(unsigned int size, 
			   void *commonAttr, 
			   unsigned int (*hashFunction)(HashTable *hash_table, void *value),
			   boolean (*equalFunction)(HashTable *hash_table, void *entryA, void *entryB))
{  
  static const unsigned int 
    initTable[] = {64, 128, 256, 512, 1024, 2048, 4096, 
		   8192, 16384, 32768, 65536, 131072, 
		   262144, 524288, 1048576, 2097152,
		   4194304, 8388608, 16777216, 33554432, 
		   67108864, 134217728, 268435456, 
		   536870912, 1073741824, 2147483648U};
  
  HashTable 
    *hashTable = CALLOC(1, sizeof(HashTable));
  
  unsigned int
    tableSize,
    i,
#ifdef DEBUG
    maxSize = (unsigned int)-1,    
#endif
    primeTableLength = sizeof(initTable)/sizeof(initTable[0]);

  hashTable->hashFunction = hashFunction;
  hashTable->equalFunction = equalFunction;
  hashTable->commonAttributes = commonAttr;

#ifdef DEBUG 
  assert(size <= maxSize);  
#endif
  for(i = 0; initTable[i] < size && i < primeTableLength; ++i);
  assert(i < primeTableLength);

  tableSize = initTable[i];

#ifdef PARALLEL
  hashTable->cntLock = CALLOC(1,sizeof(pthread_mutex_t));
  pthread_mutex_init(hashTable->cntLock, (pthread_mutexattr_t *)NULL);
  
  hashTable->lockPerSlot = (pthread_mutex_t**)CALLOC(tableSize,sizeof(pthread_mutex_t*));  
  FOR_0_LIMIT(i, tableSize)
    {
      hashTable->lockPerSlot[i] = (pthread_mutex_t*) CALLOC(1,sizeof(pthread_mutex_t));
      pthread_mutex_init(hashTable->lockPerSlot[i], (pthread_mutexattr_t *)NULL);
    }
#endif

  hashTable->table = CALLOC(tableSize, sizeof(HashElem*));
  hashTable->tableSize = tableSize;  
  hashTable->entryCount = 0;  

  return hashTable;
}


/* remove without destroying value of hashelem */
boolean removeElementFromHash(HashTable *hashtable, void *value)
{
  unsigned int
    hashValue = hashtable->hashFunction(hashtable, value),
    position = hashValue % hashtable->tableSize;
  
  HashElem *elem = hashtable->table[position];  

  if( NOT elem )
    {
      assert(0);
      return FALSE;
    }

  if(
     elem->fullKey == hashValue &&
     hashtable->equalFunction(hashtable, elem->value, value))
    {
      hashtable->table[position] = elem->next ;
      free(elem);
      hashtable->entryCount-- ; 

      return TRUE; 
    }
  
  while(elem->next)
    {
      if(elem->next->fullKey == hashValue && hashtable->equalFunction(hashtable, elem->next->value, value))
	{
	  void *nextOne = elem->next->next; 
	  free(elem->next);
	  elem->next = nextOne;
	  hashtable->entryCount-- ; 
	  return TRUE; 
	}
      elem = elem->next;
    }

  return FALSE;
}


void *searchHashTableWithInt(HashTable *hashtable, unsigned int hashValue)
{
  unsigned int 
    position = hashValue % hashtable->tableSize;
  
  HashElem 
    *elem;
  
  for(elem = hashtable->table[position]; 
      elem; 
      elem = elem->next)
    if(elem->fullKey  == hashValue)
      return  elem->value; 

  return NULL;
}


void *searchHashTable(HashTable *hashtable, void *value, unsigned int hashValue)
{
  unsigned int 
    position = hashValue % hashtable->tableSize;
  
  HashElem 
    *elem;
  
  for(elem = hashtable->table[position]; 
      elem; 
      elem = elem->next)
    if(elem->fullKey  == hashValue && 
       hashtable->equalFunction(hashtable, elem->value, value))
      return  elem->value; 

  return NULL;
}


void insertIntoHashTable(HashTable *hashTable, void *value, unsigned int index)
{  
  /* just copied this */
  HashElem 
    *hashElem = CALLOC(1, sizeof(HashElem));
  
  hashElem->fullKey = index;
  
  index =  hashElem->fullKey % hashTable->tableSize;
  
  hashElem->value = value;
  hashElem->next = hashTable->table[index];  
  hashTable->table[index] = hashElem;
#ifdef PARALLEL
  pthread_mutex_lock(hashTable->cntLock);
  hashTable->entryCount++;
  pthread_mutex_unlock(hashTable->cntLock);
#else
  hashTable->entryCount++;
#endif
}


void destroyHashTable(HashTable *hashTable, void (*freeValue)(void *value))
{
  unsigned 
    int i; 
  
  HashElem 
    *elemA, 
    *elemB,
    **table = hashTable->table;
  
  for(i = 0; i < hashTable->tableSize; ++i)
    {
      elemA = table[i];
      while(elemA != NULL)
	{
	  elemB = elemA; 
	  elemA = elemA->next; 
	  if(freeValue)
	    freeValue(elemB->value);
	  free(elemB);
	}
    }

#ifdef PARALLEL
  pthread_mutex_destroy(hashTable->cntLock);
  free(hashTable->cntLock);
  
  FOR_0_LIMIT(i,hashTable->tableSize)
    {
      pthread_mutex_destroy(hashTable->lockPerSlot[i]);
      free(hashTable->lockPerSlot[i]);
    }
  free(hashTable->lockPerSlot);
#endif

  free(hashTable->commonAttributes);  
  free(hashTable->table);
  free(hashTable);
}


void updateEntryCount(HashTable *hashTable)
{
  unsigned int 
    i, 
    result = 0;

  for(i = 0; i < hashTable->tableSize; ++i)
    {
      HashElem 
	*elem = ((HashElem**)hashTable->table)[i];
      
      while(elem)
	{
	  result++;
	  elem = elem->next;
	}
    }

  hashTable->entryCount = result;
}


HashTableIterator *createHashTableIterator(HashTable *hashTable) 
{
  unsigned 
    int i; 
  
  HashTableIterator 
    *hashTableIterator = CALLOC(1, sizeof(HashTableIterator));
  
  hashTableIterator->hashTable = hashTable;
  hashTableIterator->hashElem = NULL;
  hashTableIterator->index = hashTable->tableSize;
  
  if( NOT hashTable->entryCount)
    return hashTableIterator;

  FOR_0_LIMIT(i,hashTable->tableSize)
    {
      if(hashTable->table[i])
	{
	  hashTableIterator->hashElem = hashTable->table[i];
	  hashTableIterator->index = i;
	  break;
	}
    }
  
  return hashTableIterator;
}

boolean hashTableIteratorNext(HashTableIterator *hashTableIterator)
{
  unsigned int 
    i, 
    tableSize = hashTableIterator->hashTable->tableSize;
  
  HashElem 
    *next = hashTableIterator->hashElem->next;
  
  if(next)
    {
      hashTableIterator->hashElem = next;
      return TRUE;
    }

  i = hashTableIterator->index + 1;
  
  if(i >= tableSize)
    {
      hashTableIterator->index = i;
      return FALSE;
    }
  
  while( NOT (next = hashTableIterator->hashTable->table[i]))
    {
      if( ++i >= tableSize) 
	{
	  hashTableIterator->index = i; 
	  return FALSE;
	}
    }
  
  hashTableIterator->index = i;
  hashTableIterator->hashElem = next;

  return next != NULL ;
}


void *getCurrentValueFromHashTableIterator(HashTableIterator *hashTableIterator)
{  
  return ((hashTableIterator->hashElem) 
	  ?  hashTableIterator->hashElem->value
	  : NULL);
}

