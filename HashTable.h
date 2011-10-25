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

#ifndef HASHTABLE_H
#define HASHTABLE_H

#include <assert.h>
#include <stdlib.h>
#include <string.h>

#include "common.h"


typedef struct hash_el
{
  unsigned int fullKey;
  void *value;
  struct hash_el *next;
} HashElem;

typedef struct hash_table
{
  unsigned int tableSize;
  unsigned int entryCount;
  void *commonAttributes;
  unsigned int (*hashFunction)(struct hash_table *h, void *value);
  boolean (*equalFunction)(struct hash_table *hashtable, void *entrA, void *entryB);
  HashElem **table;
#ifdef PARALLEL
  pthread_mutex_t **lockPerSlot;
  pthread_mutex_t *cntLock; 
#endif
} HashTable;

typedef struct 
{
  HashTable *hashTable;
  HashElem *hashElem;
  unsigned int index;
} HashTableIterator;

#define FOR_HASH(htIter, hashtable)  boolean hasNext = TRUE; for(htIter = createHashTableIterator(hashtable); htIter && hasNext ; hasNext = hashTableIteratorNext(htIter))
#define FOR_HASH_2(htIter, hashtable)  hasNext = TRUE; for(htIter = createHashTableIterator(hashtable); hasNext ; hasNext = hashTableIteratorNext(htIter))

HashTable *createHashTable(unsigned int size, void *commonAttr, unsigned int (*hashFunction)(HashTable *hash_table, void *value), boolean (*equalFunction)(HashTable *hash_table, void *entryA, void *entryB));
HashTableIterator *createHashTableIterator(HashTable *hashTable); 
void *getCurrentValueFromHashTableIterator(HashTableIterator *hashTableIterator);
boolean hashTableIteratorNext(HashTableIterator *hashTableIterator);
void *searchHashTableWithInt(HashTable *hashtable, unsigned int hashValue);
void *searchHashTable(HashTable *hashtable, void *value, unsigned int hashValue);
void insertIntoHashTable(HashTable *hashTable, void *value, unsigned int index);
boolean removeElementFromHash(HashTable *hashtable, void *value);
void destroyHashTable(HashTable *hashTable, void (*freeValue)(void *value));

#endif
