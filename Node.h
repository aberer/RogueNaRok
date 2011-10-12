#ifndef NODE_H
#define NODE_H

#include "common.h"
#include "List.h"
#include "HashTable.h"

typedef struct _node
{
  int id; 
  boolean visited; 
  IndexList *edges;
} Node; 

boolean nodeEqual(HashTable *hashTable, void *entryA, void *entryB);
unsigned int nodeHashValue(HashTable *hashTable, void *value);
void freeNode(void *value);
IndexList *findAnIndependentComponent(HashTable *allNodes, Node *thisNode);


#endif
