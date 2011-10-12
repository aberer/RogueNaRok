#include "Node.h"


IndexList *findAnIndependentComponent(HashTable *allNodes, Node *thisNode)
{
  if(thisNode->visited)
    return NULL; 

  IndexList *iter  = thisNode->edges;   
  thisNode->visited = TRUE;
  IndexList *result = NULL; 
  APPEND_INT(thisNode->id, result);

  FOR_LIST(iter)
  {
    Node *found = searchHashTableWithInt(allNodes, iter->index);    
    
    if(  NOT found->visited)
      {
	IndexList *list = findAnIndependentComponent(allNodes, found);
	result = concatenateIndexList(list, result);
      }
  }
  
  return result; 
}


void freeNode(void *value)
{
  Node *node = value; 
  freeIndexList(node->edges);
  free(node);
}



boolean nodeEqual(HashTable *hashTable, void *entryA, void *entryB)
{
  return ((Node*)entryA)->id  == ((Node*)entryB)->id;  
}


unsigned int nodeHashValue(HashTable *hashTable, void *value)
{
  return ((Node*)value)->id; 
}
