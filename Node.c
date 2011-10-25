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
