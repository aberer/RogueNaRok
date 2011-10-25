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

#ifndef LIST_H
#define LIST_H

#include <stdlib.h>
#include "common.h"


typedef struct List_{
  void *value; 			
  struct List_ *next; 
} List;

typedef struct _indexList
{
  struct _indexList *next;
  int index;
} IndexList;


void freeList(List* list);
IndexList *doubleAppendToIndexList(int valueA, int valueB, IndexList *list);
boolean isSubsetOf(IndexList *subset, IndexList *set);
List* appendToList(void *value, List *list);
IndexList *appendToIndexListIfNotThere(int elem, IndexList *list);
IndexList *appendToIndexListIfNotThere2(int elem, IndexList *list);
void freeIndexList(IndexList *list);
IndexList* appendToIndexList(int value, IndexList *list) ;
boolean indexListEqual(IndexList *aList, IndexList *bList);
boolean isInIndexList(int index, IndexList *list);
int isInIndexListSpecial(int a, int b, IndexList *list);
void freeIndexList(IndexList *list);
IndexList* appendToIndexList(int value, IndexList *list) ;
void freeIndexList(IndexList *list);
IndexList *concatenateIndexList(IndexList *listA, IndexList *listB) ;
List *concatenateLists(List *listA, List *listB);
void printIndexList(IndexList *list);
IndexList *findFirstCommonElem(IndexList *listA, IndexList *listB);
boolean elemIsInIndexList(int index, IndexList *list);
List *joinLists(List* a, List*b) ;
IndexList *joinIndexListsNonRedundant(IndexList *listA, IndexList *listB);
void freeListFlat(List* list);
boolean indexListContainsIndexListUnordered(IndexList *list, IndexList *subList);
boolean haveIntersection(IndexList *listA, IndexList *listB);
int length_indexList(IndexList *list);
IndexList *setMinusOf(IndexList *list, IndexList *subtract);
int lengthIndexList(IndexList *list);
void printIndexListToFile(FILE *file, IndexList *list);
int lengthOfList(List *list);
boolean isSubsetOfReverseOrdered(IndexList *subset, IndexList *set); 

#define FOR_LIST(iter) for(;iter;iter = iter->next)
#define APPEND(elem,list) list = appendToList(elem,list)
#define APPEND_INT(elem,list) list = appendToIndexList(elem,list)
#define CONCAT_LISTS(lista,listb)  lista = concatenateLists(lista,listb)

#endif 
