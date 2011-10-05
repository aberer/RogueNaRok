#include "List.h"



boolean isInIndexList(int index, IndexList *list)
{
  IndexList *iter = list;
  FOR_LIST(iter)
    if(iter->index == index)
      return TRUE;

  return FALSE;
}



IndexList *setMinusOf(IndexList *list, IndexList *subtract)
{
  IndexList
    *iter = list,
    *result = NULL; 
  
  FOR_LIST(iter)
  {
    IndexList
      *iter2 = subtract; 
   
    boolean
      isThere = FALSE; 

    FOR_LIST(iter2)
      isThere |= iter->index == iter2->index;

    if(NOT isThere)
      APPEND_INT(iter->index, result);
  }
  freeIndexList(list);

  return result;
}


/* TODO-opt */
#ifdef NEW_SUBSET
boolean isSubsetOfReverseOrdered(IndexList *subset, IndexList *set)
{
  IndexList
    *subsetIter = subset;
  boolean found = FALSE;

  FOR_LIST(subsetIter)
    {
      IndexList
	*setIter = set; 
      int ind = subsetIter->index; 
      FOR_LIST(setIter)
      {
	if(ind > setIter->index)
	  return FALSE;
	else if((found = (ind == setIter->index)))
	  break; 
      }
      
      if(NOT found)
	return FALSE;
    }

  return TRUE; 
}
#else
boolean isSubsetOf(IndexList *subset, IndexList *set)
{
  boolean result = TRUE; 

  FOR_LIST(subset)
  {
    result &= isInIndexList(subset->index, set);
    if( NOT result)
      return result; 
  }

  return result;
}
#endif







/* TODO-opt?  */
IndexList *appendToIndexListIfNotThere(int elem, IndexList *list)
{
  IndexList *iter = list;

  boolean contains = FALSE;
  FOR_LIST(iter)  
  {
    contains |= elem == iter->index;
    if(contains)
      break;
  }

  if( NOT contains)
    return appendToIndexList(elem, list);
  else
    return list;
}



IndexList *appendToIndexListIfNotThere2(int elem, IndexList *list)
{
  IndexList *iter = list;

  boolean contains = FALSE;
  FOR_LIST(iter)  
  {
    contains |= elem == iter->index;
    if(contains)
      break;
  }

  if( NOT contains)
    {
      PR("FAIL: it was not there\n");
      return appendToIndexList(elem, list);
    }
  else
    return list;
}



boolean haveIntersection2(IndexList *listA, IndexList *listB)
{
  FOR_LIST(listA)
  {
    IndexList *iter = listB; 
    FOR_LIST(iter)
      if(listA->index == iter->index)
	return TRUE;
  }

  return FALSE;
}

boolean haveIntersection(IndexList *listA, IndexList *listB)
{  
  FOR_LIST(listA)
  {
    IndexList *iter = listB; 
    FOR_LIST(iter)
      if(listA->index == iter->index)
	return TRUE;
  }

  return FALSE;
}


void freeIndexList(IndexList *list)
{
  IndexList
    *iter;
  
  for(iter = list; iter;)
    {
      list = iter->next;
      free(iter);
      iter = list;
    }
}


boolean indexListContainsIndexListUnordered(IndexList *list, IndexList *subList)
{  
  FOR_LIST(subList)
  {
    boolean contained = FALSE; 
    IndexList
      *iter = list;
    FOR_LIST(iter)
      if(iter->index == subList->index)
	{
	  contained = TRUE;
	  break;
	}
    
    if( NOT contained)
      return FALSE;
  }
  
  return TRUE;
}


List *joinLists(List* a, List*b) 
{
  List *iter; 
  
  if( NOT a)
    return b;

  for(iter = a; iter ; iter = iter->next)
    {
      if( NOT iter->next)
	{
	  iter->next = b;
	  return a;
	}
    }

  assert(0);
  return NULL; 
}


int getMax_indexList(IndexList *list)
{
  IndexList
    *iter = list;
  int
    result = 0;

  FOR_LIST(iter)
    if(result < iter->index)
      result = iter->index;
  
  return result;
}

int lengthIndexList(IndexList *list)
{
  IndexList
    *iter = list;
  
  int result = 0;
  FOR_LIST(iter)
    result++;

  return result;
}


void printIndexList(IndexList *list)
{  
  PR(">");
  FOR_LIST(list)
    PR("%d,", list->index);
  PR("<");
}

void printIndexListToFile(FILE *file, IndexList *list)
{
  boolean isFirst = TRUE;
  FOR_LIST(list)
    if(isFirst)
      {
	fprintf(file, "%d", list->index);
	isFirst = FALSE;
      }
    else
      fprintf(file, ",%d", list->index);
}


IndexList *findFirstCommonElem(IndexList *listA, IndexList *listB)
{
  IndexList 
    *iterA = listA, 
    *iterB = listB;

  FOR_LIST(iterA)
  {
    iterB = listB; 
    FOR_LIST(iterB)
      if(iterA->index == iterB->index)
	return iterA;
  }

  return NULL;
}


IndexList *concatenateIndexList(IndexList *listA, IndexList *listB) 
{
  IndexList
    *iterA = listA; 
  
  if( NOT iterA)
    return listB;
  
  FOR_LIST(iterA)
    if( NOT iterA->next)
      {
	iterA->next = listB;
	return listA;
      }

  return NULL;
}



IndexList *doubleAppendToIndexList(int valueA, int valueB, IndexList *list)
{
  IndexList
    *listElemB = CALLOC(1, sizeof(IndexList)),
    *listElemA =  CALLOC(1, sizeof(IndexList));

  listElemA->index = valueA; 
  listElemB->index = valueB; 
  
  listElemA->next = list; 
  listElemB->next = listElemA; 

  return listElemB; 
}

IndexList* appendToIndexList(int value, IndexList *list) 
{  
  IndexList 
    *listElem = CALLOC(1, sizeof(IndexList));
  
  listElem->index = value;
  listElem->next = list;
  
  return listElem;
}


boolean indexListEqual(IndexList *aList, IndexList *bList)
{
  IndexList *iterA = aList;

  FOR_LIST(iterA)
  {
    boolean indexFound = FALSE;
    IndexList *iterB = bList; 
    FOR_LIST(iterB)
    {
      indexFound = iterA->index == iterB->index;
      if(indexFound)
	break;
    }

    if( NOT indexFound)
      return FALSE;
  }

  return lengthIndexList(aList) == lengthIndexList(bList);
}


int lengthOfList(List *list)
{
  int length = 0 ; 
  List
    *iter = list; 
  FOR_LIST(iter)
    length++;
  return length; 
}


List* appendToList(void *value, List *list)
{  
  List 
    *listElem = CALLOC(1, sizeof(List));
  
  listElem->value = value;
  listElem->next = list;
  
  return listElem;
}

List *concatenateLists(List *listA, List *listB)
{
  List
    *iter = listA ; 

  if( NOT iter)
    return listB;

  FOR_LIST(iter)
  {
    if( NOT iter->next)
      {
	iter->next = listB;
	return listA;
      }    
  }

  return NULL;
}


void freeList(List* list)
{
  List 
    *iter;
  for(iter = list; iter;)
    {
      list = iter->next;
      free(iter->value);
      free(iter);
      iter = list;
    }
}


void freeListFlat(List* list)
{
  List 
    *iter;
  for(iter = list; iter;)
    {
      list = iter->next;
      free(iter);
      iter = list;
    }
}


boolean elemIsInIndexList(int index, IndexList *list)
{
  FOR_LIST(list)
  {
    if(index == list->index)
      return TRUE;
  }
  return FALSE;
}


IndexList *joinIndexListsNonRedundant(IndexList *listA, IndexList *listB)
{
  IndexList *result = NULL;

  FOR_LIST(listA)
    APPEND_INT(listA->index, result);
  freeIndexList(listA);

  FOR_LIST(listB)
    if( NOT elemIsInIndexList(listB->index, result))
      APPEND_INT(listB->index, result);

  freeIndexList(listB);

  return result;
}
