#include "newFunctions.h"

IndexList *parseToDrop(All *tr, FILE *toDrop) 
{
  IndexList 
    *result = NULL;

  char
    line[1024];
  
  while(fgets(line, 1024, toDrop) != NULL)
    {
      char bla[1024];
      sscanf(line, "%s\n", bla);
      PR("will drop %d\n", treeFindTipByLabelString(bla, tr));
      result = appendToIndexList(treeFindTipByLabelString(bla, tr), result);
    }  

  return result;
}

void pruneTaxon(All *tr, unsigned int k) 
{
  assert(k > 0 && k <= ((unsigned int)(tr->mxtips)));

  nodeptr 
    p = tr->nodep[k],    
    q = p->back,
    q1 = q->next->back,
    q2 = q->next->next->back;

  hookupDefault(q1, q2, tr->numBranches);
  
  tr->start = findAnyTip(q1, tr->mxtips);
  
  assert(p != tr->start && q != tr->start);
}


BitVector *neglectThoseTaxa(All *tr, char *toDrop)  
{
  int
    i = 0; 
  
  BitVector
    *result = CALLOC(tr->bitVectorLength, sizeof(BitVector));  

  for(i = 0; i < tr->mxtips; i++)
    FLIP_NTH_BIT(result,i);
  
  if(strlen(toDrop) == 0)
    return result;
  
  FILE
    *toDropFile = myfopen(toDrop,"r");
  
  IndexList 
    *iter, 
    *toConvert = parseToDrop(tr, toDropFile);

  for(iter = toConvert; iter; iter = iter->next) 
    {
      UNFLIP_NTH_BIT(result, (iter->index -1) );
      assert( NOT NTH_BIT_IS_SET(result, iter->index -1 ));
    }

  freeIndexList(toConvert);

  return result;
}


Array *getOriginalBipArray(All *tr, FILE *bestTree, FILE *treeFile) 
{
  Array *result = CALLOC(1, sizeof(Array));

  int 
    i,j,bCount = 0,
    treeVectorLength = GET_BITVECTOR_LENGTH((tr->numberOfTrees+1));
  unsigned int 
    vectorLength = 0,
    **setBitVectors = initBitVector(tr, &vectorLength);
  hashtable
    *setHtable =  initHashTable(tr->mxtips * FC_INIT * 10);
  nodeptr 
    commonStart = NULL;  

  BitVector 
    lastByte = 0;  

  for(i = tr->mxtips; i < MASK_LENGTH * vectorLength; ++i)
    lastByte |= mask32[i % MASK_LENGTH];

  unsigned int 
    *randForTaxa = CALLOC(tr->mxtips, sizeof(unsigned int));
  
  for(i = 0; i < tr->mxtips; ++i)  
    randForTaxa[i] = rand();

  rewind(treeFile);
  
  if(bestTree)
    rewind(bestTree);

  /* get bipartitions of bootstrap set */
  for( i = 1; i <= tr->numberOfTrees; ++i)
    {      
      readBootstrapTree(tr, treeFile);     
      
      if( NOT commonStart)
	commonStart = tr->start;
      bCount = 0; 
      bitVectorInitravSpecial(setBitVectors, commonStart->back, tr->mxtips, vectorLength, setHtable, (i - 1), BIPARTITIONS_BOOTSTOP, (branchInfo *)NULL, &bCount, treeVectorLength, FALSE, FALSE);
    }
  
  if(bestTree)
    {
      readBestTree(tr,bestTree);      
      
      bCount = 0;
      bitVectorInitravSpecial(setBitVectors, commonStart->back, tr->mxtips, vectorLength, setHtable, tr->numberOfTrees, BIPARTITIONS_BOOTSTOP, (branchInfo *)NULL, &bCount, treeVectorLength, FALSE, FALSE);
      assert(bCount == tr->mxtips - 3);
    }
  
  result->commonAttributes = CALLOC(1,sizeof(ProfileElemAttr));
  ((ProfileElemAttr*)result->commonAttributes)->bitVectorLength = vectorLength;
  ((ProfileElemAttr*)result->commonAttributes)->treeVectorLength = treeVectorLength;
  ((ProfileElemAttr*)result->commonAttributes)->lastByte = lastByte;
  ((ProfileElemAttr*)result->commonAttributes)->randForTaxa = randForTaxa;
  result->length = setHtable->entryCount;
  result->arrayTable = CALLOC(result->length, sizeof(ProfileElem*));
  
  j = 0; 
  for(i = 0; i < setHtable->tableSize; ++i)
    {
      entry
	*elem = setHtable->table[i];
      
      while(elem)
	{
	  ((ProfileElem**)result->arrayTable)[j] = addProfileElem(elem, vectorLength, treeVectorLength, tr->numberOfTrees);
	  j++;
	  elem = elem->next;
	}      
    }
  assert(j == setHtable->entryCount);
  
  int cnt= 0; 
  for(i = 0; i < result->length; ++i)
    if(((ProfileElem**)result->arrayTable)[i]->isInMLTree)
      cnt++;
  if(bestTree)
    assert(cnt == tr->mxtips - 3);

  freeHashTable(setHtable);
  freeBitVectors(setBitVectors, 2 * tr->mxtips);
  free(setBitVectors);
  free(setHtable);
  free(randForTaxa);

  /* TEST */
  for(i = 0; i < result->length; ++i )
    ((ProfileElem**)result->arrayTable)[i]->id  = i; 
  /* END */

  return result;
}
