#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include <unistd.h>

#include "common.h"
#include "List.h"
#include "Tree.h"
#include "BitVector.h"
#include "sharedVariables.h"
#include "legacy.h"


#define PROG_NAME "RnR-mast"
#define PROG_VERSION "0.1"
#define PROG_RELEASE_DATE "not yet"

/* 
   NOTICE: as the numbering of the taxa is 1-based, throughout this
   file all array indices (that use the taxon number as access code)
   are also 1-based. This means that we waste a bit of memory, but
   hopefully avoids index wars...
*/


/* SWITCHES */
/* #define PRINT_VERY_VERBOSE */
#define FAST_BV_COMPARISON

boolean areSameBitVectors(BitVector *a, BitVector *b, int bitVectorLength);

typedef struct _BitVectorList
{
  struct _BitVectorList *next;
  BitVector *taxaInMast;
} MastList; 


IndexList *traverseForTriples(All *tr, nodeptr node, boolean isStart, BitVector ***rootedTriples, BitVector **bitVectors, int bitVectorLength) 
{ 
  int 
    i;
  IndexList 
    *leftList, *rightList, *iterA, *iterB;
  BitVector
    *leftBv, *rightBv; 
  
  if( NOT isTip(node->number, tr->mxtips) || isStart)
    {
      nodeptr 
	leftNode = (isStart) ? node->back->next->back : node->next->back, 
	rightNode = (isStart) ? node->back->next->next->back : node->next->next->back;
      leftBv = bitVectors[leftNode->number];
      rightBv = bitVectors[rightNode->number];
      leftList = traverseForTriples(tr, leftNode, FALSE, rootedTriples, bitVectors, bitVectorLength);
      rightList = traverseForTriples(tr, rightNode, FALSE, rootedTriples, bitVectors, bitVectorLength);
    }
  else 
    {
      IndexList
	*elem = malloc(sizeof(IndexList)); 
      elem->index = node->number;
      elem->next = NULL;
      FLIP_NTH_BIT(bitVectors[node->number], node->number);
      assert(NTH_BIT_IS_SET(bitVectors[node->number], node->number));
      return elem;
    }

  /* enter entries in the triple structure */
  for(iterA = leftList; iterA; iterA = iterA->next)
    {
      int
	indexA = iterA->index;
      
      assert(NTH_BIT_IS_SET(leftBv, indexA));
      for(iterB = rightList; iterB; iterB = iterB->next)
	{
	  int 
	    indexB = iterB->index;

	  assert(NTH_BIT_IS_SET(rightBv, indexB));
	  for(i = 0; i < bitVectorLength; ++i)
	    {
	      rootedTriples[indexA][indexB][i] |= leftBv[i];
	      rootedTriples[indexB][indexA][i] |= rightBv[i];
	    }

	  UNFLIP_NTH_BIT(rootedTriples[indexA][indexB], indexA);
	  UNFLIP_NTH_BIT(rootedTriples[indexB][indexA], indexB);
	}
    }

  /* prepare bitvector and list for this root */
  for(i = 0; i < bitVectorLength; ++i)
    bitVectors[node->number][i] = (leftBv[i] | rightBv[i]);

  for(iterA = leftList; iterA; iterA = iterA->next)
    if( NOT iterA->next)
      {
	iterA->next = rightList;
	
	for(iterA = leftList; iterA; iterA = iterA->next)
	  assert(NTH_BIT_IS_SET(bitVectors[node->number], iterA->index));
	break;
      }
  
  /* insert the triples containing the root -- TODO we could omit
     that, as the root always must be part of the MAST */
  i = 0;
  if(isStart)
    {
      for(iterA = leftList; iterA; iterA = iterA->next)
  	{
  	  for(iterB = leftList; iterB; iterB = iterB->next)
  	    if(iterA->index != iterB->index)
  	      FLIP_NTH_BIT(rootedTriples[iterA->index][node->number], iterB->index);
    	  i++;
    	}
      /* assert(i == tr->mxtips-1); */
    }

  return leftList; 
}


BitVector*** initializeTriplesStructure(All *tr, int bitVectorLength)
{
  int 
    i, j;
  BitVector
    ***result = malloc((tr->mxtips + 1) * sizeof(BitVector**));

  for(i = 0; i < tr->mxtips+1; ++i)
    {
      result[i] = malloc((tr->mxtips+1) * sizeof(BitVector*));
      for(j = 0; j < tr->mxtips+1; ++j)	
	result[i][j] = calloc(bitVectorLength, sizeof(BitVector));
    }

  return result; 
}

void freeTriplesStructure(All *tr, BitVector ***triples)  
{
  int 
    i,j;
  
  for(i = 0; i < tr->mxtips+1; ++i)
    {
      for(j = 0; j < tr->mxtips+1;++j)
	free(triples[i][j]);
      free(triples[i]);
    }
  free(triples);
}


BitVector*** getIntersectionOfRootedTriples(FILE *bootstrapFile, All *tr, int startingNodeIndex, BitVector *neglectThose)
{
  int 
    i,j,k,treeNum=0, 
    bitVectorLength = GET_BITVECTOR_LENGTH((tr->mxtips+1));
  IndexList *iter; 
  BitVector 			
    **bitVectors = malloc((2 * tr->mxtips - 1) * sizeof(BitVector*)),
    ***result = initializeTriplesStructure(tr, bitVectorLength), 
    ***theseTriples = initializeTriplesStructure(tr, bitVectorLength);
  for(i = 0; i < 2 * tr->mxtips - 1; ++i)
    bitVectors[i] = malloc(bitVectorLength * sizeof(BitVector)); 

  rewind(bootstrapFile);
  
  while(treeNum++ < tr->numberOfTrees)
    {
      readBootstrapTree(tr,bootstrapFile);
      
      FOR_0_LIMIT(i, tr->mxtips)
	if( NOT NTH_BIT_IS_SET(neglectThose, i))
	  pruneTaxon(tr, i+1);
      
      /* reset bit-vectors */
      FOR_0_LIMIT(i, (tr->mxtips+1))
      	memset(bitVectors[i], 0, sizeof(BitVector) * bitVectorLength);
      
      /* reset triples */
      for(i = 0; i < tr->mxtips+1; ++i)
	for(j = 0; j < tr->mxtips+1; ++j)
	  memset(theseTriples[i][j], 0, sizeof(BitVector) * bitVectorLength);

      /* get all triples in the current tree and clean up */
      IndexList *il = traverseForTriples(tr, tr->nodep[startingNodeIndex], TRUE, theseTriples, bitVectors, bitVectorLength);
      for(iter = il; iter; )
	{
	  il = iter->next;
	  free(iter);
	  iter = il; 
	}      
      
      /* intersect the triples */
      for(i = 0; i < tr->mxtips+1; ++i)
	for(j = 0; j < tr->mxtips+1; ++j)
	  for(k = 0; k < bitVectorLength; ++k)
	    if(treeNum == 1)
	      result[i][j][k] |= theseTriples[i][j][k];
	    else
	      result[i][j][k] &= theseTriples[i][j][k];
    }

  freeTriplesStructure(tr, theseTriples);
  
  for(i = 0; i < 2 * tr->mxtips - 1; ++i)
    free(bitVectors[i]);
  free(bitVectors);

  return result; 
}


typedef struct _agreementMatrix
{
  int score;
  IndexList* valuesA;
  IndexList* valuesB;
} AgreementMatrix;


/* we mean ax|b by that */
IndexList* findBestXinAxB(All *tr, int A, int B, AgreementMatrix **agreementMatrix, BitVector ***triples, int *resultScore) 
{
  int 
    i = 0,
    maximum = 0;
  BitVector 
    *currentBv = triples[A][B];
  IndexList
    *result = NULL;
  
  while(i < tr->mxtips)
    {    
#ifdef FAST_BV_COMPARISON
      if( ((i & 31) == 0) && NOT currentBv[i]) /* NOTE portability issue!!! */
      	{
	  /* assert((i % 32) == i & (MASK_LENGTH - 1)); */
      	  i += MASK_LENGTH;
      	  continue;
      	}
#endif

      if(NTH_BIT_IS_SET(currentBv, i))
	{
	  int score = ( A < i ) ? agreementMatrix[A][i].score : agreementMatrix[i][A].score;
	  assert(score);
	  if(score > maximum)
	    maximum = score;
	}
      ++i;
    }
  
  i = 0;
  while(i < tr->mxtips)
    {
#ifdef FAST_BV_COMPARISON
      if( ((i & 31) == 0) && NOT currentBv[i]) /* NOTE portability issue!!! */
      	{
	  /* assert((i % 32) == i & (MASK_LENGTH - 1)); */
      	  i += MASK_LENGTH;
      	  continue;
      	}
#endif

      if(NTH_BIT_IS_SET(currentBv,i))
	{
	  int score = ( A < i ) ? agreementMatrix[A][i].score : agreementMatrix[i][A].score;
	  assert(score); 
	  if(score == maximum)
	    result = appendToIndexList(i,result);
	}
      i++;
    }

  *resultScore = maximum ? maximum : 1;

  return result;
}


IndexList *traverseForMastTable(All *tr, nodeptr node, BitVector ***triples, AgreementMatrix **agreementMatrix, boolean isStart)
{
  if(isTip(node->number, tr->mxtips) && NOT isStart)
    {
      IndexList *elem = malloc(sizeof(IndexList));
      elem->next = NULL;
      elem->index = node->number;
      return elem;
    }
  else
    {
      nodeptr 
	leftNode = (isStart) ? node->back->next->back : node->next->back, 
	rightNode = (isStart) ? node->back->next->next->back : node->next->next->back;
      IndexList
	*nodesOnLeft = traverseForMastTable(tr, leftNode, triples, agreementMatrix, FALSE),
	*nodesOnRight = traverseForMastTable(tr, rightNode, triples, agreementMatrix, FALSE), 
	*iterA, *iterB;
      
      for(iterA = nodesOnLeft; iterA; iterA = iterA->next)
	{
	  int indexA = iterA->index;
	  for(iterB = nodesOnRight; iterB; iterB = iterB->next)
	    {
	      int indexB = iterB->index;
	      
	      int
		scoreA = 0, scoreB = 0; 
	      
	      IndexList
		*xs = findBestXinAxB(tr, indexA, indexB, agreementMatrix, triples, &scoreA),
		*ys = findBestXinAxB(tr, indexB, indexA, agreementMatrix, triples, &scoreB);
#ifdef PRINT_VERY_VERBOSE
	      printf("%s,x|%s => x = %s\t%s,y|%s => y = %s \n", tr->nameList[indexA], tr->nameList[indexB],
		     (x) ? tr->nameList[x] : "0", tr->nameList[indexB], tr->nameList[indexA], (y) ? tr->nameList[y] : "0");
#endif
	      
	      AgreementMatrix
		*currentElem = (indexA > indexB) ? &(agreementMatrix[indexB][indexA]) : &(agreementMatrix[indexA][indexB]);
	      currentElem->score = scoreA + scoreB;
#ifdef PRINT_VERY_VERBOSE
	      printf("MAST(%s,%s) = %d\n", tr->nameList[indexA], tr->nameList[indexB], currentElem->score);
#endif
	      currentElem->valuesA = (indexA > indexB) ? ys : xs;
	      currentElem->valuesB = (indexA > indexB) ? xs : ys;
	    } 
	}
      
      /* join the lists */
      for(iterA = nodesOnLeft; iterA; iterA = iterA->next)
	if( NOT iterA->next)
	  {
	    iterA->next = nodesOnRight;
	    break;
	  }
      
      if(isStart)
	{
	  int indexA = node->number; 
	  for(iterA = nodesOnLeft; iterA; iterA = iterA->next)
	    {
	      int indexB = iterA->index;
	      int
		scoreA = 0,  scoreB = 0; 
	      
	      IndexList *xs = findBestXinAxB(tr, indexA, indexB, agreementMatrix, triples, &scoreA);
	      IndexList *ys = findBestXinAxB(tr, indexB, indexA, agreementMatrix, triples, &scoreB);
	      
	      AgreementMatrix *currentElem = 
		(indexA > indexB) ? &(agreementMatrix[indexB][indexA]) : &(agreementMatrix[indexA][indexB]);
	      currentElem->score = scoreA + scoreB;
	      currentElem->valuesA = (indexA > indexB) ? ys : xs; 
	      currentElem->valuesB = (indexA > indexB) ? xs : ys; 
	    }
	}
      
      return nodesOnLeft;
    }
}


void printAgreementTable(All *tr, AgreementMatrix **matrix)
{
  int 
    i,j;

  printf("\n\n  ");
  
  for(i = 1; i <= tr->mxtips; ++i)
    {
      printf("%s,", tr->nameList[i]);
    }
  printf("\n");

  for(i = 1; i <= tr->mxtips; ++i)
    {
      printf("%s,", tr->nameList[i]);
      for(j = 1; j <= tr->mxtips; ++j)
	{
	  if( NOT matrix[i][j].score)
	    printf(" ,");
	  else
	    printf("%d,", matrix[i][j].score); 
	}
      printf("\n");
    }
  printf("\n");
}


List *filterDuplicates(List *masts, int bitVectorLength)
{
  List
    *iter, 
    *result = NULL;
  int
    i = 0,j,
    numElem = 0; 
  boolean *isDuplicate; 
  BitVector **bitVectors;
  
  for(iter = masts; iter; iter = iter->next)
    numElem++;
  
  bitVectors = calloc(numElem, sizeof(BitVector*));  
  i = 0; 
  for(iter = masts; iter;)
    {
      bitVectors[i] = (BitVector*)iter->value;
      i++;
      masts = iter->next;
      free(iter);
      iter = masts;
    }
  
  isDuplicate = calloc(numElem, sizeof(boolean));
  for(i = 0; i < numElem; i++)
    {
      if(isDuplicate[i])
	continue; 
      
      for(j = i+1; j < numElem; ++j)
	{
	  if(isDuplicate[j])
	    continue;

	  isDuplicate[j] =  areSameBitVectors(bitVectors[i], bitVectors[j], bitVectorLength);
	}
    }
  
  for(i = 0; i < numElem; ++i)
    if(isDuplicate[i])
      free(bitVectors[i]);
    else
      result = appendToList(bitVectors[i], result);
  

  free(isDuplicate);
  free(bitVectors);

  return result; 
}


List *backtrace(All *tr, AgreementMatrix **matrix, int cntI, int cntJ, int bitVectorLength, boolean allMasts) 
{
  List
    *iIter;
  IndexList 
    *iter;
  
  int 
    i;
  
  List
    *result = NULL, *fromLeft = NULL, *fromRight = NULL; 
  
  if(cntI > cntJ)
    SWAP(cntI,cntJ);
  
  AgreementMatrix
    *matrixElem = &(matrix[cntI][cntJ]);

  /* int a = 0, b = 0;  */

  for(iter = matrixElem->valuesA; iter; iter = iter->next)    
    {
      List *tmp = backtrace(tr, matrix, cntI,iter->index, bitVectorLength, allMasts);
      fromLeft = joinLists(fromLeft, tmp);
      if( NOT allMasts)
	break;
      
      /* a++; */
    }

  for(iter = matrixElem->valuesB; iter; iter = iter->next)
    {
      List *tmp = backtrace(tr, matrix, cntJ,iter->index, bitVectorLength, allMasts); 
      fromRight = joinLists(fromRight, tmp);
      if( NOT allMasts)
	break;
      /* b++; */
    }

  /* printf("%d\t%d\n", a,b); */

  if( NOT (fromLeft  || fromRight))
    result = appendToList(calloc(bitVectorLength, sizeof(BitVector)),result);
  else if( NOT (fromLeft && fromRight))
    result = fromLeft ? fromLeft : fromRight;
  else
    {
      List
	*iterA, *iterB; 
      for(iterA = fromLeft; iterA; iterA = iterA->next)
	{
	  for(iterB = fromRight; iterB; iterB = iterB->next)
	    {
	      BitVector
		*new = calloc(bitVectorLength, sizeof(BitVector));
	      for(i = 0; i < bitVectorLength; ++i)
		new[i] = ((BitVector*)iterA->value)[i] | ((BitVector*)iterB->value)[i];
	      result = appendToList(new, result); 
	    }
	}
      freeList(fromRight);
      freeList(fromLeft);

      if(allMasts)
	result = filterDuplicates(result, bitVectorLength);
    }

  

  /* set the bits of the current entry */
  for(iIter = result; iIter; iIter = iIter->next)
    {
      FLIP_NTH_BIT(((BitVector*)iIter->value), cntI);
      FLIP_NTH_BIT(((BitVector*)iIter->value), cntJ);
    }
  
  return result;
}


List *backtraceMasts(All *tr, AgreementMatrix **matrix, boolean allMasts)
{
  int
    bitVectorLength = GET_BITVECTOR_LENGTH(tr->mxtips),
    maxScore = 0,
    i, j;
  
  List 
    *result = NULL;

  /* find maximum in agreement matrix */
  for(i = 1; i <= tr->mxtips; ++i)
    for(j = 1; j <= tr->mxtips; ++j)
      if(maxScore < matrix[i][j].score)
	maxScore  = matrix[i][j].score;
  assert(maxScore);

  for(i = 1; i <= tr->mxtips; ++i)
    for(j = 1; j <= tr->mxtips; ++j)
      if(matrix[i][j].score == maxScore)
	{
	  List
	    *elems = backtrace(tr, matrix, i,j, bitVectorLength, allMasts);
	  result = joinLists(result, elems);
	}

  return result; 
}


AgreementMatrix** computeAgreementTable(All *tr, BitVector ***triples, int startingNodeIndex)
{
  int 
    i; 
  IndexList *iter;
  AgreementMatrix
    **agreementMatrix = malloc((tr->mxtips+1) * sizeof(AgreementMatrix*));
  for(i = 1; i <= tr->mxtips; ++i)
    agreementMatrix[i] = calloc(tr->mxtips+1, sizeof(AgreementMatrix));

  IndexList
    *list = traverseForMastTable(tr, tr->nodep[startingNodeIndex], triples, agreementMatrix, TRUE);

  for(iter = list; iter;)
    {
      list = iter->next;
      free(iter);
      iter = list;
    }
  
#ifdef PRINT_VERY_VERBOSE
  printAgreementTable(tr, agreementMatrix);
#endif

  return agreementMatrix;
}


void printBipartition(All *tr, BitVector *bv, int contentLength)
{
  int i; 
  for(i = 0; i < contentLength;++i)
    if(NTH_BIT_IS_SET(bv, i))
      printf("%s,", tr->nameList[i+1]);
}


void verifyMasts(All *tr, FILE *bootstrapFile, BitVector *taxaToKeep)
{
  int
    bCount = 0,
    treeVectorLength = GET_BITVECTOR_LENGTH(tr->numberOfTrees),
    i, j;
  unsigned int vectorLength = 0;
  hashtable
    *htable =  initHashTable(tr->mxtips * FC_INIT * 10);
  BitVector
    **bitVectors = initBitVector(tr, &vectorLength);
  entry *e;
  nodeptr commonStart = NULL;
  
  rewind(bootstrapFile);  

  for(i = 1; i <= tr->numberOfTrees; ++i)
    {
      /* treeReadLen(bootstrapFile, tr, FALSE, FALSE, TRUE, adef, TRUE); */
      readBootstrapTree(tr, bootstrapFile);
      
      /* prune taxa */
      for(j = 1; j <= tr->mxtips; ++j)
	if( NOT (NTH_BIT_IS_SET(taxaToKeep, j)) )
	  pruneTaxon(tr, j+1);
      
      if(i == 1 )
	commonStart = tr->start;
      bitVectorInitravSpecial(bitVectors, commonStart->back, tr->mxtips, vectorLength, htable, (i - 1), BIPARTITIONS_BOOTSTOP, (branchInfo *)NULL, &bCount, treeVectorLength, FALSE, FALSE);
    }

  for(i = 0; i < htable->tableSize; ++i)
    for(e = htable->table[i]; e; e = e->next)
      {
	if(genericBitCount(e->treeVector, treeVectorLength) != tr->numberOfTrees)
	  {
	    printf("error %d/%d\nbv: ", genericBitCount(e->treeVector, treeVectorLength), tr->numberOfTrees);
	    printBipartition(tr,e->bitVector, tr->mxtips);
	    printf("\nrepresented in: ");
	    printBitVector(e->treeVector, tr->numberOfTrees);
	    printf("\n");
	  }
	assert(genericBitCount(e->treeVector, treeVectorLength) == tr->numberOfTrees);
      }

  freeBitVectors(bitVectors, 2 * tr->mxtips);
  freeHashTable(htable);
}

int countTips(nodeptr p, int numsp)
{
  if(isTip(p->number, numsp))  
    return 1;    
  {
    nodeptr q;
    int tips = 0;

    q = p->next;
    while(q != p)
      { 
	tips += countTips(q->back, numsp);
	q = q->next;
      } 
    
    return tips;
  }
}


void printMastToFile(All *tr, FILE *bootstrapFile,  BitVector *mast, FILE *result)
{
  int
    droppedTaxaNum = 0, 
    tips, i; 

  rewind(bootstrapFile);  
  readBootstrapTree(tr,bootstrapFile);
    
  for(i = 1; i <= tr->mxtips; ++i)
    if( NOT NTH_BIT_IS_SET(mast,i) )
      {
	pruneTaxon(tr,i+1);
	droppedTaxaNum++;
      }
  
  tips = countTips(tr->start, tr->mxtips) + countTips(tr->start->back, tr->mxtips);
  assert((unsigned)tips == ((unsigned)tr->mxtips - droppedTaxaNum));
  char *tmp = writeTreeToString(tr, FALSE);
  fprintf(result, "%s", tmp);
  printBothOpen("MAST tree: %s", tmp);
}


void printMastsToFile(All *tr, FILE *bootstrapFile, List *masts)
{
  FILE
      *result = getOutputFileFromString("MaximumAgreementSubtree");

  List
    *iter;

  for(iter = masts; iter ; iter = iter->next)
    printMastToFile( tr, bootstrapFile, iter->value, result);

  fclose(result);
}


void freeAgreementMatrix(All *tr, AgreementMatrix **agm)
{
  int
    i,j;

  for(i = 1; i <= tr->mxtips; ++i)
    {
      for(j = 1; j <= tr->mxtips; ++j)
	{
	  freeIndexList(agm[i][j].valuesA);
	  freeIndexList(agm[i][j].valuesB);
	}
      free(agm[i]);
    }
  free(agm);
}


void freeAmatList(All *tr,List *aml)
{
  List
    *iter;

  for(iter = aml; iter;)
    {
      aml = iter->next;
      freeAgreementMatrix(tr, iter->value);
      free(iter);
      iter = aml;
    }
}


void calculateMast(char *bootStrapFileName, All *tr, char *excludeFileName, boolean allMasts) 
{
  int 
    bitVectorLength = GET_BITVECTOR_LENGTH((tr->mxtips+1)),
    mast = 0,
    i,j,k;

  FILE 
    *bootstrapFile = getNumberOfTrees(tr, bootStrapFileName);

  BitVector
    *taxaToNeglect = neglectThoseTaxa(tr, excludeFileName);

  List
    *iter,
    *amastList = NULL;

  for(i = 1; i <= tr->mxtips; ++i)
    {
      int
	currentMast = 0; 

      printBothOpen("rooting %d/%d\t%s\n", i, tr->mxtips, tr->nameList[i]);
      
      BitVector
	***commonRootedTriples = getIntersectionOfRootedTriples(bootstrapFile, tr, i, taxaToNeglect);
      AgreementMatrix
	**amat = computeAgreementTable(tr, commonRootedTriples, i);
            
      for(j = 1; j <= tr->mxtips; ++j)
	for(k = 1 ; k <= tr->mxtips; ++k)
	  if(currentMast < amat[j][k].score)
	    currentMast = amat[j][k].score;
     
      if(currentMast > mast)
	{
	  mast = currentMast;
	  freeAmatList(tr,amastList);
	  amastList = NULL;
	  amastList = appendToList(amat, amastList);
	}
      else if( allMasts && currentMast == mast )
	amastList = appendToList(amat, amastList);
      else
	freeAgreementMatrix(tr, amat);
      
      freeTriplesStructure(tr, commonRootedTriples);
    }
  
  List
    *accMasts = NULL;

  for(iter = amastList; iter; iter = iter->next)
    {
      List
	*Mast = backtraceMasts(tr, iter->value, allMasts);
      accMasts = joinLists(Mast, accMasts);
    }

  accMasts = filterDuplicates(accMasts, bitVectorLength);

  /* this works, commenting it out, as it only costs additional time */
  for(iter = accMasts; iter; iter = iter->next)
    verifyMasts(tr, bootstrapFile, ((BitVector*)iter->value));

  int cnt = 0; 
  for(iter = accMasts; iter; iter = iter->next)
    cnt++;
  printBothOpen("number of alternative MASTs (this is not the number of all possible MASTs, if you did not use \"ALL_MAST\"): %d\n", cnt);
  cnt = genericBitCount((BitVector*)accMasts->value, GET_BITVECTOR_LENGTH(tr->mxtips));
  printBothOpen("MAST size is: %d\n", cnt);

  /* print */
  printMastsToFile(tr,bootstrapFile, accMasts);

  freeList(accMasts);
  freeAmatList(tr, amastList);
  fclose(bootstrapFile);
}


void printRootedTriples(BitVector ***rootedTriples, All *tr)
{
  int 
    i, j, k;
    
  for(i = 1; i <= tr->mxtips; ++i)
    for(j = 1; j <= tr->mxtips;++j)
      for(k = 1; k <= tr->mxtips; ++k)
	if(NTH_BIT_IS_SET(rootedTriples[i][j], k))
	  printf("%s,%s|%s\n", tr->nameList[i], tr->nameList[k], tr->nameList[j]);
  
}


void printHelpFile()
{
  printVersionInfo(FALSE);
  printf("This program computes maximum agreement trees for unrooted input sets.\n\nSYNTAX: ./%s -i <bootTrees> -n <runId> [-w <workingDir>] [-h] [-a] [-x <excludeFile>]\n", lowerTheString(programName));
  printf("\nOBLIGATORY:\n");
  printf("-i <bootTrees>\n\tA collection of bootstrap trees.\n");
  printf("-n <runId>\n\tAn identifier for this run.\n");
  printf("\nOPTIONAL\n");
  printf("-a\n\tCompute all possible MAST trees. Without this flag, you will\n\t\
 only get a few MASTs that are easy to compute. As there may be an\n\t\
 exponential number of MASTs, use this option with care.\n");
  printf("-w <workDir>\n\tA working directory where output files are created.\n");
  printf("-x <excludeFile>\n\tExclude the taxa in this file (one taxon per line)\n\t\
 prior to computing the MAST. If you compute all MASTs anyway, this\n\t\
 option is option will not be useful. However, you can use this option\n\t\
 to speed up things.\n");
  printf("-h\n\tThis help file.\n");
}


int main(int argc, char *argv[])
{
  programName = PROG_NAME; 
  programVersion = PROG_VERSION;
  programReleaseDate = PROG_RELEASE_DATE; 

  int
    c;

  boolean
    computeAllMasts = FALSE;

  char
    *excludeFile = "",
    *bootTrees = "";

   while ((c = getopt (argc, argv, "hi:n:aw:x:")) != -1)
    {
      switch(c)
	{
	case 'i':
	  bootTrees = optarg;
	  break;
	case 'n':
	  strcpy(run_id, optarg);
	  break; 
	case 'a': 
	  computeAllMasts = TRUE;
	  break; 
	case 'x':
	  excludeFile = optarg;
	  break;
	case 'w':	  
	  strcpy(workdir, optarg);
	  break;
	case 'h':
	default:	
	  {
	    printHelpFile();
	    abort ();
	  }
	}
    }


  if( NOT strcmp(bootTrees, ""))
    {
      printf("Please specify a file containing bootstrap trees via -i.\n");
      printHelpFile();
      exit(-1);
    }  

  if( NOT strcmp(run_id, ""))
    {
      printf("Please specify a run-id via -n\n");
      printHelpFile();
      exit(-1);
    }


  /* maybe? */
  compute_bits_in_16bits();
  initializeMask();


  All 
    *tr = CALLOC(1,sizeof(All));  
  setupInfoFile();
  if  (NOT setupTree(tr, bootTrees))
    {
      PR("Something went wrong during tree initialisation. Sorry.\n");
      exit(-1);
    }   

  tr->tree_string = CALLOC(getTreeStringLength(bootTrees), sizeof(char));
  calculateMast(bootTrees, tr, excludeFile, computeAllMasts);

  return 0;
}
