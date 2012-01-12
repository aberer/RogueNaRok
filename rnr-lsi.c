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


#include <math.h>
#include <time.h> 
#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <assert.h>
#include <unistd.h>

#include "common.h"
#include "sharedVariables.h"
#include "newFunctions.h"
#include "legacy.h"
#include "Tree.h"
#include "BitVector.h"


#define PROG_NAME "RnR-lsi"
#define PROG_VERSION "1.0"
#define PROG_RELEASE_DATE "2011-10-25"

#define FLIP_NTH_BIT(bitVector,n) (bitVector[(n) / MASK_LENGTH] |= mask32[ (n) % MASK_LENGTH ])
#define NTH_BIT_IS_SET(bitVector,n) (bitVector[(n) / MASK_LENGTH] & mask32[(n) % MASK_LENGTH])

#define NUMBER_TO_PRUNE(x) ((x) / 2)
#define SHORT_UNSIGNED_MAX 65535


typedef unsigned short int smallNumber;
extern void freeIndexList(IndexList *list);

typedef struct _scorePerTaxon
{
  int taxonId;

  double mrConsensus;
  double strictConsensus; 
  double bestTreeScore;

  double leafStab;
  double taxonInstability;
} ScorePerTaxon;

void printQuadruples(All *tr, smallNumber ***quadruples)
{
  int
    i,j,k;

  for(i = 0; i < tr->mxtips; ++i)
    {
      for(j = 0; j < tr->mxtips; ++j)
	{
	  for(k = 0; k < tr->mxtips; ++k)
	    {
	      if(quadruples[i][j][k])
		printf("%s|%s,%s\t%d\n", tr->nameList[i+1], tr->nameList[j+1], tr->nameList[k+1], quadruples[i][j][k]);
	    }
	}
    }
}


IndexList *extractQuadruplesRecursively(All *tr, nodeptr node, smallNumber ***quadruples, BitVector *neglectThose, boolean isStart) 
{
  IndexList 
    *leftList, *rightList, *iterA, *iterB, *iterC;

  if( NOT isTip(node->number, tr->mxtips) || isStart)
    {
      nodeptr 
	leftNode = (isStart) ? node->back->next->back : node->next->back, 
	rightNode = (isStart) ? node->back->next->next->back : node->next->next->back;

      leftList = extractQuadruplesRecursively(tr, leftNode,  quadruples, neglectThose, FALSE);
      rightList = extractQuadruplesRecursively(tr, rightNode, quadruples, neglectThose,FALSE);
    }
  else 
    {
      IndexList
	*elem = CALLOC(1,sizeof(IndexList)); 
      elem->index = node->number;
      elem->next = NULL;
      return elem;
    }

  /* insert the information */
  for(iterA = leftList; iterA; iterA = iterA->next)
    {
      for(iterB = rightList; iterB; iterB = iterB->next)
	{
	  for(iterC = iterB->next; iterC; iterC = iterC->next)
	    {
	      int 
		indexB = iterB->index,
		indexC = iterC->index;

	      if(indexB > indexC)
		SWAP(indexB,indexC);

	      USE_UPPER_TRIANGLE_LSI(quadruples, (iterA->index-1), (indexB-1), (indexC-1))++;
	    }
	}
    }

  for(iterA = rightList; iterA; iterA = iterA->next)
    {
      for(iterB = leftList; iterB; iterB = iterB->next)
	{
	  for(iterC = iterB->next; iterC; iterC = iterC->next)
	    {
	      int
		indexB = iterB->index,
		indexC = iterC->index;
	      
	      if(indexB > indexC)
		SWAP(indexB,indexC);
	      
	      USE_UPPER_TRIANGLE_LSI(quadruples, (iterA->index-1), (indexB-1), (indexC-1))++;
	    }
	}
    }

  for(iterA = leftList; iterA; iterA = iterA->next)
    if( NOT iterA->next)
      {
	iterA->next = rightList; 
	break;
      }

  return leftList;
}

void freeQuads(All *tr, smallNumber ***quads)
{
  int i,j;
  for(i = 0; i < tr->mxtips; ++i)
    {
      for(j = 0; j < tr->mxtips; ++j)
	free(quads[i][j]);
      free(quads[i]);
    }
  free(quads);
}


smallNumber ***initQuads(All *tr) 
{
  int i, j; 
  smallNumber ***result =  CALLOC(tr->mxtips, sizeof(smallNumber**)); 


  for(i = 0; i < tr->mxtips; ++i)
    {
      result[i] = CALLOC(tr->mxtips, sizeof(smallNumber*));
      for(j = 0 ; j < tr->mxtips; ++j)
	result[i][j] = CALLOC(tr->mxtips - j, sizeof(smallNumber));
    } 
  
  return result;
}


int intcmp(const void *a, const void *b)
{
    return *(int *)b - *(int *)a;
}


void calculateLeafStability(All *tr, char *bootstrapFileName, char *excludeFileName)
{  
  FILE 
    *outf = getOutputFileFromString("leafStabilityIndices"),
    *bootstrapFile =  getNumberOfTrees(tr,bootstrapFileName);
 
  BitVector
    *neglectThose = neglectThoseTaxa(tr, excludeFileName);

  int 
    i, j, k, l;

  if(tr->numberOfTrees >= SHORT_UNSIGNED_MAX)
    {
      PR("Sorry, %s is not  capable of handling more than %d trees. You may want to adjust the code, if you have sufficient memory at disposal.\n");
      exit(-1);
    }

  fprintf(outf,"taxon\tlsDif\tlsEnt\tlsMax\n");
  PR("taxon\tlsDif\tlsEnt\tlsMax\n");

  /* calculate leaf stability for n-th taxon */
  for(i = 0; i < tr->mxtips; i++) 
    { 
      if( NOT NTH_BIT_IS_SET(neglectThose,i))
	{
	  PR("%s\tNA\tNA\tNA\n", tr->nameList[i+1]);
	  fprintf(outf,"%s\tNA\tNA\tNA\n", tr->nameList[i+1]);
	  continue;
	}

      smallNumber
	***quadruples = initQuads(tr);

      double
	lsDif = 0.0, 
	lsEnt = 0.0,
	lsMax = 0.0;

      rewind(bootstrapFile);
      
      for(j = 1; j <= tr->numberOfTrees; ++j)
	{
	  int k; 
	  readBootstrapTree(tr,bootstrapFile);
	  FOR_0_LIMIT(k,tr->mxtips)
	    if( NOT NTH_BIT_IS_SET(neglectThose, k))
	      pruneTaxon(tr,k+1);
	  extractQuadruplesRecursively(tr, tr->nodep[i+1], quadruples, neglectThose, TRUE);
	}
 
      /* calculate leaf stability of this leaf */
      int
	sum = 0;
      double 
	cnt = 0.0;
      for(j = 0; j  < tr->mxtips; ++j)
	{
	  if(NOT NTH_BIT_IS_SET(neglectThose, j))
	    continue;
	  
	  for(k = j+1; k < tr->mxtips; ++k) 
	    {
	      if(NOT NTH_BIT_IS_SET(neglectThose, k))
		continue;

	      for(l = k+1; l < tr->mxtips; ++l) 
		{
		  if(NOT NTH_BIT_IS_SET(neglectThose, l))
		    continue;

		  int rels[3] = { 

		    GET_FROM_UPPER_TRIANGLE(quadruples,j,k,l),
		    GET_FROM_UPPER_TRIANGLE(quadruples,k,j,l),
		    GET_FROM_UPPER_TRIANGLE(quadruples,l,j,k)
		  };
		  
		  sum = rels[0] + rels[1] + rels[2];

		  qsort(rels, 3, sizeof(int), intcmp);

		  assert( sum == 0  || sum == tr->numberOfTrees);		  
		  assert(rels[0] >= 0 && rels[1] >= 0 && rels[0] >= rels[1]);
		  assert(   (NTH_BIT_IS_SET(neglectThose, j) &&  NTH_BIT_IS_SET(neglectThose, k) && NTH_BIT_IS_SET(neglectThose, l) && NTH_BIT_IS_SET(neglectThose, i)) );

		  if(sum)
		    {		  
		      lsDif += (double)(rels[0] - rels[1]) / (double)tr->numberOfTrees;
		      lsMax += (double) rels[0]    / (double)  tr->numberOfTrees ;
		      /* PR(">> %f\n", ((double) rels[0]    / (double)  tr->numberOfTrees )); */
		      
		      double tmp = (double)rels[0] / (double)tr->numberOfTrees;
		      lsEnt -= tmp * log(tmp);
		      if(rels[1] > 0)
			{
			  double tmp = (double) rels[1]  / (double) tr->numberOfTrees;
			  lsEnt -= tmp * log(tmp);
			}
		      if(rels[2] > 0)
			{
			  double tmp = (double) rels[2] / (double) tr->numberOfTrees;
			  lsEnt -= tmp * log(tmp);
			}
		      cnt++;
		    }
		}
	    }
	}
      lsDif /= cnt;
      
      /* normalize between 0 and 1 */
      lsEnt /= cnt;
      lsEnt /= 3. * log(1./3.)  * 1./3.;
      lsEnt = lsEnt + 1;
      
      /* also normalize */
      lsMax /= cnt  ;
      lsMax = (lsMax - (1. / 3.)) * 3. / 2. ;
      
      PR("%s\t%f\t%f\t%f\n", tr->nameList[i+1], lsDif, lsEnt, lsMax);
      fprintf(outf, "%s\t%f\t%f\t%f\n", tr->nameList[i+1], lsDif, lsEnt, lsMax);
      freeQuads(tr,quadruples);
    }
}


void printHelpFile()
{
  printVersionInfo(FALSE);
  printf("This program computes three flavors of leaf stability index for each taxon.\n\nSYNTAX: ./%s -i <bootTrees> -n <runId> [-w <workingDir>] [-h]\n", lowerTheString(programName));
  printf("\n\tOBLIGATORY:\n");
  printf("-i <bootTrees>\n\tA collection of bootstrap trees.\n");
  printf("-n <runId>\n\tAn identifier for this run.\n");
  printf("\n\tOPTIONAL:\n");
  printf("-x <excludeFile>\n\tPrune the taxa in the file first (one taxon per line), before computing the lsi.\n");
  printf("-w <workDir>\n\tA working directory where output files are created.\n");
  printf("-h\n\tThis help file.\n");
}


int main(int argc, char *argv[])
{
  programName = PROG_NAME; 
  programVersion = PROG_VERSION;
  programReleaseDate = PROG_RELEASE_DATE; 

  int
    c;

  char
    *excludeFile = "",
    *bootTrees = "";

   while ((c = getopt (argc, argv, "hi:n:w:m:x:")) != -1)
    {
      switch(c)
	{
	case 'i':
	  bootTrees = optarg;
	  break;
	case 'n':
	  strcpy(run_id, optarg);
	  break; 
	case 'w':
	  /* printf("optarg : %s\n", optarg);  */
	  strcpy(workdir, optarg);
	  break;
	case 'x':
	  excludeFile = optarg;
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
      printHelpFile(FALSE);
      exit(-1);
    }  

  if( NOT strcmp(run_id, ""))
    {
      printf("Please specify a run-id via -n\n");
      printHelpFile(FALSE);
      exit(-1);
    }

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
 
  calculateLeafStability(tr, bootTrees, excludeFile); 

  return 0;
}

