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


#ifndef WIN32
#include <unistd.h>
#endif

#include <math.h>
#include <time.h> 
#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>

#include "common.h"
#include "BitVector.h"
#include "Tree.h"
#include "List.h"
#include "sharedVariables.h"

#define PROG_NAME  "RnR-tii"
#define PROG_VERSION "1.0"
#define PROG_RELEASE_DATE "2011-10-25"


typedef unsigned short int nodeDistance_t;

double tiiZ = 2.; 

typedef struct _distanceElem 
{
  nodeptr node;
  nodeDistance_t distance; 
} DistanceElem;


List *gatherDistances(All *tr, nodeptr node, nodeDistance_t **resultBranchDistances, boolean isStart) 
{
  List
    *distancesFromLeft,
    *distancesFromRight,
    *iterA, *iterB,
    *end = NULL;

  if( NOT isTip(node->number, tr->mxtips))
    {
      /* go down to the tips (DFS) */
      distancesFromLeft = gatherDistances(tr, node->next->back, resultBranchDistances, FALSE);
      distancesFromRight = gatherDistances(tr, node->next->next->back, resultBranchDistances, FALSE);
    }
  else if(isStart)
    {
      distancesFromLeft = gatherDistances(tr, node->back->next->back, resultBranchDistances, FALSE);
      distancesFromRight = gatherDistances(tr, node->back->next->next->back, resultBranchDistances, FALSE);
    }
  else
    {
      /* recursion ends */
      List *distanceList = CALLOC(1,sizeof(List));
      distanceList->value = CALLOC(1,sizeof(DistanceElem));
      ((DistanceElem*)distanceList->value)->node = node; 
      ((DistanceElem*)distanceList->value)->distance =  1 ; 
      return distanceList;
    }

  /* compute distances of branches from the right to branches from the left side  */
  iterA = distancesFromLeft;
  FOR_LIST(iterA)
    {
      for(iterB = distancesFromRight; iterB; iterB = iterB->next)
	{
	  int
	    a = ((DistanceElem*)iterA->value)->node->number - 1,
	    b = ((DistanceElem*)iterB->value)->node->number - 1;

	  if(b < a)
	    SWAP(a,b);

	  assert(a != b );
	  USE_UPPER_TRIANGLE_TII(resultBranchDistances,a,b) = ((DistanceElem*)iterA->value)->distance  + ((DistanceElem*)iterB->value)->distance;
	}
      
      if( NOT iterA->next)
	end = iterA;
    }

  assert(end);
  end->next = distancesFromRight;

  /* adding own value */
  for(iterA = distancesFromLeft; iterA; iterA = iterA->next)
    ((DistanceElem*)iterA->value)->distance += 1;

  if(isStart)
    {
      int a = node->number - 1;
      for(iterA = distancesFromLeft; iterA; iterA = iterA->next)
	{
	  int b = ((DistanceElem*)iterA->value)->node->number -1;
	  if(b < a)
	    SWAP(a,b);
	  assert(b != a);
	  USE_UPPER_TRIANGLE_TII(resultBranchDistances,a,b) = ((DistanceElem*)iterA->value)->distance;
	}

      /* for(i = 0; i < tr->mxtips; ++i) */
      /* 	for(j = i+1; j < tr->mxtips; ++j) */
      /* 	  printf("%d\t%d\t%d\n", i,j,USE_UPPER_TRIANGLE(resultBranchDistances,i,j)); */
    }

  return distancesFromLeft;
}


double getOneTaxonomicInstability(All *tr, int i, nodeDistance_t ***distances, unsigned int *remainingTaxa)
{
  int 
    j,k,l;
  double
    taxInstab = 0.0;

  FOR_0_LIMIT(j,tr->mxtips)
    {
      if(j == i || NOT  NTH_BIT_IS_SET(remainingTaxa, j))
	continue;
	  
      FOR_0_LIMIT(k, tr->numberOfTrees)
	{
	  for(l = k+1; l < tr->numberOfTrees; ++l)
	    {
	      double
		d1 = (i < j) ? (double)(USE_UPPER_TRIANGLE_TII(distances[k], i, j)) : (double)(USE_UPPER_TRIANGLE_TII(distances[k], j, i)),
		d2 = (i < j) ? (double)(USE_UPPER_TRIANGLE_TII(distances[l], i, j)) : (double)(USE_UPPER_TRIANGLE_TII(distances[l], j, i)),
		diff = (d1 > d2) ? (d1 - d2) : (d2 - d1), 
		sum = d1 + d2;
		  
	      assert( d1 > 0 && d2 > 0);
	      
	      taxInstab += diff / pow(sum, tiiZ); /* NOTE: could also be another exponent, but the mesquite 
						  guys use this one (emphasises close relationship)  */
	    }
	}
    }
  
  return taxInstab;
}


void getTaxonomicInstability(All *tr, char *treesFileName, char *excludeFile)
{
  FILE
    *outf = getOutputFileFromString("taxonomicInstabilityIndex"),
    *treesFile = getNumberOfTrees(tr, treesFileName);

  BitVector
    *neglectThose = neglectThoseTaxa(tr, excludeFile);
  
  int 
    i,j,k;
  nodeDistance_t 
    ***distances = CALLOC(tr->numberOfTrees, sizeof(nodeDistance_t**));

  List
    *iter;

  rewind(treesFile);
  
  FOR_0_LIMIT(i,tr->numberOfTrees)
    {
      distances[i] = CALLOC(tr->mxtips,  sizeof(nodeDistance_t*));
      for(j = 0; j < tr->mxtips; ++j)
	distances[i][j] = CALLOC(tr->mxtips - j, sizeof(nodeDistance_t));
    }
  
  /* calculate distances */
  FOR_0_LIMIT(i,tr->numberOfTrees)    
    {
      readBootstrapTree(tr, treesFile);
      
      FOR_0_LIMIT(j,tr->mxtips)
	if( NOT NTH_BIT_IS_SET(neglectThose, j))
	  pruneTaxon(tr,j+1);

      /* assert that we get the right labels */
      iter = gatherDistances(tr, tr->nodep[1], distances[i], TRUE);
      freeList(iter);
      
      FOR_0_LIMIT(j,tr->mxtips)	
	for(k = j+1; k < tr->mxtips; ++k)
	  {
	    if( NOT USE_UPPER_TRIANGLE_TII(distances[i],j,k))
	      printf("%s\t%s\t%d\t%d\n", tr->nameList[j+1], tr->nameList[k+1], j,k);
	    assert(USE_UPPER_TRIANGLE_TII(distances[i],j,k));
	  }
    }

  /* calculate taxonomic instability */
  FOR_0_LIMIT(i,tr->mxtips)
    {
      if(NTH_BIT_IS_SET(neglectThose,i))
	{
	  double tii = getOneTaxonomicInstability(tr, i, distances, neglectThose);
	  PR("%s\t%f\n", tr->nameList[i+1], tii);
	  fprintf(outf, "%s\t%f\n", tr->nameList[i+1], tii);
	}
      else
	{
	  PR("%s\tNA\n", tr->nameList[i+1]);
	  fprintf(outf, "%s\tNA\n", tr->nameList[i+1]);
	}
    }  
  fclose(outf);
}


void printHelpFile()
{ 
  printVersionInfo(FALSE);
  printf("This program computes the taxonomic instability index.\n\nSYNTAX: ./%s -i <bootTrees> -n <runId> [-w <workingDir>] [-h] [-x <excludeFile>]\n", lowerTheString(programName));
  printf("\nOBLIGATORY:\n");
  printf("-i <bootTrees>\n\tA collection of bootstrap trees.\n");
  printf("-n <runId>\n\tAn identifier for this run.\n");
  printf("\nOPTIONAL\n");
  printf("-z <z>\n\tThe exponent used in the TII formula. Use small values to emphasize close relationships and vice versa. DEFAULT: 2\n");
  printf("-w <workDir>\n\tA working directory where output files are created.\n");
  printf("-x <excludeFile>\n\tExclude the taxa in this file (one taxon per line) prior to computing the TII.\n");
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

   while ((c = getopt (argc, argv, "hi:n:x:w:z:")) != -1)
    {
      switch(c)
	{
	case 'i':
	  bootTrees = optarg;
	  break;
	case 'n':
	  strcpy(run_id, optarg);
	  break; 
	case 'x':
	  excludeFile =  optarg;
	  break;
	case 'w':
	  strcpy(workdir, optarg);
	  break;
	case 'z':
	  tiiZ = wrapStrToL(optarg); 
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

  compute_bits_in_16bits();
  initializeMask();

  All 
    *tr = CALLOC(1,sizeof(All));  
  setupInfoFile();
  if (NOT setupTree(tr, bootTrees))
    {
      PR("Something went wrong during tree initialisation. Sorry.\n");
      exit(-1);
    }     

  getTaxonomicInstability(tr, bootTrees, excludeFile);

  return 0;
}
