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
#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <assert.h>
#include <string.h>

#include "Tree.h"
#include "List.h"
#include "BitVector.h"
#include "sharedVariables.h"
#include "newFunctions.h"

#define PROG_NAME "RnR-prune"
#define PROG_VERSION "1.0"
#define PROG_RELEASE_DATE "2011-10-25"


extern char run_id[128];
extern double masterTime;
extern int tid;
extern int NumberOfThreads; 
extern volatile int NumberOfJobs;

void pruneTaxaFromTreeset(char *bootstrapFileName, char *bestTreeFile, char *toDropFileName, All *tr)
{
  int 
    i;

  IndexList
    *iter;

  FILE 
    *toDrop = myfopen(toDropFileName, "r");
  
  /* drop taxa from bootstrap trees  */
  if( strcmp(bootstrapFileName, ""))
    {
      tr->tree_string = CALLOC(10 * getTreeStringLength(bootstrapFileName), sizeof(char));
      if  (NOT setupTree(tr, bootstrapFileName))
	{
	  PR("Something went wrong during tree initialisation. Sorry.\n");
	  exit(-1);
	}   
      FILE *TMP = getNumberOfTrees(tr, bootstrapFileName);
      fclose(TMP);
    }

  /* drop taxa from best-known tree */
  if( strcmp(bestTreeFile, ""))
    {
      tr->tree_string = CALLOC(10 * getTreeStringLength(bestTreeFile), sizeof(char));
      if  (NOT setupTree(tr, bestTreeFile))
	{
	  PR("Something went wrong during tree initialisation. Sorry.\n");
	  exit(-1);
	}   
      FILE *tmp = getNumberOfTrees(tr, bestTreeFile);
      fclose(tmp);
    }

  IndexList
    *indicesToDrop = parseToDrop(tr, toDrop);

  if( strcmp(bootstrapFileName, ""))
    {   
      FILE
	*bootstrapFile = getNumberOfTrees(tr, bootstrapFileName),
	*outf = getOutputFileFromString("prunedBootstraps");
      
      FOR_0_LIMIT(i,tr->numberOfTrees)
	{
	  readBootstrapTree(tr, bootstrapFile);

	  iter = indicesToDrop;
	  FOR_LIST(iter)
	    pruneTaxon(tr, iter->index, FALSE );

	  char *tmp = writeTreeToString(tr , FALSE);	
	  fprintf(outf, "%s", tmp);
	}  
  
      fclose(outf);
    }

  if( strcmp(bestTreeFile, ""))
    {
      FILE 
	*bestTree = myfopen(bestTreeFile, "r"),
	*outf = getOutputFileFromString("prunedBestTree");
      
      readBestTree(tr, bestTree);
      iter = indicesToDrop;
      FOR_LIST(iter)
	pruneTaxon(tr, iter->index, TRUE );
      
      char *tmp = writeTreeToString(tr, TRUE); 
      fprintf(outf, "%s", tmp);	  
    }
  
  freeIndexList(indicesToDrop);
  exit(EXIT_SUCCESS);
}

static void printHelpFile()
{
  printVersionInfo(FALSE);
  printf("This program prunes a list of taxa from a bootstrap tree or a single tree with branch lengths (such as a best-known ML/MP-tree).\n\nSYNTAX: ./%s [-i <bootTrees> | -t <treeFile>] -x <excludeFile> -n <runId> [-w <workingDir>] [-h]\n", lowerTheString(programName));
  printf("\n\tOBLIGATORY:\n");
  printf("-x <excludeFile>\n\tA list of taxa (one taxon per line) to prune from either the bootstrap trees or the single best-known tree.\n");
  printf("-i <bootTrees>\n\tA collection of bootstrap trees.\n");
  printf("-t <treeFile>\n\tA single tree with branch lengths. Use either this flag or the -i flag.\n");
  printf("-n <runId>\n\tAn identifier for this run.\n");  
  printf("\n\tOPTIONAL:\n");
  printf("-w <workDir>\n\tA working directory where output files are created.\n");
  printf("-h\n\tThis help file.\n");
}


int main(int argc, char *argv[])
{
  int
    c; 
  
  char
    *bestTreeFileName = "",
    *excludeFileName = "",
    *bootTreesFileName = "";

  programName = PROG_NAME; 
  programVersion = PROG_VERSION;
  programReleaseDate = PROG_RELEASE_DATE; 

  while((c = getopt(argc,argv, "hi:t:x:n:w:")) != -1)
    {
      switch(c)
	{
	case 'i':
	  bootTreesFileName = optarg;
	  break;
	case 'w':
	  strcpy(workdir,optarg);
	  break;
	case 'n':
	  strcpy(run_id,optarg);
	  break;
	case 't':
	  bestTreeFileName = optarg;
	  break;
	case 'x':	  
	  excludeFileName = optarg;
	  break;
	case 'h': 
	default:
	  {
	    printHelpFile();
	    abort();
	  }
	}
    }

  if( NOT (strcmp(bootTreesFileName, "") ||  strcmp(bestTreeFileName, "")))
    {
      printf("Please specify either a set of bootstrap trees (-i) or a single best-known tree (-t) from which you want to prune taxa.\n");
      printHelpFile();
      exit(-1);
    }

  if( NOT strcmp(excludeFileName, "") )
    {
      printf("Please specify a file containing taxa to prune (one taxon a line) via -x\n");
      printHelpFile();
      exit(-1);
    }

  if ( NOT strcmp(run_id, ""))
    {
      printf("Please specify a run id via -n.\n");
      printHelpFile();
      exit(-1);
    }

  compute_bits_in_16bits();
  initializeMask();	
  
  setupInfoFile();
   
  All *tr = CALLOC(1,sizeof(All));
  tr->numBranches = 1;
  pruneTaxaFromTreeset(bootTreesFileName, bestTreeFileName, excludeFileName, tr);

  return 0;
}
