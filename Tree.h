#ifndef TREES_H
#define TREES_H

#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>

#include "common.h"
#include "legacy.h"

extern unsigned int *mask32;

boolean isTip(int number, int maxTips);
/* char *writeTreeToString(All *tr, boolean printBranchLengths, int whichBranch); */
char *writeTreeToString(All *tr, boolean printBranchLengths);
void readTree(char *fileName);
boolean setupTree (All *tr, char *bootstrapTrees);
void readBestTree(All *tr, FILE *file);  
void readBootstrapTree(All *tr, FILE *file);
void hookupDefault (nodeptr p, nodeptr q, int numBranches);
nodeptr findAnyTip(nodeptr p, int numsp);
int treeFindTipByLabelString(char  *str, All *tr);
int getTreeStringLength(char *fileName);
FILE *getNumberOfTrees(All *tr, char *fileName);
void freeTree(All *tr);
#endif
