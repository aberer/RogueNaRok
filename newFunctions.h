#ifndef NEW_FUNCTIONS_H
#define NEW_FUNCTIONS_H

#include "List.h"
#include "common.h"
#include "legacy.h"
#include "Tree.h"
#include "ProfileElem.h"

IndexList *parseToDrop(All *tr, FILE *toDrop);
void pruneTaxon(All *tr, unsigned int k);
BitVector *neglectThoseTaxa(All *tr, char *toDrop);
Array *getOriginalBipArray(All *tr, FILE *bestTree, FILE *treeFile) ;


#endif
