#ifndef LEGACY_H
#define LEGACY_H

#include "common.h"
#include "BitVector.h"
#include "List.h"
#include "Array.h"
#include "ProfileElem.h"

#define NUM_BRANCHES   128

typedef struct ent
{
  unsigned int *bitVector;
  unsigned int *treeVector;
  unsigned int amountTips;
  int *supportVector;
  unsigned int bipNumber;
  unsigned int bipNumber2;
  unsigned int supportFromTreeset[2]; 
  struct ent *next;
} entry;

typedef struct
{
  unsigned int *vector; 
  int support;   
  struct noderec *oP;
  struct noderec *oQ;
} branchInfo;

typedef  struct noderec
{
  unsigned int    isPresent[NUM_BRANCHES / MASK_LENGTH];
  struct noderec  *backs[NUM_BRANCHES];
  char            xs[NUM_BRANCHES];
  branchInfo      *bInf;
  double           z[NUM_BRANCHES];
  struct noderec  *next;
  struct noderec  *back;
  unsigned int   hash;
  int              support;
  int              number;
  char             x;
  double **insertionLWs;
} node, *nodeptr;

typedef struct stringEnt
{
  int nodeNumber;
  char *word;
  struct stringEnt *next;
} stringEntry ;

typedef struct
{
  unsigned int tableSize;
  stringEntry **table;
}  stringHashtable;


typedef struct _All
{
  nodeptr start;
  int mxtips;
  int numberOfTrees;
  int bitVectorLength;
  nodeptr *nodep;
  int ntips;
  int nextnode; 
  int numBranches;
  boolean partitionSmoothed[NUM_BRANCHES];
  boolean rooted;
  stringHashtable *nameHash;
  boolean grouped;
  int *constraintVector;
  char **nameList;
  char *tree_string;
  int treeStringLength;
  double fracchange;
} All;

typedef struct
{
  unsigned int tableSize;
  entry **table;
  unsigned int entryCount;
}  hashtable;


#define FC_INIT               20
#define zmin       1.0E-15  /* max branch prop. to -log(zmin) (= 34) */
#define NO_BRANCHES      -1

#define BIPARTITIONS_ALL       0
#define GET_BIPARTITIONS_BEST  1
#define DRAW_BIPARTITIONS_BEST 2
#define BIPARTITIONS_BOOTSTOP  3
#define BIPARTITIONS_RF  4
#define defaultz       0.9         /* value of z assigned as starting point */
#define nmlngth        1024         /* number of characters in species name */
#define VECTOR_LENGTH (NUM_BRANCHES / MASK_LENGTH)

void bitVectorInitravSpecial(unsigned int **bitVectors, nodeptr p, int numsp, unsigned int vectorLength, hashtable *h, int treeNumber, int function, branchInfo *bInf, int *countBranches, int treeVectorLength, boolean traverseOnly, boolean computeWRF);
hashtable *initHashTable(unsigned int n);
void freeHashTable(hashtable *h);
ProfileElem *addProfileElem(entry *helem, int vectorLength, int treeVectorLength, int numberOfTrees) ;


BitVector *neglectThoseTaxa(All *tr, char *toDrop);
void pruneTaxon(All *tr, unsigned int k);
BitVector **initBitVector(All *tr, BitVector *vectorLength);
#endif
