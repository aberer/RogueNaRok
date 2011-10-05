#include "Tree.h"

static int treeGetCh (FILE *fp) ;
static void insertHashBootstop(unsigned int *bitVector, hashtable *h, unsigned int vectorLength, int treeNumber, int treeVectorLength, unsigned int position);
static void  treeEchoContext (FILE *fp1, FILE *fp2, int n);
static double getBranchLength(All *tr, int perGene, nodeptr p);
boolean isTip(int number, int maxTips);
void getxnode (nodeptr p);
static void insertHashAll(unsigned int *bitVector, hashtable *h, unsigned int vectorLength, int treeNumber,  unsigned int position);
static void insertHash(unsigned int *bitVector, hashtable *h, unsigned int vectorLength, int bipNumber, unsigned int position);
static int countHash(unsigned int *bitVector, hashtable *h, unsigned int vectorLength, unsigned int position);


static unsigned int KISS32(void)
{
  static unsigned int 
    x = 123456789, 
    y = 362436069,
    z = 21288629,
    w = 14921776,
    c = 0;

  unsigned int t;

  x += 545925293;
  y ^= (y<<13); 
  y ^= (y>>17); 
  y ^= (y<<5);
  t = z + w + c; 
  z = w; 
  c = (t>>31); 
  w = t & 2147483647;

  return (x+y+w);
}


stringHashtable *initStringHashTable(unsigned int n)
{
  static const unsigned int initTable[] = {53, 97, 193, 389, 769, 1543, 3079, 6151, 12289, 24593, 49157, 98317,
					     196613, 393241, 786433, 1572869, 3145739, 6291469, 12582917, 25165843,
					     50331653, 100663319, 201326611, 402653189, 805306457, 1610612741};
  
  stringHashtable *h = (stringHashtable*)malloc(sizeof(stringHashtable));
  
  unsigned int
    tableSize,
    i,
    primeTableLength = sizeof(initTable)/sizeof(initTable[0]);



  i = 0;

  while(initTable[i] < n && i < primeTableLength)
    i++;

  assert(i < primeTableLength);

  tableSize = initTable[i];  

  h->table = (stringEntry**)calloc(tableSize, sizeof(stringEntry*));
  h->tableSize = tableSize;    

  return h;
}


static unsigned int  hashString(char *p, unsigned int tableSize)
{
  unsigned int h = 0;
  
  for(; *p; p++)
    h = 31 * h + *p;
  
  return (h % tableSize);
}


void addword(char *s, stringHashtable *h, int nodeNumber)
{
  unsigned int position = hashString(s, h->tableSize);
  stringEntry *p = h->table[position];
  
  for(; p!= NULL; p = p->next)
    {
      if(strcmp(s, p->word) == 0)		 
	return;	  	
    }

  p = (stringEntry *)malloc(sizeof(stringEntry));

  assert(p);
  
  p->nodeNumber = nodeNumber;
  p->word = (char *)malloc((strlen(s) + 1) * sizeof(char));

  strcpy(p->word, s);
  /* free(s); */
  
  p->next =  h->table[position];
  
  h->table[position] = p;
}


int getNumberOfTaxa(All *tr, char *bootStrapFile)
{
  FILE *f = myfopen(bootStrapFile, "rb");

  char 
    **nameList,
    buffer[nmlngth + 2]; 

  int
    i = 0,
    c,
    taxaSize = 1024,
    taxaCount = 0;
   
  nameList = (char**)malloc(sizeof(char*) * taxaSize);  

  while((c = fgetc(f)) != ';')
    {
      if(c == '(' || c == ',')
	{
	  c = fgetc(f);
	  if(c ==  '(' || c == ',')
	    ungetc(c, f);
	  else
	    {	      
	      i = 0;	      	     
	      
	      do
		{
		  buffer[i++] = c;
		  c = fgetc(f);
		}
	      while(c != ':' && c != ')' && c != ',');
	      buffer[i] = '\0';	    

	      for(i = 0; i < taxaCount; i++)
		{
		  if(strcmp(buffer, nameList[i]) == 0)
		    {
		      printf("A taxon labelled by %s appears twice in the first tree of tree collection %s, exiting ...\n", buffer, bootStrapFile);
		      exit(-1);
		    }
		}	     
	     
	      if(taxaCount == taxaSize)
		{		  
		  taxaSize *= 2;
		  nameList = (char **)realloc(nameList, sizeof(char*) * taxaSize);		 
		}
	      
	      nameList[taxaCount] = (char*)malloc(sizeof(char) * (strlen(buffer) + 1));
	      strcpy(nameList[taxaCount], buffer);
	     
	      taxaCount++;
			    
	      ungetc(c, f);
	    }
	}   
    }
  
  printf("Found a total of %d taxa in first tree of tree collection %s\n", taxaCount, bootStrapFile);
  printf("Expecting all remaining trees in collection to have the same taxon set\n\n");

  tr->nameList = (char **)malloc(sizeof(char *) * (taxaCount + 1));  
  for(i = 1; i <= taxaCount; i++)
    tr->nameList[i] = nameList[i - 1];
  
  free(nameList);

  tr->nameHash = initStringHashTable(10 * taxaCount);
  for(i = 1; i <= taxaCount; i++)
    addword(tr->nameList[i], tr->nameHash, i);

  fclose(f);

  return taxaCount;
}


boolean setupTree (All *tr, char *bootstrapFile)
{
  nodeptr  p0, p, q;
  int
    i,
    j,
    k,
    tips,
    inter; 

  tips = getNumberOfTaxa(tr, bootstrapFile);
  tr->mxtips = tips;
  
  tips  = tr->mxtips;
  inter = tr->mxtips - 1;
  tr->numberOfTrees = -1;

  if (NOT(p0 = (nodeptr) malloc((tips + 3*inter) * sizeof(node))))
    {
      printf("ERROR: Unable to obtain sufficient tree memory\n");
      return  FALSE;
    }

  if (NOT(tr->nodep = (nodeptr *) malloc((2*tr->mxtips) * sizeof(nodeptr))))
    {
      printf("ERROR: Unable to obtain sufficient tree memory, too\n");
      return  FALSE;
    }

  tr->nodep[0] = (node *) NULL;    /* Use as 1-based array */

  for (i = 1; i <= tips; i++)
    {
      p = p0++;

      p->hash   =  KISS32(); /* hast table stuff */
      p->x      =  0;
      p->number =  i;
      p->next   =  p;
      p->back   = (node *)NULL;
      p->bInf   = (branchInfo *)NULL;

      
      for(k = 0; k < NUM_BRANCHES; k++)
	{
	  p->xs[k]    = 0;
	  p->backs[k] = (nodeptr)NULL;
	}

      for(k = 0; k < VECTOR_LENGTH; k++)
	p->isPresent[k] = 0;

      tr->nodep[i] = p;
    }

  for (i = tips + 1; i <= tips + inter; i++)
    {
      q = (node *) NULL;
      for (j = 1; j <= 3; j++)
	{	 
	  p = p0++;
	  if(j == 1)
	    p->x = 1;
	  else
	    p->x =  0;
	  p->number = i;
	  p->next   = q;
	  p->bInf   = (branchInfo *)NULL;
	  p->back   = (node *) NULL;
	  p->hash   = 0;

	  if(j == 1)
	    for(k = 0; k < NUM_BRANCHES; k++)
	      {
		p->xs[k]    = 1;
		p->backs[k] = (nodeptr)NULL;
	      }
	  else
	    for(k = 0; k < NUM_BRANCHES; k++)
	      {
		p->xs[k]    = 0;
		p->backs[k] = (nodeptr)NULL;
	      }

	  for(k = 0; k < VECTOR_LENGTH; k++)
	    p->isPresent[k] = 0;


	  q = p;
	}
      p->next->next->next = p;
      tr->nodep[i] = p;
    }

  tr->start       = (node *) NULL;

  tr->ntips       = 0;
  tr->nextnode    = 0;

  for(i = 0; i < tr->numBranches; i++)
    tr->partitionSmoothed[i] = FALSE;

  return TRUE;
}


nodeptr findAnyTip(nodeptr p, int numsp)
{   
  return  isTip(p->number, numsp) ? p : findAnyTip(p->next->back, numsp);
} 


void hookup (nodeptr p, nodeptr q, double *z, int numBranches)
{
  int i;

  p->back = q;
  q->back = p;

  for(i = 0; i < numBranches; i++)
    p->z[i] = q->z[i] = z[i];
}


static void newviewBipartitions(BitVector **bitVectors, nodeptr p, int numsp, unsigned int vectorLength)
{
  if(isTip(p->number, numsp))
    return;
  {
    nodeptr 
      q = p->next->back, 
      r = p->next->next->back;
    unsigned int       
      *vector = bitVectors[p->number],
      *left  = bitVectors[q->number],
      *right = bitVectors[r->number];
    unsigned 
      int i;           

    while(NOT p->x)
      {	
	if(NOT p->x)
	  getxnode(p);
      }

    p->hash = q->hash ^ r->hash;

    if(isTip(q->number, numsp) && isTip(r->number, numsp))
      {		
	for(i = 0; i < vectorLength; i++)
	  vector[i] = left[i] | right[i];	  	
      }
    else
      {	
	if(isTip(q->number, numsp) || isTip(r->number, numsp))
	  {
	    if(isTip(r->number, numsp))
	      {	
		nodeptr tmp = r;
		r = q;
		q = tmp;
	      }	   
	    	    
	    while(NOT r->x)
	      {
		if(NOT r->x)
		  newviewBipartitions(bitVectors, r, numsp, vectorLength);
	      }	   

	    for(i = 0; i < vectorLength; i++)
	      vector[i] = left[i] | right[i];	    	 
	  }
	else
	  {	    
	    while((NOT r->x) || (NOT q->x))
	      {
		if(NOT q->x)
		  newviewBipartitions(bitVectors, q, numsp, vectorLength);
		if(NOT r->x)
		  newviewBipartitions(bitVectors, r, numsp, vectorLength);
	      }	   	    	    	    	   

	    for(i = 0; i < vectorLength; i++)
	      vector[i] = left[i] | right[i];	 
	  }

      }     
  }     
}



void bitVectorInitravSpecial(unsigned int **bitVectors, nodeptr p, int numsp, unsigned int vectorLength, hashtable *h, int treeNumber, int function, branchInfo *bInf, int *countBranches, int treeVectorLength, boolean traverseOnly, boolean computeWRF)
{
  if(isTip(p->number, numsp))
    return;
  else
    {
      nodeptr q = p->next;          

      do 
	{
	  bitVectorInitravSpecial(bitVectors, q->back, numsp, vectorLength, h, treeNumber, function, bInf, countBranches, treeVectorLength, traverseOnly, computeWRF);
	  q = q->next;
	}
      while(q != p);
      
      newviewBipartitions(bitVectors, p, numsp, vectorLength);
      
      assert(p->x);

      if(traverseOnly)
	{
	  if(NOT(isTip(p->back->number, numsp)))
	    *countBranches =  *countBranches + 1;
	  return;
	}

      if(NOT(isTip(p->back->number, numsp)))
	{
	  unsigned int *toInsert  = bitVectors[p->number];
	  unsigned int position = p->hash % h->tableSize;

	  switch(function)
	    {
	    case BIPARTITIONS_ALL:	      
	      insertHashAll(toInsert, h, vectorLength, treeNumber, position);
	      *countBranches =  *countBranches + 1;	
	      break;
	    case GET_BIPARTITIONS_BEST:	   	     
	      insertHash(toInsert, h, vectorLength, *countBranches, position);	     
	      
	      p->bInf            = &bInf[*countBranches];
	      p->back->bInf      = &bInf[*countBranches];        
	      p->bInf->support   = 0;	  	 
	      p->bInf->oP = p;
	      p->bInf->oQ = p->back;
	      
	      *countBranches =  *countBranches + 1;		
	      break;
	    case DRAW_BIPARTITIONS_BEST:
	      {
		int found = countHash(toInsert, h, vectorLength, position);
		if(found >= 0)
		  bInf[found].support =  bInf[found].support + 1;
		*countBranches =  *countBranches + 1;
	      }	      
	      break;
	    case BIPARTITIONS_BOOTSTOP:	      
	      insertHashBootstop(toInsert, h, vectorLength, treeNumber, treeVectorLength, position);
	      *countBranches =  *countBranches + 1;
	      break;
	    default:
	      assert(0);
	    }	  	  
	}
      
    }
}




void hookupDefault (nodeptr p, nodeptr q, int numBranches)
{
  int i;

  p->back = q;
  q->back = p;

  for(i = 0; i < numBranches; i++)
    p->z[i] = q->z[i] = defaultz;
}


static boolean treeLabelEnd (int ch)
{
  switch (ch) 
    {
    case EOF:  
    case '\0':  
    case '\t':  
    case '\n':  
    case '\r': 
    case ' ':
    case ':':  
    case ',':   
    case '(':   
    case ')':  
    case ';':
      return TRUE;
    default:
      break;
    }
  return FALSE;
}

static boolean  treeGetLabel (FILE *fp, char *lblPtr, int maxlen)
{
  int      ch;
  boolean  done, quoted, lblfound;

  if (--maxlen < 0) 
    lblPtr = (char *) NULL; 
  else 
    if (lblPtr == NULL) 
      maxlen = 0;

  ch = getc(fp);
  done = treeLabelEnd(ch);

  lblfound = NOT done;
  quoted = (ch == '\'');
  if (quoted && NOT done) 
    {
      ch = getc(fp); 
      done = (ch == EOF);
    }

  while (NOT done) 
    {
      if (quoted) 
	{
	  if (ch == '\'') 
	    {
	      ch = getc(fp); 
	      if (ch != '\'') 
		break;
	    }
        }
      else 
	if (treeLabelEnd(ch)) break;     

      if (--maxlen >= 0) *lblPtr++ = ch;
      ch = getc(fp);
      if (ch == EOF) break;
    }

  if (ch != EOF)  (void) ungetc(ch, fp);

  if (lblPtr != NULL) *lblPtr = '\0';

  return lblfound;
}

static boolean  treeFlushLabel (FILE *fp)
{ 
  return  treeGetLabel(fp, (char *) NULL, (int) 0);
} 

int lookupWord(char *s, stringHashtable *h)
{
  unsigned int position = hashString(s, h->tableSize);
  stringEntry *p = h->table[position];
  
  for(; p!= NULL; p = p->next)
    {
      if(strcmp(s, p->word) == 0)		 
	return p->nodeNumber;	  	
    }

  return -1;
}


int treeFindTipByLabelString(char  *str, All *tr)
{
  int lookup = lookupWord(str, tr->nameHash);

  if(lookup > 0)
    {
      /* assert(! tr->nodep[lookup]->back); */
      return lookup;
    }
  else
    { 
      printf("ERROR: Cannot find tree species: %s\n", str);
      return  0;
    }
}


int treeFindTipName(FILE *fp, All *tr)
{
  char    str[nmlngth+2];
  int      n;

  if(treeGetLabel(fp, str, nmlngth+2))
    n = treeFindTipByLabelString(str, tr);
  else
    n = 0;
   

  return  n;
}


static boolean treeProcessLength (FILE *fp, double *dptr)
{
  int  ch;
  
  if ((ch = treeGetCh(fp)) == EOF)  return FALSE;    /*  Skip comments */
  (void) ungetc(ch, fp);
  
  if (fscanf(fp, "%lf", dptr) != 1) {
    printf("ERROR: treeProcessLength: Problem reading branch length\n");
    treeEchoContext(fp, stdout, 40);
    printf("\n");
    return  FALSE;
  }
  
  return  TRUE;
}


static int treeFlushLen (FILE  *fp)
{
  double  dummy;  
  int     ch;
  
  ch = treeGetCh(fp);
  
  if (ch == ':') 
    {
      ch = treeGetCh(fp);
      
      ungetc(ch, fp);
      if(NOT treeProcessLength(fp, & dummy)) return 0;
      return 1;	  
    }
  
  
  
  if (ch != EOF) (void) ungetc(ch, fp);
  return 1;
}


boolean isTip(int number, int maxTips)
{
  assert(number > 0);

  if(number <= maxTips)
    return TRUE;
  else
    return FALSE;
}


static char *Tree2StringREC(char *treestr, All *tr, nodeptr p, boolean printBranchLengths, boolean printNames, boolean printLikelihood, boolean rellTree, boolean finalPrint, int perGene, boolean branchLabelSupport, boolean printSHSupport)
{
  char  *nameptr;            
      
  if(isTip(p->number, tr->mxtips)) 
    {	       	  
      if(printNames)
	{
	  nameptr = tr->nameList[p->number];     
	  sprintf(treestr, "%s", nameptr);
	}
      else
	sprintf(treestr, "%d", p->number);    
	
      while (*treestr) treestr++;
    }
  else 
    {                 	 
      *treestr++ = '(';
      treestr = Tree2StringREC(treestr, tr, p->next->back, printBranchLengths, printNames, printLikelihood, rellTree, 
			       finalPrint, perGene, branchLabelSupport, printSHSupport);
      *treestr++ = ',';
      treestr = Tree2StringREC(treestr, tr, p->next->next->back, printBranchLengths, printNames, printLikelihood, rellTree, 
			       finalPrint, perGene, branchLabelSupport, printSHSupport);
      if(p == tr->start->back) 
	{
	  *treestr++ = ',';
	  treestr = Tree2StringREC(treestr, tr, p->back, printBranchLengths, printNames, printLikelihood, rellTree, 
				   finalPrint, perGene, branchLabelSupport, printSHSupport);
	}
      *treestr++ = ')';                    
    }

  if(p == tr->start->back) 
    {	      	 
      if(printBranchLengths && NOT rellTree)
	sprintf(treestr, ":0.0;\n");
      else
	sprintf(treestr, ";\n");	 	  	
    }
  else 
    {                   
      if(rellTree || branchLabelSupport || printSHSupport)
	{	 	 
	  if(( NOT isTip(p->number, tr->mxtips)) && 
	     ( NOT isTip(p->back->number, tr->mxtips)))
	    {			      
	      assert(p->bInf != (branchInfo *)NULL);
	      
	      if(rellTree)
		sprintf(treestr, "%d:%8.20f", p->bInf->support, p->z[0]);
	      if(branchLabelSupport)
		sprintf(treestr, ":%8.20f[%d]", p->z[0], p->bInf->support);
	      if(printSHSupport)
		sprintf(treestr, ":%8.20f[%d]", getBranchLength(tr, perGene, p), p->bInf->support);
	      
	    }
	  else		
	    {
	      if(rellTree || branchLabelSupport)
		sprintf(treestr, ":%8.20f", p->z[0]);	
	      if(printSHSupport)
		sprintf(treestr, ":%8.20f", getBranchLength(tr, perGene, p));
	    }
	}
      else
	{
	  if(printBranchLengths)	    
	    /* sprintf(treestr, ":%8.20f", getBranchLength(tr, perGene, p));	      	    */
	    sprintf(treestr, ":%8.20f", p->z[0]);
	  else	    
	    sprintf(treestr, "%s", "\0");	    
	}      
    }
  
  while (*treestr) treestr++;
  return  treestr;
}


static double getBranchLength(All *tr, int perGene, nodeptr p)
{
  double 
    z = 0.0,
    x = 0.0;

  assert(perGene != NO_BRANCHES);
  assert(tr->fracchange != -1.0);
  z = p->z[0];
  if (z < zmin) 
    z = zmin;      	 
  
  x = -log(z) * tr->fracchange;           

  return x;			/* here! */
}


static nodeptr uprootTree (All *tr, nodeptr p, boolean readBranchLengths, boolean readConstraint)
{
  nodeptr  q, r, s, start;
  int      n, i;              

  for(i = tr->mxtips + 1; i < 2 * tr->mxtips - 1; i++)
    assert(i == tr->nodep[i]->number);

  
  if(isTip(p->number, tr->mxtips) || p->back) 
    {
      printf("ERROR: Unable to uproot tree.\n");
      printf("       Inappropriate node marked for removal.\n");
      assert(0);
    }
  
  assert(p->back == (nodeptr)NULL);
  
  tr->nextnode = tr->nextnode - 1;

  assert(tr->nextnode < 2 * tr->mxtips);
  
  n = tr->nextnode;               
  
  assert(tr->nodep[tr->nextnode]);

  if (n != tr->mxtips + tr->ntips - 1) 
    {
      printf("ERROR: Unable to uproot tree.  Inconsistent\n");
      printf("       number of tips and nodes for rooted tree.\n");
      assert(0);
    }

  q = p->next->back;                  /* remove p from tree */
  r = p->next->next->back;
  assert(p->back == (nodeptr)NULL);
    
  if(readBranchLengths)
    {
      double b[NUM_BRANCHES];
      int i;
      for(i = 0; i < tr->numBranches; i++)
	b[i] = (r->z[i] + q->z[i]);
      hookup (q, r, b, tr->numBranches);
    }
  else    
    hookupDefault(q, r, tr->numBranches);    

  if(readConstraint && tr->grouped)
    {    
      if(tr->constraintVector[p->number] != 0)
	{
	  printf("Root node to remove should have top-level grouping of 0\n");
	  assert(0);
	}
    }  
 
  assert(NOT(isTip(r->number, tr->mxtips) && isTip(q->number, tr->mxtips))); 

  assert(p->number > tr->mxtips);

  if(tr->ntips > 2 && p->number != n) 
    {    	
      q = tr->nodep[n];            /* transfer last node's conections to p */
      r = q->next;
      s = q->next->next;
      
      if(readConstraint && tr->grouped)	
	tr->constraintVector[p->number] = tr->constraintVector[q->number];       
      
      hookup(p,             q->back, q->z, tr->numBranches);   /* move connections to p */
      hookup(p->next,       r->back, r->z, tr->numBranches);
      hookup(p->next->next, s->back, s->z, tr->numBranches);           
      
      q->back = q->next->back = q->next->next->back = (nodeptr) NULL;
    }
  else    
    p->back = p->next->back = p->next->next->back = (nodeptr) NULL;
  
  assert(tr->ntips > 2);
  
  start = findAnyTip(tr->nodep[tr->mxtips + 1], tr->mxtips);
  
  assert(isTip(start->number, tr->mxtips));
  tr->rooted = FALSE;
  return  start;
}


char *Tree2String(char *treestr, All *tr, nodeptr p, 
		  boolean printBranchLengths, boolean printNames, 
		  boolean printLikelihood, boolean rellTree, 
		  boolean finalPrint, int perGene, 
		  boolean branchLabelSupport, boolean printSHSupport)
{ 
  Tree2StringREC(treestr, tr, p, printBranchLengths, printNames, printLikelihood, rellTree, 
		 finalPrint, perGene, branchLabelSupport, printSHSupport);      
  
  while (*treestr) treestr++;
  
  return treestr;
}


boolean whitechar (int ch)
{
  return (ch == ' ' || ch == '\n' || ch == '\t' || ch == '\r');
}


static void  treeEchoContext (FILE *fp1, FILE *fp2, int n)
{ /* treeEchoContext */
  int      ch;
  boolean  waswhite;
  
  waswhite = TRUE;
  
  while (n > 0 && ((ch = getc(fp1)) != EOF)) {
    if (whitechar(ch)) {
      ch = waswhite ? '\0' : ' ';
      waswhite = TRUE;
    }
    else {
      waswhite = FALSE;
    }
    
    if (ch > '\0') {putc(ch, fp2); n--;}
  }
}


int treeFinishCom (FILE *fp, char **strp)
{
  int  ch;
  
  while ((ch = getc(fp)) != EOF && ch != ']') {
    if (strp != NULL) *(*strp)++ = ch;    /* save character  */
    if (ch == '[') {                      /* nested comment; find its end */
      if ((ch = treeFinishCom(fp, strp)) == EOF)  break;
      if (strp != NULL) *(*strp)++ = ch;  /* save closing ]  */
    }
  }
  
  if (strp != NULL) **strp = '\0';        /* terminate string  */
  return  ch;
}


static int treeGetCh (FILE *fp) 
{
  int  ch;

  while ((ch = getc(fp)) != EOF) {
    if (whitechar(ch)) ;
    else if (ch == '[') {                   /* comment; find its end */
      if ((ch = treeFinishCom(fp, (char **) NULL)) == EOF)  break;
    }
    else  break;
  }
  
  return  ch;
}


static boolean treeNeedCh (FILE *fp, int c1, char *where)
{
  int  c2;
  
  if ((c2 = treeGetCh(fp)) == c1)  return TRUE;
  
  printf("ERROR: Expecting '%c' %s tree; found:", c1, where);
  if (c2 == EOF) 
    {
      printf("End-of-File");
    }
  else 
    {      	
      ungetc(c2, fp);
      treeEchoContext(fp, stdout, 40);
    }
  putchar('\n');

  if(c1 == ':')    
    printf("RAxML may be expecting to read a tree that contains branch lengths\n");

  return FALSE;
}


static boolean addElementLen (FILE *fp, All *tr, nodeptr p, boolean readBranchLengths, boolean readNodeLabels, int *lcount)
{   
  nodeptr  q;
  int      n, ch, fres;
  
  if ((ch = treeGetCh(fp)) == '(') 
    { 
      n = (tr->nextnode)++;
      if (n > 2*(tr->mxtips) - 2) 
	{
	  if (tr->rooted || n > 2*(tr->mxtips) - 1) 
	    {
	      printf("ERROR: Too many internal nodes.  Is tree rooted?\n");
	      printf("       Deepest splitting should be a trifurcation.\n");
	      return FALSE;
	    }
	  else 
	    {
	      assert(NOT readNodeLabels);
	      tr->rooted = TRUE;
	    }
	}
      
      q = tr->nodep[n];

      if (NOT addElementLen(fp, tr, q->next, readBranchLengths, readNodeLabels, lcount))        return FALSE;
      if (NOT treeNeedCh(fp, ',', "in"))             return FALSE;
      if (NOT addElementLen(fp, tr, q->next->next, readBranchLengths, readNodeLabels, lcount))  return FALSE;
      if (NOT treeNeedCh(fp, ')', "in"))             return FALSE;
      
      if(readNodeLabels)
	{
	  char label[64];
	  int support;

	  if(treeGetLabel (fp, label, 10))
	    {	
	      int val = sscanf(label, "%d", &support);      
	      assert(val == 1);
	      if(val != 1 )
		PR("error\n");

	      /*printf("LABEL %s Number %d\n", label, support);*/
	      p->support = q->support = support;
	      /*printf("%d %d %d %d\n", p->support, q->support, p->number, q->number);*/
	      assert(p->number > tr->mxtips && q->number > tr->mxtips);
	      *lcount = *lcount + 1;
	    }
	}
      else	
	(void) treeFlushLabel(fp);
    }
  else 
    {   
      ungetc(ch, fp);
      if ((n = treeFindTipName(fp, tr)) <= 0)          return FALSE;
      q = tr->nodep[n];
      if (tr->start->number > n)  tr->start = q;
      (tr->ntips)++;
    }
  
  if(readBranchLengths)
    {
      double branch;
      if (NOT treeNeedCh(fp, ':', "in"))                 return FALSE;
      if (NOT treeProcessLength(fp, &branch))            return FALSE;
      
      /* printf("Branch %8.20f %d\n", branch, tr->numBranches); */
      hookup(p, q, &branch, tr->numBranches);
      /* PR(">>> %d\n", tr->numBranches); */
    }
  else
    {
      fres = treeFlushLen(fp);
      if(NOT fres) return FALSE;
      
      hookupDefault(p, q, tr->numBranches);
    }
  return TRUE;          
}


int getTreeStringLength(char *fileName)
{
  int
    i = 0;
  char
    cc;
  FILE
    *f = myfopen(fileName,"r"); 

  while((cc=getc(f))!='\n')  i++;
  fclose(f);

  return i; 
}


int treeReadLen (FILE *fp, All *tr, 
		 boolean readBranches, boolean readNodeLabels, 
		 boolean topologyOnly, boolean completeTree)
{
  nodeptr  
    p;
  
  int      
    i, 
    ch, 
    lcount = 0; 

  for (i = 1; i <= tr->mxtips; i++) 
    {
      tr->nodep[i]->back = (nodeptr) NULL; 
      if(topologyOnly)
	tr->nodep[i]->support = -1;
    }

  for(i = tr->mxtips + 1; i < 2 * tr->mxtips; i++)
    {
      tr->nodep[i]->back = (nodeptr)NULL;
      tr->nodep[i]->next->back = (nodeptr)NULL;
      tr->nodep[i]->next->next->back = (nodeptr)NULL;
      tr->nodep[i]->number = i;
      tr->nodep[i]->next->number = i;
      tr->nodep[i]->next->next->number = i;

      if(topologyOnly)
	{
	  tr->nodep[i]->support = -2;
	  tr->nodep[i]->next->support = -2;
	  tr->nodep[i]->next->next->support = -2;
	}
    }

  if(topologyOnly)
    tr->start       = tr->nodep[tr->mxtips];
  else
    tr->start       = tr->nodep[1];

  tr->ntips       = 0;
  tr->nextnode    = tr->mxtips + 1;      
 
  for(i = 0; i < tr->numBranches; i++)
    tr->partitionSmoothed[i] = FALSE;
  
  tr->rooted      = FALSE;     

  p = tr->nodep[(tr->nextnode)++]; 
  
  while((ch = treeGetCh(fp)) != '(');
  
  if(NOT topologyOnly)
    assert(readBranches == FALSE && readNodeLabels == FALSE);
  
       
  if (NOT addElementLen(fp, tr, p, readBranches, readNodeLabels, &lcount))                 
    assert(0);
  if (NOT treeNeedCh(fp, ',', "in"))                
    assert(0);
  if (NOT addElementLen(fp, tr, p->next, readBranches, readNodeLabels, &lcount))
    assert(0);
  if (NOT tr->rooted) 
    {
      if ((ch = treeGetCh(fp)) == ',') 
	{ 
	  if (NOT addElementLen(fp, tr, p->next->next, readBranches, readNodeLabels, &lcount))
	    assert(0);	    
	}
      else 
	{                                    /*  A rooted format */
	  tr->rooted = TRUE;
	  if (ch != EOF)  (void) ungetc(ch, fp);
	}	
    }
  else 
    {      
      p->next->next->back = (nodeptr) NULL;
    }
  if (NOT treeNeedCh(fp, ')', "in"))                
    assert(0);

  if(topologyOnly)
    assert(NOT(tr->rooted && readNodeLabels));

  (void) treeFlushLabel(fp);
  
  if (NOT treeFlushLen(fp))                         
    assert(0);
 
  if (NOT treeNeedCh(fp, ';', "at end of"))       
    assert(0);
  
  if (tr->rooted) 
    {     
      assert(NOT readNodeLabels);

      p->next->next->back = (nodeptr) NULL;      
      tr->start = uprootTree(tr, p->next->next, FALSE, FALSE);      
      if (NOT tr->start)                              
	{
	  printf("FATAL ERROR UPROOTING TREE\n");
	  assert(0);
	}    
    }
  else    
    tr->start = findAnyTip(p, tr->mxtips);    
  
   if(NOT topologyOnly)
    {
      if(tr->ntips < tr->mxtips)
	{
	  if(completeTree)
	    {
	      printBothOpen("Hello this is your friendly RAxML tree parsing routine\n");
	      printBothOpen("The RAxML option you are uisng requires to read in only complete trees\n");
	      printBothOpen("with %d taxa, there is at least one tree with %d taxa though ... exiting\n", tr->mxtips, tr->ntips);
	      exit(-1);
	    }
	}
    }
 
  
  return lcount;
}


void getxnode (nodeptr p)
{
  nodeptr  s;

  if ((s = p->next)->x || (s = s->next)->x)
    {
      p->x = s->x;
      s->x = 0;
    }

  assert(p->x);
}


entry *initEntry(void)
{
  entry *e = (entry*)CALLOC(1,sizeof(entry));

  e->bitVector     = (unsigned int*)NULL;
  e->treeVector    = (unsigned int*)NULL;
  e->supportVector = (int*)NULL;
  e->bipNumber  = 0;
  e->bipNumber2 = 0;
  e->supportFromTreeset[0] = 0;
  e->supportFromTreeset[1] = 0;
  e->next       = (entry*)NULL;

  return e;
}


static void insertHashAll(unsigned int *bitVector, hashtable *h, unsigned int vectorLength, int treeNumber,  unsigned int position)
{    
  if(h->table[position] != NULL)
    {
      entry *e = h->table[position];     

      do
	{	 
	  unsigned int i;
	  
	  for(i = 0; i < vectorLength; i++)
	    if(bitVector[i] != e->bitVector[i])
	      break;
	  
	  if(i == vectorLength)
	    {
	      if(treeNumber == 0)
		e->bipNumber = 	e->bipNumber  + 1;
	      else
		e->bipNumber2 = e->bipNumber2 + 1;
	      return;
	    }
	  
	  e = e->next;	 
	}
      while(e != (entry*)NULL); 

      e = initEntry(); 
  
      e->bitVector  = (unsigned int*)CALLOC(vectorLength, sizeof(unsigned int));
      /* e->bitVector = (unsigned int*)malloc_aligned(vectorLength * sizeof(unsigned int)); */
      memset(e->bitVector, 0, vectorLength * sizeof(unsigned int));


      memcpy(e->bitVector, bitVector, sizeof(unsigned int) * vectorLength);

      if(treeNumber == 0)	
	e->bipNumber  = 1;       	
      else		 
	e->bipNumber2 = 1;
	
      e->next = h->table[position];
      h->table[position] = e;              
    }
  else
    {
      entry *e = initEntry(); 
  
      e->bitVector  = (unsigned int*)CALLOC(vectorLength, sizeof(unsigned int));
      /* e->bitVector = (unsigned int*)malloc_aligned(vectorLength * sizeof(unsigned int)); */
      memset(e->bitVector, 0, vectorLength * sizeof(unsigned int));

      memcpy(e->bitVector, bitVector, sizeof(unsigned int) * vectorLength);

      if(treeNumber == 0)	
	e->bipNumber  = 1;	  	
      else    
	e->bipNumber2 = 1;	

      h->table[position] = e;
    }

  h->entryCount =  h->entryCount + 1;
}


static void insertHash(unsigned int *bitVector, hashtable *h, unsigned int vectorLength, int bipNumber, unsigned int position)
{
  entry *e = initEntry();

  e->bipNumber = bipNumber; 
  /*e->bitVector = (unsigned int*)CALLOC(vectorLength, sizeof(unsigned int)); */

  e->bitVector = (unsigned int*)CALLOC(vectorLength , sizeof(unsigned int));
  memset(e->bitVector, 0, vectorLength * sizeof(unsigned int));
 
  memcpy(e->bitVector, bitVector, sizeof(unsigned int) * vectorLength);
  
  if(h->table[position] != NULL)
    {
      e->next = h->table[position];
      h->table[position] = e;           
    }
  else
    h->table[position] = e;

  h->entryCount =  h->entryCount + 1;
}


static int countHash(unsigned int *bitVector, hashtable *h, unsigned int vectorLength, unsigned int position)
{ 
  if(h->table[position] == NULL)         
    return -1;
  {
    entry *e = h->table[position];     

    do
      {	 
	unsigned int i;

	for(i = 0; i < vectorLength; i++)
	  if(bitVector[i] != e->bitVector[i])
	    goto NEXT;
	   
	return (e->bipNumber);	 
      NEXT:
	e = e->next;
      }
    while(e != (entry*)NULL); 
     
    return -1;   
  }
}


static void insertHashBootstop(unsigned int *bitVector, hashtable *h, unsigned int vectorLength, int treeNumber, int treeVectorLength, unsigned int position)
{    
  if(h->table[position] != NULL)
    {
      entry *e = h->table[position];     

      do
	{	 
	  unsigned int i;
	  
	  for(i = 0; i < vectorLength; i++)
	    if(bitVector[i] != e->bitVector[i])
	      break;
	  
	  if(i == vectorLength)
	    {
	      e->treeVector[treeNumber / MASK_LENGTH] |= mask32[treeNumber % MASK_LENGTH];
	      return;
	    }
	  
	  e = e->next;
	}
      while(e != (entry*)NULL); 

      e = initEntry(); 

      e->bipNumber = h->entryCount;
       
      /*e->bitVector  = (unsigned int*)CALLOC(vectorLength, sizeof(unsigned int));*/
      e->bitVector = (unsigned int*)CALLOC(vectorLength, sizeof(unsigned int));
      memset(e->bitVector, 0, vectorLength * sizeof(unsigned int));


      e->treeVector = (unsigned int*)CALLOC(treeVectorLength, sizeof(unsigned int));
      
      e->treeVector[treeNumber / MASK_LENGTH] |= mask32[treeNumber % MASK_LENGTH];
      memcpy(e->bitVector, bitVector, sizeof(unsigned int) * vectorLength);
     
      e->next = h->table[position];
      h->table[position] = e;          
    }
  else
    {
      entry *e = initEntry(); 

      e->bipNumber = h->entryCount;

      e->bitVector = (unsigned int*)CALLOC(vectorLength , sizeof(unsigned int));
      memset(e->bitVector, 0, vectorLength * sizeof(unsigned int));

      e->treeVector = (unsigned int*)CALLOC(treeVectorLength, sizeof(unsigned int));

      e->treeVector[treeNumber / MASK_LENGTH] |= mask32[treeNumber % MASK_LENGTH];
      memcpy(e->bitVector, bitVector, sizeof(unsigned int) * vectorLength);     

      h->table[position] = e;
    }

  h->entryCount =  h->entryCount + 1;
}


FILE *getNumberOfTrees(All *tr, char *fileName)
{
  FILE 
    *f = myfopen(fileName, "r");

  int 
    trees = 0,
    ch;

  while((ch = fgetc(f)) != EOF)
    if(ch == ';')
      trees++;

  assert(trees > 0);

  tr->numberOfTrees = trees;

  rewind(f);

  return f;
}


/* INTERFACE TO OUTSIDE WOLRD */
void readBestTree(All *tr, FILE *file)
{
  treeReadLen(file, tr, TRUE, FALSE, TRUE,  TRUE);
}


void readBootstrapTree(All *tr, FILE *file)
{
  treeReadLen(file, tr, FALSE, FALSE, TRUE, TRUE);
}


char *writeTreeToString(All *tr, boolean printBranchLengths)
{
  Tree2String(tr->tree_string, tr, tr->start->back, /*  */
	      printBranchLengths, TRUE, 
	      FALSE, FALSE, 
	      TRUE, 0,
	      FALSE, FALSE);
  return tr->tree_string;
}

void freeTree(All *tr)
{
  int i; 
  for(i = 1; i <= tr->mxtips; ++i)
    free(tr->nameList[i]);
  free(tr->nameList);

  for(i = 0; i < tr->nameHash->tableSize; ++i)
    {
      stringEntry *elem = tr->nameHash->table[i];
      while(elem)
	{
	  stringEntry *iter = elem->next;
	  free(elem->word);
	  free(elem);	  
	  elem = iter;
	}      
    }
  free(tr->nameHash->table);
  free(tr->nameHash);  
  /* FOR_0_LIMIT(i,((tr->mxtips-1) )) */
  /*   free(tr->nodep[i]); */
  free(tr->nodep);

  free(tr);
}
/* END OF INTERFACE */


