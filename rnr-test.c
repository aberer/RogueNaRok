#include <CUnit/Basic.h>
#include <CUnit/CUnit.h>
#include "Dropset.h"
#include "List.h"

int mxtips = 32; 


Dropset *dropset = NULL;
HashTable *mergingHash;
int bitVectorLength; 
BitVector *droppedTaxa, *paddingBits; 

extern int maxDropsetSize;


int setup()
{
  bitVectorLength = 1;
  droppedTaxa = CALLOC(1,sizeof(BitVector));
  paddingBits = CALLOC(1,sizeof(BitVector));

  dropset = CALLOC(1,sizeof(Dropset));

  IndexList
    *tmp3 = NULL, 
    *tmp2 = NULL,
    *tmp = NULL; 
  APPEND_INT(1,tmp);
  APPEND_INT(2,tmp);  
  
  APPEND_INT(3,tmp2);
  APPEND_INT(4,tmp2);  

  APPEND_INT(2,tmp3);
  APPEND_INT(3,tmp3);

  addEventToDropset(dropset, tmp);
  addEventToDropset(dropset, tmp2);
  addEventToDropset(dropset, tmp3);

  mergingHash = createHashTable(20, NULL, dropsetHashValue, dropsetEqual);
  
  

  
  return 0; 
}


int tearDown()
{
  freeDropsetDeep(dropset);
  return 0;
}


void test_addEventToDropset()
{
  MergingEvent
    *me = dropset->mergingEvents;
  IndexList    
    *il = me->mergingBipartitions,
    *list = il ;
  
  CU_ASSERT(NOT me->next);
  printf("INFO: %d\n", lengthIndexList(il)) ; 
  puts("\n");
  FOR_LIST(list)
    printf("%d,", list->index);
  puts("\n");
  CU_ASSERT(lengthIndexList(il) == 4);
  
  IndexList *tmp = NULL;
  APPEND_INT(1,tmp);
  APPEND_INT(2,tmp);
  APPEND_INT(3,tmp);
  
  CU_ASSERT(indexListContainsIndexListUnordered(il, tmp));
}

void test_addEventToDropset2()
{
  
}



int main()
{
  CU_initialize_registry();
  CU_pSuite st1 = CU_add_suite("suite1", &setup, &tearDown);
  CU_add_test(st1, "test adding to dropsets", &test_addEventToDropset);
  CU_basic_run_tests();
  CU_cleanup_registry();	
  return 0; 
}


