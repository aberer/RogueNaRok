#include "ProfileElem.h"


CLONE_ARRAY_FLAT(cloneProfileArrayFlat, ProfileElem*, ProfileElemAttr*)


void addElemToArray(ProfileElem *elem, Array *array)
{
  GET_PROFILE_ELEM(array, array->length) = elem;
  array->length++;
}


void freeProfileElem(ProfileElem *elem)
{
  free(elem->treeVector);
  free(elem->bitVector);  
  free(elem);
}

Array* profileToArray(HashTable *profile, boolean updateFrequencyCount, boolean assignIds)
{
  HashTableIterator* 
    hashTableIterator = createHashTableIterator(profile);
  
  Array 
    *result = CALLOC(1, sizeof(Array));
  
  unsigned int 
    count = 0;

  /* remember to always copy s.t. free() runs w/o problems */
  
  ProfileElemAttr 
    *profileElemAttr;
  
  result->commonAttributes = CALLOC(1, sizeof(ProfileElemAttr));
  result->commonAttributes = memcpy(result->commonAttributes, profile->commonAttributes, sizeof(ProfileElemAttr));
  profileElemAttr = result->commonAttributes;

  result->length = profile->entryCount;
  result->arrayTable = CALLOC(profile->entryCount, sizeof(ProfileElem*));

  if( NOT hashTableIterator)
    return result;
  
  do
    {
      ProfileElem *profileElem = getCurrentValueFromHashTableIterator(hashTableIterator);

      assert(profileElem);

      if(updateFrequencyCount)
	profileElem->treeVectorSupport = genericBitCount(profileElem->treeVector, profileElemAttr->treeVectorLength);

      if(assignIds)
	profileElem->id = count;
      
      ((ProfileElem**)result->arrayTable)[count] = profileElem;
      assert(profileElem->bitVector && profileElem->treeVector);
      count++;
    }
  while(hashTableIteratorNext(hashTableIterator));
  
  assert(count == profile->entryCount);
  
  free(hashTableIterator);
  
  return result;
}

int sortBipProfile(const void *a, const void *b)
{
  if(*((ProfileElem**)a) == NULL)
    return 1; 
  if(*((ProfileElem**)b) == NULL)
    return -1;
  
  int
    as = (*(ProfileElem**)a)->numberOfBitsSet,
    bs = (*(ProfileElem**)b)->numberOfBitsSet;

  if(as == bs)
    return 0; 

  return as < bs ? -1 : 1; 
}


int sortById(const void *a, const void *b)
{
  int
    as = (*(ProfileElem**)a)->id,
    bs = (*(ProfileElem**)b)->id;

  if(as == bs)
    return 0; 

  return as < bs ? -1 : 1; 
}


int sortBySupport(const void *a, const void *b)
{
  unsigned int
    as = (*(ProfileElem**)a)->treeVectorSupport,
    bs = (*(ProfileElem**)b)->treeVectorSupport;

  if(as == bs)
    return 0; 
  return as > bs ? -1 : 1; 
}


/* what is the index in the (ordered) profile of the first element to
   have at least i bits set?  */
int *createNumBitIndex(Array *bipartitionProfile, int mxtips)
{
  int *result  = CALLOC(mxtips, sizeof(int));   
  memset(result, -1, mxtips * sizeof(int));
  qsort(bipartitionProfile->arrayTable, bipartitionProfile->length, sizeof(ProfileElem**), sortBipProfile); 
  
  int
    i,    
    max = 0,
    current = 0; 
  
  FOR_0_LIMIT(i,bipartitionProfile->length)
    {
      ProfileElem
	*elem = GET_PROFILE_ELEM(bipartitionProfile, i);

      if(NOT elem )
	break;

      if(elem->numberOfBitsSet != current)
	{
	  current = elem->numberOfBitsSet;
	  result[current] = i;
	  max = i; 
	}
    }
  
  current = max;
  for(i = mxtips-1; i >= 0; --i)
    {
      if(result[i] == -1)
	result[i] = current;
      else
  	current = result[i];
    }
  
  return result;
}
