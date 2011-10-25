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


#include "common.h"

#define ALLOW_OVERWRITE_INFO_FILE TRUE

extern void printHelpFile();

/* GLOBAL VARIABLES =(  */
char *infoFileName = "",
  *programName = "",
  *programVersion = "", 
  *programReleaseDate ="",
  run_id[128] = "",
  workdir[1024] = "";


double updateTime(double* time)
{
  double tmp = *time;
  return (*time = gettime()) - tmp;
}


int wrapStrToL(char *string)
{
  int errno = 0; 
  char *p; 
  int result = strtol(string, &p, 10);
  if (errno != 0 || *p != 0 || p == string)
    {
      printf("Something went wrong parsing a number.");
      exit(-1);
    }
  return result; 
}


double wrapStrToDouble(char *string)
{
  int errno = 0; 
  char *p; 
  double result = strtod(string, &p);
  if (errno != 0 || *p != 0 || p == string)
    {
      printf("Something went wrong parsing a number.");
      exit(-1);
    }
  return result; 
}


int filexists(char *filename)
{
  FILE *fp;
  int res;
  fp = fopen(filename,"rb");

  if(fp)
    {
      res = 1;
      fclose(fp);
    }
  else
    res = 0;

  return res;
}

void  printVersionInfo(boolean toInfoFile)
{
  if(toInfoFile)
    PR("\nThis is %s version %s released by Andre J. Aberer in %s.\n\n", programName, programVersion, programReleaseDate);
  else
    printf("\nThis is %s version %s released by Andre J. Aberer in %s.\n\n", programName, programVersion, programReleaseDate);
}


FILE *myfopen(const char *path, const char *mode)
{
  FILE *fp = fopen(path, mode);

  if(strcmp(mode,"r") == 0 || strcmp(mode,"rb") == 0)
    {
      if(fp)
	return fp;
      else
	{
	  if(processID == 0)
	    printf("The file %s you want to open for reading does not exist, exiting ...\n", path);
	  assert(0);
	  exit(-1);
	  return (FILE *)NULL;
	}
    }
  else
    {
      if(fp)
	return fp;
      else
	{
	  if(processID == 0)
	    printf("The file %s RAxML wants to open for writing or appending can not be opened [mode: %s], exiting ...\n",
		   path, mode);
	  exit(-1);
	  return (FILE *)NULL;
	}
    }
}


void setupInfoFile()
{
  char *result = CALLOC(1024, sizeof(char));
  strcpy(result,         workdir);
  strcat(result,         programName);
  strcat(result,         "_info");
  strcat(result,         ".");
  strcat(result,         run_id);
  
  if( NOT ALLOW_OVERWRITE_INFO_FILE && filexists(result))
    {
      printf("The run-id >%s< you specified already exists. Aborting in order to protect the output of this previous run.\n", run_id);
      exit(-1);
    }

  FILE *tmp = myfopen(result, "w");

  fclose(tmp);
  infoFileName = result;
  printVersionInfo(TRUE);
}


char *lowerTheString(char *string)
{
  int
    i,
    stringLength = strlen(string);
  char
    *result = CALLOC(stringLength, sizeof(char));
  
  FOR_0_LIMIT(i,stringLength)
     result[i] = tolower(string[i]);
  
  
  return result;
}


FILE *getOutputFileFromString(char *fileName)
{

  char result[1024];
  strcpy(result,         workdir);
  strcat(result,         programName);
  strcat(result,         "_");
  strcat(result,         fileName);
  strcat(result,         ".");
  strcat(result,         run_id);

  return  myfopen(result, "w");
}


double gettime(void)
{
#ifdef WIN32
  time_t tp;
  struct tm localtm;
  tp = time(NULL);
  localtm = *localtime(&tp);
  return 60.0*localtm.tm_min + localtm.tm_sec;
#else
  struct timeval ttime;
  gettimeofday(&ttime , NULL);
  return ttime.tv_sec + ttime.tv_usec * 0.000001;
#endif
}


void printBothOpen(const char* format, ... )
{
  FILE *f = myfopen(infoFileName, "ab");

  va_list args;
  va_start(args, format);
  vfprintf(f, format, args );
  va_end(args);

  va_start(args, format);
  vprintf(format, args );
  va_end(args);

  fclose(f);
}
