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


#ifndef COMMON_H
#define COMMON_H

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdarg.h>
#include <assert.h>
#include <math.h>
#include <ctype.h>

#ifdef WIN32
#include <direct.h>
#endif

#ifndef WIN32
#include <sys/times.h>
#include <sys/types.h>
#include <sys/time.h>
#include <sys/mman.h>
#include <unistd.h>
#endif

typedef char boolean;

/* GENERAL  */
#define NOT ! 
#define TRUE             1
#define FALSE            0
#define CALLOC(num, size) calloc(num, size)
#define ABS(x)    (((x)<0)   ?  (-(x)) : (x))
#define FOR_0_LIMIT(iter,limit) for(iter=0;iter < (limit); iter++)
#define FOR_N_LIMIT(iter,n,limit) for(iter=(n); iter < (limit); iter++)
#define PR printBothOpen
#define SWAP(a, b) (((a) ^= (b)), ((b) ^= (a)), ((a) ^= (b)))
#define USE_UPPER_TRIANGLE_LSI(matrix,a,b,c) (matrix[a][b][(c)-(b)])
#define USE_UPPER_TRIANGLE_TII(matrix,x,y) (matrix[(x)][(y-x)]) /* assumes that x < y */
#define MIN(a,b) (((a) > (b)) ? (b) : (a))
#define MAX(a,b) (((a) < (b)) ? (b) : (a))

#define GET_FROM_UPPER_TRIANGLE(matrix,a,b,c) ((b<c) ? matrix[a][b][(c)-(b)] : matrix[a][c][(b)-(c)])

int processID;
void  printVersionInfo(boolean toInfoFile);
int wrapStrToL(char *string);
void printBothOpen(const char* format, ... );
double wrapStrToDouble(char *string);
char *lowerTheString(char *string);
FILE *getOutputFileFromString(char *fileName);
double gettime(void);
void setupInfoFile();
double updateTime(double* time);
FILE *myfopen(const char *path, const char *mode);

#endif
