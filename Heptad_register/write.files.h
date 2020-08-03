/***********************************************
************************************************
**                                            **
** Marcoil                                    **
** a program for predicting coiled-coil       **
** domains in protein sequences               **     
** Copyright (C) 2001 Mauro C. Delorenzi      **
**                                            **
************************************************
************************************************/

#include <cstdlib> // for access to the Unix-system
#include <stdio.h>
#include <stdlib.h> 

#include <iostream> // for  cin, cout, <<  etc
#include <fstream> // for file streams like fout
#include <iomanip>  // for  endl;

#include "globals.h"

#include <math.h>
#include <float.h>
#include <string.h>

/***********************************
* function-declarations  */

void  WriteWarning (int seqNb, int  n);
void  ParseIntoDomains (int totalNumber,  float BgrProbProfile[], int thresholdNb, const
		double thresholds[]);

void  WriteSeq (int seqNb,  int seqLen, const int seq[], const char   SeqName[150]); 
void  WriteSeqP (int seqNb,  int seqLen, const int seq[], const char   SeqName[150]) ;
void  WriteSeqId (int seqNb,  const char   SeqName[150]) ;
void  WriteSeqIdP (int seqNb,  const char   SeqName[150]); 

/* function-declarations 
************************************/
