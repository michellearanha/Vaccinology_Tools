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

#include <math.h>
#include <float.h>
#include <string.h>

#include "PlatformMain.h"

#include "globals.h"
 

/* global vars 
************************************/
int Decode[22];
int Code[27];
bool  verbose;
bool writeControls;
double CutoffForWriting, DomainThreshold;
int CutoffForConnecting;
bool FBDdetailsS, FBDdetailsC, FBDdetailsL, FBDdetailsD;
bool bWriteSeq;
double thresholds[6];
int    thresholdNb;
bool optMTIK, optCmatrix;

int 	   seqNb, seqLen;
int      seq[kMaxSeqLen+2];
bool     completedFile;
char	  SeqName[kMaxSeqName+2];
FILE *cFpWritePP, *cFpWritePPD,*cFpWriteCP; // must be open for all seqs
FILE *cFpWriteCoP, *cFpWriteCoL;
FILE *cFpWriteDom; // must be open for all seqs
FILE *cFpReadSeq; // must be open for all seqs
