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

#include "PlatformMain.h"

extern int 	seqNb, seqLen;
extern int  seq[kMaxSeqLen+1];
extern bool  completedFile;
extern char	 SeqName[kMaxSeqName+2];
extern FILE *cFpWritePP, *cFpWritePPD,*cFpWriteCP; 
extern FILE *cFpWriteDom; 
extern FILE *cFpReadSeq; 
extern bool verbose;
extern bool FBDdetailsD;
extern bool FBDdetailsS, FBDdetailsC, FBDdetailsL;


// ---------------------------------------------------------------------------
//		 PlatformMain 					
// ---------------------------------------------------------------------------

void PlatformMain(const char *transProbFile, const char *emissProbFile, char *seqProbFile, int modus )
{
if (verbose) cout << " PlatformMain, matrixType = " << modus << "  start initCodes" << endl;
initCodes();

if (modus == 0)  FBmain(transProbFile,  emissProbFile,  seqProbFile );
 else
 { 
 if (modus == 1)  CoilsMain(transProbFile,  emissProbFile, seqProbFile, modus); // MTK 
 else
  {
  if (modus == 2)  CoilsMain(transProbFile,  emissProbFile, seqProbFile, modus); // MTIDK 
 else
  {
  if (modus == 3)  CoilsMain(transProbFile,  emissProbFile, seqProbFile, modus); // HMM 
  else { cout << "Inconsistent Modus in PlatformMain, bailing out" ;   exit (1);}
 }}}
}
/* PlatformMain
************************************/
