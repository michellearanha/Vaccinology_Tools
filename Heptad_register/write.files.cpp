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

#include "write.files.h"

extern int Decode[22];
extern bool verbose;
extern bool FBDdetailsD;
extern bool FBDdetailsS, FBDdetailsC, FBDdetailsL;
extern FILE *cFpWriteDom; 
extern FILE  *cFpWritePP, *cFpWritePPD, *cFpWriteCP;  
extern FILE  *cFpWriteCoP, *cFpWriteCoL; 

// ---------------------------------------------------------------------------
//		 WriteWarning 					
// ---------------------------------------------------------------------------
void  WriteWarning ( int seqNb,  int  n)
{
int i;
std::cerr << "\n SEQUENCE NB " << seqNb << " IS LONGER THAN " << kMaxSeqLen << ", WAS TRUNCATED" << std::endl; 
if (verbose)  std::cout  <<  "\n WriteWarning  concluded " <<  std::endl ; 	
}
/* WriteWarning
************************************/
// ---------------------------------------------------------------------------
//		 WriteSeq 					
// ---------------------------------------------------------------------------
void  WriteSeq (int seqNb,  int seqLen, const int seq[], const char   SeqName[150]) 
{
int i, j=0, k=0;
char c;

if (verbose)
	{
	std::cout << "\n>SEQUENCE NB " << seqNb; 
	std::cout << "\n" << SeqName << "\n";
	for (i=1; i<= seqLen; i++)
	     {
	     j++; k++; c = Decode[seq[i]];
	     std::cout << c;
	     if (j == 60) {j=0; k=0; std::cout << "\n";}
	     else { if (k == 10){k=0; std::cout << " ";}}
	     }
	}		

if (FBDdetailsL)
	 {
	  fprintf(  cFpWritePP, "\n%s   ## %d \n", SeqName , seqNb );
	  for (i=1; i<= seqLen; i++)
 		  {
		  j++; k++; 
		  fprintf(  cFpWritePP, "%c", Decode[seq[i]] );
		  if (j == 60) {j=0; k=0; fprintf(  cFpWritePP,"\n"); }
		  else { if (k == 10){k=0; fprintf(  cFpWritePP," "); }}
		  }
	  fprintf(  cFpWritePP, "*\n\n");j=0; k=0;
	  }

if (FBDdetailsS)
	 {
	  fprintf(  cFpWritePPD, "\n%s   ## %d \n", SeqName , seqNb );
	  for (i=1; i<= seqLen; i++)
 		  {
		  j++; k++; 
		  fprintf(  cFpWritePPD, "%c", Decode[seq[i]] );
		  if (j == 60) {j=0; k=0; fprintf(  cFpWritePPD,"\n"); }
		  else { if (k == 10){k=0; fprintf(  cFpWritePPD," "); }}
		  } // type
	  fprintf(  cFpWritePPD, "*\n");j=0; k=0;
	  }

if (FBDdetailsC)
	 {
	  fprintf(  cFpWriteCP, "\n%s   ## %d \n", SeqName , seqNb );
	  for (i=1; i<= seqLen; i++)
 		  {
		  j++; k++; 
		  fprintf(  cFpWriteCP, "%c", Decode[seq[i]] );
		  if (j == 60) {j=0; k=0; fprintf(  cFpWriteCP,"\n"); }
		  else { if (k == 10){k=0; fprintf(  cFpWriteCP," "); }}
		  } // type
	  fprintf(  cFpWriteCP, "*\n");j=0; k=0;
	  }

if (FBDdetailsD)
	 {
	  fprintf(  cFpWriteDom, "\n%s   ## %d \n", SeqName , seqNb );
	  for (i=1; i<= seqLen; i++)
 		  {
		  j++; k++; 
		  fprintf(  cFpWriteDom, "%c", Decode[seq[i]] );
		  if (j == 60) {j=0; k=0; fprintf(  cFpWriteDom,"\n"); }
		  else { if (k == 10){k=0; fprintf(  cFpWriteDom," "); }}
		  } // type
	  fprintf(  cFpWriteDom, "*\n");
	 }
}

/* WriteSeq
************************************/
// ---------------------------------------------------------------------------
//		 WriteSeqId 					
// ---------------------------------------------------------------------------
void  WriteSeqId (int seqNb,  const char   SeqName[150]) 
{
int i, j=0, k=0;
char c;

if (FBDdetailsL)  fprintf(  cFpWritePP, "\n%s  ## %d \n", SeqName , seqNb );
if (FBDdetailsS)  fprintf(  cFpWritePPD,"\n%s  ## %d \n", SeqName , seqNb );
if (FBDdetailsC)  fprintf(  cFpWriteCP, "\n%s  ## %d \n", SeqName , seqNb );
if (FBDdetailsD)  fprintf(  cFpWriteDom,"\n%s  ## %d \n", SeqName , seqNb );
}

/* WriteSeqId
************************************/
// ---------------------------------------------------------------------------
//		 WriteSeqIdP					
// ---------------------------------------------------------------------------
void  WriteSeqIdP (int seqNb,  const char   SeqName[150]) 
{
int i, j=0, k=0;
char c;

if (FBDdetailsS || FBDdetailsC)  fprintf(  cFpWriteCoP, "\n%s  ## %d \n", SeqName , seqNb );
if (FBDdetailsL)  fprintf(  cFpWriteCoL, "\n%s  ## %d \n", SeqName , seqNb );
if (FBDdetailsD)  fprintf(  cFpWriteDom, "\n%s  ## %d \n", SeqName , seqNb );
}

/* WriteSeqIdP
************************************/
// ---------------------------------------------------------------------------
//		 WriteSeqP					
// ---------------------------------------------------------------------------
void  WriteSeqP (int seqNb,  int seqLen, const int seq[], const char   SeqName[150]) 
{//make second one with only seqname
int i, j=0, k=0;
char c;

if (verbose)
	{
	std::cout << "\n>SEQUENCE NB " << seqNb; 
	std::cout << "\n" << SeqName << "\n";
	for (i=1; i<= seqLen; i++)
	     {
	     j++; k++; c = Decode[seq[i]];
	     std::cout << c;
	     if (j == 60) {j=0; k=0; std::cout << "\n";}
	     else { if (k == 10){k=0; std::cout << " ";}}
	     }
	}		

if ((FBDdetailsS) || (FBDdetailsC))
	 {
	  fprintf(  cFpWriteCoP, "\n%s   ## %d \n", SeqName , seqNb );
	  for (i=1; i<= seqLen; i++)
 		  {
		  j++; k++; 
		  fprintf(  cFpWriteCoP, "%c", Decode[seq[i]] );
		  if (j == 60) {j=0; k=0; fprintf(  cFpWriteCoP,"\n"); }
		  else { if (k == 10){k=0; fprintf(  cFpWriteCoP," "); }}
		  }
	  fprintf(  cFpWriteCoP, "*\n\n");j=0; k=0;
	  }

if (FBDdetailsL )
	 {
	  fprintf(  cFpWriteCoL, "\n%s   ## %d \n", SeqName , seqNb );
	  for (i=1; i<= seqLen; i++)
 		  {
		  j++; k++; 
		  fprintf(  cFpWriteCoL, "%c", Decode[seq[i]] );
		  if (j == 60) {j=0; k=0; fprintf(  cFpWriteCoL,"\n"); }
		  else { if (k == 10){k=0; fprintf(  cFpWriteCoL," "); }}
		  }
	  fprintf(  cFpWriteCoL, "*\n\n");j=0; k=0;
	  }

if (FBDdetailsD)
	 {
	  fprintf(  cFpWriteDom, "\n%s   ## %d \n", SeqName , seqNb );
	  for (i=1; i<= seqLen; i++)
 		  {
		  j++; k++; 
		  fprintf(  cFpWriteDom, "%c", Decode[seq[i]] );
		  if (j == 60) {j=0; k=0; fprintf(  cFpWriteDom,"\n"); }
		  else { if (k == 10){k=0; fprintf(  cFpWriteDom," "); }}
		  } // type
	  fprintf(  cFpWriteDom, "*\n");
	 }
}

/* WriteSeqP
************************************/
// ---------------------------------------------------------------------------
//		 ParseIntoDomains 					
// ---------------------------------------------------------------------------
void  ParseIntoDomains (int seqLen,  float BgrProbProfile[], int thresholdNb, const double thresholds[])
{
int i, k, j, domNb, first[101], last[101];  
float  threshold, Prob,  max, maximum[101];

//bool CCD;

fprintf(cFpWriteDom, ". . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .\n");
BgrProbProfile[seqLen+1] = 1;
for (k=1; k<=thresholdNb; k++)
	{
	domNb=0;threshold = thresholds[k-1];
	for (i=1; i<=seqLen; i++)
		{
		Prob = 1 - BgrProbProfile[i];
		if ( Prob > threshold)
			{
			domNb++; if (domNb>100)  domNb = 100;
			first[domNb] =i; max=Prob;
			while ((Prob > threshold) & (i <= seqLen))
				{i++; Prob=(1 - BgrProbProfile[i]);
				if (Prob > max)   max= Prob;}
			last[domNb] =i-1; 
			maximum[domNb] = max;
		}	}	
	fprintf(cFpWriteDom, "\nNUMBER PREDICTED COILED-COIL DOMAINS WITH THRESHOLD %3.1f : %d",  100 * threshold,domNb );
	for (j=1; j<=domNb; j++)
		{
		fprintf(cFpWriteDom, "\n %2d. from %d to %d (length = %d) with max = %3.1f ", j, first[j], last[j], 
			(last[j] + 2 - first[j]), 100 * maximum[j]);
	}	}
fprintf(cFpWriteDom,"\n\n******************************************************************************"); 
fprintf(cFpWriteDom,"**********************************************************************************\n"); 
}
/* ParseIntoDomains
************************************/

