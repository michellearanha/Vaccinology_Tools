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

#include "read.seqfileC.h"

extern FILE * cFpReadSeq; // must be open for all seqs
extern int	   Code[27];
extern bool	   verbose;	

// ---------------------------------------------------------------------------
//		 ReadSeq 					
// ---------------------------------------------------------------------------
bool  ReadSeq (int seqNb, int &seqLen, int Seq[], bool &completedFile, char SeqName[])			
{
int	 	a, i=0, k, t;
static int		c; // char
bool	initiated,stopProtein, overlong;
stopProtein = false; initiated = false;  overlong = false; 
if (seqNb == 1)  c = (int) '.';

while (( !completedFile) && (i < kMaxSeqLen) && (!stopProtein ))   
	{	
	if (strchr (">", (char) c )!= NULL )  
		{
		if (initiated)	{stopProtein = true;}
		else {
			k = 0;
			while ((c != (int)'\n') && (c  != EOF))  
			{
			if (k < kMaxSeqName) {SeqName[k] = (char) c; t=k;}
			c = (int) fgetc( cFpReadSeq ); k++; }  SeqName[t+1] = '\0'; 
    		}    } 
	if (strchr ("abcdefghijklmnopqrstuvwxyz", c )!= NULL ) 
		{
		c = c -32; // conversion to capital letter
		}
	if (strchr ("ABCDEFGHIJKLMNOPQRSTUVWXYZ", c )!= NULL )
	     {
	     i++;
	     if (i <= kMaxSeqLen)
		     {
		     a = c - kAtoOne;
		     if ((Code[a]> 0) && (Code[a]< 21))
			     {initiated = true; Seq[i]= Code[a];}
		     else
			     {
			     std::cout << "Warning: NonStandard amino acid {" << c << "=" << (char) c << "} at Position " << i << 
			               " in Sequence " << seqNb << std::endl; 
			     
			     initiated = true; Seq[i]= 0;
			}     }
	     else
		     {
		     std::cout << "\nSequence too long, sequence: " << seqNb << std::endl; 
		     overlong = true;
		     while ((c != (int)'>') && (c  != EOF)) { c = (int) fgetc( cFpReadSeq );}// go to end of overlong sequence
	     }	}
	if ((!stopProtein ) && (!overlong) && (c  != EOF))   c = (int) fgetc( cFpReadSeq );
	if  (c == EOF)	{completedFile  = true; 
				 stopProtein = true;}
	}
seqLen =  i;
return (stopProtein);
}
// ---------------------------------------------------------------------------

 

