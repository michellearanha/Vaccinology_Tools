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

#include "init.h"

extern int Decode[22];
extern int Code[27];
extern bool  verbose;

// ---------------------------------------------------------------------------
//		 init.codes 					
// ---------------------------------------------------------------------------

void initCodes(void)
{
int i, j;

if (verbose)   cout << "\n start in initCodes"<< endl;
 
  Code[0] = 0;Code[1] = 1;// 1 for Ala A
  Code[2] = 0;Code[3] = 2;
  Code[4] = 3;Code[5] = 4;
  Code[6] = 5;Code[7] = 6;
  Code[8] = 7;Code[9] = 8;
  Code[10] = 0;Code[11] = 9;
  Code[12] = 10;Code[13] = 11;
  Code[14] = 12;Code[15] = 0;
  Code[16] = 13;Code[17] = 14;
  Code[18] = 15;Code[19] = 16;
  Code[20] = 17;Code[21] = 0;
  Code[22] = 18;Code[23] = 19;
  Code[24] = 0;Code[25] = 20;// 20 for Tyr Y
  Code[26] = 0;

for (i=1; i < 21 ; i++)	
	{
	for (j=1; j <= 26 ; j++)	
		{
		if (Code[j] == i)  Decode[i] = j + kAtoOne;
		}
	}
Decode[0] = 24 + kAtoOne;  Decode[21] = 24 + kAtoOne;

if (verbose)   cout << "\n finished in initCodes"<< endl;
}
/* init.codes
************************************/
