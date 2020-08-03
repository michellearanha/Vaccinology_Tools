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
#include "CoilsAlgo.h"

extern  int Decode[22];
extern bool	verbose;	
extern bool FBDdetailsD;
extern bool writePSSMProfile;
extern double CutoffForWriting;
extern int   matrixNb;
extern bool FBDdetailsD;
extern bool FBDdetailsS, FBDdetailsC, FBDdetailsL;

int   Shiftdown ( int a );
int   Shiftup ( int a );
double Normalpdf  ( double m, double s, double x );
double CCProb ( double score );

void  WriteProbProfilePSSM (int seqLen,  const int seq[] , const double PosStateProb[][kStates], const int posWin[], 
const int posPh[], FILE * cFpWritePC ,FILE * cFpWriteCoL, const double RawScore[]  );

/********************************************************************************************************** */
// ---------------------------------------------------------------------------
//	 CoilsAlgo
// --------------------------------------------------------- ------------------
void	CoilsAlgo (int totalNumber, const int Seq[] , const double score [][26], double PosStateProb[][kStates],
		 FILE * cFpWriteCoP,   FILE * cFpWrite, FILE * cFpWriteCoL, int matrixType ) 
{
int	  i,j,k,m,t, ph[8], maxph, mode, maxWin, maxt; 
double  max, sco[8][29], RawScore[totalNumber+2] ;
static int	  posPh[kLength+2], posWin[kLength+2];
static float   posScore[kLength+2] ; 
bool debugging = false;
bool showRawScores = false;

// "LUPAS-ALGORITHMUS";  1.  BUILD-UP,  POS UP TO 28, FIRST COMPLETE ARRAY
for (k=1; k<= 7; k++)	 {sco[k][28] = 0;  ph[k]=k;}
for (i=0; i<= 27; i++)
	{
	for (k=1; k<= 7; k++)
		 {
		 sco[k][28] += score[ ph[k] ][Seq[28-i]];  // k refers to phase at pos 28
	       ph[k] = Shiftdown (ph[k]);
	}	}

max = -100;
for (k=1; k<= 7; k++)	 
	{
	if (sco[k][28]> max)   {max= sco[k][28];  maxph=k;}
	ph[k]=k;
	}
posScore[1] = max ;    posPh[1]= Shiftup(maxph); posWin[1]= 1; // phase at 1 = phase at 29
maxt= 0;
for (j=1; j<= 27; j++)
	{
	for (k=1; k<=7; k++)	
		{
		ph[k]=Shiftup(ph[k]);
		sco[k][28-j] = sco[k][29-j] + score[ ph[k] ][Seq[28+j]] - score[ ph[k] ][Seq[j]] ; 
		if (sco[k][28-j] > max)   {max= sco[k][28-j];  maxph=k;  maxt= j;}
		}
	posScore[1+j] = max;    posPh[1+j]=  ((maxph+j+1) % 7);    posWin[1+j]= 1+j-maxt; 
	if  ( posPh[1+j] == 0) 	posPh[1+j] = 7;
	}
	
/**/        //  2. SHIFTS TO NEXT POSITIONS
for (i=29; i<= totalNumber; i++)
	{
	for (k=2; k<=7; k++)	 {ph[k]=k-1;}
	ph[1]=7;
	max = -10000000; maxWin = -1;
	for (j=28; j>= 2; j--)  //  27x7 scores remain equal but change phase and type
		{
		for (k=1; k<=7; k++)	
			{
			sco[k][j] = sco[ph[k]][j-1] ; 
			if (sco[k][j] > max)   {max= sco[k][j];  maxph=k;  maxWin = j;}
		}	}
	if (i+27 <= totalNumber) // generate new scores for window over this + next 27 aa
		{
		for (k=1; k<= 7; k++)	
			{
			sco[k][1] = 0 ;  ph[k] = k; 
			for (j=0; j<= 27; j++)
				{
				sco[k][1] += score[ ph[k] ][Seq[i+j]] ;
				ph[k]=Shiftup(ph[k]);// for next aa
				} 
			}
		for (k=1; k<=7; k++)  
			{
			if (sco[k][1] > max)   {max= sco[k][1];  maxph=k; maxWin = 1;}
			}
		}
	else
		{ 
		for (k=1; k<=7; k++)	{sco[k][1] = -100 ; }
		}
	posScore[i] = max ;    posPh[i]= maxph;  posWin[i]=  maxWin;
	}

for (i=1; i<= totalNumber; i++)
	{
	RawScore[i] = exp ( (double) posScore[i] / 28) ;
	PosStateProb[i][0] =  1 - CCProb(exp ( (double) posScore[i] / 28)  );
	if (showRawScores)
		{
		fprintf(  cFpWrite, "\n  log-Score[i=%d] : %f ", i, ( (double) posScore[i] / 28)   );
		fprintf(  cFpWrite, "\t  Score : %f ",  RawScore[i]  );
		fprintf(  cFpWrite, "\t CCProb= : %f", 1 - PosStateProb[i][0]  );
		}
	}

if (verbose) std::cout << "\n Coils computation concluded"  << std::endl;
if (writePSSMProfile)  WriteProbProfilePSSM (totalNumber, Seq , PosStateProb, posWin, posPh, cFpWriteCoP, cFpWriteCoL, RawScore);
}
/* CoilsAlgo
************************************/
// ---------------------------------------------------------------------------
//		¥ 	Shiftdown\\<*/> 					
// ---------------------------------------------------------------------------
int   Shiftdown ( int a )
{
int b;
if (a == 1)  b = 7;
else		 b = a - 1;
return b;
}
/* Shiftdown
************************************/
// ---------------------------------------------------------------------------
//		¥ 	Shiftup\\<*/> 					
// ---------------------------------------------------------------------------
int   Shiftup ( int a )
{
int b;
if (a == 7)  b = 1;
else		 b = a + 1;
return b;
}
/* Shiftup
************************************/
// ---------------------------------------------------------------------------
//		¥ 	CCProb\\<*/> 					
// ---------------------------------------------------------------------------
double   CCProb ( double score )
{
double b, c, d, ratio;

//  different implementations of COILS give slightly different values
// because of different normal distributions average/var
// and maybe rounding error => I use double precision 
//paper  :   cc   1.63  0.24;    bgr  0.77  0.20. 
//pascal pgm : 1.628,0.243,0.770,0.202,
// unix / C version of coils:
// uw 28 1.63 0.22 0.77 0.20 30 in the file old.mat
// uw 28 1.69 0.18 0.80 0.18 30 in the file new.mat
//  MacStripe ??: 1.63, 0.22, ; 0.772, 0.2005, 
//  or 1.63, 0.24, 0.78, 0.20 ??
// unix version of coils:
ratio = 30;
if (matrixNb == 1)
	{
 	 c = Normalpdf(1.63, 0.22, score);  d = Normalpdf(0.77, 0.20, score);  
 	 }
else  // MTIDK and 'pseudoprobabilities"
	{
	 c = Normalpdf(1.69, 0.18, score);   d = Normalpdf(0.80, 0.18, score);
 	 }
b =  c / (ratio * d  +  c );
return b;
}
/* CCProb
************************************/
// ---------------------------------------------------------------------------
//		¥ 	Normalpdf\\<*/> 					
// ---------------------------------------------------------------------------
double   Normalpdf  ( double m, double s, double x )
#define 	kpi 	3.141592654
{
double  a, u;
u =  (x - m) / s;
a = exp ( (-0.5* u*u) )  / sqrt (2*s* kpi);
return  a;
}
#undef kpi
/* Normalpdf
************************************/
// ---------------------------------------------------------------------------
//		 WriteProbProfilePSSM 					
// ---------------------------------------------------------------------------
void  WriteProbProfilePSSM (int seqLen,  const int seq[] , const double PosStateProb[][kStates], const int posWin[], 
const int posPh[], FILE * cFpWritePC ,FILE * cFpWriteCoL, const double RawScore[] )
{
int    i, mode;
float  ProbForWritingPSSM;
bool   writing;
double Prob;
ProbForWritingPSSM = CutoffForWriting;

if (verbose) std::cout << "\n WriteProbProfilePSSM"   << ";   limit ="   << ProbForWritingPSSM  <<  std::endl;

if (FBDdetailsL)
	{ 
	fprintf(  cFpWriteCoL, "    raw scores, cc-probability in percent and best heptad phase\n");
	for (i=1; i<= seqLen; i++)
		{
		if ((  PosStateProb[i][0] <= (1- ProbForWritingPSSM))  || (i==1) || (i==seqLen))
			{
			fprintf(  cFpWriteCoL, "%4d %c  ",  i , Decode[seq[i]] ); 
			fprintf(  cFpWriteCoL, "%6.3f  ",  RawScore[i] ); 
			fprintf(  cFpWriteCoL, "%5.1f %c \n", 100.0 * ((double) 1.0 - PosStateProb[i][0]), (char) (posPh[i]+96)); 
		}	}
	fprintf(cFpWriteCoL,"\n\n******************************************************************************"); 
	fprintf(cFpWriteCoL,"**********************************************************************************\n"); 
	}

if  (FBDdetailsS ||  FBDdetailsC)
	{
	 fprintf( cFpWritePC, "\n coiled-coil probability in percent and heptad position with highest probability\n");
	 mode = 0;
	 for (i=1; i<= seqLen; i++)
		 {
		 Prob = 1 - PosStateProb[i][0]; 
		 if  (( Prob >=  ProbForWritingPSSM ) || (i==1) || (i==seqLen)) 
			 {
			 writing = true;
			 mode ++;
			 fprintf(  cFpWritePC, "%4d%c", i, Decode[seq[i]]);
			 if (Prob >= 0.9995) fprintf(  cFpWritePC, " 100 "); 
			 else fprintf(  cFpWritePC, " %4.1f", 100 * Prob ); 
			 fprintf(  cFpWritePC, " %2d", (posWin[i]) );
			 fprintf(  cFpWritePC, "%c ", (char) (posPh[i]+96) );
			 if (i==seqLen) fprintf(cFpWritePC, "  **");
			 else
				 {
				 if  (mode ==7)  {mode = 0; fprintf(cFpWritePC, "\n");}
				 else fprintf(cFpWritePC, "");
			 }	}
		 else  {
			 if (writing)  {writing = false; fprintf(cFpWritePC, "\n[..]\n"); mode=0;}
		 }	}
	//		 }
	//	 else  
	//		 {
	//		 if (writing)  {writing = false; fprintf(cFpWritePC, "\n[..]\n"); mode = 0; }
	//		 }
	//	 if  (mode ==7)  {mode = 0; fprintf(cFpWritePC, "\n");}
	//	 if (i==seqLen)   fprintf(cFpWritePC, " *");
	//	 }
	fprintf(cFpWritePC,"\n\n******************************************************************************"); 
	fprintf(cFpWritePC,"**********************************************************************************\n"); 
	}
}
/* WriteProbProfilePSSM
************************************/
// ---------------------------------------------------------------------------

