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
#include "FBalgorithm.h"

extern bool	   verbose;	
extern  int	   Code[27], Decode[22];
extern bool    writeControls;
extern  int    emProb[kStates][21], trProb[kStates][kStates] ;
extern  int    CutoffForConnecting;

static int ForwConnFirst[kStates], BackConnFirst[kStates];   
static int ForwConnLast[kStates], BackConnLast[kStates];   
static int ForwConnList[kStates * kStates], BackConnList[kStates * kStates];

void    SetUpConnections(void);
void    SetUpScaledParameters( long double scEmProb[ kStates][ 21], long  double scTrProb[ kStates][ kStates]);
void    InitialiseFBvars(void);

/********************************************************************************************************** */
// ---------------------------------------------------------------------------
//	 Forward-backward algorithm
// --------------------------------------------------------- ------------------
void	FBalgorithm (int seqNb, int totalNumber, const int Seq[] , double PosStateProb[][kStates] ) 
{
static int	st, pos, s, i, n, scalePeriod, FBFactor[kLength+2];
static long double probSeq,nb, factor, Rfactor; 
static bool FBdebugging;
static FILE *   cFpWrite;
static long  double  scForw[kStates][kLength+2];
static  long  double  scBackw[kStates][kLength+2] ;
static long  double  scEmProb[ kStates][ 21], scTrProb[ kStates][ kStates];
static long double   intermediate[kStates]; 

if (verbose)  std::cout << "\n SetUpConnections "  << std::endl;
 if (seqNb == 1)   SetUpConnections();
if (verbose)  std::cout << "\n SetUpScaledParameters "  << std::endl;
 if (seqNb == 1)   SetUpScaledParameters(scEmProb, scTrProb);

FBdebugging = false;
factor = exp(100);scalePeriod = 50;
for (i=0; i<kLength+2; i++)		 {FBFactor[i] = 0;}
if (writeControls || FBdebugging) 
	{
	if  ( (cFpWrite = fopen("Outputs/Checks/FBalgo.check","w" ))  == NULL )
		{ std::cerr << "\n could not open the file Outputs/Checks/FBalgo.check for writing \n";}
	}
		
if (verbose)  std::cout << "\n FB computation  kStates = " <<  kStates  << std::endl;
if (writeControls) fprintf(cFpWrite,"\n\n\n  FB computation  kStates = %d  \n\n\n", kStates );
								
scForw[0][0] = 1.0;
for (st=1; st<kStates; st++)		scForw[st][0] = 0.0;
if (writeControls)
	{
	for (st=0; st<kStates; st++)	fprintf(cFpWrite,"\n scForw[st %d][pos 0] = %6.15f; ", st,  (double) scForw[st][0] );
	}

for (pos=1; pos<=totalNumber; pos++)		
	{
	for (st=0; st<kStates; st++)	
		{
		scForw[st][pos] = 0;
		for (n=BackConnFirst[st]; n <= BackConnLast[st] ; n++)
			{
			s = BackConnList[n];
			scForw[st][pos] += (scForw[s][pos-1] *  scTrProb[s][st] ) ;
			}
		scForw[st][pos] = scForw[st][pos] * scEmProb[st][Seq[pos]] ; 
		if (FBdebugging  )
			{
			fprintf(  cFpWrite, "\n  scForw[st %d][pos %d ]   = %6.15f;  ", st,  pos,  (double) scForw[st][pos] );
		}	}
		//   modified Rabiner-Type overflow-protection / SCALING
	if (writeControls)  fprintf(  cFpWrite, "\n scForw[0][pos] = %f  at pos %d  \n" , (double) scForw[0][pos] , pos );
	if (   (pos % scalePeriod)  == 0)	//  done only one in scalePeriod aa  
		{
		nb = scForw[0][pos];  
		for (st=1; st<kStates; st++)	{ if  ( scForw[st][pos] >   nb )    nb = scForw[st][pos];}
		if ( log (nb) > 100 )
			{for (st=0; st<kStates; st++)	scForw[st][pos]  = scForw[st][pos]  /  factor ;
			FBFactor[pos] = 1;	
			if (writeControls)  fprintf(  cFpWrite, "\n scFBFactor at pos %d  \n", pos );
	}	}	}

for (st=0; st<kStates; st++)  scForw[st][totalNumber+1] = 0.0;
for (n=BackConnFirst[0]; n <= BackConnLast[0] ; n++)
		{
		s = BackConnList[n];
		scForw[0][totalNumber+1] += (scForw[s][totalNumber] *  scTrProb[s][0] ) ;
		}
for (st=1; st<kStates; st++)	scForw[st][totalNumber+1] = 0.0 ;
if (FBdebugging)
	{
	for (st=0; st<kStates; st++)	
	fprintf(  cFpWrite, "\n  scForw[st %d][totalNumber+1 = %d]   = %6.15f;  ", st, totalNumber+1,  (double) scForw[st][totalNumber+1] );
	fprintf(  cFpWrite, "\n\n");
	}
					// BACKWARD
Rfactor = 1 / factor ;
for (s=1; s<kStates; s++)	scBackw[s][totalNumber+1] = 0.0 ;
scBackw[0][totalNumber+1] = 1.0;
if (FBdebugging)
	{
	for (st=0; st<kStates; st++)	
	fprintf(  cFpWrite, "\n  scBackw[st %d][totalNumber+1 = %d]   = %6.15f;  ", st, totalNumber+1,  (double) scBackw[st][totalNumber+1] );
	}

for (s=0; s<kStates; s++) {scBackw[s][totalNumber] = scTrProb[s][0]; 
					intermediate[s] = scBackw[s][totalNumber] *scEmProb[s][Seq[totalNumber]];}
if (FBdebugging)
	{
	for (st=0; st<kStates; st++)	
	fprintf(  cFpWrite, "\n  scBackw[st %d][%d ]   = %6.15f;  ", st,totalNumber , (double)  scBackw[st][totalNumber] );
	}
for (pos=totalNumber-1; pos>=1; pos--)		
	{
	
	if (  FBFactor[pos+1] == 1  )	
		{
		for (st=0; st<kStates; st++)	
			{
			scBackw[st][pos] = 0.0 ;
			for (n=ForwConnFirst[st]; n <= ForwConnLast[st] ; n++)
				{
				s = ForwConnList[n];
				scBackw[st][pos] += scTrProb[st][s] * intermediate[s] * Rfactor ;
		}	}	}
	else	
		{
		for (st=0; st<kStates; st++)	
			{
			scBackw[st][pos] = 0.0 ;
			for (n=ForwConnFirst[st]; n <= ForwConnLast[st] ; n++)
				{
				s = ForwConnList[n];
				scBackw[st][pos] += scTrProb[st][s] * intermediate[s];
				}
		}	}	
	for (st=0; st<kStates; st++)   intermediate[st] = scBackw[st][pos] * scEmProb[st][Seq[pos]];
	}			    

scBackw[0][0] = scTrProb[0][0] * intermediate[0]; 
for (s=1; s<kStates; s++)	
	{scBackw[s][0] = 0.0;  scBackw[0][0] += scTrProb[0][s] * intermediate[s];}  

if (writeControls)
      {
	fprintf(  cFpWrite, "\n  scBackw[0][0]  = %6.15f;  ;" ,   (double) scBackw[0][0]) ;
	fprintf(  cFpWrite, "\n  scForw[0][totalNumber+1=%d] = %6.15f;  . \n\n\n", totalNumber+1,  (double) scForw[0][totalNumber+1] ); 
	}
nb = 	(scBackw[0][0] - scForw[0][totalNumber+1]) / scBackw[0][0];
if (nb < 0.0)  nb = - nb; 
if (nb > 00001)  std::cout  << "\n Failed ForwBackw-PrecisionCheck" << std::endl;  	
	
for (pos=0; pos<=totalNumber+1; pos++)
     {
     probSeq = 0 ;
     for (st=0; st<kStates; st++)
	     {
	     PosStateProb[pos][st] = (double)  (scForw[st][pos] * scBackw[st][pos] / scBackw[0][0] );
 	     probSeq += PosStateProb[pos][st] ;
	     if (FBdebugging)    
	    		{
 			fprintf(  cFpWrite, "\n  PosStateProb[st][pos] st = %d pos %d = %6.15f;  ",st,  pos,  (double) PosStateProb[st][pos] );
	     }	}
     if (FBdebugging)   fprintf(  cFpWrite, "\n  probSeq at pos %d = %6.20f;  ", pos,  (double) probSeq);
     nb = probSeq - 1.0;
     if (nb < 0.0)  nb = - nb; 
     if (nb > 0.00001)   std::cout  << "\n Failed probSeq=1 PrecisionCheck" << std::endl;  	
     }
if (writeControls || FBdebugging)     fclose( cFpWrite);
if (verbose)  std::cout << "\n FB computation concluded"  << std::endl;
}
/* FBalgorithm
************************************/
// ---------------------------------------------------------------------------
//	 SetUpConnections for Marcoil			
// ---------------------------------------------------------------------------
void    SetUpConnections(void)
{
int		t, s, n;
bool	      first, checkconnections;

checkconnections = false;
if (CutoffForConnecting > 0)  checkconnections = true;

if (verbose)  std::cout << "\n start SetUpConnections " << std::endl;
n=-1;
for (s=0; s < kStates; s++)
	{
	first=true;
	for (t=0; t < kStates ; t++)
		{
		if (trProb[s][t] > kMinProb)
			{
			ForwConnList[++n] = t;
			ForwConnLast[s] = n;
			if (first==true) {first=false; ForwConnFirst[s] = n;}
	}	}	}

n=-1;
for (s=0; s < kStates; s++)
	{
	first=true;
	for (t=0; t < kStates ; t++)
		{
		if (trProb[t][s] > kMinProb)
			{
			BackConnList[++n] = t;
			BackConnLast[s] = n;
			if (first==true) {first=false; BackConnFirst[s] = n;}
	}	}	}
	
if (checkconnections)
	{
	for (s=0; s < kStates; s++)
		{
		for (n=ForwConnFirst[s]; n <= ForwConnLast[s] ; n++)
			{
			t = ForwConnList[n];
	      	//std::cout  << n << ": " << t << "  |  ";
		}	}	
	for (s=0; s < kStates; s++)
		{
		for (n=BackConnFirst[s]; n <= BackConnLast[s] ; n++)
			{
			t = BackConnList[n];
	      	std::cout  << n << ": " << t << "  |  ";
			if ((n % 10) == 0 )   std::cout << "\n"  <<  std::endl;
	}	}	}	
}

/* SetUpConnections
************************************/
// ---------------------------------------------------------------------------
//	 SetUpScaledParameters for Marcoil			
// ---------------------------------------------------------------------------
void  SetUpScaledParameters(long  double scEmProb[kStates][21], long  double scTrProb[kStates][kStates])
{
int	i,j, aa;	
long double factor;
FILE *cFpWrite;


factor = 0.0;
for (aa=1; aa <=20 ; aa++)
	{
	 factor +=  (( (long double ) emProb[0][aa]) / 1000000) * (( ( long double) emProb[0][aa]) / 1000000);
	}
factor = ( long double) trProb[0][0] * factor / 100000000 ;
factor = 1.15 *  ( 1 /factor ) ; //1.15: empirical factor to avoid underflow, only overflow is catched
if (writeControls)
	{
	if  ( (cFpWrite = fopen("Outputs/Checks/SCALE.check","w" ))  == NULL )
	{ std::cerr << "\n could not open the file Outputs/Checks/SCALE.check for writing \n";}
	else	{fprintf(  cFpWrite, "\n factor  %f ", (double) factor  );}
	}

for (i=0; i<kStates; i++)
	{
	for (j=0; j<=20; j++)
		{
		scEmProb[i][j] =  factor *  emProb[i][j] / 1000000 ;
		} 
	for (j=0; j<kStates; j++)
		{
		scTrProb[i][j] = ( ( long double) 1.0 *  trProb[i][j] )  / 100000000;
		}
	}
}
/* SetUpScaledParameters
************************************/
// ---------------------------------------------------------------------------

