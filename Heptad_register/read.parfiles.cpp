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

#include "read.parfiles.h"
/***********************************
* global vars   */

extern  bool verbose;
extern  bool optMTIK;
extern  int  emProb[kStates][21], trProb[kStates][kStates] ;
//extern  long double	   emissions[ kStates][ 21], transitions[ kStates][ kStates] ;
extern  int ForwConnFirst[kStates], BackConnFirst[kStates];   
extern  int ForwConnLast[kStates], BackConnLast[kStates];   
extern  int ForwConnList[kStates * kStates], BackConnList[kStates * kStates];

/* global vars 
************************************/

// ---------------------------------------------------------------------------
//		ReadTransProb 					
// ---------------------------------------------------------------------------
bool	ReadTransProb (FILE *cFrTrans)
{
int		s, t, number;
bool	      end;
if (verbose)  std::cout << "\n  starting ReadTransProb \n" ;

end = false;
for (s=0; s < kStates; s++)
	{
	if ( fscanf(cFrTrans, "%d\n ", &number ) == EOF )   	 end = true;//read state-nb

	for (t=0; t < kStates ; t++)
		{
		 if ( end || (fscanf(cFrTrans, "%d\n ",&number) == EOF)) 
			    { trProb[s][t] = 0 ; end = true;}
		 else	
		    {trProb[s][t] = number; 
	}	}   }

return (!end);
}
// ---------------------------------------------------------------------------
//		ReadEmissProb 					
// ---------------------------------------------------------------------------
bool	ReadEmissProb(FILE *cFrPar)
{
int		state, aa, number, upTo;
bool		end, tying;

if (verbose)  std::cout << "\n  starting ReadEmissProb \n" ;
end = false;tying = true;
upTo = kStates;
if (tying)  upTo = 8;
emProb[0][0]= 50000;// just to have a number (is 1/20), equalised to other states	
for (state=0; state < upTo; state++)
     {
     
     if ( fscanf(cFrPar, "%d\n ", &number ) == EOF )  end = true; //read state-nb
     for (aa=1; aa <=20 ; aa++)
	     {
	     emProb[state][0] = emProb[0][0];  //( equiv. to 1/20, aa 0 for nonstd aa)
	     if ( end || (fscanf(cFrPar, "%d\n ",&number) == EOF))  end = true; //read aa-nb
	     if ( end || (fscanf(cFrPar, "%d\n ",&number) == EOF)) 
		     {emProb[state][aa] = 0; end = true; std::cout << "\n end" ;}
	     else  {emProb[state][aa] = number;} 
     }	}
  			
if (tying)  
	{
	for (state=8; state < kStates;state++)
		{
		for (aa=0; aa <=20 ; aa++)  emProb[state][aa] = emProb[state-7][aa] ;
	}	}	

return (!end);
}
// ---------------------------------------------------------------------------
//		 ReadProp 					
// ---------------------------------------------------------------------------
bool	ReadProp ()
{
FILE * cFrPar;
int  score[8][26];
int  state, aa, number;
double sum, factor;
bool	end;
char	*tok;
char	buffer[90], * filename;

if   (optMTIK)  filename = "Inputs/R5.MTIDK";
else		    filename = "Inputs/R5.MTK";

if ( (cFrPar = fopen( filename,"r" ))  == NULL )
	{ std::cout << "\n could not open the file " << filename << "\n"; exit(1); }
else	{
	//if (verbose)  cout << "\n reading propensities from %s" << filename << endl;
	}
std::cout << "\n reading propensities from %s" << filename << std::endl;

end = false;
for (state=0; state < 8;state++)//1. read propensities
	{
	if ( fscanf(cFrPar, "%d\n ", &number ) == EOF )
		{ 
		 end = true;
		}
	if (  end || (fgets ( buffer, 80, cFrPar ) == NULL ) )
		{
		 tok  = "name"; 
		 end = true;
		 }
	else
		{
		tok = strtok   (buffer, " \t  ") ;
		}
	for (aa=1; aa <=20 ; aa++)
		{
		if ( end || (fscanf(cFrPar, "%d\n ",&number) == EOF)) 
			{
			 end = true; 
			 }
		if ( end || (fscanf(cFrPar, "%d\n ",&number) == EOF)) 
			{
			score[state][aa] = 0 ;
			 end = true;
			 }
		else
			{
			score[state][aa] = number ;
	}	}	}
			//2, derive em probs from propensities and bgr em probs, aa 0 remains as above
for (state=1; state < 8;state++)
	{
	sum = 0;
	for (aa=1; aa <=20 ; aa++)
		{
		sum += ((double) score[state][aa] / 1000) *  ((double) emProb[0][aa]) ;
		}
	factor = 1000000.0 / sum;//renormalise
	for (aa=1; aa <=20 ; aa++)
		{
		emProb[state][aa] = (int) floor( 0.5 + factor * ((double) score[state][aa] / 1000) *  ((double) emProb[0][aa]));
		}
	}
	
for (state=8; state < kStates;state++)	for (aa=0; aa <=20 ; aa++)  emProb[state][aa] = emProb[state-7][aa] ;
	
fclose( cFrPar);
return (! end);
} 
/* ReadProp
************************************/
// ---------------------------------------------------------------------------

// ---------------------------------------------------------------------------
//      ¥ >    WriteParametersToFile: EMISSION PROB			
// ---------------------------------------------------------------------------
void WriteEmissProb(FILE *cFwPar) 
{
int	state, aa;

fprintf(  cFwPar, "\n");
for (state=0; state < 8 ; state++)   // up to 7 for info, up to 65 to keep usual format
		{
		fprintf(cFwPar, "%3d\t\n", state );
		for (aa=1; aa <21 ; aa++)
			{
			 fprintf(cFwPar, "%2d  %6d\t",aa,emProb[state][aa]);
			 if (aa%5 ==0)  fprintf(cFwPar,"\n\n");
			}
		fprintf(cFwPar,"\n\n");
		}
fprintf(cFwPar,"\n");
}
// ---------------------------------------------------------------------------
//      ¥ >    WriteTransProbToFile								
// ---------------------------------------------------------------------------
void WriteTransProb (FILE *cFwTrans) 
{
int i,j;

fprintf(cFwTrans,"\n");
for (i=0; i < kStates ; i++)
	{
	fprintf(cFwTrans,"%6d \t",i);
	fprintf(cFwTrans,"\n");
	for (j=0; j < kStates ; j++)
		{
		fprintf (cFwTrans, "%10d",trProb[i][j]) ;
		if (  (j==0) || ((j % 7) == 0)    )
			fprintf(cFwTrans,"\n\n");
		}
	fprintf(cFwTrans,"\n\n");
	}
fprintf(cFwTrans,"\n");
}
// ---------------------------------------------------------------------------

