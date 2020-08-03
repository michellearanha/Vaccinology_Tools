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

#include "CoilsMain.h"

extern int	   Code[27], Decode[22];
extern bool	   verbose;	
extern bool FBDdetailsD;
extern bool FBDdetailsS, FBDdetailsC, FBDdetailsL;
extern bool bWriteSeq;

extern int 	 seqNb, seqLen;
extern int   seq[kMaxSeqLen+1];
extern bool  completedFile;
extern char	  SeqName[kMaxSeqName+2];
extern FILE *cFpWriteDom; 
extern FILE *cFpWriteCoP, *cFpWriteCoL; 
extern FILE *cFpReadSeq; 
extern double thresholds[6];
extern int    thresholdNb;

FILE *cFpWriteAlgo; 
bool writePSSMProfile;
int matrixNb;
/***********************************
* function-declarations  */

bool	ReadEmissProbForCoils (FILE *cFrPar,double logScore[][26]);
void  DoCoils (int seqLen, const int Seq[] ,const double logScore[][26], int matrixType) ;
bool	ReadScoresCoils (FILE *cFrTrans, double logScore[][26]);

// ---------------------------------------------------------------------------
//		 CoilsMain 					
// ---------------------------------------------------------------------------
void CoilsMain(const char *transProbFile, const char *emissProbFile, char *seqProbFile, int matrixType)
{
char *scoreFile;
double logScore[8][26];
int i, j, n, l;
FILE	*ReadFilePointer[4]; 
char 	*ReadFileName[11], *WriteFileName[11];
if (verbose) std::cout << " CoilsMain, matrixType = " << matrixType << std::endl;
if (verbose)   std::cout << " CoilsMain, matrixType = " << matrixType << std::endl;
if  (matrixType==1)   scoreFile = "Inputs/R5.MTK";
if  (matrixType==2)   scoreFile = "Inputs/R5.MTIDK";
if  (matrixType<3) 
	{
	if ( (ReadFilePointer[2] = fopen( scoreFile,"r" ))  == NULL )
		{ std::cout << "\n could not open the file " << scoreFile << "\n"; exit(1); }
	else {ReadScoresCoils( ReadFilePointer[2], logScore); 
		}  	
	if (verbose)  std::cout << "finished ReadScoresCoils" << " using the file " << scoreFile << std::endl;
	fclose( ReadFilePointer[2] );
	}
else
	{
	if  ((ReadFilePointer[2] = fopen(emissProbFile,"r" ))  == NULL )
		{ std::cerr << "\n could not open the file " << emissProbFile << "\n"; exit(1); }
	else {ReadEmissProbForCoils( ReadFilePointer[2], logScore);}  	
	if (verbose) std::cout << "\n finished ReadEmissProb \n";
	fclose( ReadFilePointer[2] );
		}

if  ( (cFpWriteAlgo = fopen("Outputs/Checks/PSSM","w" ))  == NULL )
	{ std::cerr << "\n could not open the file Outputs/Checks/PSSM for writing \n";}

if (FBDdetailsD)
	{
	 if  ( (cFpWriteDom = fopen("Outputs/DomainsPSSM","w" ))  == NULL )
		 { std::cerr << "\n could not open the file Outputs/DomainsPSSM for writing \n";}
	else   fprintf(  cFpWriteDom, " PREDICTED COILED-COIL DOMAINS: OVERVIEW\n\n");
	}

if (FBDdetailsS || FBDdetailsC)
	{writePSSMProfile = true;}
else {writePSSMProfile = false;}

if (writePSSMProfile)
	{
	if  ( (cFpWriteCoP = fopen("Outputs/CompactProfilePSSM","w" ))  == NULL )
		{ std::cerr << "\n could not open the file Outputs/ProbProfilePSSM for writing \n";}
	else   
		{
		fprintf(cFpWriteCoP, " COILED-COIL PROBABILITY PER RESIDUE, COMPACT REPRESENTATION\n\n");
		 if (matrixType > 2)  fprintf(  cFpWriteCoP, "This matrix is not calibrated, will use pseudo-probabilities (strictly monotone to raw score) \n\n"); 
		 else
	 		{
	 		fprintf(  cFpWriteCoP, "\nUsing a calibration from one of the implementations of the Lupas et al. method (approximately)");
	 		if (matrixType == 1)  fprintf(cFpWriteCoP, "\nUsing: cc: Normalpdf(1.63, 0.22, score)   others: Normalpdf(0.77, 0.20, score)");
	 		else  fprintf(cFpWriteCoP, "\nUsing: cc: Normalpdf((1.69, 0.18)   others: Normalpdf(0.80, 0.18)"); 
	       fprintf( cFpWriteCoP, "\nFor raw scores consult the file ProbListPSSM");
	      fprintf( cFpWriteCoP, "\nCutoff is %8.4f \n\n", kdefCutoff );
			}
		}
	}

if (FBDdetailsL)
	{
	if  ( (cFpWriteCoL = fopen("Outputs/ProbListPSSM","w" ))  == NULL )
		{ std::cerr << "\n could not open the file Outputs/ProbListPSSM for writing \n";}
	else   
		{
      	 fprintf(  cFpWriteCoL, " COILED-COIL PROBABILITY LIST PER RESIDUE\n\n");
		 if (matrixType > 2)  fprintf(  cFpWriteCoL, "These matrix is not calibrated, will use pseudo-probabilities (strictly monotone to raw score) \n\n"); 
		 else
	 		{
	 		fprintf(  cFpWriteCoL, "\nUsing a calibration from one of the implementations of the Lupas et al. method (approximately)");
	 		if (matrixType == 1)  fprintf(cFpWriteCoL, "\nUsing: cc: Normalpdf(1.63, 0.22, score)   others: Normalpdf(0.77, 0.20, score)");
	 		else  fprintf(cFpWriteCoL, "\nUsing: cc: Normalpdf((1.69, 0.18)   others: Normalpdf(0.80, 0.18) "); 
			fprintf(cFpWriteCoL, "\nFor raw scores consult the file ProbListPSSM"); 
	      	fprintf( cFpWriteCoL, "\nCutoff is %8.4f \n\n", kdefCutoff );
			}
		}
	}

if  ((cFpReadSeq= fopen( seqProbFile, "r" ))  == NULL )
	{ std::cerr << "\n could not open the file  " << seqProbFile << "\n"; 
		exit(1); }
	
if (verbose)  std::cout << "\n opened all Files" <<  std::endl;

completedFile = false; seqNb = 0;
while (! completedFile)
	{
	seqNb++;seqLen = 0;
	if ( ! (ReadSeq ( seqNb, seqLen, seq, completedFile, SeqName) )   )// if not completed 
		{
		if (verbose)  std::cout << "\n ReadSeq done" <<  std::endl;
		if (seqLen > 0)  { WriteWarning(seqNb,1);}	
		else	{seqNb--;}
		} 
	if (seqLen > 0)
		{ 
		if (verbose)   std::cout << "\n read seq-Nb " << seqNb << "  name = " << SeqName << std::endl;
		if (bWriteSeq) WriteSeqP (seqNb, seqLen, seq, SeqName);
		else WriteSeqIdP(seqNb, SeqName);
		if (verbose)   std::cout << "\n processing seq-Nb " << seqNb << "  name = " << SeqName << std::endl;
		if (kMinSeqCoils > 0)
			{DoCoils (seqLen, seq , logScore, matrixType); }
		if ( (seqNb % kperiodStdout1) == 0)	 
			{
			if (seqNb == kperiodStdout1) printf ("One point %d sequences, one line %d\n", kperiodStdout1,kperiodStdout2); 
			printf (".");  fflush (stdout);
			}
		if ( (seqNb % kperiodStdout2) == 0)  printf ("$\n");  	 
		if ( (seqNb % kperiodStdout3) == 0)  printf ( "\n processed seq-nb %d \n", seqNb);
		}
	}
	
if (verbose)    std::cout << "Coils file completed"  << std::endl;

if (writePSSMProfile)	fclose( cFpWriteCoP);
if (FBDdetailsD)  fclose( cFpWriteDom);
fclose( cFpWriteAlgo);
fclose( cFpReadSeq);

std::cout  <<  "\n processed  " << seqNb << "sequences" <<  std::endl;; 
}
/* 
************************************/
// ---------------------------------------------------------------------------
//		 DoCoils 					
// ---------------------------------------------------------------------------
void  DoCoils (int seqLen, const int Seq[] ,const double logScore [][26],int matrixType ) 
{
int 		i;
float    BgrProbProfile[seqLen+2];
double   PosStateProb[seqLen+2][kStates];

matrixNb = matrixType; 


if (verbose)   std::cout << "\n starting Coils-algo "  << std::endl;
CoilsAlgo (seqLen, Seq , logScore, PosStateProb , cFpWriteCoP, cFpWriteAlgo, cFpWriteCoL, matrixType ); 
if (verbose)   std::cout << "\n DoCoils after Coils-algo  " << seqNb << std::endl;
if (verbose)   std::cout << "\n WRITING PROB PROFILE   "  << std::endl;
		
for (i=1; i<= seqLen; i++)	BgrProbProfile[i] = (float) PosStateProb[i][0];
if (FBDdetailsD)   ParseIntoDomains( seqLen, BgrProbProfile, thresholdNb, thresholds );
}
/* DoCoils
************************************/
// ---------------------------------------------------------------------------
//		 ReadScoresCoils 					
// ---------------------------------------------------------------------------
bool	ReadScoresCoils (FILE *cFrPar, double  logScore[][26])
{
int  score[8][26];
int		state, aa, number;
bool	end;
char			 *tok;
char			buffer[90];
if (verbose) std::cout <<  "\n  READING COILS-MATRIX" << std::endl;

end = false;
for (state=0; state < 8;state++)
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
			if (number > 0)   logScore[state][aa] = log(((double) number / 1000)) ;
			else  		logScore[state][aa] = -2000 ;
			}
		}
		logScore[state][0] = 0 ;
	}
return (! end);
} 
/* ReadScoresCoils
************************************/
// ---------------------------------------------------------------------------
//		 ReadEmissProbForCoils 					
// ---------------------------------------------------------------------------

bool	ReadEmissProbForCoils (FILE *cFrPar, double logScore[][26])
{
int		k, state, aa, number, emProb[8][21];
bool		end;

if (verbose)  std::cout << "\n  ReadEmissProbForCoils ReadEmissProb \n" ;
end = false;

for (state=0; state < 8; state++)
     {
     if ( fscanf(cFrPar, "%d\n ", &number ) == EOF )  end = true; //read state-nb
     for (aa=1; aa <=20 ; aa++)
	     {
	     emProb[state][0] = 50000; //( equiv. to 1/20, aa 0 for nonstd aa)
	     if ( end || (fscanf(cFrPar, "%d\n ",&number) == EOF))  end = true; //read aa-nb
	     if ( end || (fscanf(cFrPar, "%d\n ",&number) == EOF)) 
		     {emProb[state][aa] = 0; end = true; std::cout << "\n end" ;}
	     else  {emProb[state][aa] = number;} 
     }	}
  			
if (!end)  
	{
	 if (verbose)  std::cout << "\n coils with parameters from emission matrix \n" ;
	for (k=1; k<=7; k++)
		{
		logScore[k][0] = 0; // nonstd aa equipared to background
		}
	 for (aa=1; aa<= 20; aa++)
		 {
		 for (k=1; k<=7; k++)
			 {
			 if  (log(emProb[k][aa]) > 0)    logScore[k][aa] = log(emProb[k][aa])  -  log(emProb[0][aa])  ;
			 else  logScore[state][aa] = -2000 ;
			 }
	}	 }
return (!end);
}
/* ReadEmissProbForCoils
************************************/
// ---------------------------------------------------------------------------

