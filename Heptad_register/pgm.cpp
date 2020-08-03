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

#include "pgm.h"  



// ---------------------------------------------------------------------------
//		 main 					
// ---------------------------------------------------------------------------

int main(int argc, char *argv[])
{
 char  * P;
int i, j, n, let, len;
char 	*ReadFileName[11], *WriteFileName[11];
char 	*thisarg, optSeqFile[101], optTransFile[101], optEmissFile[101];
char  thisletter;
bool  help;
bool  gotSeqFile, gotTransFile, gotEmissFile, optPSSM;
bool  optTransLow, optTransHigh;
 
 cout  <<  "\n ********** starting  at:  " << endl ; 
 system ("date"); 
 verbose = false;writeControls=false;
 gotSeqFile=false;gotTransFile=false;gotEmissFile=false;optCmatrix = false;
 optPSSM=false;optMTIK=false;optTransLow=false;optTransHigh=false;
 CutoffForWriting = -1.0;DomainThreshold=-1.0;
 FBDdetailsS=false; FBDdetailsC=false; FBDdetailsL=false; FBDdetailsD=false;

for (i = 1; i < argc; i++)  // index 0 forp gmNAme
  {
   char *optionFlags = "HLPimv";
   switch  (argv [i][0])
       {
       case  '-' :
          if (strpbrk ( argv [i] , optionFlags)!= NULL )
		{
		len = strlen(argv [i]);
		for (let = 1; let < len; let++) 
			{ 
			thisletter = argv [i][let]; // cout << "\n thisletter = " << thisletter;
			switch  (thisletter)
             		{	
                  	case  'v' :  verbose = true;      break;
                        case  'H' :  optTransHigh = true; break;
                        case  'L' :  optTransLow = true;  break;
                	      case  'P' :  optPSSM = true;     break;
                	      case  'm' :  optCmatrix = true;  cout << "\noptCmatrix = true ";  break;
                	      case  'i' :  optMTIK = true;    cout << "\nMTIDK matrix = true ";  break;  
                 	      default :    cout <<  "Unrecognized option " <<  argv [i] << endl; break;
		}	}	}
	   else {
	   	switch  (argv [i][1])
             	{
                   case  'E' : 
				     i++;
				     if(i >= argc)     {cout << "\nError:missing filename after option -E \n";exit(1);} 
				     else  {n = strlen(argv[i]);
            			     if (n > 100) {cout << "\n Warning: filename too long max = 100 char \n"; exit(1);}
            			     strncpy(optEmissFile,argv[i], 100);
            			     gotEmissFile = true;
					     }
                  break;
                  case  'T' : 
				     i++;
				     if(i >= argc)   {cout << "\nError:missing filename after option -T \n";exit(1);} 
				     else  {n = strlen(argv[i]);
            			     if (n > 100) {cout << "\n Warning: filename too long max = 100 char \n"; exit(1);}
            			     strncpy(optTransFile,argv[i], 100);
            			     gotTransFile = true;
					     }
                  break;
                  case  'c' :
					if  (argv [i][2] != '\0'){  P = argv [i] + 2; 
									} // if no space
                  		else	{
						i++;
				            if(i >= argc)   {cout << "\nError:missing number after option -c \n";exit(1);} 
  						else  P = argv [i];  
									} // if space
			 		CutoffForWriting = (double) atof (P);
                  break;
                  case  't' :
					if  (argv [i][2] != '\0'){  P = argv [i] + 2; 
									} // if no space
                  		else	{
						i++;
				            if(i >= argc)   {cout << "\nError:missing number after option -t \n";exit(1);} 
  						else  P = argv [i];  
									} // if space
			 		DomainThreshold = (double) atof (P);
                  break;
        	     default : cout <<  "Unrecognized option %s\n" <<  argv [i];
 		}	} 
	 break;
       case  '+' :
		    len = strlen(argv [i]);
		    for (let = 1; let < len; let++) 
			    { 
			    thisletter = argv [i][let];
			    switch  (thisletter)
             		    {	
 	            	    case  'd' :  FBDdetailsD = true;  cout << "\n +d => set write domain parsing" << endl; break;
                  	    case  'l' :  FBDdetailsL = true;  cout << "\n +l => set write list of coiled-coil probabilities" << endl; break;
                  	    case  'm' :  FBDdetailsC = true;  cout << "\n +m => set write compact list of probabilities" << endl; break;
                  	    case  's' :  FBDdetailsS = true;  cout << "\n +s => set write probabilities for each state (not for PSSM)" << endl; break;
                  	    case  'S' :  bWriteSeq = true;    cout << "\n +S => set write protein sequence  " << endl; break;
                  	    case  'c' : 
						if  (argv [i][2] != '\0'){  P = argv [i] + 2; 
									} // if no space
                  			else	{
							i++;
				            	if(i >= argc)   {cout << "\nError:missing number after option +c \n";exit(1);} 
  							else  P = argv [i];  
									} // if space
			 				CutoffForConnecting =  atoi (P);
                		  default :  cout <<  "Unrecognized option +" <<  thisletter << endl; break;
			    }	    }
        break;
	default :
      	  {
      	  if(!gotSeqFile) 
			{
			n = strlen(argv[i]);
       		if (n > 100) {cout << "\n Warning: filename too long max = 100 char \n"; exit(1);}
      		strncpy(optSeqFile,argv[i], 100);
      		gotSeqFile = true;
      		} 
	         else { cerr <<  "Unrecognized option: " <<  argv [i]  << endl ;}
   }	}	  }          


if ((CutoffForWriting >= 0) && (CutoffForWriting <= 1)  )   
	{ if (verbose) cout << "\n Set CutoffForWriting  to " <<  CutoffForWriting ; 
	}
else {CutoffForWriting = kdefCutoff; 
	if (verbose) cout << "\n Using default CutoffForWriting = " << CutoffForWriting ; 
	}

if(!gotSeqFile) 
  	{
	if (verbose)	 cout << "\nNo sequence file was specified, will use default: SEQUENCES/seqFile \n";
    ReadFileName[1] = "SEQUENCES/seqFile";
  	}
else {ReadFileName[1] = optSeqFile;  }	

if(!gotTransFile) 
  	{
	if (optTransHigh)  ReadFileName[3] = "Inputs/R3.transProbHigh";  
	else  
		{
		if (optTransLow)  ReadFileName[3] = "Inputs/R3.transProbLow";
		else {  
			if (verbose)	 cout << "\nNo file for transition Probabilities  was specified, will use default:Inputs/R3.transProb \n";
 			ReadFileName[3] = "Inputs/R3.transProb"; 
			} 
  	}	}
else {ReadFileName[3] = optTransFile;  }	

if(!gotEmissFile) 
  	{
		 if (verbose)	 
cout << "\nNo file for emission Probabilities was specified, will use default: Inputs/R2.emissProb \n" <<  endl;
		  ReadFileName[2] = "Inputs/R2.emissProb";  
  	}	
else {ReadFileName[2] = optEmissFile;  }		

if (DomainThreshold != -1.0) {thresholdNb = 1; thresholds[0]= DomainThreshold;}
else
	{
	 thresholdNb = 5;
	 thresholds[0]= kParsingthreshold1;
	 thresholds[1]= kParsingthreshold2;
	 thresholds[2]= kParsingthreshold3;
	 thresholds[3]= kParsingthreshold4;
	 thresholds[4]= kParsingthreshold5;
	}

if (! optPSSM)
	{
	PlatformMain(ReadFileName[3] , ReadFileName[2]  , ReadFileName[1], 0 );
	}
else {
      if (gotEmissFile)
		PlatformMain(ReadFileName[3] , ReadFileName[2]  , ReadFileName[1], 3 );
	else
		{	
 		if (optMTIK)
			PlatformMain(ReadFileName[3] , ReadFileName[2]  , ReadFileName[1], 2 );
		else 
			PlatformMain(ReadFileName[3] , ReadFileName[2]  , ReadFileName[1], 1 );//default:MTK or MTIK, otherwise matrix as entered with -E
	}	}

cout  <<  "\n       *********  finishing  at:  " << endl ; 
system ("date"); 
cout  <<  "************8********** " << endl << endl; 
return 0;
}
/* MAIN
************************************/
