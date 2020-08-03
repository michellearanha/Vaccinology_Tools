 /***********************************************************************
************************************************************************
**                                                                    **
** Marcoil                                                            **
** a program for predicting coiled-coil domains in protein sequences  **     
** Copyright (C) 2001 Mauro C. Delorenzi                              **
**                                                                    **
** This program is free software; you can redistribute it and/or      **
** modify it under the terms of the GNU General Public License        **
** as published by the Free Software Foundation; either version 2     **
** of the License, or (at your option) any later version.             **
**                                                                    **
** This program is distributed in the hope that it will be useful,    **
** but WITHOUT ANY WARRANTY; without even the implied warranty of     **
** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the      **
** GNU General Public License for more details.                       **
**                                                                    **
** You should have received a copy of the GNU General Public License  **
** along with this program; if not, write to the Free Software        **
** Foundation, Inc.,                                                  **
** 59 Temple Place - Suite 330,                                       **
** Boston, MA  02111-1307, USA.                                       **
**                                                                    **
************************************************************************
***********************************************************************/

#include "Marcoil.h"

// #include <cstlib> cstlib: No such file or directory
// -lc   in  "g++ -o inverse pgm.o utility.o -lm -lc"  helps out


/***********************************
* function-declarations  */

 void  AnalyseRes (int  n);
 void  SwapData (void);
 void  AnalyseResDiff (int  n);
 
 
/* function-declarations 
************************************/
/***********************************
* functions-definitions  */


//Option processArgs(int argc, char *argv[]);
//***********************************
main(int argc, char *argv[])
{

int gotInFile=0;

// READ OPTIONS, FILENAMES, SUBPROGRAM THAT IS REQUIRED
// SET FLAGS, VALUES, CHECK-OPEN FILES, SHOW MENU, RUN SUBPROGRAMS

for (i = 1; i < argc; i++)
  {
   if (!strcmp(argv[i], "-penalty")) 
  	{
        if (!sscanf(argv[++i],"%d", &parameters->penalty)) errorMessage();
        penalty_flag = 1;
        }
        
  thisarg = argv[i];
  if(thisarg[0] == '-')
  	{
        if(0==strcmp(thisarg+1,"acc")) 
        	{
        	i++;
        	if(i >= argc) 
        		{
          		fprintf(stderr,"Error: must provide acceptor data file with option -acc\n");
          		exit(1);
          		}
          	else
          	strcpy(accFile,argv[i]);
          	}
        }
   else if (0==strcmp(thisarg+1,"geo")) 
   	{
        opt.geometric = 1;
    	}
    else 
    	{
        fprintf(stderr,"Error: unknown option %s\n",thisarg);
        exit(1);
      }       
    } else {
      if(!gotInFile) {
        strcpy(opt.inFile,thisarg);
        gotInFile = 1;
      } else {
        fprintf(stderr,"Warning: ignoring command line option %s\n",thisarg);
      }
  }

 /* make sure we got a file to process */
  if(!gotInFile) 
  	{
    fprintf(stderr,"Error: no sequence file was specified\n");
    exit(1);
  	}



bool  DiffAnalysis = true;
bool  SimpleAnalysis = true;

int i;
char filename[20] = "file";
char c;
ifstream xx;

 cout  <<  "\n ********** starting  at:  " << endl ; 
 system ("date"); // calling a UNIX cmd
 
//cout   <<   "\n main start test-insert" <<  endl ; 
//OpenReadFile (xx );  // REFUSES TO PASS IFSTREAM OBJECT AS ARGUMENT ????!!!!!!  
//for (i=0; i<+10; i++)
//	{
//	xx.get(c);
//	cout << c << endl;
//	}
//cout   <<   "\n test-insert end" <<  endl ; 
  
  Readcvtest(); // READ INPUT: "parameters" from the crossvalidation 
  			// and the sequence sets  ,  fct decl in  sublist.h
  
  ReadFile(1);  // READ INPUT: RAW DATA FROM MARCOIL  ,  fct decl in  read.h
  
  if (DiffAnalysis)
  	{
	AnalyseRes(1);
	SwapData(); //keep into DomResM, PosResM, PosDomM, what needed for comparison case by case
 	ReadFile(2); // second set of results (file "W13.Res2")  IN STD VARS, FIRST SET (file "W13.Res1")  COPIED INTO "M-VARS"
	AnalyseRes(2);  // AnalyseRes2(2);
	AnalyseResDiff(2);
 	ReadFile(3); // thirdset of results (file "W13.Res3")  IN STD VARS, FIRST SET (file "W13.Res1")  still in "M-VARS"
	AnalyseRes(3);  // AnalyseRes2(2);
	AnalyseResDiff(3);
	}
  else 
  	{
  	if (SimpleAnalysis)
 		 {
		 cout   <<   "\n calling  AnalyseRes(1)" <<  endl; 
		 AnalyseRes(1);
		 }
	}
 
cout  <<  "\n\n       *********  finishing  at:  " << endl ; 
system ("date"); 
cout  <<  "************8********** " << endl << endl; 
 
return 0;
}


/* MAIN
************************************/
/***********************************
*   AnalyseRes  */


void  AnalyseRes (int  n)
{
  int repet, level, type, k, gr, dom, llimit, hlimit, cl;
  int nbd[12], nbGroups, Lengthclass;
  int domsInTest[kRepetNb+2][kNbGroupsIndexes];
  double diff;
  double diffAa;
  double TypeSum [kPosTypes][kNbGroupsIndexes];  
  double TypeSumAa [4];  
  double TypeSumSq [kPosTypes][kNbGroupsIndexes]; // ii index 6 for types  but  12 for groups  AND length-classes
  double TypeSumSqAa[4] ; // ii index 6 for types  but  12 for groups  AND length-classes
  double TypeSumDiff [kPosTypes][kNbGroupsIndexes];
  double TypeSumDiffAa[4]   ;
  double TypeSumF [kPosTypes][kNbGroupsIndexes];
  double TypeSumFSq [kPosTypes][kNbGroupsIndexes]; // ii index 6 for types  but  12 for groups   AND length-classes
  double TypeSumFDiff [kPosTypes][kNbGroupsIndexes];
  double NegtypeAvrg [kPosTypes][4];
  double NegtypeDev [kPosTypes][4];
  double inverse;
  int    DomRight[kNbDom+2][kNbGroupsIndexes];
  int    SurvClass[52][kNbLev][kNbGroupsIndexes];


  ofstream fout;
  ofstream foutL;
  ofstream foutdr;
  ofstream foutsp;
  ofstream foutaaR;
  
 
  
  if (n == 1) 
 	{
	fout.open ("ResFiles/Out1");
 	if (! fout) 
 		{
		cerr << "\n could not open Out1 \n" ; exit(1);
		}
  cout   <<   "\n starting  AnalyseRes, writing to Out1" <<  endl; 
	foutL.open ("ResFiles/OutLen1");
 	if (! foutL) 
 		{
		cerr << "\n could not open OutLen1 \n" ; exit(1);
		}
  cout   <<   "starting  AnalyseRes, writing also to OutLen1" <<  endl; 
	}
	
 if (n == 2) 
 	{
	fout.open ("ResFiles/Out2");
 	if (! fout) 
 		{
		cerr << "\n could not open Out2 \n" ; exit(1);
		}
  cout   <<   "\n starting  AnalyseRes, writing to Out2" <<  endl; 
	foutL.open ("ResFiles/OutLen2");
 	if (! foutL) 
 		{
		cerr << "\n could not open OutLen2 \n" ; exit(1);
		}
  cout   <<   "starting  AnalyseRes, writing also to OutLen2" <<  endl; 
	}
	
	
 if (n == 3) 
 	{
	fout.open ("ResFiles/Out3");
 	if (! fout) 
 		{
		cerr << "\n could not open Out3 \n" ; exit(1);
		}
  cout   <<   "\n starting  AnalyseRes, writing to Out3" <<  endl; 
	foutL.open ("ResFiles/OutLen3");
 	if (! foutL) 
 		{
		cerr << "\n could not open OutLen3 \n" ; exit(1);
		}
  cout   <<   "starting  AnalyseRes, writing also to OutLen3" <<  endl; 
	}

  ofstream foutR;
  
  if (n == 1) 
 	{
	foutR.open ("Rfiles/forR1");
 	if (! foutR) 
 		{
		cerr << "\n could not open forR1 \n" ; exit(1);
		}
  cout  << "writing matrix to forR1  " <<  endl; 
	}
 if (n == 2) 
 	{
	foutR.open ("Rfiles/forR2");
 	if (! foutR) 
 		{
		cerr << "\n could not open forR2 \n" ; exit(1);
		}
  cout << "writing matrix to forR2 " <<  endl; 
	}

 if (n == 3) 
 	{
	foutR.open ("Rfiles/forR3");
 	if (! foutR) 
 		{
		cerr << "\n could not open forR3 \n" ; exit(1);
		}
  cout << "writing matrix to forR3 " <<  endl; 
	}


	
	
//foutdr  <<   " DomRight-RESULTS  for  SURVIVAL-PLOTS  \n"  ; 

 nbGroups = 9; // sequences per group is the triple of nb
  
   	 // doms-END per group is nbd, DATA CREATED WITH MARCOIL, SEE FILE "seqDom" 
//cv1:  
/* nbd[0]  = 0;
  nbd[1]  = 78;	
  nbd[2]  = 128;	 //41*3=123 SEQDOM[123] = 128 
  nbd[3]  = 136;	 //44*3=132 SEQDOM[132] = 137
  nbd[4]  = 327;   //66*3=198 SEQDOM[ ] = 327
  nbd[5]  = 472;	 //73*3=219 SEQDOM[ ] =  472
  nbd[6]  = 580;	 //95*3=285 
  nbd[7]  = 644 ;	 //102*3=306 
  nbd[8]  = 734;    //112*3=336 
  nbd[9]  = 795;	 //125*3=375 
  nbd[10]  = 930;	 //170*3=510 
  nbd[11]  = 1047;	 //193*3=579 
  */

//cv2: 
 nbd[0]  = 0;
 nbd[1]  = 42;	
  nbd[2]  = 83;	 
  nbd[3]  = 285;	 
  nbd[4]  = 416;   
  nbd[5]  = 506;	 
  nbd[6]  = 565;	 
  nbd[7]  = 611 ;	 
  nbd[8]  = 689;    
  nbd[9]  = 820;	 

 	
  inverse = 1 / 150.0;
  fout.precision(6);
  fout.setf(ios::showpoint);
 
 fout   <<   "\n RESULTS ON THE NTS "; 
  for (type=0; type<= 3; type++)
 	{
	for (level=0; level<= 4; level++)
		{
		TypeSum [level][type] = 0;
		TypeSumSq [level][type] = 0;
		for (repet=1; repet<= kRepetNb; repet++)
			{
			TypeSum[level][type] += NegRes[repet][level][type] ;
			TypeSumSq[level][type] += NegRes[repet][level][type] * NegRes[repet][level][type] ;
			} // repet 
		} // level
	} // type
 for (type=0; type<= 3; type++)
 	{
	fout   <<   "\n NegtypeAvrg values for type = "   <<  type  << " "  ; 
	fout   <<   "\n level:  \t0\t\t1\t\t2\t\t3\t\t4\n\t"  ; 
	for (level=0; level<= 4; level++)
		{
		NegtypeAvrg[level][type] =  inverse * TypeSum[level][type]; 
		fout   << setw(10) << NegtypeAvrg [level][type] << "\t" ; 
		TypeSumDiff [level][type] = 0;  // SECOND METHOID FOR stddev since first gets sqrt(<0) sometimes
		for (repet=1; repet<= kRepetNb; repet++)
			{
			diff = NegRes[repet][level][type] - NegtypeAvrg[level][type]  ;
			TypeSumDiff [level][type] += diff * diff;  // sum square deviations
			} // repet 
		} // level
	fout   <<  "\n\t"; 
	for (level=0; level<= 4; level++)
		{
		NegtypeDev[level][type] = sqrt ( inverse * ( TypeSumSq[level][type] -
			(inverse * TypeSum[level][type]) * TypeSum[level][type] ) ); 
		fout  << setw(10) <<  NegtypeDev [level][type] << "\t" ; 
		} // level
	fout   <<  "\n\t"; 
	for (level=0; level<= 4; level++) // SECOND METHOID FOR stddev
		{
		NegtypeDev[level][type] = sqrt ( inverse * TypeSumDiff [level][type] ); 
		fout   << setw(10) <<  NegtypeDev [level][type] << "\t" ; 
		} // level
	fout   << "\n\n\n"; 
	} // type
	
	
 // cout   <<   "\n pos2 \n "; 
fout   <<   "\n RESULTS ON THE PTS, all sequences (LS+TS) "; 
 fout  <<   "\n\n  ADD.RES - STATISTICS \n"  ; 
 for (type=0; type<= 3; type++)
 	{
	for (level=0; level<= 4; level++)
		{
		TypeSum [level][type] = 0;
		TypeSumSq [level][type] = 0;
			{
			for (repet=1; repet<= kRepetNb; repet++)
				{
				TypeSum[level][type] += AddRes[repet][level][type] ;
				TypeSumSq[level][type] += AddRes[repet][level][type] * AddRes[repet][level][type] ;
				} // repet 
			} // else
		} // level
	} // type

  for (type=0; type<= 3; type++)
 	{
	fout   <<   "\n PostypeAvrg values for type = "   <<  type  << " "  ; 
	fout   <<   "\n level:   \t0\t\t1\t\t2\t\t3\t\t4\n\t"  ; 
	for (level=0; level<= 4; level++)
		{
		PostypeAvrg[level][type] =  inverse * TypeSum[level][type]; 
		fout << setw(10)   <<  PostypeAvrg [level][type] << "\t" ; 
		TypeSumDiff [level][type] = 0;
			{
			for (repet=1; repet<= kRepetNb; repet++)
				{
				diff = AddRes[repet][level][type] - PostypeAvrg[level][type]  ;
				TypeSumDiff [level][type] += diff * diff   ;
				} // repet 
			}//else
		} // level
	fout   <<   "\n\t"  ; 
	for (level=0; level<= 4; level++)
		{
		PostypeDev[level][type] =  sqrt  ( inverse * ( TypeSumSq[level][type] -
			(inverse * TypeSum[level][type]) * TypeSum[level][type] ) );
		fout  << setw(10)   <<  PostypeDev [level][type] << "\t" ; 
		} // level
	fout   <<  "\n\t"; 
	for (level=0; level<= 4; level++)
		{
		PostypeDev[level][type] = sqrt ( inverse * TypeSumDiff [level][type] ); 
		fout  << setw(10)   <<  PostypeDev [level][type] << "\t" ; 
		} // level
	fout   <<   "\n\n"  ; 
	} // type
	



 if (n == 1) 
 	{
	foutaaR.open ("Rfiles/forRaa1");
 	if (! fout) 
 		{
		cerr << "\n could not open forRaa1 \n" ; exit(1);
		}
  cout   <<   "AnalyseRes, writing to forRaa1" <<  endl; 
	}
 if (n == 2) 
 	{
	foutaaR.open ("Rfiles/forRaa2");
 	if (! fout) 
 		{
		cerr << "\n could not open forRaa2 \n" ; exit(1);
		}
  cout   << "AnalyseRes, writing to forRaa2" <<  endl; 
	}
 if (n == 3) 
 	{
	foutaaR.open ("Rfiles/forRaa3");
 	if (! fout) 
 		{
		cerr << "\n could not open forforRaa3 \n" ; exit(1);
		}
  cout   << "AnalyseRes, writing to forforRaa3" <<  endl; 
	}

// header for R  150x5  matrix (kRepetNb x kLevelNb)

 for (level=0; level<= 4; level++)
	{
	if (n ==1)	foutaaR << " CL" << level ;
	if (n ==2)	foutaaR << " HL" << level ;
	if (n ==3)	foutaaR << " LL" << level ;		 
	}
foutaaR  <<  "\n"  ; 

	
fout  <<   "\n\n  NUMBER OF TRUE POSITIVES \n"  ; 
// cout   <<   "\n pos3 \n "; 

  for (type=0; type < kPosTypes; type++)
 	{
	for (level=0; level<= 4; level++)
		{
		TypeSum [level][type] = 0;
		TypeSumSq [level][type] = 0;
		if (type == 5)
			{
			TypeSumAa[level] = 0;
			TypeSumSqAa [level] = 0;
			for (repet=1; repet<= kRepetNb; repet++)
				{
				PosResAa[repet][level] = ((double) PosRes[repet][level][4] ) / ( (double) PosRes[repet][level][5] );
				TypeSumAa[level]+= PosResAa[repet][level] ;
				TypeSumSqAa[level] += PosResAa[repet][level] * PosResAa[repet][level] ;
				} // changes values read from file, type 5 becomes the TP-aa-rate  
			}//if
		if (type == 2)
			{
			for (repet=1; repet<= kRepetNb; repet++)
				{
				TypeSum[level][type] += 0.001 * PosRes[repet][level][type] ;
				TypeSumSq[level][type] += 0.000001 * PosRes[repet][level][type] * PosRes[repet][level][type] ;
				} // repet 
			}//if
		else
			{
			for (repet=1; repet<= kRepetNb; repet++)
				{
				TypeSum[level][type] += PosRes[repet][level][type] ;
				TypeSumSq[level][type] += PosRes[repet][level][type] * PosRes[repet][level][type] ;
				} // repet 
			} // else
		} // level
	} // type
// cout   <<   "\n pos4  "  <<  endl ; 


foutaaR << "this is foutaaR file \n" ;

foutaaR.precision(4);
for (repet=1; repet<= kRepetNb; repet++)
	{
	for (level=0; level<= 4; level++)
		{
		foutaaR << PosResAa[repet][level]  << " " ;
		}
	foutaaR << "\n" ;
	}




  for (type=0; type < kPosTypes; type++)
 	{
	fout   <<   "\n PostypeAvrg values for type = "   <<  type  << " "  ; 
	fout   <<   "\n level:   \t0\t\t1\t\t2\t\t3\t\t4\n\t"  ; 
	for (level=0; level<= 4; level++)
		{
		PostypeAvrg[level][type] =  inverse * TypeSum[level][type] ; 
		fout << setw(10)   <<  PostypeAvrg [level][type] << "\t" ; 
		TypeSumDiff [level][type] = 0;
			{
			for (repet=1; repet<= kRepetNb; repet++)
				{
				if (type ==2)  {diff = 0.001 * PosRes[repet][level][type] - PostypeAvrg[level][type] ; }
				else {diff =  PosRes[repet][level][type] - PostypeAvrg[level][type]  ;}
				TypeSumDiff [level][type] += diff * diff   ;
				} // repet 
			}//else
		} // level
	fout   <<   "\n\t"  ; 
	for (level=0; level<= 4; level++)
		{
		PostypeDev[level][type] =  sqrt  ( inverse * ( TypeSumSq[level][type] -
			(inverse * TypeSum[level][type]) * TypeSum[level][type] ) );
		fout  << setw(10)   <<  PostypeDev [level][type] << "\t" ; 
		} // level
	fout   <<  "\n\t"; 
	for (level=0; level<= 4; level++)
		{
		PostypeDev[level][type] = sqrt ( inverse * TypeSumDiff [level][type] ); 
		fout  << setw(10)   <<  PostypeDev [level][type] << "\t" ; 
		} // level
	fout   <<   "\n\n"  ; 
	} // type

 type=6;
 	{
	fout   <<   "\n PostypeAvrg values for Aa TS , type = "   <<  type  << " "  ; 
	fout   <<   "\n level:   \t0\t\t1\t\t2\t\t3\t\t4\n\t"  ; 
	for (level=0; level<= 4; level++)
		{
		PostypeAvrgAa[level] =  inverse * TypeSumAa[level] ; 
		fout << setw(10)   <<  PostypeAvrgAa [level] << "\t" ; 
		TypeSumDiffAa [level] = 0;
			{
			for (repet=1; repet<= kRepetNb; repet++)
				{
				diffAa =  PosResAa[repet][level] - PostypeAvrgAa[level] ;
				TypeSumDiffAa [level] += diffAa * diffAa   ;
				} // repet 
			}//else
		} // level
	fout   <<   "\n\t"  ; 
	for (level=0; level<= 4; level++)
		{
		PostypeDevAa[level] =  sqrt  ( inverse * ( TypeSumSqAa[level] -
			(inverse * TypeSumAa[level]) * TypeSumAa[level] ) );
		fout  << setw(10)   <<  PostypeDevAa [level] << "\t" ; 
		} // level
	fout   <<  "\n\t"; 
	for (level=0; level<= 4; level++)
		{
		PostypeDevAa[level] = sqrt ( inverse * TypeSumDiffAa [level] ); 
		fout  << setw(10)   <<  PostypeDevAa [level] << "\t" ; 
		} // level
	fout   <<   "\n\n"  ; 
	} // type

	
fout  <<   "\n\n"  ; 
// cout   <<   "\n pos5  "  <<  endl ; 

fout  <<   " PosDom-RESULTS \n"  ; 

fout.precision(4);

fout  <<   "  **********   PER LENGTH CLASSES  ***********  \n"  ; 

		//  initialize
for (level=0; level<= 4; level++)
     {
     for (cl= 0 ; cl <= kNbLenCla; cl++)
	     {
		TypeSum[level][cl]  = 0;
		TypeSumSq[level][cl] = 0;
		for (repet=1; repet<= kRepetNb; repet++)  {PosLength [repet][level][cl] = 0;}
		}
	}
		// count true positives per length class   per level  at each repet in LS+TS
for (level=0; level<= 4; level++)
     {
     for (repet=1; repet<= kRepetNb; repet++)
	     {
	     for (dom= 1 ; dom <= kNbDom; dom++) // all doms 
		     {
		     Lengthclass = domLenCat[dom] ;
		     PosLength[repet][level][Lengthclass] += DomRes[repet][level][dom];
		     }
	     for (cl= 0 ; cl <= kNbLenCla; cl++) // all length classes
	     		{
	     		TypeSum[level][cl] += PosLength[repet][level][cl];
	     		TypeSumSq[level][cl] += PosLength[repet][level][cl] * PosLength[repet][level][cl];
			}
	     }// repet
     }// level

foutL << "this is foutL file \n" ;


for (level=0; level<= 4; level++)
	{
	for (cl= 0 ; cl <= kNbLenCla-1; cl++)
 		{
		if (n ==1)	foutL << " C" << level << "C" << cl;
		if (n ==2)	foutL << " H" << level << "C" << cl;
		if (n ==3)	foutL << " L" << level << "C" << cl;		 
		}
	}
foutL  <<  "\n"  ; 
foutL << "this is foutaaR file \n" ;


foutL.precision(4);
for (repet=1; repet<= kRepetNb; repet++)
	{
	for (level=0; level<= 4; level++)
		{
		for (cl= 0 ; cl <= kNbLenCla-1; cl++)
			{
			PosLengthR[repet][level][cl] = ( ( (float) PosLength[repet][level][cl] ) / ( (float)NbLenCat[cl] ) ) ;
			foutL << PosLengthR[repet][level][cl]  << " " ;
			}
		}
	foutL << "\n" ;
	}





		//  results :  averages and stddev
	fout   <<   "\nclass:\t   0\t   1\t   2\t   3\t   4\t   5\t   6\t   7\t   8\t   9\t   \nlevels\n"  ; 
for (level=0; level<= 4; level++)
	{
	fout.precision(4);
	
	fout  << "le"  <<  level  << "      " ; 
	for (cl= 0 ; cl <= kNbLenCla-1; cl++)
	    	{
		PosLengthAvg [level][cl] =  inverse * TypeSum[level][cl];
		fout  << setw(5)  <<  PosLengthAvg [level][cl] << "\t" ; 
		TypeSumDiff [level][cl] = 0;  // SECOND METHOID FOR stddev since first gets sqrt(<0) sometimes
		for (repet=1; repet<= kRepetNb; repet++)
			{
			diff = PosLength[repet][level][cl] - PosLengthAvg[level][cl]  ;
			TypeSumDiff [level][cl] += diff * diff;  // sum square deviations
			} // repet 
		} // classes
	fout   <<   "\n      "  ; 
	
	fout.precision(2);

	for (cl= 0 ; cl <= kNbLenCla-1; cl++)
	    	{
		PosLengthDev[level][cl] =  sqrt (inverse * TypeSumSq[level][cl] -
			PosLengthAvg[level][cl] * PosLengthAvg[level][cl] ); 
//		fout   << setw(5)   <<  PosLengthDev[level][cl] << "\t" ; 
		} 
//	fout   <<   "\n\t"  ; 
	for (cl= 0 ; cl <= kNbLenCla-1; cl++)
		{
		PosLengthDev[level][cl] = sqrt ( inverse * TypeSumDiff [level][cl] ); 
		fout   << setw(5) <<  PosLengthDev [level][cl] << "\t" ; 
		} 
	fout   <<   "\n\n"  ;
	} // levels

fout   <<   "Nbdoms:"  ;
for (cl= 0 ; cl <= kNbLenCla-1; cl++)
		{
		fout   <<  NbLenCat[cl]  << "\t" ; 
		} 


fout   <<   "\n\n\n\n"  ;
fout.precision(4);
fout.setf(ios::showpoint);


fout  <<   " PosDom-RESULTS \n"  ; 


// per doms:  per group and later restricted to those not in LS 
  for (gr=1; gr<= nbGroups; gr++)   // groups 1-11
 	{
	hlimit = nbd[gr]; 
	llimit = nbd[gr-1] + 1;
	for (level=0; level<= 4; level++)
		{
		TypeSum[level][gr]  = 0;
		TypeSumSq[level][gr] = 0;
		for (repet=1; repet<= 150; repet++)
			{
			PosDom [repet][level][gr] = 0;
			for (dom= llimit ; dom <= hlimit; dom++) // doms in that group
				{
				PosDom [repet][level][gr] += DomRes[repet][level][dom];
				}
			TypeSum[level][gr] += PosDom [repet][level][gr];
			TypeSumSq[level][gr] += PosDom [repet][level][gr] * PosDom [repet][level][gr];
			}// repet
		}// level
	}// group
  
  
  		//   GROUP 0  :  TOTALS 
 for (level=0; level<= 4; level++)// sums as gr 0 should give same as above but different when outside LS
	{
	TypeSum[level][0]  = 0;
	TypeSumSq[level][0] = 0;
	for (repet=1; repet<= 150; repet++)
		{
		PosDom [repet][level][0] = 0;
		for (gr= 1 ; gr <= nbGroups; gr++) // doms in that group
			{
			PosDom [repet][level][0] += PosDom [repet][level][gr];  // VALUES FOR SUM = GROUP 0
			}
		TypeSum[level][0] += PosDom [repet][level][0];
		TypeSumSq[level][0] += PosDom [repet][level][0] * PosDom [repet][level][0];
		}// repet
 	}// level

 for (gr=0; gr<= nbGroups; gr++)
 	{
	if (gr == 0)
		{
 		fout   <<   "\n PosDomAvrg values for group = "   <<  gr  << " with nb doms =  " << nbd[nbGroups] ; 
		}
	else
		{
 		fout   <<   "\n PosDomAvrg values for group = "   <<  gr  << " with nb doms =  " << (nbd[gr] - nbd[gr-1]) ; 
		}
	fout   <<   "\n level:   \t0\t\t1\t\t2\t\t3\t\t4\n\t"  ; 
	 for (level=0; level<= 4; level++)
		{
		PosDomAvrg[level][gr] =  inverse * TypeSum[level][gr]; 
		fout  << setw(10)  <<  PosDomAvrg [level][gr] << "\t" ; 

		TypeSumDiff [level][gr] = 0;  // SECOND METHOID FOR stddev since first gets sqrt(<0) sometimes
		for (repet=1; repet<= kRepetNb; repet++)
			{
			diff = PosDom[repet][level][gr] - PosDomAvrg[level][gr]  ;
			TypeSumDiff [level][gr] += diff * diff;  // sum square deviations
			} // repet 
		} // level
	fout   <<   "\n\t"  ; 
	for (level=0; level<= 4; level++)
		{
		PosDomDev[level][gr] =  sqrt (inverse * TypeSumSq[level][gr] -
			PosDomAvrg[level][gr] * PosDomAvrg[level][gr] ); 
		fout   << setw(10)   <<  PosDomDev [level][gr] << "\t" ; 
		} // level
	fout   <<   "\n\t"  ; 
	for (level=0; level<= 4; level++) // SECOND METHOD FOR stddev
		{
		PosDomDev[level][gr] = sqrt ( inverse * TypeSumDiff [level][gr] ); 
		fout   << setw(10) <<  PosDomDev [level][gr] << "\t" ; 
		} // level
	fout   <<   "\n\n"  ;
	} // group
	




fout  <<   "\n\n"  ; 
fout  <<   " PosDom-RESULTS limited to doms excluded from the LS \n"  ; 

// per doms: LIMITED TO THOSE NOT IN THE LEARNING SET,  have cvtest = 1
//  cvTest[kNbDom+2][kRepetNb+2];


  for (gr=1; gr<= nbGroups; gr++)
 	{
	hlimit = nbd[gr]; 
	llimit = nbd[gr-1] + 1;
	for (level=0; level<= 4; level++)
		{
		TypeSum[level][gr]  = 0;
		TypeSumSq[level][gr] = 0;
		TypeSumF[level][gr]  = 0;
		TypeSumFSq[level][gr] = 0;
		for (repet=1; repet<= 150; repet++)
			{
			if (level == 0)  
				{
				domsInTest[repet][gr] = 0;
				//if (gr == 1) domsInTest[repet][0] = 0;
				}

			PosDom [repet][level][gr] = 0;
			for (dom= llimit ; dom <= hlimit; dom++) // doms in that group
				{
				if (cvTest[dom][repet] == 1) // doms in testset 
					{
					PosDom [repet][level][gr] += DomRes[repet][level][dom];
					if (level == 0) 
						{
						domsInTest[repet][gr] ++;
						//domsInTest[repet][0] ++;
						}
					}
				}
			TypeSum[level][gr] += PosDom [repet][level][gr];
			TypeSumSq[level][gr] += PosDom [repet][level][gr] * PosDom [repet][level][gr];
			PosDomFrac[repet][level][gr]  =  (double) PosDom [repet][level][gr]  / ((double)  domsInTest[repet][gr]);
			TypeSumF[level][gr] += PosDomFrac [repet][level][gr];
			TypeSumFSq[level][gr] += PosDomFrac [repet][level][gr] * PosDomFrac [repet][level][gr];
			}// repet
		}// level
	}// group

 for (level=0; level<= 4; level++)
	{
	TypeSumF[level][0]  = 0;
	TypeSumFSq[level][0] = 0;
	TypeSum[level][0]  = 0;
	TypeSumSq[level][0] = 0;
	for (repet=1; repet<= 150; repet++)
		{
		if (level == 0) {domsInTest[repet][0] = 0;}
		PosDom [repet][level][0] = 0;
		for (gr= 1 ; gr <= nbGroups; gr++) // doms in that group
			{
			PosDom [repet][level][0] += PosDom [repet][level][gr];  // VALUES FOR SUM = GROUP 0
			if (level == 0) domsInTest[repet][0] += domsInTest[repet][gr];
			}
		TypeSum[level][0] += PosDom [repet][level][0];
		TypeSumSq[level][0] += PosDom [repet][level][0] * PosDom [repet][level][0];
		PosDomFrac [repet][level][0] = (double) PosDom [repet][level][0] / ((double)  domsInTest[repet][0]);  
		TypeSumF[level][0] += PosDomFrac [repet][level][0];
		TypeSumFSq[level][0] += PosDomFrac [repet][level][0] * PosDomFrac [repet][level][0];
		}// repet
 	}// level
 
 
 
 
 // write out for import into R, matrix 150x65
 
 
 for (level=0; level<= 4; level++)
	{
	if (n == 1)
 	for (gr=0; gr<= nbGroups; gr++) 	foutR << "l" << level << "g" << gr << "C  "; 
	if (n == 2)
 	for (gr=0; gr<= nbGroups; gr++) 	foutR << "l" << level << "g" << gr << "H  "; 
	if (n == 3)
 	for (gr=0; gr<= nbGroups; gr++) 	foutR << "l" << level << "g" << gr << "L  "; 
	}

foutR <<  "\n  " ;

 for (repet=1; repet<= 150; repet++)
	{
	for (level=0; level<= 4; level++)
		{
 		for (gr=0; gr<= nbGroups; gr++) 	foutR << PosDomFrac [repet][level][gr] <<  "  " ; 
		}
 	foutR << "\n  "; 
	}
 
 
 for (gr=0; gr<= nbGroups; gr++)
 	{
	if (gr == 0)
		{
 		fout   <<   "\n PosDomAvrgTEST values for group = "   <<  gr  << " with nb doms =  " << nbd[nbGroups] ; 
		}
	else
		{
 		fout   <<   "\n PosDomAvrgTEST values for group = "   <<  gr  << " with nb doms =  " << (nbd[gr] - nbd[gr-1]) ; 
		}
	fout   <<   "\n level:   \t0\t\t1\t\t2\t\t3\t\t4\n\t"  ; 
	 for (level=0; level<= 4; level++)
		{
		PosDomAvrg[level][gr] =  inverse * TypeSum[level][gr]; 
		fout  << setw(10)  <<  PosDomAvrg [level][gr] << "\t" ; 

		TypeSumDiff [level][gr] = 0;  // SECOND METHOD FOR stddev since first gets sqrt(<0) sometimes
		for (repet=1; repet<= kRepetNb; repet++)
			{
			diff = PosDom[repet][level][gr] - PosDomAvrg[level][gr]  ;
			TypeSumDiff [level][gr] += diff * diff;  // sum square deviations
			} // repet 
		} // level
	fout   <<   "\n\t"  ; 
	for (level=0; level<= 4; level++)
		{
		PosDomDev[level][gr] =  sqrt (inverse * TypeSumSq[level][gr] -
			PosDomAvrg[level][gr] * PosDomAvrg[level][gr] ); 
		fout   << setw(10)   <<  PosDomDev [level][gr] << "\t" ; 
		} // level
	fout   <<   "\n\t"  ; 
	for (level=0; level<= 4; level++) // SECOND METHOD FOR stddev
		{
		PosDomDev[level][gr] = sqrt ( inverse * TypeSumDiff [level][gr] ); 
		fout   << setw(10) <<  PosDomDev [level][gr] << "\t" ; 
		} // level
	fout   <<   "\n\n\t"  ;
	
				//  USINF FRACTIONAL NUMBER (SENSITIVITY)
	for (level=0; level<= 4; level++)
		{
		PosDomFracAvrg[level][gr] =  inverse * TypeSumF[level][gr]; 
		fout  << setw(10)  <<  (100 * PosDomFracAvrg [level][gr]) << "%\t" ; 

		TypeSumDiff [level][gr] = 0;  // SECOND METHOD FOR stddev since first gets sqrt(<0) sometimes
		for (repet=1; repet<= kRepetNb; repet++)
			{
			diff = PosDomFrac[repet][level][gr] - PosDomFracAvrg[level][gr]  ;
			TypeSumDiff [level][gr] += diff * diff;  // sum square deviations
			} // repet 
		} // level
	fout   <<   "\n\t"  ; 
	for (level=0; level<= 4; level++)
		{
		PosDomFracDev[level][gr] =  sqrt (inverse * TypeSumFSq[level][gr] -
			PosDomFracAvrg[level][gr] * PosDomFracAvrg[level][gr] ); 
		fout   << setw(10)   <<  (100 * PosDomFracDev [level][gr]) << "%\t" ; 
		} // level
	fout   <<   "\n\t"  ; 
	for (level=0; level<= 4; level++) // SECOND METHOD FOR stddev
		{
		PosDomFracDev[level][gr] = sqrt ( inverse * TypeSumDiff [level][gr] ); 
		fout   << setw(10) <<  (100 * PosDomFracDev [level][gr]) << "%\t" ; 
		} // level
	fout   <<   "\n\n"  ;
	} // group



// each dom over repets, for SURVIVAL-PLOTS and SINGLE-DOM profile

 if (n == 1) 
 	{
	foutdr.open ("Rfiles/forRdr1");
 	if (! fout) 
 		{
		cerr << "\n could not open forRdr1 \n" ; exit(1);
		}
  cout   <<   "AnalyseRes, writing to forRdr1" <<  endl; 
	}
 if (n == 2) 
 	{
	foutdr.open ("Rfiles/forRdr2");
 	if (! fout) 
 		{
		cerr << "\n could not open forRdr2 \n" ; exit(1);
		}
  cout   << "AnalyseRes, writing to forRdr2" <<  endl; 
	}
 if (n == 3) 
 	{
	foutdr.open ("Rfiles/forRdr3");
 	if (! fout) 
 		{
		cerr << "\n could not open forRdr3 \n" ; exit(1);
		}
  cout   << "AnalyseRes, writing to forRdr3" <<  endl; 
	}
//foutdr  <<   " DomRight-RESULTS  for  SURVIVAL-PLOTS  \n"  ; 
for (dom= 1 ; dom <= kNbDom ; dom++)  
		{
		if (n ==1)	foutdr << " dC" << dom;
		if (n ==2)	foutdr << " dH" << dom;
		if (n ==3)	foutdr << " dL" << dom;		 
		}
foutdr  <<  "\n"  ; 


for (k=0; k<= 50; k++)
 	{
 	for (level=0; level<= 4; level++)
		{
		for (gr=0; gr<= nbGroups; gr++)
 			{
			SurvClass[k][level][gr] = 0;
			}
		}
	} 

for (level=0; level<= 4; level++)
	{
	gr = 1;	
	for (dom= 1 ; dom <= kNbDom ; dom++) 
		{
		if   (dom ==  nbd[gr]+1 )		{gr ++;}
 		DomRight [dom][level] = 0;
 		for (repet=1; repet<= kRepetNb; repet++)
			{
			if (cvTest[dom][repet] == 1) // doms in testset 
				{
				DomRight [dom][level] += DomRes[repet][level][dom];
				}
			}// repet
	//	foutdr  <<   " DomRight[dom=" << dom  << "][" << level << "] = "  <<  DomRight [dom][level] << "\n"  ; 
		k = (int) DomRight [dom][level] ; // 0..50
		SurvClass[k][level][0] ++;  // collects the frequency of cases with k TP
		SurvClass[k][level][gr] ++;
		}// DOM
	}// level
	


 for (level=0; level<= 4; level++) 
	{
	for  (dom= 1 ; dom <= kNbDom ; dom++)
 		{
		//foutdr  <<   " DomRight[dom=" << dom  << "][" << level << "] = "  <<  DomRight [dom][level] << "\n"  ; 
		foutdr << DomRight [dom][level] << " "  ; 
		}
	foutdr  <<  "\n"  ; 
	}

if (n == 1) 
 	{
	foutsp.open ("Rfiles/forRsp1");
 	if (! fout) 
 		{
		cerr << "\n could not open forRsp1 \n" ; exit(1);
		}
  cout   << "AnalyseRes, writing to forRsp1" <<  endl; 
	}
 if (n == 2) 
 	{
	foutsp.open ("Rfiles/forRsp2");
 	if (! fout) 
 		{
		cerr << "\n could not open forRsp2 \n" ; exit(1);
		}
  cout   <<  "AnalyseRes, writing to forRsp2" << endl; 
	}
 if (n == 3) 
 	{
	foutsp.open ("Rfiles/forRsp3");
 	if (! fout) 
 		{
		cerr << "\n could not open forRsp3 \n" ; exit(1);
		}
  cout   << "AnalyseRes, writing to forRsp3" << endl; 
	}

// SurvClass[k][level][gr] go over to cdf
//foutsp  <<   " SurvClass-RESULTS  for  SURVIVAL-PLOTS  \n"  ; 
 for (level=0; level<= 4; level++)
	{
	for (gr=0; gr<= nbGroups; gr++)
 		{
		if (n ==1)	foutsp << " CL" << level << "G" << gr;
		if (n ==2)	foutsp << " HL" << level << "G" << gr;
		if (n ==3)	foutsp << " LL" << level << "G" << gr;		 
		}
	}
foutsp  <<  "\n"  ; 

k=50;
 	for (level=0; level<= 4; level++)
		{
		for (gr=0; gr<= nbGroups; gr++)
 			{
			foutsp  << SurvClass[k][level][gr] << " ";
			}
		}
	foutsp  <<  "\n"  ; 

for (k=49; k>= 0; k--)
 	{
 	for (level=0; level<= 4; level++)
		{
		for (gr=0; gr<= nbGroups; gr++)
 			{
			SurvClass[k][level][gr] += SurvClass[k+1][level][gr];  // cdf: k=49: cases with 49 or 50 TP
			foutsp  << SurvClass[k][level][gr] << " ";
			}
		}
	foutsp  <<  "\n"  ; 
	} 
		
cout  <<  "AnalyseRes  concluded " <<  endl ; 	
fout   <<  "AnalyseRes  concluded " <<  endl ;    
fout.close();
foutR.close();
foutL.close();
foutsp.close();
foutdr.close();
}



/* AnalyseRes
************************************/

 /***********************************
*   SwapData  */

void  SwapData (void)
{
int i,j,k, repet, level, cl;
 
for (j=0;j< kNbLev; j++)  
	{
	for (i=0;i<kRepetNb+2; i++)
		{
		for (k=0; k<kNbDom+2; k++)
			{
			DomResM[i][j][k] = DomRes[i][j][k];
			}
		for (k=0; k< kNbGroups; k++)
			{
			PosDomM[i][j][k] = PosDom[i][j][k];
			PosDomFracM[i][j][k] = PosDomFrac[i][j][k];
			}
		for (k=0; k<kPosTypes; k++)
			{
			PosResM[i][j][k] = PosRes[i][j][k];
			}
		}
	}

for (repet=1; repet<= kRepetNb; repet++)
	{
	for (level=0; level<= 4; level++)
		{
		for (cl= 0; cl<= kNbLenCla-1; cl++)
			{
			PosLengthM[repet][level][cl] = PosLengthR[repet][level][cl] ;
			}
		}
	}
}


/* SwapData
************************************/
/***********************************
*   AnalyseResDiff  */


//  DomResM[i][j][k] :  O OR 1 PER EACH  REPET - LEVEL - DOM , compare at single prediction res level
// can summarise at ggoup or overall level,  TS   or  TS+LS

//  PosDomM[i][j][k] :  TP-dom-number PER EACH  REPET - LEVEL - GROUP (TS ONLY)  , compare PAIRWISE DIFFERENCES AT GROUP level
//  create PosDomDiff[kRepetNb+2][5][12] and do stat as above

//  PosResM[i][j][k] :  summary stats (TP) PER EACH  REPET - LEVEL - 4 TYPES (TS AND LS)
// create PosResDiff[kRepetNb+2][5][4]   and do stat as above

 void  AnalyseResDiff (int  n)
 {
int 	    DRDoldTot[6][5]; 
int 	    DRDnewTot[6][5]; 
int       i,j,k,dom,level,repet, type, gr, d, sum, cl;
int       nbd[12], nbGroups;
int	    FDfreq[5][43];
float		code;

double TypeSumSq[kNbLev][12]; // ii index 6 for types  but  12 for groups
double TypeSum [kNbLev][12]; // ii index 6 for types  but  12 for groups
double TypeSumDiff [kNbLev][12];
double PosResDiffAvrg [kNbLev][kPosTypes]; //5 = nb levels
double PosResDiffDev [kNbLev][kPosTypes];
double PosDomDiffAvrg [kNbLev][12];
double PosDomDiffDev [kNbLev][12];
double diff;
double inverse;


bool   twodevmeth, writeLabel;
ofstream fout;

twodevmeth = false;  

if (n ==2)
	{
	fout.open ("ResFiles/OutDiffMH");
	if (! fout) 
 		{
		cerr << "\n could not open OutDiffMH \n" ; exit(1);
		}
	cout << "\nstarting  AnalyseResDiff, writing to OutDiffMH" <<  endl; 
	}
if (n ==3)
	{
	fout.open ("ResFiles/OutDiffML");
	if (! fout) 
 		{
		cerr << "\n could not open OutDiff \n" ; exit(1);
		}
	cout << "\nstarting  AnalyseResDiff, writing to OutDiffML" ; 
	}
//ofstream fout ("OutDiff");


fout.precision(4);
fout.setf(ios::showpoint);
cout.precision(3);
cout.setf(ios::showpoint);

for (level=0; level<5; level++)   { for (j=0; j<43; j++)        FDfreq[level][j] = 0;}

for (j=0; j<5; j++)
	{
	for (i=0; i<=kRepetNb; i++)
		{
		for (k=0; k<12; k++)
			{
			PosDomDiff[i][j][k] = PosDomM[i][j][k]  -  PosDom[i][j][k];
			PosDomFracDiff[i][j][k] = PosDomFracM[i][j][k]  -  PosDomFrac[i][j][k];
			if ( (k==0)   &&  (i!=0)  )
				{// frequency classes for PosDomDiff
				d = (int)  ( 200 *  PosDomFracDiff[i][j][0] +  5.5 );
				FDfreq[j][d] ++;
				//cout << setw(4) <<  d << j << k << i  << FDfreq[d] << endl;
				}
			}
		for (k=0; k<kPosTypes; k++)
			{
			PosResDiff[i][j][k] = PosResM[i][j][k] -  PosRes[i][j][k];
			}
		}
	}


for (level=0; level<5; level++)
	{
	sum = 0;
	fout  << "\n\n level = "  << level  << " \n    " ;
	fout  << "\n abs. Freq of PosDomFracDiff and  cumul. abs. Freq   in \n class \n    " ;
	for (j=0; j<26; j++)
		{
		fout  << setw(7)  << ( 0.5* (j - 5) )  << "%       " <<  FDfreq[level][j] ;
		sum += FDfreq[level][j];
		fout << setw(7)   << "          " <<  sum << "\n    ";
		}
	}


// PosResDiff[kRepetNb+2][5][4]   do stat as above

nbGroups = 11;

inverse = 1 / 150.0;
fout.precision(6);
fout.setf(ios::showpoint);
fout  <<   "\n\n  DIFFERENCES IN THE NUMBER OF TRUE POSITIVES \n"  ; 

  for (type=0; type< kPosTypes; type++)
 	{
	for (level=0; level<= 4; level++)
		{
		TypeSum [level][type] = 0;
		TypeSumSq [level][type] = 0;
		if (type == 2)
			{
			for (repet=1; repet<= kRepetNb; repet++)
				{
				TypeSum[level][type] += 0.001 * PosResDiff[repet][level][type] ;
				TypeSumSq[level][type] += 0.000001 * PosResDiff[repet][level][type] * PosResDiff[repet][level][type] ;
				} // repet 
			}//if
		else
			{
			for (repet=1; repet<= kRepetNb; repet++)
				{
				TypeSum[level][type] += PosResDiff[repet][level][type] ;
				TypeSumSq[level][type] += PosResDiff[repet][level][type] * PosResDiff[repet][level][type] ;
				} // repet 
			} // else
		} // level
	} // type

  for (type=0; type< kPosTypes; type++)
 	{
	fout   <<   "\n PostypeAvrg values for type = "   <<  type  << " "  ; 
	fout   <<   "\n level:   \t0\t\t1\t\t2\t\t3\t\t4\n\t"  ; 
	for (level=0; level<= 4; level++)
		{
		PosResDiffAvrg[level][type] =  inverse * TypeSum[level][type]; 
		fout << setw(10)   <<  PosResDiffAvrg [level][type] << "\t" ; 
		TypeSumDiff [level][type] = 0;
			{
			for (repet=1; repet<= kRepetNb; repet++)
				{
				diff = PosResDiff[repet][level][type] - PosResDiffAvrg[level][type]  ;
				if (type ==2)  diff = 0.001 * PosResDiff[repet][level][type] - PosResDiffAvrg[level][type]  ;
				TypeSumDiff [level][type] += diff * diff   ;
				} // repet 
			}//else
		} // level
	if (twodevmeth)    fout   <<   "\n\t"  ; 
	for (level=0; level<= 4; level++)
		{
		PosResDiffDev[level][type] =  sqrt  ( inverse * ( TypeSumSq[level][type] -
			(inverse * TypeSum[level][type]) * TypeSum[level][type] ) );
		if (twodevmeth)    fout  << setw(10)   <<  PosResDiffDev [level][type] << "\t" ; 
		} // level
	fout   <<  "\n\t"; 
	for (level=0; level<= 4; level++)
		{
		PosResDiffDev[level][type] = sqrt ( inverse * TypeSumDiff [level][type] ); 
		fout  << setw(10)   <<  ( 2.0  * PosResDiffDev [level][type]) << "\t" ; 
		} // level
	fout   <<   "\n\n"  ; 
	} // type
	
fout  <<   "\n\n\n"  ; 


//   PosDomDiff[kRepetNb+2][5][12] and do stat as above

 for (gr=0; gr<= nbGroups; gr++)
 	{
	for (level=0; level<= 4; level++)
		{
		TypeSum [level][gr] = 0;
		TypeSumSq [level][gr] = 0;
		for (repet=1; repet<= kRepetNb; repet++)
			{
			TypeSum[level][gr] += PosDomDiff[repet][level][gr] ;
			TypeSumSq[level][gr] += PosDomDiff[repet][level][gr] * PosDomDiff[repet][level][gr] ;
			} // repet 
		} // level
	} // type

 for (gr=0; gr<= nbGroups; gr++)
 	{
 	fout   <<   "\n PosDomDiffTEST values for group = "   <<  gr  << " "  ; 
	fout   <<   "\n level:   \t0\t\t1\t\t2\t\t3\t\t4\n\t"  ; 
	 for (level=0; level<= 4; level++)
		{
		PosDomDiffAvrg[level][gr] =  inverse * TypeSum[level][gr]; 
		fout  << setw(10)  <<  PosDomDiffAvrg [level][gr] << "\t" ; 

		TypeSumDiff [level][gr] = 0;  // SECOND METHOD FOR stddev since first gets sqrt(<0) sometimes
		for (repet=1; repet<= kRepetNb; repet++)
			{
			diff = PosDomDiff[repet][level][gr] - PosDomDiffAvrg[level][gr]  ;
			TypeSumDiff [level][gr] += diff * diff;  // sum square deviations
			} // repet 
		} // level
	if (twodevmeth)    fout   <<   "\n\t"  ; 
	for (level=0; level<= 4; level++)
		{
		PosDomDiffDev[level][gr] =  sqrt (inverse * TypeSumSq[level][gr] -
			PosDomDiffAvrg[level][gr] * PosDomDiffAvrg[level][gr] ); 
		if (twodevmeth)    fout   << setw(10)   <<  PosDomDiffDev [level][gr] << "\t" ; 
		} // level
	fout   <<   "\n\t"  ; 
	for (level=0; level<= 4; level++) // SECOND METHOD FOR stddev
		{
		PosDomDiffDev[level][gr] = sqrt ( inverse * TypeSumDiff [level][gr] ); 
		fout   << setw(10) <<  ( 2.0  * PosDomDiffDev [level][gr]) << "\t" ; 
		} // level
	fout   << "\n\n"; 
	} // group
fout   << "\n\n\n  SPECIFICITY  PLUS-MINUS FIXED, COMPARE SENSITIVITY IN % OR  A MEASURE OF EFFICIENCY "; 

//   use difference   PosDomFracDiff  =  PosDomFracM  -   PosDomFrac;

 for (gr=0; gr<= nbGroups; gr++)
 	{
	for (level=0; level<= 4; level++)
		{
		TypeSum [level][gr] = 0;
		TypeSumSq [level][gr] = 0;
		for (repet=1; repet<= kRepetNb; repet++)
			{
			TypeSum[level][gr] += PosDomFracDiff[repet][level][gr] ;
			TypeSumSq[level][gr] += PosDomFracDiff[repet][level][gr] * PosDomFracDiff[repet][level][gr] ;
			} // repet 
		} // level
	} // type

 for (gr=0; gr<= nbGroups; gr++)
 	{
 	fout   <<   "\n PosDomDiffTEST values for group = "   <<  gr  << " "  ; 
	fout   <<   "\n level:   \t0\t\t1\t\t2\t\t3\t\t4\n\t"  ; 
	 for (level=0; level<= 4; level++)
		{
		PosDomFracDiffAvrg[level][gr] =  inverse * TypeSum[level][gr]; 
		fout  << setw(10)  << (100.0 * PosDomFracDiffAvrg [level][gr]) << "%\t" ; // mmmmmmmmm

		TypeSumDiff [level][gr] = 0;  // SECOND METHOD FOR stddev since first gets sqrt(<0) sometimes
		for (repet=1; repet<= kRepetNb; repet++)
			{
			diff = PosDomFracDiff[repet][level][gr] - PosDomFracDiffAvrg[level][gr]  ;
			TypeSumDiff [level][gr] += diff * diff;  // sum square deviations
			} // repet 
		} // level
	if (twodevmeth)    fout   <<   "\n\t"  ; 
	for (level=0; level<= 4; level++)
		{
		PosDomFracDiffDev[level][gr] =  sqrt (inverse * TypeSumSq[level][gr] -
			PosDomFracDiffAvrg[level][gr] * PosDomFracDiffAvrg[level][gr] ); 
		if (twodevmeth)    fout   << setw(10)   <<  (100.0 *PosDomFracDiffDev [level][gr]) << "%\t" ; 
		} // level
	fout   <<   "\n\t"  ; 
	for (level=0; level<= 4; level++) // SECOND METHOD FOR stddev
		{
		PosDomFracDiffDev[level][gr] = sqrt ( inverse * TypeSumDiff [level][gr] ); 
		fout   << setw(10) <<  ( 2.0 * 100.0 * PosDomFracDiffDev [level][gr]) << "%\t" ; 
		} // level
	fout   << "\n\n"; 
	} // group
fout   << "\n\n\n"; 




//  DIFFERENCES AT THE SINGLE DOM - SINGLE REPET-PRED LEVEL   coded as +1 (pred by M only)  and  -1 (pred by C only)

 for (dom=1; dom<= kNbDom; dom++)
 	{
	for (level=0; level<= 4; level++)
		{
		for (k=0; k<6; k++)	{DomResDiff [k][level][dom] = 0;}//  I index: 0,3 diff=+1   1,4 diff=0    2,5 diff=-1  low:TS high: all doms
		for (repet=1; repet<= kRepetNb; repet++)
			{
			i = DomResM[repet][level][dom] - DomRes[repet][level][dom];
			DomResDiff [4-i][level][dom] ++;
			if (cvTest[dom][repet] == 1)  {DomResDiff [1-i][level][dom] ++;}
			} // repet 
		} // level
	} // dom
// print results


fout.precision(3);
fout.setf(ios::showpoint);

   
  	 // doms-END per group is nbd, DATA CREATED WITH MARCOIL, SEE FILE "seqDom" 
//cv1:  
/* nbd[0]  = 0;
  nbd[1]  = 78;	
  nbd[2]  = 128;	 //41*3=123 SEQDOM[123] = 128 
  nbd[3]  = 136;	 //44*3=132 SEQDOM[132] = 137
  nbd[4]  = 327;   //66*3=198 SEQDOM[ ] = 327
  nbd[5]  = 472;	 //73*3=219 SEQDOM[ ] =  472
  nbd[6]  = 580;	 //95*3=285 
  nbd[7]  = 644 ;	 //102*3=306 
  nbd[8]  = 734;    //112*3=336 
  nbd[9]  = 795;	 //125*3=375 
  nbd[10]  = 930;	 //170*3=510 
  nbd[11]  = 1047;	 //193*3=579 
  */

//cv2: 
 nbd[0]  = 0;
 nbd[1]  = 42;	
  nbd[2]  = 83;	 
  nbd[3]  = 285;	 
  nbd[4]  = 416;   
  nbd[5]  = 506;	 
  nbd[6]  = 565;	 
  nbd[7]  = 611 ;	 
  nbd[8]  = 689;    
  nbd[9]  = 820;	 
 
j = 1;
//fout  <<  "dom le0  +   =   -   +   =   -   le1  +   =   -   +   =   -   le2  +   =   -   +   =   -   le3";
//fout  <<  "  +   =   -   +  =   -   le4  +   =   -   +   =   - " << endl ;
fout  <<  "  dom    le0  +   -   le1  +   -   le2  +   -   le3";
fout  <<  "  +   -   le4  +   - " << endl ;

for (level=0; level<= 4; level++)
	{
	for (k=0; k<6; k++)
		{
		DRDnewTot[k][level] = 0;
		DRDoldTot[k][level] = 0;
		}
	}


 for (dom=1; dom<= kNbDom; dom++)
 	{
	d = 0; writeLabel = false;
	for (level=0; level<= 4; level++)   d += DomResDiff [0][level][dom] + DomResDiff [2][level][dom];
	if (d >= 20)  writeLabel = true;
	if (writeLabel)  fout  <<  setw(4)  <<  dom << "       ";
	for (level=0; level<= 4; level++)
		{
		for (k=0; k<6; k++)	//  (k=0; k<6; k++)
			{
			if ((k < 3) &&  (k != 1)  &&  writeLabel )  {fout  << setw(4)  <<  DomResDiff [k][level][dom];}
			DRDnewTot[k][level] += DomResDiff [k][level][dom] ;
			}//  I index: 0,3 diff=+1   1,4 diff=0    2,5 diff=-1  low:TS high: all doms
		if (writeLabel)   fout  <<   "     ";
		}
	if (writeLabel)    fout  << endl;
		// group j  summaries:
	if  ( dom  ==  nbd[j]  )      // end of group j  
		{
		fout  <<  setw(4)  <<  "\n gr "  << j << "      ";
		for (level=0; level<= 4; level++)
			{
			for (k=0; k<6; k++)	
				{
				if ((k < 3) &&  (k != 1))  fout  << setw(4)  <<  ( DRDnewTot[k][level] - DRDoldTot[k][level]  ) ;
			 	DRDoldTot[k][level] = DRDnewTot[k][level] ;
				}
			fout  <<   "     ";
			}
		fout  << "\n\n\n";
		j ++;
		}
	} // for dom

fout  <<   "\n tot       "; ;
for (level=0; level<= 4; level++)
	{
	for (k=0; k<3; k++)	
		{
		if ((k < 3) &&  (k != 1))    fout  << setw(5)  <<  DRDnewTot[k][level] ;
		}
	fout  <<   "   ";
	}
fout  << endl;

fout  << "\n\n\n          tot and all  \n\n"; ;
for (level=0; level<= 4; level++)
	{
	fout  << "level " << level << ":  " ;
	for (k=0; k<6; k++)	
		{
		fout  << setw(7)  <<  DRDnewTot[k][level]  <<  "  ";
		if (k == 2)   fout  << "         ";
		}
	fout  << endl;
	}
fout  << endl;

fout.close();


// PRODUCING TABLE OF DIFF IN LENGTH FOR R
if (n ==2) fout.open ("ResFiles/outLenDiffMH");
if (n ==3) fout.open ("ResFiles/outLenDiffML");
fout.precision(7);
fout.setf(ios::showpoint);


for (level=0; level<= 4; level++)
	{
	for (cl= 0 ; cl <= kNbLenCla-1; cl++)
 		{
		if (n ==2)	fout << " CH" << level << "C" << cl;
		if (n ==3)	fout << " CL" << level << "C" << cl;		 
		}
	}
fout  <<  "\n"  ; 


for (repet=1; repet<= kRepetNb; repet++)
	{
	for (level=0; level<= 4; level++)
		{
		for (cl= 0; cl<= kNbLenCla-1; cl++)
			{
			fout << (PosLengthM[repet][level][cl] -  PosLengthR[repet][level][cl])  << " " ;
			}
		}
	fout << "\n" ;	
	}	



fout.close();


// PRODUCING STRINGS TO BE IMPORTED INTO R AS LABELS
fout.open ("Rfiles/r.strings");
fout.precision(7);
fout.setf(ios::showpoint);

k = 1;
for (dom=1; dom<= kNbDom; dom++)
 	{
	j = domSeq[dom];
	i = dom - seqDom[ j - 1 ]  ;
	fout  << "dom." << dom << "~"  << j << "." << i << "/"  << k  << endl ;
	if  (dom == nbd[k]   )  k ++;
	//code = dom + ( (float) domSeq[dom] / 1000 );
	//fout << setw(10)  <<  code <<  endl ;
	}

fout.close();


//cout  <<   "\n AnalyseResDiff  concluded " <<  endl ; 	
//fout   <<   "\n AnalyseResDiff  concluded " <<  endl ;    

}










/*******************************************************************************

indeces for types of data from NTS:
0  Nb FP Proteins in the whole NTS
1  Nb FP Domains in the whole NTS
2  Nb FP Aa  in the whole NTS
3  Nb FP Protein Fragments in the whole NTS

indeces for types of data "Add::
0  Nb predicted
1  length error per end (1 dom 2 ends)
2  absolute error per end 
3  average length of pred  dom 

indeces for types of data from PTS:
0  Nb TP Proteins in the whole PTS
1  Nb TP Domains in the whole PTS
2  Nb TP Aa  in the whole PTS
3  Nb TP Protein Fragments in the whole PTS

4  Nb TP Aa over Testset
5  Nb cc Aa in Testset
//
 changes values read from file, type 5 becomes the TP-aa-rate									       *
cv2 protein families:
cv2: all have weight 1/9
GROUP1:  TROPOMYOSIN-GROUP	42	42			 
GROUP2:MYOSIN/PARAMYOSIN-GROUP  36 (78)				 
GROUP3: IF (KERATIN/LAMIN/DESMIN/VIMENTIN)	69 (147)		 
GROUP4:  DYNEIN-GROUP		18 (165)			 
GROUP5:  KINESIN-GROUP		54 (219)			 
GROUP6:  LAMININ(LIKE)-GROUP	18 (237)			 
GROUP7:  SNARE-GROUP		30 (267)			 			 
GROUP8:  TF (B-ZIP & MYC)	78 (345)			 
GROUP9:  MIXED POOL-GROUP	75 (420)			 

*
*******************************************************************************/










