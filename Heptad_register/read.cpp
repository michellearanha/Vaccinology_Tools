/*******************************************************************************
*									       *
*  file read.cpp  				       
*									       *
*									       *
* Written by Mauro Delorenzi							       
*									       *
*									       *
* WEHI    						       
* FEBRUARY  2001								       
*									       *
* Reading Output from Marcoil and Computing Statistics on the Results
*
*******************************************************************************/


#include "read.h"



		
void  ReadFile (int n)
{
using namespace std;
char S[220];
char a, c;
int value;
int level, repet, dom, k, sum;
float fvalue;
bool unfinished;
ifstream fin;
 
 
ofstream fout ("checkerFiles/OutRead");
if (! fout) 
 	{
	cerr << "\n could not open OutRead \n" ; exit(1);
	}
	

  if (n == 1) 
 	{
	fin.open ("Helpers/W13.Res1");
 	if (! fin) 
 		{
		cerr << "\n could not open Helpers/W13.Res1 \n" ; exit(1);
		}
	cout   <<   "\n ReadFile Helpers/W13.Res1" <<  endl ; 
	}
 if (n == 2) 
 	{
	fin.open ("Helpers/W13.Res2");
 	if (! fin) 
 		{
		cerr << "\n could not open Helpers/W13.Res2 \n" ; exit(1);
		}
	cout   <<   "\n ReadFile Helpers/W13.Res2" <<  endl ; 
	}
 if (n == 3) 
 	{
	fin.open ("Helpers/W13.Res3");
 	if (! fin) 
 		{
		cerr << "\n could not open Helpers/W13.Res3 \n" ; exit(1);
		}
	cout   <<   "\n ReadFile Helpers/W13.Res3" <<  endl ; 
	}

fout << "  File for checking of reading of input data \n\n"; // =>=>=>=>=>=>.
	
 repet = 1;
 level = -1;
 unfinished = true;

while ( unfinished  && ( fin.get(c) ) ) 
	{
	if (c == '<')  
		{
		fin.get(c);
		
		level ++;
		
		fin.get(a);
		fin.getline(S , 210);
		if (a == 'N')
			{
			if (level == 5) {level = 0; 
						repet ++;}
			for (k=0; k<kNbNegTypes; k++)
				{
				fin >> value ;
				NegRes[repet][level][k] = value;
				fout   <<  " " <<  NegRes[repet][level][k] ; // =>=>=>=>=>=>.
				}//for k
			fout << " neg  \n"; // =>=>=>=>=>=>.
			}//  a == N
		if (a == 'P')
			{
			if (level == 5) {level = 0;}

			for (k=0; k<kPosTypes; k++)
				{
				fin >> value ;
				PosRes[repet][level][k] = value;
				fout  <<  " "  <<  PosRes[repet][level][k] ; // =>=>=>=>=>=>.
				}//for k
			fout <<  endl; // =>=>=>=>=>=>.
			sum = 0;
			for (dom=1; dom<= kNbDom; dom++)
				{
				fin >> value ;
				DomRes[repet][level][dom] = value;
				if ( dom  <= 10) {fout  <<  DomRes[repet][level][dom] ;}
				sum += DomRes[repet][level][dom] ;
				//if ( (dom % 10)  == 0) {fout << "\n"; }
				}//for dom
			fout << " \n    sum DomRes  = " <<  sum << "\n"; // =>=>=>=>=>=>.
			for (k=0; k<kNbAddTypes; k++)
				{
				fin >> fvalue ;
				AddRes[repet][level][k] = fvalue;
				fout   <<  " "  <<  AddRes[repet][level][k] ;
				}//for k
			fout << " pos \n repet " << repet <<  "\n\n"; // =>=>=>=>=>=>.
			}//  a == P
		}// c == <

	if (repet > kRepetNb)  unfinished = false;
	} // while 

 fout << "\n end of read-loop"  << endl;  // =>=>=>=>=>=>.


//	while (fin.getline(S , 210) ) // stops if line is longer!!!
//		{
//		fout << S << endl; 
//		}
//while (fin >> value )
	//{
	//fout << value << endl; 
	//}

 
 fin.close();
 fout.close();
  
}



 

