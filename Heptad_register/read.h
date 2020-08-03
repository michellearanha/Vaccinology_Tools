/*******************************************************************************
*									       *
*  file read.h  				       
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


#include <stdio.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include "basics.h"  
// #define	kNbDom	
using namespace std;


void  ReadFile(int n);


extern  long int  NegRes[][5][4]; // open dim only first possible, is kRepetNb+1
extern  long int  PosRes[][5][kPosTypes]; 
extern 	 int  DomRes[][5][kNbDom+2];
extern    double  AddRes[][5][4];



