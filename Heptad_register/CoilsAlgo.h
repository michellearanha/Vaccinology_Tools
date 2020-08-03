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

#include <cstdlib> // for access to the Unix-system
#include <stdio.h>
#include <stdlib.h> 

#include <iostream> // for  cin, cout, <<  etc
#include <fstream> // for file streams like fout
#include <iomanip>  // for  endl;

#include <math.h>
#include <float.h>
#include <string.h>

#include "globals.h"

void	CoilsAlgo (int totalNumber, const int Seq[] , const double score [][26], double PosStateProb[][kStates],
		 FILE * cFpWriteCoP,   FILE * cFpWrite, FILE * cFpWriteCoL, int matrixType );

