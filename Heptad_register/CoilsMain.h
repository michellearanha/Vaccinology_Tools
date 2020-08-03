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
#include "CoilsAlgo.h"
#include "read.seqfileC.h" //////////CCCCCCCCCC
#include "read.parfiles.h"
#include "write.files.h"
 
#define   kMinSeqCoils  29

/***********************************
* function-declarations  */

void CoilsMain(const char *scoreFile, const char *emissProbFile, char *seqProbFile, int matrixType);

/* function-declarations 
************************************/
