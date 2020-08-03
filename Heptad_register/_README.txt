
***************  ********************************
**  README   **  ***   MARCOIL  version 1.2   ***
***************  ********************************
	

THIS IS NOT YET A PUBLIC RELEASE.

THE PROGRAM IS MADE DOWNLOADABLE HERE FOR GIVING 
REFEREES OF OUR MANUSCRIPT THE POSSIBILITY TO USE / 
CHECK THE PROGRAM /ALGORITHMS USED IN OUR STUDIES.
 
_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _

MARCOIL 
  is a hidden MARkov model for predicting existence and 
  localisation of coiled-coil domains.

	
  Questions and comments are welcome, email 
        delorenzi@wehi.edu.au
  Please report bugs.
_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _

INSTALLATION

  After downloading a copy of the file Marcoilcode.tar.gz in an
  appropriate directory, 
  extract with the commands
      gunzip Marcoilcode.tar.gz 
      tar Marcoilcode.tar  -xvf 
      
  This command creates a directory MARCOIL containing all the files;
  cd into the new directory and compile the code by typing make.
  The makefile calls the gcc compiler and the needed libraries.
  Various files should appear, including the program marcoil.
  With make clean you can remove the .o files, which are not further 
  needed after compilation.

_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _

TESTING

  From the MARCOIL directory, run the program on the file 
  SEQUENCES/seqfile2 by typing:
   
	marcoil +cdlsS SEQUENCES/seqfile
	
  and check out the appearance of four result files in the 
  subdirectory Outputs and their content.
   
  The 5 options  +cdlsS  in the above command do the following:
  c 	the output file   CompactProfile   	is produced
  d 	the output file   Domains   		is produced
  l 	the output file   ProbList	   	is produced
  s 	the output file   ProbPerState   	is produced 
  S 	the protein Sequence is also written to result files
  	
  Marcoil always rewrites the same output files, so to keep results 
  they have to be renamed or moved. 
  When the program has finished, a message is written to stdout, 
  as are some type of problems or errors, if encountered. 
_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _

RUNNING THE PROGRAM

  Marcoil takes its parameters from a bunch of files all located
  in the same subdirectory Inputs. This should not be modified.

  If you run the PSSM algorithm by using the -C option, then only
  the files of the options +s is inactive and the output 
  files are called CompactProfilePSSM, DomainsPSSM and ProbListPSSM.

  Running the program under UNIX or LINUX, follow the usual rules
  for adresses. For example,start the program from its directory with:
  
  ./ marcoil  <options>  <SequenceFileAdress>

  (you can omit  ./  if the directory was added to your path)

  If the name of the file for the sequences is omitted, the program
  looks for SEQUENCES/seqFile and tries to use this file. 
  The output should be self-explanatory. 
_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _

INPUT FILES

  Please use only fasta or multifasta format, as only these 
  were tested. The format is not checked by the program and you
  could get false results.
  The number of sequences is unlimited.
  Only the 20 standard amino acids are represented in the model,
  other letters are interpreted as unknown amino acids. 
  Lowercase is converted to uppercase, numbers, spaces,
  newlines, tabs and other special signs in the sequence are ignored.
  
  The program should be able to handle an input with the folling 
  schematic structure:

  > name of first sequence (name is optional)
  FIRST SEQUENCE  perhaps with digits
  61 continuation of sequence
  > name of following sequence
  FIRST SEQUENCE  perhaps with  s p a c e s  and newlines

  61 continuation of sequence

  (EOF)
_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _

USAGE & COMMAND-LINE OPTIONS

marcoil <-evwHLC<i>> <+cdlsS> <-E filename> <-T filename> 
	<-c number> <+c number> <filename.seq>
The order of the <...> blocks is irrelevant.
       
OPTIONS

  a) algorithmic / file / parameter  options

 [ -H ]  
 Read transition parameters from the file R3.transProbHigh 
 This is the default. We call it the MARCOIL-H predictor 
 [ -L ]  
 Read transition parameters from the file R3.transProbLow
 This is what we call the MARCOIL-L predictor. 
 [ -E emfileName]  
 Read emission parameters from the user's file emfileName
 [ -T trfileName]  
 Read transition parameters from the user's file trfileName
 [ -m trfileName]  
 Base emission parameters on the file R5.MTK or R5.MTIDK
 (if -i), except state 0, based on R2.emissProb (or the -E
 file)
 [ -i ]  
 Use the MTIDK matrix instead of MTK (active only under -P
 or -m)
 
 As file name the adress relative to the current directory
 has to be used.
 
 [ -P ]  
 Run the PSSM28 algorithm instead of MARCOIL, with the MTK matrix 
 (default, file R5.MTK) or the matrix specified with the -E option 
 The first is equivalent to COILS28 / MTK (up to rounding errors)

 
  b) numeric 
  
 [ -c probabilityCutoff]  
 Results are written only for positions with coiled-coil
 probability equal or superior to probabilityCutoff
 [ -t probabilityThreshold]  
 For the parsing, for defininf domain borders, use the
 threshold probabilityThreshold instead of the default
 values kParsingthreshold, defined in globals.h
 	 
  c) output options
  
 [ +d ] [ +l ] [ +m ] [ +s ] [ +S ] 
 Explained above.
 [ -v ]  
 Run in verbose mode. A fair bit of output will be printed to the 
 screen to follow execution. Useful for tracing down bugs.
  
  
************************************************************************
