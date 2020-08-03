
marcoil:  pgm.o wf.o init.o rparf.o  rseqf.o  FBmain.o FBa.o Coim.o Coia.o Pla.o
	g++ -O2 -o marcoil pgm.o wf.o init.o rparf.o  rseqf.o  FBmain.o FBa.o Coim.o Coia.o Pla.o  -lm -lc -lstdc++
	
pgm.o:  pgm.cpp pgm.h PlatformMain.h globals.h 
	g++  -O2 -c pgm.cpp 
	 
Pla.o: PlatformMain.h  PlatformMain.cpp  globals.h init.h  CoilsMain.h FBmain.h 
	g++  -O2 -c -o Pla.o PlatformMain.cpp

init.o: init.h  init.cpp globals.h 
	g++  -O2 -c  init.cpp

FBmain.o: FBmain.h  FBmain.cpp  FBalgorithm.h globals.h read.seqfileC.h read.parfiles.h write.files.h 
	g++  -O2 -c  FBmain.cpp

FBa.o: FBalgorithm.h  FBalgorithm.cpp  globals.h 
	g++  -O2 -c -o FBa.o FBalgorithm.cpp

rseqf.o: read.seqfileC.h  read.seqfileC.cpp globals.h 
	g++  -O2 -c -o rseqf.o read.seqfileC.cpp

rparf.o: read.parfiles.h  read.parfiles.cpp  globals.h
	g++  -O2 -c -o rparf.o read.parfiles.cpp

Coim.o: CoilsMain.h  CoilsMain.cpp  CoilsAlgo.h globals.h read.seqfileC.h read.parfiles.h write.files.h 
	g++  -O2 -c -o Coim.o CoilsMain.cpp

Coia.o: CoilsAlgo.h  CoilsAlgo.cpp  globals.h 
	g++  -O2 -c -o Coia.o CoilsAlgo.cpp

wf.o: write.files.h  write.files.cpp  globals.h 
	g++  -O2 -c -o wf.o write.files.cpp

clean:
	rm *.o



