MCC = icpc
MGCC = icc
-include ../../../compiler.mk

MSTOOLKIT = ../../../src/MSToolkit
COMETSEARCHPATH = ../../../src/CometSearch
override CXXFLAGS += -O3 -g -Wall -Wextra -Wno-write-strings -static -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64 -DGCC -I. -I$(MSTOOLKIT)/include

COMETSEARCH = Threading.o CometInterfaces.o CometSearch.o CometPreprocess.o CometPostAnalysis.o CometMassSpecUtils.o CometWriteOut.o CometWriteSqt.o CometWritePepXML.o CometWritePercolator.o CometWriteTxt.o CometSearchManager.o

all:  $(COMETSEARCH)
	ar rcs libcometsearch.a $(COMETSEARCH)
	rm -f *.o

clean:
	rm -f *.o libcometsearch.a

Threading.o:          $(COMETSEARCHPATH)/Threading.cpp $(COMETSEARCHPATH)/Threading.h
	${MCC} ${CXXFLAGS} $(COMETSEARCHPATH)/Threading.cpp -c
CometSearch.o:        $(COMETSEARCHPATH)/CometSearch.cpp $(COMETSEARCHPATH)/Common.h $(COMETSEARCHPATH)/CometData.h $(COMETSEARCHPATH)/CometDataInternal.h $(COMETSEARCHPATH)/CometSearch.h $(COMETSEARCHPATH)/CometInterfaces.h $(COMETSEARCHPATH)/ThreadPool.h
	${MCC} ${CXXFLAGS} $(COMETSEARCHPATH)/CometSearch.cpp -c
CometPreprocess.o:    $(COMETSEARCHPATH)/CometPreprocess.cpp $(COMETSEARCHPATH)/Common.h $(COMETSEARCHPATH)/CometData.h $(COMETSEARCHPATH)/CometDataInternal.h $(COMETSEARCHPATH)/CometPreprocess.h $(COMETSEARCHPATH)/CometInterfaces.h $(MSTOOLKIT)
	${MCC} ${CXXFLAGS} $(COMETSEARCHPATH)/CometPreprocess.cpp -c
CometMassSpecUtils.o: $(COMETSEARCHPATH)/CometMassSpecUtils.cpp $(COMETSEARCHPATH)/Common.h $(COMETSEARCHPATH)/CometData.h $(COMETSEARCHPATH)/CometDataInternal.h $(COMETSEARCHPATH)/CometMassSpecUtils.h $(COMETSEARCHPATH)/CometInterfaces.h
	${MCC} ${CXXFLAGS} $(COMETSEARCHPATH)/CometMassSpecUtils.cpp -c
CometPostAnalysis.o:  $(COMETSEARCHPATH)/CometPostAnalysis.cpp $(COMETSEARCHPATH)/Common.h $(COMETSEARCHPATH)/CometData.h $(COMETSEARCHPATH)/CometDataInternal.h $(COMETSEARCHPATH)/ThreadPool.h $(COMETSEARCHPATH)/CometPostAnalysis.h $(COMETSEARCHPATH)/CometMassSpecUtils.h $(COMETSEARCHPATH)/CometInterfaces.h $(COMETSEARCHPATH)/CometDecoys.h
	${MCC} ${CXXFLAGS} $(COMETSEARCHPATH)/CometPostAnalysis.cpp -c
CometWriteOut.o:      $(COMETSEARCHPATH)/CometWriteOut.cpp $(COMETSEARCHPATH)/Common.h $(COMETSEARCHPATH)/CometData.h $(COMETSEARCHPATH)/CometDataInternal.h $(COMETSEARCHPATH)/CometMassSpecUtils.h $(COMETSEARCHPATH)/CometWriteOut.h $(COMETSEARCHPATH)/CometInterfaces.h
	${MCC} ${CXXFLAGS} $(COMETSEARCHPATH)/CometWriteOut.cpp -c
CometWriteSqt.o:      $(COMETSEARCHPATH)/CometWriteSqt.cpp $(COMETSEARCHPATH)/Common.h $(COMETSEARCHPATH)/CometData.h $(COMETSEARCHPATH)/CometDataInternal.h $(COMETSEARCHPATH)/CometMassSpecUtils.h $(COMETSEARCHPATH)/CometWriteSqt.h $(COMETSEARCHPATH)/CometInterfaces.h
	${MCC} ${CXXFLAGS} $(COMETSEARCHPATH)/CometWriteSqt.cpp -c
CometWritePepXML.o:   $(COMETSEARCHPATH)/CometWritePepXML.cpp $(COMETSEARCHPATH)/Common.h $(COMETSEARCHPATH)/CometData.h $(COMETSEARCHPATH)/CometDataInternal.h $(COMETSEARCHPATH)/CometMassSpecUtils.h $(COMETSEARCHPATH)/CometWritePepXML.h $(COMETSEARCHPATH)/CometInterfaces.h
	${MCC} ${CXXFLAGS} $(COMETSEARCHPATH)/CometWritePepXML.cpp -c
CometWritePercolator.o:   $(COMETSEARCHPATH)/CometWritePercolator.cpp $(COMETSEARCHPATH)/Common.h $(COMETSEARCHPATH)/CometData.h $(COMETSEARCHPATH)/CometDataInternal.h $(COMETSEARCHPATH)/CometMassSpecUtils.h $(COMETSEARCHPATH)/CometWritePercolator.h $(COMETSEARCHPATH)/CometInterfaces.h
	${MCC} ${CXXFLAGS} $(COMETSEARCHPATH)/CometWritePercolator.cpp -c
CometWriteTxt.o:      $(COMETSEARCHPATH)/CometWriteTxt.cpp $(COMETSEARCHPATH)/Common.h $(COMETSEARCHPATH)/CometData.h $(COMETSEARCHPATH)/CometDataInternal.h $(COMETSEARCHPATH)/CometMassSpecUtils.h $(COMETSEARCHPATH)/CometWriteTxt.h $(COMETSEARCHPATH)/CometInterfaces.h
	${MCC} ${CXXFLAGS} $(COMETSEARCHPATH)/CometWriteTxt.cpp -c
CometSearchManager.o:      $(COMETSEARCHPATH)/CometSearchManager.cpp $(COMETSEARCHPATH)/Common.h $(COMETSEARCHPATH)/CometData.h $(COMETSEARCHPATH)/CometDataInternal.h $(COMETSEARCHPATH)/CometMassSpecUtils.h $(COMETSEARCHPATH)/CometSearch.h $(COMETSEARCHPATH)/CometPostAnalysis.h $(COMETSEARCHPATH)/CometWriteOut.h $(COMETSEARCHPATH)/CometWriteSqt.h $(COMETSEARCHPATH)/CometWriteTxt.h $(COMETSEARCHPATH)/CometWritePepXML.h $(COMETSEARCHPATH)/CometWritePercolator.h $(COMETSEARCHPATH)/Threading.h $(COMETSEARCHPATH)/ThreadPool.h $(COMETSEARCHPATH)/CometSearchManager.h $(COMETSEARCHPATH)/CometInterfaces.h
	${MCC} ${CXXFLAGS} $(COMETSEARCHPATH)/CometSearchManager.cpp -c
CometInterfaces.o:      $(COMETSEARCHPATH)/CometInterfaces.cpp $(COMETSEARCHPATH)/Common.h $(COMETSEARCHPATH)/CometData.h $(COMETSEARCHPATH)/CometDataInternal.h $(COMETSEARCHPATH)/CometMassSpecUtils.h $(COMETSEARCHPATH)/CometSearch.h $(COMETSEARCHPATH)/CometPostAnalysis.h $(COMETSEARCHPATH)/CometWriteOut.h $(COMETSEARCHPATH)/CometWriteSqt.h $(COMETSEARCHPATH)/CometWriteTxt.h $(COMETSEARCHPATH)/CometWritePepXML.h $(COMETSEARCHPATH)/CometWritePercolator.h $(COMETSEARCHPATH)/Threading.h $(COMETSEARCHPATH)/ThreadPool.h $(COMETSEARCHPATH)/CometSearchManager.h $(COMETSEARCHPATH)/CometInterfaces.h
	${MCC} ${CXXFLAGS} $(COMETSEARCHPATH)/CometInterfaces.cpp -c
