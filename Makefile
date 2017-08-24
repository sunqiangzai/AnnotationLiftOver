cc=g++
obj= Cds.o song.o myutil.o myfunctions.o alignNeedlemanForTranscript.o nucleotideCodeSubstitutionMatrix.o reAnnotationAndMsa.o model.o parameters.o InputParser.o controlLayer.o
CXXFLAGS = -std=c++0x -I ./ -std=c++11 -pthread
song: $(obj)
	$(cc) -o AnnotationLiftOver $(obj) $(CXXFLAGS)
	rm -rf *.o
Cds.o: Cds.cpp
	$(cc) -c Cds.cpp $(CXXFLAGS)
myutil.o:	myutil.cpp
	$(cc) -c myutil.cpp $(CXXFLAGS)
myfunctions.o:	myfunctions.cpp
	$(cc) -c myfunctions.cpp $(CXXFLAGS)
song.o:	song.cpp
	$(cc) -c song.cpp $(CXXFLAGS)
alignNeedlemanForTranscript.o: alignNeedlemanForTranscript.cpp
	$(cc) -c alignNeedlemanForTranscript.cpp $(CXXFLAGS)
alignNeedlemanWunsch.o: alignNeedlemanWunsch.cpp
	$(cc) -c alignNeedlemanWunsch.cpp $(CXXFLAGS)
nucleotideCodeSubstitutionMatrix.o: nucleotideCodeSubstitutionMatrix.cpp
	$(cc) -c nucleotideCodeSubstitutionMatrix.cpp $(CXXFLAGS)
reAnnotationAndMsa.o: reAnnotationAndMsa.cpp
	$(cc) -c reAnnotationAndMsa.cpp $(CXXFLAGS)
model.o: model.cpp
	$(cc) -c model.cpp $(CXXFLAGS)
parameters.o: parameters.cpp
	$(cc) -c parameters.cpp $(CXXFLAGS)
InputParser.o: InputParser.cpp
	$(cc) -c InputParser.cpp $(CXXFLAGS)
controlLayer.o: controlLayer.cpp
	$(cc) -c controlLayer.cpp $(CXXFLAGS)
clean:
	rm -rf *.o song
