cmake_minimum_required(VERSION 3.8)
project(geneStructureBasedRealignment)

set(CMAKE_CXX_STANDARD 11)

set(SOURCE_FILES
        alignNeedlemanForTranscript.cpp
        alignNeedlemanForTranscript.h
        alignNeedlemanWunsch.cpp
        alignNeedlemanWunsch.h
        Cds.cpp
        Cds.h
        controlLayer.cpp
        controlLayer.h
        InputParser.cpp
        InputParser.h
        model.cpp
        model.h
        myfunctions.cpp
        myfunctions.h
        myutil.cpp
        myutil.h
        nucleotideCodeSubstitutionMatrix.cpp
        nucleotideCodeSubstitutionMatrix.h
        parameters.cpp
        parameters.h
        reAnnotationAndMsa.cpp
        reAnnotationAndMsa.h
        sole.h
        song.cpp)

add_executable(AnnotationLiftOver ${SOURCE_FILES})