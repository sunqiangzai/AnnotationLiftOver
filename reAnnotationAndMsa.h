/*
 * =====================================================================================
 *
 *       Filename:  reAnnotationAndMsa.h
 *
 *    Description:
 *
 *        Version:  1.0
 *        Created:  06/03/2017 23:33:37
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Baoxing Song (songbx.me), song@mpipz.mpg.de
 *   Organization:  MPIPZ
 *
 * =====================================================================================
 */
/*************************************************************************




 ************************************************************************/


#ifndef VARIANTSANNOTATIONSOFTWARE_REANNOTATIONANDMSA_H
#define VARIANTSANNOTATIONSOFTWARE_REANNOTATIONANDMSA_H
#include <string>
#include "model.h"
#include "myutil.h"
#include "myfunctions.h"
#include <iostream>
#include <atomic>
#include <mutex>

void reAnnotationSingleLine( std::string& referenceGenomeFilePath, std::string& referenceGffFilePath, std::string& sdiFile,
                             std::map<std::string, std::vector<Gene> >& genes, int maxThread, std::string& regex, std::string& regexG, int & lengthThread, std::string & vcfFix);
void reAnnotationSingleLine( std::map<std::string, std::vector<Transcript> >& referenceTranscriptHashSet,
                             std::map<std::string, Fasta>& referenceGenome, std::map<std::string, Transcript>& targetTranscriptsHashMap,
                             std::map<std::string, Fasta>& targetGenome, NucleotideCodeSubstitutionMatrix& nucleotideCodeSubstitutionMatrix,
                             int& maxThread, int & lengthThread);
void reAnnotationAndMsa( std::string& referenceGenomeFilePath, std::string referenceGffFilePath, std::map<std::string, std::string> sdiFilePaths, int& maxThread,  std::string& regex, int & lengthThread );
//void * transcriptRealignment( void *arg );
//void transcriptRealignment( std::string& accessionId, std::string& chromosomeName, std::map<std::string, Fasta>& targetGenome, Transcript* it3, MyThreadCount& myThreadCount );

void transcriptRealignment( Transcript& tartgetTranscript, Transcript& referenceTranscript,
                            NucleotideCodeSubstitutionMatrix& nucleotideCodeSubstitutionMatrix,
                            std::map<std::string, Fasta>& targetGenome,
                            std::map<std::string, Fasta>& referenceGenome, std::string chromosomeName, std::map<std::string,
        Transcript>& targetTranscriptsHashMap, std::atomic_int & number_of_runing_threads, int & lengthThread  );
void transcriptRealignmentAndExonerate( Transcript& tartgetTranscript, Transcript& referenceTranscript,
                                        NucleotideCodeSubstitutionMatrix& nucleotideCodeSubstitutionMatrix,
                                        std::map<std::string, Fasta>& targetGenome,
                                        std::map<std::string, Fasta>& referenceGenome, std::string chromosomeName, std::map<std::string,
        Transcript>& targetTranscriptsHashMap, std::atomic_int & number_of_runing_threads , std::string & prefixUuid, int & lengthThread );


void reAnnotationAndExonerate( std::string& referenceGenomeFilePath, std::string referenceGffFilePath, std::string sdiFile,
                               std::map<std::string, std::vector<Gene> >& genes, int maxThread, std::string& regex, std::string& regexG, std::string & outputGffFile, int & lengthThreadre, std::string & vcfFix);
void reAnnotationAndExonerateAndNovo( std::string& referenceGenomeFilePath, std::string referenceGffFilePath, std::string novoGffFilePath, std::string sdiFile,
                                       std::map<std::string, std::vector<Gene> >& genes, int maxThread, std::string& regex, std::string& regexG,
                                       std::string novoRegex, std::string& novoRegexG, std::string & outputGffFile, int & lengthThread, std::string & vcfFix);
void runExonerateEst(std::string& transcriptName, std::string& cdsSequence, std::string& targetSequence,
                     NucleotideCodeSubstitutionMatrix& nucleotideCodeSubstitutionMatrix,
                     std::map<std::string, Transcript>& targetTranscriptsHashMap, int& startTarget, int& endTarget,
                     STRAND strand, std::string& tchromeSomeName, std::string& prefixUuid, std::map<std::string, Fasta>& targetGenome);
void readExonerateEstResult( std::string& transcriptName, std::string& cdsSequence, std::string& targetSequence,
                             NucleotideCodeSubstitutionMatrix& nucleotideCodeSubstitutionMatrix,
                             std::map<std::string, Transcript>& targetTranscriptsHashMap, int& startTarget, int& endTarget,
                             STRAND& strand, std::string& tchromeSomeName, std::string& fileLocation, std::map<std::string, Fasta>& targetGenome);
void runExonerateProtein(std::string& transcriptName, std::string& protenSequene, std::string& targetSequence,
                         NucleotideCodeSubstitutionMatrix& nucleotideCodeSubstitutionMatrix,
                         std::map<std::string, Transcript>& targetTranscriptsHashMap, int& startTarget, int& endTarget,
                         STRAND& strand, std::string& tchromeSomeName, std::string& fileLocation, std::map<std::string, Fasta>& targetGenome );

void readExonerateProteinResult( std::string& fileLocation, NucleotideCodeSubstitutionMatrix& nucleotideCodeSubstitutionMatrix,
                                 int& startTarget, int& endTarget, STRAND& strand, std::map<std::string, Transcript>& targetTranscriptsHashMap,
                                 std::string& transcriptName, std::string & tchromeSomeName, std::map<std::string, Fasta>& targetGenome);
void readAugustusGff( std::string& fileLocation );
void runSystemCommand( std::string command, std::atomic_int & number_of_runing_threads );
#endif //VARIANTSANNOTATIONSOFTWARE_REANNOTATIONANDMSA_H
