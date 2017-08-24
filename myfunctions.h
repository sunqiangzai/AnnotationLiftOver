/*
 * =====================================================================================
 *
 *       Filename:  myfunctions.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  05/31/2017 09:51:38
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

#ifndef _MYFUNCTIONS_H
#define _MYFUNCTIONS_H
#include <string>
#include "model.h"
#include <vector>
#include <map>
#include "reAnnotationAndMsa.h"

int getPseudoGenomeSequence(std::string& referenceGenomeFastaFile, std::string& sdiFile, std::map<std::string, Fasta>& targetSequences);
int getPseudoGenomeSequence(std::map<std::string, Fasta>& referenceSequences, std::string& sdiFile, std::map<std::string, Fasta>& targetSequences);
int getPseudoGenomeSequence( std::string& referenceGenomeFastaFile, std::string& sdiFile, std::string& outputFile);
int getPseudoGenomeSequence(std::map<std::string, Fasta>& referenceSequences, std::map<std::string, std::vector<Variant> >& variantsMap,
                            std::map<std::string, Fasta>& targetSequences);
int getPseudoGenomeSequence(std::map<std::string, Fasta>& referenceSequences, std::map<std::string, std::vector<Variant> >& variantsMap,
                            std::map<std::string, Fasta>& targetSequences, std::string& chromosome);
int getPseudoGenomeSequence(std::map<std::string, Fasta>& referenceSequences, std::map<std::string, std::vector<Variant> >& variantsMap,
                            std::map<std::string, Fasta>& targetSequences, std::set<std::string>& chromosomes);

void writeFasta(std::ostream& out, std::string& seqname, std::string& sequence, int& linewidth);
int myCoordinateLiftOver( std::string& sdiFile, std::string& chromosome, int& position );
void myGffCoordinateLiftOver( std::string& sdiFile, std::string& gffFile, std::string& outputFile );

void getSequences(std::string& gffFile, std::string& genome, std::string& outputProteinSequences,
                  std::string& outputCdsSequences, std::string& outputGenomeSequences);
void getSequences(std::string& gffFile, std::string& genome, std::string& outputProteinSequences,
                  std::string& outputCdsSequences, std::string& outputGenomeSequences, std::string& regex);
void myAnnotationLiftOver( std::string& gffFile, std::string& referenceGenomeFile,
                           std::string& sdiFile, std::string& outputGffFile,
                           std::string& regex, std::string& regexG);
void myReAnnotationLiftoverSingleLine( std::string& referenceGenomeFile, std::string& inputGffFile,
                                       std::string& variantsFile, std::string& outputGffFile, std::string& regex, int &maxThread, std::string& regexG, int & lengthThread );
void myReAnnotationLiftoverAndOrthologous( std::string& referenceGenomeFile, std::string& inputGffFile,
                                           std::string& variantsFile, std::string& outputGffFile, std::string& regex, int &maxThread, std::string& regexG, int & lengthThread );
void myReAnnotationAndExonerateAndNovo( std::string& referenceGenomeFile, std::string& inputGffFile,std::string novoGffFilePath,
                                         std::string& variantsFile, std::string& outputGffFile, std::string& regex, int &maxThread,
                                         std::string& regexG,std::string novoRegex, std::string& novoRegexG, int & lengthThread );
void countNumberOfTwoneighborSNP( std::string& sdiFile, std::string & outputPrefix, int & rangeLength);
void countNumberSNPAndIndel( std::string& sdiFile, std::string & outputPrefix, int & rangeLength);
void generateRandomSdi( std::string& sdiFile, std::string& referenceGenomeFile, std::string & outputPrefix );
#endif
