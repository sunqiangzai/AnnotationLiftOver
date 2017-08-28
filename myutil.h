// =====================================================================================
// 
//       Filename:  myutil.h
// 
//    Description:  
// 
//        Version:  1.0
//        Created:  04/12/2017 04:26:31 PM
//       Revision:  none
//       Compiler:  g++
// 
//         Author:  Baoxing Song (songbx.me), song@mpipz.mpg.de
//        Company:  MPIPZ
// 
// =====================================================================================
/*************************************************************************




 ************************************************************************/

#ifndef _MYUTIL_H
#define _MYUTIL_H
#include <algorithm>
#include "model.h"
#include <string>
#include <vector>
#include <map>
#include "nucleotideCodeSubstitutionMatrix.h"

std::vector<std::string> &split(const std::string &s, char delim,std::vector<std::string> &elems);
std::vector<std::string> split(const std::string &s, char delim);


//input a string, the string could be trasformed to uppercas
void songToUpCase( std::string& str );
/*  input a string, the transformed would also be returned */
void songToLowCase( std::string& str);
std::string& songStrRemoveBlank(std::string& str);
std::string& songStrReplaceAll(std::string& str, std::string& pattern, std::string& pattern2);

std::string getReverseComplementary(std::string& sequence);

//take the path of a fasta file as input and return a vector of fastaRecords
void readFastaFile( const std::string& filePath, std::map<std::string, Fasta>& sequences);
std::string getSubsequence(std::map<std::string, Fasta>& sequences, std::string seqName, int start, int end);
std::string getSubsequence(std::map<std::string, Fasta>& sequences, std::string seqName, int start, int end, STRAND strand);
//take the path of a sdi file as input and return a HashMap of sdi records
void readSdiFile (const std::string& filePath, std::map<std::string, std::vector<Variant> >& variantsMap, std::string & vcfFix);
void readSdiFile (const std::string& filePath, std::map<std::string, std::vector<Variant> >& variantsMap, const std::string& chromosome, std::string & vcfFix);
void readGffFile (const std::string& filePath, std::map<std::string, std::vector<Transcript> >& transcriptHashSet, std::string& cdsParentRegex);
int getChangedFromBasement(std::string chromosomeName, int basement, std::map<std::string, std::vector<Variant> >& variantsMap);

//take the path of a sdi file as input and return a HashMap of transcripts
void readGffFile (const std::string& filePath, std::map<std::string, std::vector<Transcript> >& transcriptHashSet);


void annotationLiftOver(std::map<std::string, std::vector<Transcript> >& refTranscriptHashSet,
                        std::map<std::string, Transcript >& targetTranscriptHashMap,
                        std::map<std::string, std::vector<Variant> >& variantsMap,
                        std::map<std::string, Fasta>& targetGenome, std::map<std::string, Fasta>& referenceGenome,
                        std::string chromosome);

void annotationLiftOver(std::map<std::string, std::vector<Transcript> >& refTranscriptHashSet,
                        std::map<std::string, Transcript >& targetTranscriptHashMap,
                        std::map<std::string, std::vector<Variant> >& variantsMap,
                        std::map<std::string, Fasta>& targetGenome, std::map<std::string, Fasta>& referenceGenome);

void annotationLiftOver(std::map<std::string, std::vector<Transcript> >& refTranscriptHashSet,
                        std::map<std::string, Transcript >& targetTranscriptHashMap,
                        std::map<std::string, std::vector<Variant> >& variantsMap,
                        std::map<std::string, Fasta>& targetGenome, std::map<std::string, Fasta>& referenceGenome,
                        std::set<std::string>& chromosomes);
void checkOrfState( Transcript& targetTranscript, Transcript &referenceTranscript,
                    std::map<std::string, Fasta>& targetGenome,
                    std::map<std::string, Fasta>& referenceGenome,
                    NucleotideCodeSubstitutionMatrix& nucleotideCodeSubstitutionMatrix);
void checkOrfState( Transcript& targetTranscript,
                    std::map<std::string, Fasta>& targetGenome,
                    NucleotideCodeSubstitutionMatrix& nucleotideCodeSubstitutionMatrix);

bool ifSpliceSitesOk(Transcript& targetTranscript, Transcript &referenceTranscript,
                     std::map<std::string, Fasta>& targetGenome,
                     std::map<std::string, Fasta>& referenceGenome, NucleotideCodeSubstitutionMatrix& nucleotideCodeSubstitutionMatrix);
bool checkSpliceSites( std::string& s1t, std::string& s2t, NucleotideCodeSubstitutionMatrix& nucleotideCodeSubstitutionMatrix);
bool ifSpliceSitesOk(Transcript& targetTranscript,
                     std::map<std::string, Fasta>& targetGenome,
                     NucleotideCodeSubstitutionMatrix& nucleotideCodeSubstitutionMatrix);
bool checkSpliceSites( std::string& s1t, std::string& s2t, NucleotideCodeSubstitutionMatrix& nucleotideCodeSubstitutionMatrix);

bool checkSpliceSites(std::string& s1, std::string& s2, std::string& s1t, std::string& s2t, NucleotideCodeSubstitutionMatrix& nucleotideCodeSubstitutionMatrix);
std::string agIUPACcodesTranslation(std::string& ag);
std::string gtIUPACcodesTranslation(std::string& gt);
bool ifLengthDivisibleByThree(std::string& sequence);
bool ifNewStopCodon(std::string& cdsSequence, NucleotideCodeSubstitutionMatrix& nucleotideCodeSubstitutionMatrix);
bool ifEndWithStopCodon(std::string& cdsSequence, NucleotideCodeSubstitutionMatrix& nucleotideCodeSubstitutionMatrix);
bool ifStartWithStartCodon(std::string& cdsSequence, NucleotideCodeSubstitutionMatrix& nucleotideCodeSubstitutionMatrix);



//translation related begin
std::string nA2AA(std::string& seq);
std::string nA2AANoDeletedIndel(std::string& seq);
std::string nA2AA(std::string& seq, NucleotideCodeSubstitutionMatrix& nucleotideCodeSubstitutionMatrix);
std::string nA2AANoDeletedIndel(std::string& seq, NucleotideCodeSubstitutionMatrix& nucleotideCodeSubstitutionMatrix);
//translation related end

bool if_file_exists (const std::string& name);
long GetFileSize(std::string& filename);
#endif

