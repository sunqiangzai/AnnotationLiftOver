/*
 * =====================================================================================
 *
 *       Filename:  myMsa.cpp
 *
 *    Description:  :
 *
 *        Version:  1.0
 *        Created:  06/14/2017 09:49:09
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (), 
 *   Organization:  
 *
 * =====================================================================================
 */
#include <stdlib.h>
/*************************************************************************




 ************************************************************************/

#include <iostream>

#include <seqan/align.h>
#include <seqan/graph_align.h>
#include <seqan/graph_msa.h>

#include "myMsa.h"

#include <uuid/uuid.h>
using namespace seqan;
std::vector<std::string>& MyMsa::getNames(){
    return _names;
}
std::vector<std::string>& MyMsa::getSeqs(){
    return _seqs;
}
std::vector<int>& MyMsa::getStarts(){
    return _starts;
}
std::vector<int>& MyMsa::getEnds(){
    return _ends;
}
std::vector<std::vector<Variant> >& MyMsa::getVariants(){
    return _variants;
}

void MyMsa::setNames( std::vector<std::string>& names){
    this->_names = names;
}
void MyMsa::setSeqs( std::vector<std::string>& seqs){
    this->_seqs = seqs;
}
void MyMsa::setStarts( std::vector<int>& starts){
    this->_starts=starts;
}
void MyMsa::setEnds( std::vector<int>& ends){
    this->_ends=ends;
}
void MyMsa::setVariants(std::vector<std::vector<Variant> >& variants){
    this->_variants=variants;
}
//
//
//// a interface for both alignment methods should be set and should set windows and overlap
//std::vector<std::string>& myMsaAlignment(MyMsa& myMsa, std::string method, int windowSize, int overlapSize){
/////TO DO
//}
//
//
//
//std::vector<std::string>& myMsaAlignment(std::vector<std::string>& seqs){
//    Align<DnaString> align;
//    resize(rows(align), seqs.size() );
//    for (int i = 0; i < seqs.size(); ++i){
//        assignSource(row(align, i), seqs[i]);
//    }
//    globalMsaAlignment(align, SimpleScore(getMyMsaMatchP(), getMyMsaMisMatchP(), getMyMsaGapP(), getMyMsaOpenGapP() )); //match, mismatch, gap, gapOpen
//    std::vector<std::string> aligns;
//    for (int i = 0; i < myMsa.getSeqs().size(); ++i){
//        aligns.push_back(row(align, i));
//    }
//    return aligns;
//}
//
//std::vector<std::string>& mafftMsaAlignment(MyMsa& myMsa){
//    std::string tempFolder = createdTempFloder();
//
//}
