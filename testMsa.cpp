/*
 * =====================================================================================
 *
 *       Filename:  testMsa.cpp
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  06/04/2017 18:39:06
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

#include<iostream>

#include <seqan/align.h>
#include <seqan/graph_align.h>
#include <seqan/graph_msa.h>

using namespace seqan;

int main()
{
    char const * strings[4] =
    {
        "DPKKPRGKMSSYAFFVQTSREEHKKKHPDASVNFSEFSKKCSERWKTMSAKEKGKFEDMA"
        "KADKARYEREMKTYIPPKGE",
        "RVKRPMNAFIVWSRDQRRKMALENPRMRNSEISKQLGYQWKMLTEAEKWPFFQEAQKLQA"
        "MHREKYPNYKYRPRRKAKMLPK",
        "FPKKPLTPYFRFFMEKRAKYAKLHPEMSNLDLTKILSKKYKELPEKKKMKYIQDFQREKQ"
        "EFERNLARFREDHPDLIQNAKK",
        "HIKKPLNAFMLYMKEMRANVVAESTLKESAAINQILGRRWHALSREEQAKYYELARKERQ"
        "LHMQLYPGWSARDNYGKKKKRKREK"
    };
    Align<DnaString> align;
    resize(rows(align), 4);
    for (int i = 0; i < 4; ++i){
        assignSource(row(align, i), strings[i]);
    }
    globalMsaAlignment(align, SimpleScore(5, -3, -1, -3)); //match, mismatch, gap, gapOpen
    std::cout << row(align,1) << "\n";
   
//    typedef StringSet<DnaString, Dependent<> > TStringSet;
//    typedef Graph<Alignment<TStringSet, void> > TAlignmentGraph;
//    Score<int> scoringScheme(5, -3, -1, -3);
//    TStringSet stringSet;
//    for (int i = 0; i < 4; ++i){
//        DnaString seq = strings[i];
//        appendValue(stringSet, seq);
//    }
//    TAlignmentGraph alignmentGraph(stringSet);
//
//    int score = globalAlignment(alignmentGraph, scoringScheme, Gotoh());
//
//    std::cout << "Score = " << score << "\nsong"
//                      << alignmentGraph << std::endl;
    return 0;
}


