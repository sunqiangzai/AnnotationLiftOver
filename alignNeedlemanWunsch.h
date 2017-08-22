/*
 * =====================================================================================
 *
 *       Filename:  alignNeedlemanWunsch.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  06/02/2017 15:11:41
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

#ifndef _ALIGNNEEDLEMANWUNSCH_H
#define _ALIGNNEEDLEMANWUNSCH_H

#include <string>
#include <map>
#include <set>
#include "nucleotideCodeSubstitutionMatrix.h"

/*
double nucleotide_substitution_matrix[17][17] =
    {
        //A,   T,  C,  G,  U,  N,  R,  Y,  S,  W,  K,  M,  B,   D,   H,  V
        { 1,  -1, -1, -1, -1,  0, 0.5,-1, -1,0.5, -1,0.5, -1, 0.3, 0.3,0.3} ,   //A
        {-1,   1, -1, -1, -1,  0, -1,0.5, -1,0.5,0.5, -1,0.3, 0.3, 0.3, -1} ,   //T
        {-1,  -1,  1, -1, -1,  0, -1,0.5,0.5, -1, -1,0.5,0.3, -1,  0.3,0.3} ,   //C
        {-1,  -1, -1,  1, -1,  0, 0.5,-1,0.5, -1,0.5, -1,0.3, 0.3,  -1,0.3} ,   //G
        {-1,  -1, -1, -1,  1,  0, -1,0.5, -1,0.5,0.5, -1,0.3, 0.3, 0.3, -1} ,   //U
        { 0,   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,   0,   0,  0} ,   //N
        {0.5, -1, -1,0.5, -1,  0,  1, -1, -1, -1, -1, -1, -1,  -1,  -1, -1} ,   //R
        {-1, 0.5,0.5, -1,0.5,  0, -1,  1, -1, -1, -1, -1, -1,  -1,  -1, -1} ,   //Y
        {-1,  -1,0.5,0.5, -1,  0, -1, -1,  1, -1, -1, -1, -1,  -1,  -1, -1} ,   //S
        {0.5,0.5, -1, -1,0.5,  0, -1, -1, -1,  1, -1, -1, -1,  -1,  -1, -1} ,   //W
        {-1, 0.5, -1,0.5,0.5,  0, -1, -1, -1, -1,  1, -1, -1,  -1,  -1, -1} ,   //K
        {0.5, -1,0.5, -1, -1,  0, -1, -1, -1, -1, -1,  1, -1,  -1,  -1, -1} ,   //M
        {-1, 0.3,0.3,0.3,0.3,  0, -1, -1, -1, -1, -1, -1,  1,  -1,  -1, -1} ,   //B
        {0.3,0.3, -1,0.3,0.3,  0, -1, -1, -1, -1, -1, -1, -1,   1,  -1, -1} ,   //D
        {0.3,0.3,0.3, -1,0.3,  0, -1, -1, -1, -1, -1, -1, -1,  -1,   1, -1} ,   //H
        {0.3, -1,0.3,0.3, -1,  0, -1, -1, -1, -1, -1, -1, -1,  -1,  -1,  1}     //V
    };
std::map<char,int> create_map()
{
    std::map<char,int> m;
    m['A'] = 0;
    m['T'] = 1;
    m['C'] = 2;
    m['G'] = 3;
    m['U'] = 4;
    m['N'] = 5;
    m['R'] = 6;
    m['Y'] = 7;
    m['S'] = 8;
    m['W'] = 9;
    m['K'] = 10;
    m['M'] = 11;
    m['B'] = 12;
    m['D'] = 13;
    m['H'] = 14;
    m['V'] = 15;
  return m;
}

std::map<char,int> dna_acid_map = create_map();
*/
/*
std::set<std::string> create_start(){
    std::set<std::string> m;
    m.insert("TTG");
    m.insert("ATG");
    m.insert("CTG");
    m.insert("YTG");
    m.insert("WTG");
    m.insert("MTG");
    m.insert("HTG");
    return m;
}

std::set<std::string> create_stop(){
    std::set<std::string> m;
    m.insert("TAA");
    m.insert("TAG");
    m.insert("TGA");
    m.insert("TRA");
    m.insert("TAR");
    return m;
}

std::set<std::string> allStartCodons = create_stop();
std::set<std::string> allStopCodons = create_stop();
*/
// enum VARIANTCATEGORY
// {
//     SNP, INSERTION, DELETION, SNPORINSERTION, SNPORDELETION, INSERTIONORDELETION, SNPORINSERTIONORDELETION
// };

class NeedlemanWunsch
{
    private:
    std::string _alignment_a;
    std::string _alignment_b;
    std::string _signs;
    std::string _dna_a, _dna_b;
    size_t _length_of_a, _length_of_b;
    double _open_gap_penalty;
    double _extend_gap_penalty;
    double **_similarity_matrix;
    VARIANTCATEGORY **_track_matrix;

    double _substitute_matrix[17][17];
    double _similarity, _gaps, _identity;
    NucleotideCodeSubstitutionMatrix nucleotideCodeSubstitutionMatrix;
    public:
    NeedlemanWunsch(std::string dna_a, std::string dna_b, double match_score, double mis_match_score, double open_gap_penalty, double extend_gap_penalty);
    ~NeedlemanWunsch();
    NeedlemanWunsch()=delete;
    
    void calculate_similarity();
    void populate_subs_matrix(std::string subs_matrix_file);
    void dna_align();
    void print_results();
};
#endif
