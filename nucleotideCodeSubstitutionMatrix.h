/*
 * =====================================================================================
 *
 *       Filename:  nucleotideCodeSubstitutionMatrix.h
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

#ifndef _NUCLEOTIDECODESUBSTITUTIONMATRIX_H
#define _NUCLEOTIDECODESUBSTITUTIONMATRIX_H

#include <map>
#include <set>
#include <string>

enum VARIANTCATEGORY
{
    SNP, INSERTION, DELETION, SNPORINSERTION, SNPORDELETION, INSERTIONORDELETION, SNPORINSERTIONORDELETION
};

enum BEGINMIDDLEEND //  FOR DNA TO PROTEIN TRANSLATION
{
    BEGIN, MIDDLE, END
};

class NucleotideCodeSubstitutionMatrix{
    private:
        std::map<char, int> dna_acid_map;
//        std::map<std::string, double> allStartCodons; // it was designed for transcript alignment, do not use it anymore
//        std::map<std::string, double> allStopCodons;
        std::set<std::string> mustStartCodons;
        std::set<std::string> mustStopCodons;
        std::set<std::string> possibleStartCodons;
        std::set<std::string> possibleStopCodons;
        std::set<std::string> basicStartCodons;
        std::set<std::string> basicStopCodons;
        std::map<std::string, std::set<std::string> > dnaIupacCode;
        std::map<std::string, char > middleStandardGeneticCode;

        double** _exon_subsitition_matrix;
        double** _intron_subsitition_matrix;
        double** _start_stop_codon_subsitition_matrix;
        double** _splice_sites_subsitition_matrix;

        std::set<std::string> legalNasString;
        std::set<char> legalNasChar;


    public:
        NucleotideCodeSubstitutionMatrix();
        ~NucleotideCodeSubstitutionMatrix();
        std::map<char, int>& get_dna_acid_map();
        //std::map<std::string, double>& get_allStartCodons();
        //std::map<std::string, double>& get_allStopCodons();
        std::set<std::string>& getMustStartCodons();
        std::set<std::string>& getMustStopCodons();
        std::set<std::string>& getPossibleStartCodons();
        std::set<std::string>& getPossibleStopCodons();
        std::set<std::string>& getBasicStartCodons();
        std::set<std::string>& getBasicStopCodons();
        //give a sequence, return all the possiable iupac combinations
        void getAllPossibleWithIupac(std::string& seq, std::set<std::string>& currentPossibleCombinations);
        std::map<std::string, std::set<std::string> >& getDnaIupacCode();
        double nucleotide_substitution_matrix[17][17];
        char getGeneticCode(std::string coding, BEGINMIDDLEEND beginmiddleend );
        char getGeneticCode(std::string coding);

        double** get_exon_subsitition_matrix();
        double** get_intron_subsitition_matrix();
        double** get_start_stop_codon_subsitition_matrix();
        double** get_splice_sites_subsitition_matrix();
        std::set<std::string>& getLegalNasString();
        std::set<char>& getLegalNasChar();
};
#endif
//std::unordered_set<std::string> stopCodons = {"ATG", "TTG", "CTG"};// = {"TAA", "TAG", "TGA"}; // TO DO, use iupac code



