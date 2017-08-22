/*
 * =====================================================================================
 *
 *       Filename:  alignNeedlemanWunsch.cpp
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  06/02/2017 15:11:39
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

#include<iostream>
#include <map>
#include "alignNeedlemanWunsch.h"
//#include "nucleotideCodeSubstitutionMatrix.h"
#include "nucleotideCodeSubstitutionMatrix.h"
#include <string>
#include <cstring>
#include <fstream>
#include <cstdlib>
#include <limits>
#include <sstream>
#include <iomanip>

NeedlemanWunsch::NeedlemanWunsch(std::string dna_a, std::string dna_b, double match_score, double mis_match_score, double open_gap_penalty, double extend_gap_penalty)
{
    this->_dna_a = dna_a;
    this->_dna_b = dna_b;
    this->_open_gap_penalty = open_gap_penalty;
    this->_extend_gap_penalty = extend_gap_penalty;
//    std::cout << "row 44" << std::endl;
    this->_length_of_a = dna_a.length();
    this->_length_of_b = dna_b.length();

    this->_similarity = 0;
   // std::cout << "row 49" << std::endl;
    this->_similarity_matrix = new double*[this->_length_of_a + 1];
    for (size_t i = 0; i < (this->_length_of_a + 1); i++)
    {
        this->_similarity_matrix[i] = new double[this->_length_of_b + 1];
    }

    for (size_t i = 0; i<= _length_of_a; i++)
    {
        for (size_t j = 0; j<= _length_of_b; j++)
        {
            _similarity_matrix[i][j] = 0;
        }
    }
    //std::cout << "row 63" << std::endl;
    // this matrix is for set different penalty for open gap and extend gap begin
    // 0 for match, 1 for deletion, 2 for insertation
    // and the track also changed to use this matrix
    this->_track_matrix = new VARIANTCATEGORY*[this->_length_of_a + 1];
    for (size_t i = 0; i < (this->_length_of_a + 1); i++)
    {
        _track_matrix[i] = new VARIANTCATEGORY[this->_length_of_b + 1];
    }
    for (size_t i = 0; i <= _length_of_a; i++)
    {
        for (size_t j = 0; j <= _length_of_b; j++)
        {
            _track_matrix[i][j] = SNP;
        }
    }
    for( int i=0; i<17; i++){
        for( int j=0; j<17; j++){
            if( this->nucleotideCodeSubstitutionMatrix.nucleotide_substitution_matrix[i][j] != 
                    this->nucleotideCodeSubstitutionMatrix.nucleotide_substitution_matrix[j][i] ){
                std::cerr << "there is something wrong with matirx, please check it " << i << " " << j << std::endl;
            }
            if( this->nucleotideCodeSubstitutionMatrix.nucleotide_substitution_matrix[i][j] > 0 ){
                this->_substitute_matrix[i][j] = this->nucleotideCodeSubstitutionMatrix.nucleotide_substitution_matrix[i][j] * match_score;
            }else{
                this->_substitute_matrix[i][j] = -(this->nucleotideCodeSubstitutionMatrix.nucleotide_substitution_matrix[i][j] * mis_match_score);
            }
            //std::cout << i << " " << j << " " << _subsitition_matrix[i][j] << std::endl;
        }
    }
}

NeedlemanWunsch::~NeedlemanWunsch() {
    for (size_t i = 0; i < (_length_of_a + 1); i++) {
        delete[] this->_similarity_matrix[i];
        delete[] this->_track_matrix[i];
    }
    delete[] this->_similarity_matrix;
    delete[] this->_track_matrix;
//    for (int i = 0; i < 17; i++) {
//        delete[] this->_substitute_matrix[i];
//    }
//    delete[] this->_substitute_matrix;
}

// Calculating similarity matrix
void NeedlemanWunsch::calculate_similarity()
{

    // I don't think the following score matrix is correct
    // for (int i = 0; i < (_length_of_a + 1); i++)
    // {
    //     _similarity_matrix[i][0] = i * gap_penalty;
    // }

    // for (int j = 0; j < (_length_of_a + 1); j++)
    // {
    //     _similarity_matrix[0][j] = j * gap_penalty;
    // }


    int match = 0, insert = 0, del = 0 , selected = 0;
    for (size_t i = 1; i < _length_of_a + 1; i++)
    {
        for (size_t j = 1; j < _length_of_b + 1; j++)
        {
            match = _similarity_matrix[i - 1][j - 1] + 
                _substitute_matrix[this->nucleotideCodeSubstitutionMatrix.get_dna_acid_map()[_dna_a[i - 1] ] ][this->nucleotideCodeSubstitutionMatrix.get_dna_acid_map()[_dna_b[j - 1] ] ];
     //       std::cout << i << " " <<  j << " match " << match << std::endl;    
            del = _similarity_matrix[i - 1][j] + this->_open_gap_penalty;
            insert = _similarity_matrix[i][j - 1] + this->_open_gap_penalty;

            if( _track_matrix[i - 1][j] == DELETION || _track_matrix[i - 1][j] == SNPORDELETION 
             || _track_matrix[i - 1][j] == SNPORINSERTIONORDELETION ){ //deletion
                del = _similarity_matrix[i - 1][j] + this->_extend_gap_penalty;
            }
            if( _track_matrix[i][j - 1] == INSERTION || _track_matrix[i][j - 1] == SNPORINSERTION ||
             _track_matrix[i][j - 1] == SNPORINSERTIONORDELETION ){ //insertion
                insert = _similarity_matrix[i][j - 1] + this->_extend_gap_penalty;
            }
     //       std::cout << i << " " <<  j << " del " << del << std::endl;
      //      std::cout << i << " " <<  j << " insert " << insert << std::endl;
            if( del==match && del == insert  ){
                selected = del;
                _track_matrix[i][j] = SNPORINSERTIONORDELETION;
            }else if( del >insert && del==match  ){
                selected = del;
                _track_matrix[i][j] = SNPORDELETION;
            }else if( del >match && insert==del  ){
                selected = del;
                _track_matrix[i][j] = INSERTIONORDELETION;
            }else if( insert >del && insert == match  ){
                selected = match;
                _track_matrix[i][j] = SNPORINSERTION;
            } else if ( del > match && del > insert){// prefer deletion  
                if( _track_matrix[i - 1][j] == SNPORDELETION || _track_matrix[i - 1][j] == SNPORINSERTIONORDELETION  ){
                    _track_matrix[i - 1][j] = DELETION;
                }
                selected = del;
                _track_matrix[i][j] = DELETION;
            }else if( insert > match ){//prefer insertion, so that the INDELs could be put together
                if( _track_matrix[i][j - 1] == SNPORINSERTION || _track_matrix[i][j - 1] == SNPORINSERTIONORDELETION ){
                   _track_matrix[i][j - 1] = INSERTION;
                }
                selected = insert;
                _track_matrix[i][j] = INSERTION;
            }else{
                selected = match;
                _track_matrix[i][j] = SNP;
            }
            _similarity_matrix[i][j] = selected;
 //           std::cout << i << " " << j << " " << _track_matrix[i][j] << std::endl;
        }
    }
}


// Trace back step.
void NeedlemanWunsch::dna_align()
{

    _alignment_a = "";
    _alignment_b = "";
    _signs = "";

    int i = static_cast<int>(this->_length_of_a);
    int j = static_cast<int>(this->_length_of_b);

    for( int a=0; a<=i; a++ ){
        for( int b=0; b<=j; b++ ){
            std::cout << "\t" << _similarity_matrix[a][b];
        }
        std::cout << std::endl;
    }
    while (i > 0 || j > 0)
    {
        // Going to S(i-1, j-1) //match
        if (i > 0 && j > 0 && (_track_matrix[i][j]==SNP || _track_matrix[i][j]==SNPORDELETION ||
                  _track_matrix[i][j]==SNPORINSERTION || _track_matrix[i][j]==SNPORINSERTIONORDELETION  ))
        {
            _alignment_a = _dna_a[i - 1] + _alignment_a;
            _alignment_b = _dna_b[j - 1] + _alignment_b;

            if (_substitute_matrix[this->nucleotideCodeSubstitutionMatrix.get_dna_acid_map()[_dna_a[i - 1]]][this->nucleotideCodeSubstitutionMatrix.get_dna_acid_map()[_dna_b[j - 1]]] > 0)
            {
                if (_dna_a[i - 1] != _dna_b[j - 1])
                {
                    _signs = ":" + _signs;
                }
                else
                {
                    _signs = "|" + _signs;
                    _identity += 1;
                }
                _similarity += 1;
            }
            else
            {
                _signs = "." + _signs;
            }

            i -= 1;
            j -= 1;
        }
        else if ( (i>0 && (_track_matrix[i][j]==INSERTIONORDELETION || _track_matrix[i][j]==DELETION )) || (j ==0 && i >0) ) //deletion
        {
            _alignment_a = _dna_a[i - 1] + _alignment_a;
            _alignment_b = '-' + _alignment_b;
            _signs = " " + _signs;
            _gaps += 1;
            i -= 1;
        }
        // Going to S(i, j-1) //insertion
        else 
        {
            _alignment_a = '-' + _alignment_a;
            _alignment_b = _dna_b[j - 1] + _alignment_b;
            _signs = " " + _signs;
            _gaps += 1;
            j -= 1;
        }
    }
}

void NeedlemanWunsch::print_results()
{
    //std::cout << _similarity_matrix[102][138];
    //std::cout << "\nAlignment: \n\n";
    for (size_t i = 0; i < _alignment_a.length(); i++)
    {
        std::cout << _alignment_a[i];
    }
    std::cout << "\n";

    for (size_t i = 0; i < _alignment_a.length(); i++)
    {
        std::cout << _signs[i];
    }
    std::cout << "\n";

    for (size_t i = 0; i < _alignment_b.length(); i++)
    {
        std::cout << _alignment_b[i];
    }
    std::cout << "\n";
    std::cout << std::setfill(' ') << std::setw(50) << "\n";

    double percentage;
    std::cout << "Score: " << _similarity_matrix[_length_of_a][_length_of_b] << "\n";

    std::cout << "Length: " << _alignment_a.length() << " (with gaps)\n";

    percentage = ((double)_identity / (double)_alignment_a.length()) * 100;
    std::cout << "Identity: " << _identity << "/" <<_alignment_a.length() << " ( %" << percentage << " ) " << "\n";

    percentage = ((double)_similarity / (double)_alignment_a.length()) * 100;
    std::cout << "Similarity: " << _similarity << "/" << _alignment_a.length() <<" ( %" << percentage << " ) " << "\n";

    percentage = ((double)_gaps / (double)_alignment_a.length()) * 100;
    std::cout << "Gaps: " << _gaps << "/" << _alignment_a.length() << " ( %" << percentage << " ) " << "\n";
}

