/*
 * =====================================================================================
 *
 *       Filename:  alignNeedlemanForTranscript.cpp
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  06/02/2017 15:11:39
 *       Revision:  none
 *       Compiler:  g++
 *
 *         Author:  Baoxing Song (songbx.me), song@mpipz.mpg.de
 *   Organization:  MPIPZ
 *
 * =====================================================================================
 */
#include <stdlib.h>
/*************************************************************************




 ************************************************************************/

#include<iostream>
#include <map>
#include "alignNeedlemanForTranscript.h"
#include <fstream>
#include <sstream>
#include <iomanip>
#include <algorithm>
#include "parameters.h"

SpliceSitePosition::SpliceSitePosition(int donorsSpliceSitePosition, int acceptorSpliceSitePosition) {
    this->_donorSpliceSitePosition=donorsSpliceSitePosition;
    this->_acceptorSpliceSitePosition=acceptorSpliceSitePosition;
}
int SpliceSitePosition::getDonorSpliceSitePosition(){
    return _donorSpliceSitePosition;
}

int SpliceSitePosition::getAcceptorSpliceSitePosition(){
    return _acceptorSpliceSitePosition;
}
NeedlemanWunschForTranscript::NeedlemanWunschForTranscript(std::string dna_a, std::string dna_b, int startCodonPosition, int stopCodonPosition, std::vector<SpliceSitePosition>& splitSitePositions){
    this->_dna_a = dna_a;
    this->_dna_b = dna_b;
    this->_exon_match_score=getAlignmentExonMatchP();
    this->_exon_mismatch_score=getAlignmentExonMismatchP();
    this->_exon_open_gap_penalty = getAlignmentExonOpenGapP();
    this->_exon_extend_gap_penalty = getAlignmentExonExtendGapP();

    this->_intron_match_score=getAlignmentIntronMatchP();
    this->_intron_mismatch_score=getAlignmentIntronMismatchP();
    this->_intron_open_gap_penalty = getAlignmentIntronOpenGapP();
    this->_intron_extend_gap_penalty = getAlignmentIntronExtendGapP();

    this->_splice_sites_match_score=getAlignmentSpliceSitesMatchP();
    this->_splice_sites_mismatch_score=getAlignmentSpliceSitesMismatchP();
    this->_splice_sites_open_gap_penalty = getAlignmentSpliceSitesOpenGapP();
    this->_splice_sites_extend_gap_penalty = getAlignmentSpliceSitesExtendGapP();

    this->_start_stop_codon_match_score=getAlignmentStartStopCodonMatchP();
    this->_start_stop_codon_mismatch_score=getAlignmentStartStopCodonMismatchP();
    this->_start_stop_codon_open_gap_penalty = getAlignmentStartStopCodonOpenGapP();
    this->_start_stop_codon_extend_gap_penalty = getAlignmentStartStopCodonExtendGapP();
    this->_startCodonPosition = startCodonPosition;
    this->_stopCodonPosition = stopCodonPosition;
    this->_spliceSitePositions=splitSitePositions;
    //std::cout << "needleman transcript initialize" << std::endl;
    this->initialize();
    //std::cout << "needleman transcript calculate_similarity" << std::endl;
    // std::cout << dna_a << std::endl;
    // std::cout << dna_b << std::endl;

    // std::cout << "start " <<  _startCodonPosition << std::endl;
    // std::cout << "stop " << _stopCodonPosition << std::endl;
    this->calculate_similarity();
    // std::cout << "needleman transcript dna_align" << std::endl;
    this->dna_align();
}
NeedlemanWunschForTranscript::NeedlemanWunschForTranscript(
    std::string dna_a, std::string dna_b, 
        double exon_match_score, double exon_mis_match_score, 
        double exon_open_gap_penalty, double exon_extend_gap_penalty,
        double intron_match_score, double intron_mis_match_score, 
        double intron_open_gap_penalty, double intron_extend_gap_penalty,
        double splice_sites_match_score, double splice_sites_mis_match_score,
        double splice_sites_open_gap_penalty, double splice_sites_extend_gap_penalty,
        double start_stop_codon_match_score, double start_stop_codon_mis_match_score, 
        double start_stop_codon_open_gap_penalty, double start_stop_codon_extend_gap_penalty,
        int startCodonPosition, int stopCodonPosition, std::vector<SpliceSitePosition>& spliceSitePositions
    ) {
    this->_dna_a = dna_a;
    this->_dna_b = dna_b;
    this->_exon_match_score=exon_match_score;
    this->_exon_mismatch_score=exon_mis_match_score;
    this->_exon_open_gap_penalty = exon_open_gap_penalty;
    this->_exon_extend_gap_penalty = exon_extend_gap_penalty;

    this->_intron_match_score=intron_match_score;
    this->_intron_mismatch_score=intron_mis_match_score;
    this->_intron_open_gap_penalty = intron_open_gap_penalty;
    this->_intron_extend_gap_penalty = intron_extend_gap_penalty;

    this->_splice_sites_match_score=splice_sites_match_score;
    this->_splice_sites_mismatch_score=splice_sites_mis_match_score;
    this->_splice_sites_open_gap_penalty = splice_sites_open_gap_penalty;
    this->_splice_sites_extend_gap_penalty = splice_sites_extend_gap_penalty;

    this->_start_stop_codon_match_score=start_stop_codon_match_score;
    this->_start_stop_codon_mismatch_score=start_stop_codon_mis_match_score;
    this->_start_stop_codon_open_gap_penalty = start_stop_codon_open_gap_penalty;
    this->_start_stop_codon_extend_gap_penalty = start_stop_codon_extend_gap_penalty;
    this->_startCodonPosition = startCodonPosition;
    this->_stopCodonPosition = stopCodonPosition;
    this->_spliceSitePositions=spliceSitePositions;
    this->initialize();
    this->calculate_similarity();
    this->dna_align();
}
void NeedlemanWunschForTranscript::initialize(){
    this->_length_of_a = this->_dna_a.length();
    this->_length_of_b = this->_dna_b.length();

    this->_similarity = 0;

    this->_similarity_matrix = new double*[this->_length_of_a + 1];
    for (int i = 0; i < (this->_length_of_a + 1); i++)
    {
        this->_similarity_matrix[i] = new double[this->_length_of_b + 1];
    }

    for (int i = 0; i<= _length_of_a; i++) {
        for (int j = 0; j<= 1; j++) {
            _similarity_matrix[i][j] = 0;
        }
    }

    for (int j = 0; j<= _length_of_b; j++){
        for (int i = 0; i<= 1; i++) {
            _similarity_matrix[i][j] = 0;
        }
    }

    // this matrix is for set different penalty for open gap and extend gap begin
    // 0 for match, 1 for deletion, 2 for insertation
    // and the track also changed to use this matrix
    this->_track_matrix = new VARIANTCATEGORY*[this->_length_of_a + 1];
    for (int i = 0; i < (this->_length_of_a + 1); i++) {
        _track_matrix[i] = new VARIANTCATEGORY[this->_length_of_b + 1];
    }

    for (int i = 0; i <= 1; i++) {
        for (itn j = 0; j <= _length_of_b; j++) {
            _track_matrix[i][j] = SNP;
        }
    }
    for (int j = 0; j <= 1; j++){
        for (int i = 0; i <= _length_of_a; i++) {
            _track_matrix[i][j] = SNP;
        }
    }

    for( int i=0; i<17; i++){
        for( int j=0; j<17; j++){
//            if( this->nucleotideCodeSubstitutionMatrix.nucleotide_substitution_matrix[i][j] !=
//                this->nucleotideCodeSubstitutionMatrix.nucleotide_substitution_matrix[j][i] ){
//                std::cerr << "there is something wrong with matirx, please check it " << i << " " << j << std::endl;
//            }
            if( this->nucleotideCodeSubstitutionMatrix.nucleotide_substitution_matrix[i][j] > 0 ){
                this->_exon_subsitition_matrix[i][j] = this->nucleotideCodeSubstitutionMatrix.nucleotide_substitution_matrix[i][j] * _exon_match_score;
                this->_intron_subsitition_matrix[i][j] = this->nucleotideCodeSubstitutionMatrix.nucleotide_substitution_matrix[i][j] * _intron_match_score;
                this->_splice_sites_subsitition_matrix[i][j] = this->nucleotideCodeSubstitutionMatrix.nucleotide_substitution_matrix[i][j] * _splice_sites_match_score;
                this->_start_stop_codon_subsitition_matrix[i][j] = this->nucleotideCodeSubstitutionMatrix.nucleotide_substitution_matrix[i][j] * _start_stop_codon_match_score;
            }else{
                this->_exon_subsitition_matrix[i][j] = -(this->nucleotideCodeSubstitutionMatrix.nucleotide_substitution_matrix[i][j] * _exon_mismatch_score);
                this->_intron_subsitition_matrix[i][j] = -(this->nucleotideCodeSubstitutionMatrix.nucleotide_substitution_matrix[i][j] * _intron_mismatch_score);
                this->_splice_sites_subsitition_matrix[i][j] = -(this->nucleotideCodeSubstitutionMatrix.nucleotide_substitution_matrix[i][j] * _splice_sites_mismatch_score);
                this->_start_stop_codon_subsitition_matrix[i][j] = -(this->nucleotideCodeSubstitutionMatrix.nucleotide_substitution_matrix[i][j] * _start_stop_codon_mismatch_score);
            }
        }
    }
}

NeedlemanWunschForTranscript::~NeedlemanWunschForTranscript()
{
    for (int i = 0; i < (_length_of_a + 1); i++)
    {
        delete[] this->_similarity_matrix[i];
        delete[] this->_track_matrix[i];
    }
    delete[] this->_similarity_matrix;
    delete[] this->_track_matrix;
}

ELEMENTS NeedlemanWunschForTranscript::checkElements( int position ){
    
    if( position < this->_startCodonPosition || position > this->_stopCodonPosition+2  ){
        return INTRON;
    }
    for(std::vector<SpliceSitePosition>::size_type i = 0; i != this->_spliceSitePositions.size(); i++ ){
    	SpliceSitePosition spliceSite = this->_spliceSitePositions[i];
        if(position > (spliceSite.getDonorSpliceSitePosition()+1) && position < (spliceSite.getAcceptorSpliceSitePosition()-1) ){
            return INTRON;
        }else if( position== spliceSite.getAcceptorSpliceSitePosition() || position == spliceSite.getAcceptorSpliceSitePosition()-1){
            return SPLICEACCEPTOR;
        }else if( position== spliceSite.getDonorSpliceSitePosition() || position == spliceSite.getDonorSpliceSitePosition()+1  ){
            return SPLICEDONOR;
        }
    }
    if (position == this->_startCodonPosition || position == this->_startCodonPosition+1 || position == this->_startCodonPosition+2 ){
        return START;
    }
    if (position == this->_stopCodonPosition || position == this->_stopCodonPosition+1 || position == this->_stopCodonPosition+2 ){
        return STOP;
    }
    return EXON;
}

// Calculating similarity matrix
void NeedlemanWunschForTranscript::calculate_similarity(){
    double match = 0, insert = 0, del = 0;
    for (int i=1; i < _length_of_a + 1; ++i){
        //std::cout << "234 " << i << std::endl;

        if( INTRON == this->checkElements(i) ) {
            for (int j = 1; j < _length_of_b + 1; ++j) {
                match = _similarity_matrix[i - 1][j - 1] + _intron_subsitition_matrix
                [this->nucleotideCodeSubstitutionMatrix.get_dna_acid_map()[_dna_a[i - 1]]]
                [this->nucleotideCodeSubstitutionMatrix.get_dna_acid_map()[_dna_b[j - 1]]];

                if (_track_matrix[i - 1][j] == DELETION || _track_matrix[i - 1][j] == SNPORDELETION
                    || _track_matrix[i - 1][j] == SNPORINSERTIONORDELETION) { //deletion
                    del = _similarity_matrix[i - 1][j] + this->_intron_extend_gap_penalty;
                } else {
                    del = _similarity_matrix[i - 1][j] + this->_intron_open_gap_penalty;
                }
                //std::cout << "185 " << std::endl;
                if (i == 1) {
                    insert = 0;
                } else if (_track_matrix[i][j - 1] == INSERTION || _track_matrix[i][j - 1] == SNPORINSERTION ||
                           _track_matrix[i][j - 1] == SNPORINSERTIONORDELETION) { //insertion
                    insert = _similarity_matrix[i][j - 1] + this->_intron_extend_gap_penalty;
                } else {
                    insert = _similarity_matrix[i][j - 1] + this->_intron_open_gap_penalty;
                }

                if (i == _length_of_a && insert < _similarity_matrix[i][j - 1] && _similarity_matrix[i][j - 1] > 0) {
                    insert = _similarity_matrix[i][j - 1];
                }
                //std::cout << "190 " << std::endl;

                setScore( match, insert, del, i, j);
                //std::cout << "226 " << std::endl;
            }
        }else if( EXON == this->checkElements(i)  ){
            for (int j = 1; j < _length_of_b + 1; ++j) {
                //std::cout << "229 " << std::endl;
                match = _similarity_matrix[i - 1][j - 1] + this->_exon_subsitition_matrix
                [this->nucleotideCodeSubstitutionMatrix.get_dna_acid_map()[_dna_a[i - 1]]]
                [this->nucleotideCodeSubstitutionMatrix.get_dna_acid_map()[_dna_b[j - 1]]];

                if (_track_matrix[i - 1][j] == DELETION || _track_matrix[i - 1][j] == SNPORDELETION
                    || _track_matrix[i - 1][j] == SNPORINSERTIONORDELETION) { //deletion
                    del = _similarity_matrix[i - 1][j] + this->_exon_extend_gap_penalty;
                } else {
                    del = _similarity_matrix[i - 1][j] + this->_exon_open_gap_penalty;
                }

                if (i == 1) {
                    insert = 0;
                } else if (_track_matrix[i][j - 1] == INSERTION || _track_matrix[i][j - 1] == SNPORINSERTION ||
                           _track_matrix[i][j - 1] == SNPORINSERTIONORDELETION) { //insertion
                    insert = _similarity_matrix[i][j - 1] + this->_exon_extend_gap_penalty;
                } else {
                    insert = _similarity_matrix[i][j - 1] + this->_exon_open_gap_penalty;
                }
                if (i == _length_of_a && insert < _similarity_matrix[i][j - 1] && _similarity_matrix[i][j - 1] > 0) {
                    insert = _similarity_matrix[i][j - 1];
                }
                //std::cout << "249 " << std::endl;
                setScore( match, insert, del, i, j);
                //std::cout << "283 " << std::endl;
            }
        }else if ( START == this->checkElements(i) || STOP == this->checkElements(i)){
            for (int j = 1; j < _length_of_b + 1; ++j) {
                //std::cout << "285 " << std::endl;
                match = _similarity_matrix[i - 1][j - 1] + this->_start_stop_codon_subsitition_matrix
                [this->nucleotideCodeSubstitutionMatrix.get_dna_acid_map()[_dna_a[i - 1]]]
                [this->nucleotideCodeSubstitutionMatrix.get_dna_acid_map()[_dna_b[j - 1]]];


                if (_track_matrix[i - 1][j] == DELETION || _track_matrix[i - 1][j] == SNPORDELETION
                    || _track_matrix[i - 1][j] == SNPORINSERTIONORDELETION) { //deletion
                    del = _similarity_matrix[i - 1][j] + this->_start_stop_codon_extend_gap_penalty;
                } else {
                    del = _similarity_matrix[i - 1][j] + this->_start_stop_codon_open_gap_penalty;
                }
                if (i == 1) {
                    insert = 0;
                } else if (_track_matrix[i][j - 1] == INSERTION || _track_matrix[i][j - 1] == SNPORINSERTION ||
                           _track_matrix[i][j - 1] == SNPORINSERTIONORDELETION) { //insertion
                    insert = _similarity_matrix[i][j - 1] + this->_start_stop_codon_extend_gap_penalty;
                } else {
                    insert = _similarity_matrix[i][j - 1] + this->_start_stop_codon_open_gap_penalty;
                }
                if (i == _length_of_a && insert < _similarity_matrix[i][j - 1] && _similarity_matrix[i][j - 1] > 0) {
                    insert = _similarity_matrix[i][j - 1];
                }

                setScore( match, insert, del, i, j);

                // test for new start and stop codon score method start
                if (i == _startCodonPosition + 2 && j > 2) {
                    std::stringstream targetthree;
                    targetthree << _dna_b[j - 3];
                    targetthree << _dna_b[j - 2];
                    targetthree << _dna_b[j - 1];
                    std::string startThree = targetthree.str();
                    if (nucleotideCodeSubstitutionMatrix.getPossibleStartCodons().find(startThree) !=
                        nucleotideCodeSubstitutionMatrix.getPossibleStartCodons().end()) {
                        _similarity_matrix[i - 2][j - 2] =
                                _similarity_matrix[i - 3][j - 3] + _start_stop_codon_subsitition_matrix[0][0];
                        _track_matrix[i - 2][j - 2] = SNP;
                        _similarity_matrix[i - 1][j - 1] =
                                _similarity_matrix[i - 2][j - 2] + _start_stop_codon_subsitition_matrix[0][0];
                        _track_matrix[i - 1][j - 1] = SNP;
                        _similarity_matrix[i][j] =
                                _similarity_matrix[i - 1][j - 1] + _start_stop_codon_subsitition_matrix[0][0];
                        _track_matrix[i][j] = SNP;
                    }
                }

                if (i == _stopCodonPosition + 2 && j > 2) {
                    std::stringstream targetthree;
                    targetthree << _dna_b[j - 3];
                    targetthree << _dna_b[j - 2];
                    targetthree << _dna_b[j - 1];
                    std::string stopThree = targetthree.str();
                    // the the query has been the end, and the snp is not higher than insertion, do not change
                    if (nucleotideCodeSubstitutionMatrix.getPossibleStopCodons().find(stopThree) !=
                        nucleotideCodeSubstitutionMatrix.getPossibleStopCodons().end()
                        && !(i == _length_of_a &&
                             (_similarity_matrix[i - 3][j - 3] + 3 * _start_stop_codon_subsitition_matrix[0][0]) <
                             _similarity_matrix[i][j])) {
                        _similarity_matrix[i - 2][j - 2] =
                                _similarity_matrix[i - 3][j - 3] + _start_stop_codon_subsitition_matrix[0][0];
                        _track_matrix[i - 2][j - 2] = SNP;
                        _similarity_matrix[i - 1][j - 1] =
                                _similarity_matrix[i - 2][j - 2] + _start_stop_codon_subsitition_matrix[0][0];
                        _track_matrix[i - 1][j - 1] = SNP;
                        _similarity_matrix[i][j] =
                                _similarity_matrix[i - 1][j - 1] + _start_stop_codon_subsitition_matrix[0][0];
                        _track_matrix[i][j] = SNP;
                    }
                }

                // test for new start and stop codon score method end
            }
        }else if ( SPLICEDONOR == this->checkElements(i) || SPLICEACCEPTOR == this->checkElements(i)  ){
            for (int j = 1; j < _length_of_b + 1; ++j) {
                match = _similarity_matrix[i - 1][j - 1] + this->_splice_sites_subsitition_matrix
                [this->nucleotideCodeSubstitutionMatrix.get_dna_acid_map()[_dna_a[i - 1] ] ]
                [this->nucleotideCodeSubstitutionMatrix.get_dna_acid_map()[_dna_b[j - 1] ] ];

                if( _track_matrix[i - 1][j] == DELETION || _track_matrix[i - 1][j] == SNPORDELETION 
             		|| _track_matrix[i - 1][j] == SNPORINSERTIONORDELETION ){ //deletion
                		del = _similarity_matrix[i - 1][j] + this->_splice_sites_extend_gap_penalty;
            	}else{
                    del = _similarity_matrix[i - 1][j] + this->_splice_sites_open_gap_penalty;
                }

            	if( _track_matrix[i][j - 1] == INSERTION || _track_matrix[i][j - 1] == SNPORINSERTION ||
             		_track_matrix[i][j - 1] == SNPORINSERTIONORDELETION ){ //insertion
                	insert = _similarity_matrix[i][j - 1] + this->_splice_sites_extend_gap_penalty;
            	}else{
                    insert = _similarity_matrix[i][j - 1] + this->_splice_sites_open_gap_penalty;
                }
                setScore( match, insert, del, i, j);
            }
        }
    }
}

void NeedlemanWunschForTranscript::setScore(double match, double insert, double del, int i, int j){
    double selected=0;
    if( del >insert && del==match  ){
        selected = del;
        _track_matrix[i][j] = SNPORDELETION;
    }else if( insert >del && insert == match  ){
        selected = match;
        _track_matrix[i][j] = SNPORINSERTION;
    } else if ( del > match && del > insert){// prefer deletion
        int t = 1;
        while( i-t >=1 && (_track_matrix[i - t][j] == SNPORDELETION || _track_matrix[i - t][j] == SNPORINSERTIONORDELETION ) ){
            _track_matrix[i - t][j] = DELETION;
        }
        selected = del;
        _track_matrix[i][j] = DELETION;
    }else if( insert > match && insert > del ){//prefer insertion, so that the INDELs could be put together
        int t = 1;
        while( j-t >=1 && (_track_matrix[i][j-t] == SNPORINSERTION || _track_matrix[i][j-t] == SNPORINSERTIONORDELETION ) ){
            _track_matrix[i][j-t] = INSERTION;
        }
        selected = insert;
        _track_matrix[i][j] = INSERTION;
    }else if (match > insert && match > del){
        selected = match;
        _track_matrix[i][j] = SNP;
    }else if ( del >match && insert==del  ){
        selected = del;
        _track_matrix[i][j] = INSERTIONORDELETION;
    } else{
        selected = del;
        _track_matrix[i][j] = SNPORINSERTIONORDELETION;
    }
    _similarity_matrix[i][j] = selected;
}

// Trace back step.
void NeedlemanWunschForTranscript::dna_align()
{
    std::string reserve;
    reserve.reserve(_length_of_b*2);
    std::stringstream _alignment_a_sst(reserve);
    std::stringstream _alignment_b_sst(reserve);
    std::stringstream _signs_sst(reserve);

    int i = this->_length_of_a;
    int j = this->_length_of_b;

    while (i > 0 || j > 0)
    {
        //std::cout << i << " " << j << " " << std::endl;
        // Going to S(i-1, j-1) //match
        if (i > 0 && j > 0 && (_track_matrix[i][j]==SNP || _track_matrix[i][j]==SNPORDELETION ||
                  _track_matrix[i][j]==SNPORINSERTION || _track_matrix[i][j]==SNPORINSERTIONORDELETION) ){
            _alignment_a_sst << _dna_a[i - 1];
            _alignment_b_sst << _dna_b[j - 1];

            if (this->nucleotideCodeSubstitutionMatrix.nucleotide_substitution_matrix[this->nucleotideCodeSubstitutionMatrix.get_dna_acid_map()[_dna_a[i - 1]]][this->nucleotideCodeSubstitutionMatrix.get_dna_acid_map()[_dna_b[j - 1]]] > 0){
                if (_dna_a[i - 1] != _dna_b[j - 1]) {
                    _signs_sst << ":";
                } else {
                    _signs_sst << "|";
                    _identity += 1;
                }
                _similarity += 1;
            } else {
                _signs_sst <<  "." ;
            }
            --i;
            --j;
        } else if ( (i>0 && (_track_matrix[i][j]==INSERTIONORDELETION || _track_matrix[i][j]==DELETION )) || (j ==0 && i >0) ) //deletion
        {
            _alignment_a_sst << _dna_a[i - 1];
            _alignment_b_sst << '-';
            _signs_sst << "-";
            ++_gaps;
            --i;
        }// Going to S(i, j-1) //insertion
        else {
            _alignment_a_sst << '-';
            _alignment_b_sst << _dna_b[j - 1];
            _signs_sst << "+";
            ++_gaps;
            --j;
        }
    }

    _alignment_a = _alignment_a_sst.str();
    _alignment_b = _alignment_b_sst.str();
    _signs = _signs_sst.str();
    reverse(_alignment_a.begin(), _alignment_a.end());
    reverse(_alignment_b.begin(), _alignment_b.end());
    reverse(_signs.begin(), _signs.end());
}

void NeedlemanWunschForTranscript::print_results()
{
    // for (size_t j = 0; j <= _length_of_b; ++j){
    //    for (size_t i = 0; i <= _length_of_a; ++i){
    //        std::cout << _similarity_matrix[i][j] << " ";
    //    }
    //    std::cout <<std::endl;
    // }
    std::cout << "\nAlignment: \n\n";
    std::stringstream elements;
    int letter = 0;
    for (ing i = 0; i < _alignment_a.length(); i++)
    {
  //      std::cout << _alignment_a[i];
        if( _alignment_a[i] != '-'  ){
            ++letter;
        }
        if( INTRON == this->checkElements(letter) ){
            elements << 'I';
        }else if( EXON == this->checkElements(letter)  ){
            elements << 'E';
        }else if( START == this->checkElements(letter)  ){
            elements << 'T';
        }else if( STOP == this->checkElements(letter)  ){
            elements << 'P';
        }else if( SPLICEDONOR == this->checkElements(letter)   ){
            elements << 'D';
        }else if( SPLICEACCEPTOR == this->checkElements(letter)   ){
            elements << 'A';
        }
    }
    std::cout << _alignment_a << std::endl;

    std::cout << _signs << std::endl;
    
    std::cout << _alignment_b << std::endl;
    std::cout << elements.str() << std::endl;
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

