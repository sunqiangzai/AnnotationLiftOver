/*
 * =====================================================================================
 *
 *       Filename:  alignNeedlemanForTranscript.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  06/02/2017 15:11:41
 *       Revision:  none
 *       Compiler:  g++
 *
 *         Author:  Baoxing Song (songbx.me), song@mpipz.mpg.de
 *   Organization:  MPIPZ
 *
 * =====================================================================================
 */
/*************************************************************************




 ************************************************************************/
#ifndef _ALIGNNEEDLEMANWUNSCHFORTRANSCRIPT_H
#define _ALIGNNEEDLEMANWUNSCHFORTRANSCRIPT_H

#include <string>
#include <map>
#include <vector>
#include "nucleotideCodeSubstitutionMatrix.h"


enum ELEMENTS
{
    START, STOP, EXON, INTRON, SPLICEDONOR, SPLICEACCEPTOR
};

class SpliceSitePosition{
    private:
        int _donorSpliceSitePosition;
        int _acceptorSpliceSitePosition;
    public:
        SpliceSitePosition(int donorsSpliceSitePosition, int acceptorSpliceSitePosition);
        int getDonorSpliceSitePosition();
        int getAcceptorSpliceSitePosition();
};

class NeedlemanWunschForTranscript {
    private:
    std::string _alignment_a;
    std::string _alignment_b;
//    std::string _signs;

    std::string _dna_a, _dna_b;
    int _length_of_a, _length_of_b;

    double _exon_match_score;
    double _exon_mismatch_score;
    double _exon_open_gap_penalty;
    double _exon_extend_gap_penalty;
    double** _exon_subsitition_matrix;
    //std::map<char, <std::map<char, double> > _exon_subsitition_matrix;

    double _intron_match_score;
    double _intron_mismatch_score;
    double _intron_open_gap_penalty;
    double _intron_extend_gap_penalty;
    double** _intron_subsitition_matrix;

    double _start_stop_codon_match_score;
    double _start_stop_codon_mismatch_score;
    double _start_stop_codon_open_gap_penalty;
    double _start_stop_codon_extend_gap_penalty;
    double** _start_stop_codon_subsitition_matrix;

    double _splice_sites_match_score;
    double _splice_sites_mismatch_score;
    double _splice_sites_open_gap_penalty;
    double _splice_sites_extend_gap_penalty;
    double** _splice_sites_subsitition_matrix;

    double **_similarity_matrix;
    VARIANTCATEGORY **_track_matrix;

//    double _similarity, _gaps, _identity;
    NucleotideCodeSubstitutionMatrix nucleotideCodeSubstitutionMatrix;
    int _startCodonPosition;
    int _stopCodonPosition;
    std::vector<SpliceSitePosition> _spliceSitePositions;

    public:
        NeedlemanWunschForTranscript(std::string& dna_a, std::string& dna_b,
        double& exon_match_score, double& exon_mis_match_score,
        double& exon_open_gap_penalty, double& exon_extend_gap_penalty,
        double& intron_match_score, double& intron_mis_match_score,
        double& intron_open_gap_penalty, double& intron_extend_gap_penalty,
        double& split_sites_match_score, double& split_sites_mis_match_score,
        double& split_sites_open_gap_penalty, double& split_sites_extend_gap_penalty,
        double& start_stop_codon_match_score, double& start_stop_codon_mis_match_score,
        double& start_stop_codon_open_gap_penalty, double& start_stop_codon_extend_gap_penalty,
        int& startCodonPosition, int& stopCodonPosition, std::vector<SpliceSitePosition>& spliceSitePositions);
        NeedlemanWunschForTranscript(std::string& dna_a, std::string& dna_b, int startCodonPosition,
                                     int stopCodonPosition, std::vector<SpliceSitePosition>& splitSitePositions);
        void initialize();
        void initialize1();
        void initialize2();
        NeedlemanWunschForTranscript()=delete;
        ~NeedlemanWunschForTranscript();
        std::string getAlignment_a(){return _alignment_a;}
        std::string getAlignment_b(){return _alignment_b;}
        void setScore(double match, double insertion, double deletion, int i, int j);
        void calculate_similarity();
    //    void populate_subs_matrix(std::string subs_matrix_file);
        void dna_align();
        void print_results();
        ELEMENTS checkElements( int position );
};

#endif

