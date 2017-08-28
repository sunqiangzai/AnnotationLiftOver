/*
 * =====================================================================================
 *
 *       Filename:  InputParser.cpp
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  06/25/2017 01:13:39
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Baoxing Song (songbx.me), song@mpipz.mpg.de
 *   Organization:  MPIPZ
 *
 * =====================================================================================
 */
#include <stdlib.h>
/*************************************************************************




 ************************************************************************/

#include <iostream>
#include "InputParser.h"

InputParser::InputParser (int &argc, char **argv){
    for (int i=1; i < argc; ++i)
        this->tokens.push_back(std::string(argv[i]));
}

std::string InputParser::getCmdOption( std::string &option) {
    std::vector<std::string>::iterator itr;
    itr =  std::find(this->tokens.begin(), this->tokens.end(), option);
    if (itr != this->tokens.end() && ++itr != this->tokens.end()){
        return *itr;
    }
    return "";
}

std::string InputParser::getCmdOption( const char* o) {
    std::string option = o;
    return getCmdOption(option);
}

bool InputParser::cmdOptionExists(std::string &option) {
    return std::find(this->tokens.begin(), this->tokens.end(), option)
           != this->tokens.end();
}

bool InputParser::cmdOptionExists( const char* o){
    std::string option = o;
    return cmdOptionExists(option);
}

void usage( ){
    std::string progName = "AnnotationLiftOver";
    std::cout << "Program " << progName << std::endl <<
    "Usage:  "<<progName<<" <command> [options]"<< std::endl <<
    "Commands:"<< std::endl <<
        "getGenomeSequence              create fasta file using reference genome sequence and variant calling records" << std::endl <<
        "coordinateLiftOver             get the corresponding coordinate of reference coordinate at another accession" << std::endl <<
        "gffCoordinateLiftOver          transform the reference GFF/GTF coordinate to the coordinate of another accession," << std::endl <<
        "                               don't care about the completeness of start/stop codon, splice sites" << std::endl<<
        "getSequences                   get the protein/CDS/gene sequence with genome sequence file and GFF/GTF file" << std::endl<<
//        "annotationLiftOver             transform the reference GFF/GTF coordinate to the coordinate of another accession," << std::endl <<
//        "                               complementing with trying to keep the ORF/splice sites and complete as possible by " << std::endl <<
//        "                               genome sequence alignment" << std::endl <<
//        "annotationLiftOverAndOrth      transform the reference GFF/GTF coordinate to the coordinate of another accession," << std::endl <<
//        "                               complementing with trying to keep the ORF/splice sites and complete as possible by " << std::endl <<
//        "                               genome sequence alignment and then complementing with aligning CDS sequence of " << std::endl <<
//        "                               reference to the genome sequence of target line and then complementing with " << std::endl <<
//        "                               aligning protein sequence of reference to the genome sequence of target line" << std::endl <<
        "annotationAndExonerateAndNovo  transform the reference GFF/GTF coordinate to the coordinate of another accession," << std::endl <<
        "                               complementing with trying to keep the ORF/splice sites and complete as possible by"  << std::endl <<
        "                               genome sequence alignment. And then complementing with aligning CDS sequence of " << std::endl <<
        "                               reference to the genome sequence of target line and then complementing with"  << std::endl <<
        "                               aligning protein sequence of reference to the genome sequence of target line" << std::endl <<
        "                               Finally the missing gene annotation would be complemented with" << std::endl <<
        "                               other annotation inputted " << std::endl <<
                                        std::endl;
}
