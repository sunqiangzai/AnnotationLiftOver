/*
 * =====================================================================================
 *
 *       Filename:  testAlignNeedlemanWunsch.cpp
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  06/02/2017 17:48:28
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Baoxing Song, song@mpipz.mpg.de
 *   Organization:  MPIPZ
 *
 * =====================================================================================
 */

/*************************************************************************




 ************************************************************************/



#include <fstream>
#include <iostream>
#include <string>
#include <cstdlib>

#include "alignNeedlemanForTranscript.h"

// Usage: ./program dna_file1 dna_file2
int main(int argc, char **argv) 
{
//	// File Parsing.
//	std::fstream compare_one, compare_two;
//	if (argc != 3)
//	{
//		std::cout << "Warning: Using default values! \n";
//		compare_one.open("mouse-hemoglobin-sequence.fasta");
//		compare_two.open("human-hemoglobin-sequence.fasta");
//	}
//	else
//	{
//		compare_one.open(argv[1]);
//		compare_two.open(argv[2]);
//		std::cout << argv[1] << std::endl;
//		std::cout << argv[2] << std::endl;
//	}
//
//
//	if (!compare_one || !compare_two)
//	{
//		std::cerr << "Error (2): Necessary files not found! \n";
//		exit(2);
//	}
//
//	std::string dna_a = "", dna_b = "", line;
//	for (int i = 0; std::getline(compare_one, line); i++)
//	{
//		if (i != 0)
//		{
//			dna_a += line;
//		}
//	}
//
//	for (int i = 0; std::getline(compare_two, line); i++)
//	{
//		if (i != 0)
//		{
//			dna_b += line;
//		}
//	}
	std::string dna_a = "ATGCTTATCAGTATCAGCCCACTGATATTTGTGATACCAGTATCATCTGACGTGGCTTCTTCTGATTGGTTACATTTGACAAAAGCAAAAAATATTATATATATTTATTAA";
	std::string dna_b = "CAGTAGTGAATGTTGCAATCAAGAAATAGTTATAAAACCATGTACATGGGTGATATTTTTAGATTATCACAAGTACAAAAATTATTTTTATATTAGTAAAGCAATTAGAATATGCTTATCAGTATCAGCCCACTGAATATTTGTGATACCAGTATCATCTGACGTGGCTTCTTCTGATTGGTTACATTTGACAAAAGCAAAAAATATTATATATATTTATTAAACAAACATTATAAATAATACTAAAAAATATTCTTGATGATATTTTGTGAGTATCACCTATGAACATGCTCTAAAGATCAATATCTTAAAACAAATTTAGTGATATAATTAA";
    //std::cout << "77" <<std::endl;
	// Application of algorithm.
    int startCodonPosition=1;
    int stopCodonPosition=109;
    std::vector<SpliceSitePosition> spliceSitePositions;
//    SpliteSitePosition spliteSitePosition(321, 464);
//    splitSitePositions.push_back(spliteSitePosition);
//    SpliteSitePosition spliteSitePosition2(752, 874);
//    splitSitePositions.push_back(spliteSitePosition2);
//	NeedlemanWunschForTranscript *nw = new NeedlemanWunschForTranscript(dna_a, dna_b,  6, -1, -4, -2, 1, -1, -2, -1, 10,
//																		-10, -10, -10, 10, -10, -10, -10, startCodonPosition, stopCodonPosition, spliceSitePositions);

    NeedlemanWunschForTranscript *nw = new NeedlemanWunschForTranscript(dna_a, dna_b, startCodonPosition, stopCodonPosition, spliceSitePositions);
	//std::cout << "80" <<std::endl;
	//nw->calculate_similarity();
	//std::cout << "82" <<std::endl;
	//nw->dna_align();
	//std::cout << "84" <<std::endl;
	nw->print_results();
	//std::cout << "86" <<std::endl;
	//getchar();
	return EXIT_SUCCESS;
}
