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
 *         Author:  YOUR NAME (), 
 *   Organization:  
 *
 * =====================================================================================
 */

/*************************************************************************




 ************************************************************************/



#include <fstream>
#include <iostream>
#include <string>
#include <cstdlib>

#include "alignNeedlemanWunsch.h"

// Usage: ./program dna_file1 dna_file2
int main(int argc, char **argv) 
{
	// File Parsing.
	std::fstream compare_one, compare_two;
	if (argc != 3) 
	{
		std::cout << "Warning: Using default values! \n";
		compare_one.open("mouse-hemoglobin-sequence.fasta");
		compare_two.open("human-hemoglobin-sequence.fasta");
	}
	else
	{
		compare_one.open(argv[1]);
		compare_two.open(argv[2]);
		std::cout << argv[1] << std::endl;
		std::cout << argv[2] << std::endl;
	}

	
	if (!compare_one || !compare_two)
	{
		std::cerr << "Error (2): Necessary files not found! \n";
		exit(2);
	}
	
	std::string dna_a = "", dna_b = "", line;
	for (int i = 0; std::getline(compare_one, line); i++)
	{
		if (i != 0) 
		{
			dna_a += line;
		}
	}

	for (int i = 0; std::getline(compare_two, line); i++)
	{
		if (i != 0)
		{
			dna_b += line;
		}
	}
	//std::cout << "77" <<std::endl;
	// Application of algorithm.
	NeedlemanWunsch *nw = new NeedlemanWunsch(dna_a, dna_b, 1, -1, -2, -1);
	//std::cout << "80" <<std::endl;
	nw->calculate_similarity();
	//std::cout << "82" <<std::endl;
	nw->dna_align();
	//std::cout << "84" <<std::endl;
	nw->print_results();
	//std::cout << "86" <<std::endl;
	//getchar();
	return EXIT_SUCCESS;
}
