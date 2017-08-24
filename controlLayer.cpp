/*
 * =====================================================================================
 *
 *       Filename:  controlLayer.cpp
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  06/25/2017 09:38:26
 *       Revision:  none
 *       Compiler:  gcc
 *
 *        Author:  Baoxing Song (songbx.me), song@mpipz.mpg.de
 *   Organization:  MPIPZ
 *
 * =====================================================================================
 */

/*************************************************************************




 ************************************************************************/

#include <iostream>
#include "InputParser.h"
#include <string>
#include <sstream>
#include "myfunctions.h"
int getGenomeSequence(int argc, char** argv){
    std::stringstream usage;
    usage <<  "Usage:    * getGenomeSequence -r reference -v variants -o output" << std::endl<<
        "Options" << std::endl <<
        "   -h        produce help message" << std::endl <<
        "   -r        reference genome in fasta format" << std::endl <<
        "   -v        variant calling result in vcf/sdi format" << std::endl <<
        "   -o        the output pseudo genome sequence in fasta format" << std::endl;
    InputParser inputParser (argc, argv);
    if(inputParser.cmdOptionExists("-h") ||inputParser.cmdOptionExists("--help")  ){
        std::cerr << usage.str();
    }else if( inputParser.cmdOptionExists("-r") && inputParser.cmdOptionExists("-v") && inputParser.cmdOptionExists("-o")  ){
        std::string reference = inputParser.getCmdOption("-r");
        std::string variants = inputParser.getCmdOption("-v");
        std::string output = inputParser.getCmdOption("-o");
        getPseudoGenomeSequence(reference, variants, output);
        return 0;
    }else{
        std::cerr << usage.str();
    }
    return 1;
}

int coordinateLiftOver( int argc, char** argv ){
    std::stringstream usage;
    usage <<  "Usage:    * coordinateLiftOver -v variants -c chromosome -p position " << std::endl<<
        "Options" << std::endl <<
        "   -h        produce help message" << std::endl <<
        "   -v        variant calling result in vcf/sdi format" << std::endl <<
        "   -c        chromosome, should be consistent with the chromosome information in sdi file (The coordinate starts from 1)" << std::endl <<
        "   -p        the position/coordinate in reference genome" <<std::endl;
    InputParser inputParser (argc, argv);
    if(inputParser.cmdOptionExists("-h") ||inputParser.cmdOptionExists("--help")  ){
        std::cerr << usage.str();
    }else if( inputParser.cmdOptionExists("-v") && inputParser.cmdOptionExists("-c") && inputParser.cmdOptionExists("-p")  ){
        std::string variantsFile = inputParser.getCmdOption("-v");
        std::string chromosome = inputParser.getCmdOption("-c");
        std::string coordinateS = inputParser.getCmdOption("-p");
        int coordinate = std::stoi(coordinateS);
        int liftCoordinate = myCoordinateLiftOver( variantsFile, chromosome, coordinate  );
        std::cout << chromosome << " " << liftCoordinate << std::endl;
        return 0;
    }else{
        std::cerr << usage.str();
    }
    return 1;
}
int gffCoordinateLiftOver( int argc, char** argv  ){
    std::stringstream usage;
    usage << "simplely coordinate liftover, gffCoordinateLiftOver"  <<std::endl <<
        "Usage:    * gffCoordinateLiftOver -v variants -i inputGffFile -o outputGffFile " << std::endl<<
        "Options" << std::endl <<
        "   -h        produce help message" << std::endl <<
        "   -v        variant calling result in vcf/sdi format" << std::endl <<
        "   -i        the input GFF/GTF file of reference line/accession" << std::endl <<
        "   -o        the output GFF/GTF file of target line/accession" << std::endl;
    InputParser inputParser (argc, argv);
    if(inputParser.cmdOptionExists("-h") ||inputParser.cmdOptionExists("--help")  ){
        std::cerr << usage.str();
    }else if( inputParser.cmdOptionExists("-v") && inputParser.cmdOptionExists("-i") && inputParser.cmdOptionExists("-o")   ){
        std::string variantsFile = inputParser.getCmdOption("-v");
        std::string inputGff = inputParser.getCmdOption("-i");
        std::string outputGff = inputParser.getCmdOption("-o");
        myGffCoordinateLiftOver(variantsFile, inputGff, outputGff);
        return 0;
    }else{
        std::cerr << usage.str();
    }
    return 1;
}
int getSequences( int argc, char** argv   ){
    std::stringstream usage;
    usage <<  "Usage:    * getSequences -i inputGffFile -r inputGenome -p outputProteinSequences -c outputCdsSequences -g outputGenomeSequences " << std::endl<<
        "   -h        produce help message" << std::endl <<
        "   -i        reference genome in GFF/GTF format" << std::endl <<
        "   -r        genome sequence in fasta format" << std::endl <<
        "   -p        output file of protein sequence in fasta format" << std::endl <<
        "   -c        output file of protein coding sequence (without intron) in fasta format" << std::endl <<
        "   -g        output file of protein coding sequence (with intron) in fasta frormat" << std::endl <<
        "   -x        default: ([\\s\\S]*)Parent=([\\s\\S]*?)[;,][\\s\\S]*$ it works for parese the TAIR10 annotation" << std::endl<<
        "                 regex to parser the structure of CDS elements and parent transcript in GFF/GTF file" << std::endl;
    InputParser inputParser (argc, argv);
    if(inputParser.cmdOptionExists("-h") ||inputParser.cmdOptionExists("--help")  ){
        std::cerr << usage.str();
    }else if( inputParser.cmdOptionExists("-i") && inputParser.cmdOptionExists("-r") && inputParser.cmdOptionExists("-p") && inputParser.cmdOptionExists("-c") && inputParser.cmdOptionExists("-g") ){
        std::string inputGffFile = inputParser.getCmdOption("-i");
        std::string genome = inputParser.getCmdOption("-r");
        std::string outputProteinSequences = inputParser.getCmdOption("-p");
        std::string outputCdsSequences = inputParser.getCmdOption("-c");
        std::string outputGenomeSequences = inputParser.getCmdOption("-g");
        std::string regex;
        if( inputParser.cmdOptionExists("-x") ){
            regex = inputParser.getCmdOption("-x");
        }else{
            regex = "([\\s\\S]*)Parent=([\\s\\S]*?)[;,][\\s\\S]*$";
        }
        getSequences( inputGffFile, genome, outputProteinSequences, outputCdsSequences, outputGenomeSequences, regex );
        return 0;
    }else{
        std::cerr << usage.str();
    }
    return 1;
}

int annotationLiftOver( int argc, char** argv ){
    std::stringstream usage;
    usage <<  "Usage:    * annotationLiftOver -i inputGffFile -r referenceGenomeSequence -v variants -o outputGffFile" << std::endl <<
        "   -h        produce help message" << std::endl <<
        "   -i        the input GFF/GTF file of reference line/accession" << std::endl <<
        "   -r        reference genome in fasta format" << std::endl <<
        "   -v        variant calling result in vcf/sdi format" << std::endl <<
        "   -o        the output GFF/GTF file of target line/accession" << std::endl <<
        "   -x        default: ([\\s\\S]*)Parent=([\\s\\S]*?)[;,][\\s\\S]*$  it works for parese the TAIR10 annotation" << std::endl<<
        "                 regex to parser the structure of CDS elements and parent transcript in GFF/GTF file" << std::endl <<
        "   -g        default: (.+?)\\\\.   it works for parese the TAIR10 annotation" << std::endl <<
        "             regex to parser the parent gene id from transcript id"  << std::endl <<
        "   -t        (int) number of threads, default: 4 "<< std::endl <<
        "   -l        (int) if sequence is longer than this threshold would be aligned, for RAM and time saving porpose  "<< std::endl;
    InputParser inputParser (argc, argv);
    if(inputParser.cmdOptionExists("-h") ||inputParser.cmdOptionExists("--help")  ){
        std::cerr << usage.str();
        return 1;
    } else if( inputParser.cmdOptionExists("-i") && inputParser.cmdOptionExists("-r") && inputParser.cmdOptionExists("-v") && inputParser.cmdOptionExists("-o") ){
        std::string inputGffFile = inputParser.getCmdOption("-i");
        std::string referenceGenomeSequence = inputParser.getCmdOption("-r");
        std::string variants = inputParser.getCmdOption("-v");
        std::string outputGffFile = inputParser.getCmdOption("-o");
        std::string regex;
        if( inputParser.cmdOptionExists("-x")  ){
            regex = inputParser.getCmdOption("-x");
        }else{
            regex = "([\\s\\S]*)Parent=([\\s\\S]*?)[;,][\\s\\S]*$";
        }
        std::string regexG;
        if( inputParser.cmdOptionExists("-g") ){
            regexG = inputParser.getCmdOption("-g");
        }else{
            regexG = "(.+?)\\.";
        }
        int threads;
        if( inputParser.cmdOptionExists("-t") ){
            threads = std::stoi( inputParser.getCmdOption("-t") );
        }else{
            threads = 4;
        }
        int lengthThread;
        if( inputParser.cmdOptionExists("-l") ){
            lengthThread = std::stoi( inputParser.getCmdOption("-l") );
        }else{
            lengthThread=100000;
        }
        myReAnnotationLiftoverSingleLine( referenceGenomeSequence, inputGffFile, variants,
                                          outputGffFile, regex, threads, regexG, lengthThread );
        return 0;
    } else{
        std::cerr << usage.str();
        return 1;
    }
    return 1;
}

int annotationLiftOverAndOrth( int argc, char** argv ){
    std::stringstream usage;
    usage <<  "Usage:    * annotationLiftOverAndOrth -i inputGffFile -r referenceGenomeSequence -v variants -o outputGffFile" << std::endl <<
          "   -h        produce help message" << std::endl <<
          "   -i        the input GFF/GTF file of reference line/accession" << std::endl <<
          "   -r        reference genome in fasta format" << std::endl <<
          "   -v        variant calling result in vcf/sdi format" << std::endl <<
          "   -o        the output GFF/GTF file of target line/accession" << std::endl <<
          "   -x        default: ([\\s\\S]*)Parent=([\\s\\S]*?)[;,][\\s\\S]*$  it works for parese the TAIR10 annotation" << std::endl<<
          "                 regex to parser the structure of CDS elements and parent transcript in GFF/GTF file" << std::endl <<
          "   -g        default: (.+?)\\\\.   it works for parese the TAIR10 annotation" << std::endl <<
          "             regex to parser the parent gene id from transcript id"  << std::endl <<
          "   -t        (int) number of threads, default: 4 "<< std::endl <<
          "   -l        (int) if sequence is longer than this threshold would be aligned, for RAM and time saving porpose  "<< std::endl;
    InputParser inputParser (argc, argv);
    if(inputParser.cmdOptionExists("-h") ||inputParser.cmdOptionExists("--help")  ){
        std::cerr << usage.str();
        return 1;
    } else if( inputParser.cmdOptionExists("-i") && inputParser.cmdOptionExists("-r") && inputParser.cmdOptionExists("-v") && inputParser.cmdOptionExists("-o") ){
        std::string inputGffFile = inputParser.getCmdOption("-i");
        std::string referenceGenomeSequence = inputParser.getCmdOption("-r");
        std::string variants = inputParser.getCmdOption("-v");
        std::string outputGffFile = inputParser.getCmdOption("-o");
        std::string regex;
        if( inputParser.cmdOptionExists("-x")  ){
            regex = inputParser.getCmdOption("-x");
        }else{
            regex = "([\\s\\S]*)Parent=([\\s\\S]*?)[;,][\\s\\S]*$";
        }
        std::string regexG;
        if( inputParser.cmdOptionExists("-g") ){
            regexG = inputParser.getCmdOption("-g");
        }else{
            regexG = "(.+?)\\.";
        }
        int threads;
        if( inputParser.cmdOptionExists("-t") ){
            threads = std::stoi( inputParser.getCmdOption("-t") );
        }else{
            threads = 4;
        }
        int lengthThread;
        if( inputParser.cmdOptionExists("-l") ){
            lengthThread = std::stoi( inputParser.getCmdOption("-l") );
        }else{
            lengthThread=30000;
        }
        myReAnnotationLiftoverAndOrthologous( referenceGenomeSequence, inputGffFile, variants,
                                          outputGffFile, regex, threads, regexG , lengthThread);
        return 0;
    } else{
        std::cerr << usage.str();
        return 1;
    }
    return 1;
}


int reAnnotationAndExonerateAndNovo( int argc, char** argv ){
    std::stringstream usage;
    usage <<  "Usage:    * AnnotationAndExonerateAndNovo -i inputGffFile -r referenceGenomeSequence -v variants -o outputGffFile" << std::endl <<
          "   -h        produce help message" << std::endl <<
          "   -i        the input GFF/GTF file of reference line/accession" << std::endl <<
          "   -n        the de novo annotation GFF of the target accession" << std::endl <<
          "   -r        reference genome in fasta format" << std::endl <<
          "   -v        variant calling result in vcf/sdi format" << std::endl <<
          "   -o        the output GFF/GTF file of target line/accession" << std::endl <<
          "   -x        default: ([\\s\\S]*)Parent=([\\s\\S]*?)[;,][\\s\\S]*$  it works for parese the TAIR10 annotation" << std::endl<<
          "                 regex to parser the structure of CDS elements and parent transcript in GFF/GTF file" << std::endl <<
          "   -g        default: ([\\s\\S]+?)\\\\.   it works for parese the TAIR10 annotation" << std::endl <<
          "             regex to parser the parent gene id from transcript id"  << std::endl <<
          "   -nx        default: ID=([\\s\\S]*?);Parent=([\\s\\S]*?)$  it works for parese the output of AUGUSTUS" << std::endl<<
          "                 regex to parser the structure of CDS elements and parent transcript in de novo GFF/GTF file" << std::endl <<
          "   -ng        default: ([\\s\\S]+?)\\\\.   it works for parese the output of AUGUSTUS" << std::endl <<
          "             regex to parser the parent gene id from transcript id"  << std::endl <<
          "   -t        (int) number of threads, default: 4 "<< std::endl <<
          "   -l        (int) if sequence is longer than this threshold would be aligned, for RAM and time saving purpose  "<< std::endl;
    InputParser inputParser (argc, argv);
    if(inputParser.cmdOptionExists("-h") ||inputParser.cmdOptionExists("--help")  ){
        std::cerr << usage.str();
        return 1;
    } else if( inputParser.cmdOptionExists("-i") && inputParser.cmdOptionExists("-r") && inputParser.cmdOptionExists("-n") && inputParser.cmdOptionExists("-v") && inputParser.cmdOptionExists("-o") ){
        std::string referenceGenomeSequence = inputParser.getCmdOption("-r");
        std::string inputGffFile = inputParser.getCmdOption("-i");
        std::string novoGffFilePath = inputParser.getCmdOption("-n");
        std::string variants = inputParser.getCmdOption("-v");
        std::string outputGffFile = inputParser.getCmdOption("-o");
        std::string regex;
        if( inputParser.cmdOptionExists("-x")  ){
            regex = inputParser.getCmdOption("-x");
        }else{
            regex = "([\\s\\S]*)Parent=([\\s\\S]*?)[;,][\\s\\S]*$";
        }
        std::string regexG;
        if( inputParser.cmdOptionExists("-g") ){
            regexG = inputParser.getCmdOption("-g");
        }else{
            regexG = "([\\s\\S]+?)\\.";
        }
        int threads;
        if( inputParser.cmdOptionExists("-t") ){
            threads = std::stoi( inputParser.getCmdOption("-t") );
        }else{
            threads = 4;
        }
        std::string novoRegex;
        if( inputParser.cmdOptionExists("-nx") ){
            novoRegex = std::stoi( inputParser.getCmdOption("-nx") );
        }else{
            novoRegex = "ID=([\\s\\S]*?);Parent=([\\s\\S]*?)$";
        }
        std::string novoRegexG;
        if( inputParser.cmdOptionExists("-ng") ){
            novoRegexG = std::stoi( inputParser.getCmdOption("-ng") );
        }else{
            novoRegexG = "([\\s\\S]+?)\\.";
        }
        int lengthThread;
        if( inputParser.cmdOptionExists("-l") ){
            lengthThread = std::stoi( inputParser.getCmdOption("-l") );
        }else{
            lengthThread=50000;
        }
        myReAnnotationAndExonerateAndNovo( referenceGenomeSequence, inputGffFile, novoGffFilePath,
                                            variants, outputGffFile, regex, threads,
                                            regexG, novoRegex, novoRegexG , lengthThread);
        return 0;
    } else{
        std::cerr << usage.str();
        return 1;
    }
    return 1;
}

int myCountNumberOfTwoneighborSNP( int argc, char** argv  ){
    std::stringstream usage;
    usage <<  "Usage:    * countNumberOfTwoneighborSNP -v variants" << std::endl <<
        "   -h        produce help message" << std::endl <<
        "   -v        variant calling result in vcf/sdi format" << std::endl <<
        "   -l        (int) number of bases in which range there should be no variants for double SNP" << std::endl <<
        "   -o        prefix of output file" <<std::endl;
    InputParser inputParser (argc, argv);
    if(inputParser.cmdOptionExists("-h") ||inputParser.cmdOptionExists("--help")  ){
        std::cerr << usage.str();
        return 1;
    } else if( inputParser.cmdOptionExists("-v") && inputParser.cmdOptionExists("-o") ) {
        std::string variants = inputParser.getCmdOption("-v");
        std::string outputPrefix = inputParser.getCmdOption("-o");
        int rangeLength;
        if( inputParser.cmdOptionExists("-l") ){
            rangeLength = std::stoi( inputParser.getCmdOption("-l") );
        }else{
            rangeLength = 3;
        }
        countNumberOfTwoneighborSNP(variants, outputPrefix, rangeLength);
        return 0;
    }else {
        std::cerr << usage.str();
        return 1;
    }
}


int mycountNumberSNPAndIndel( int argc, char** argv  ){
    std::stringstream usage;
    usage <<  "Usage:    * countNumberSNPAndIndel -v variants" << std::endl <<
          "   -h        produce help message" << std::endl <<
          "   -v        variant calling result in vcf/sdi format" << std::endl <<
          "   -l        (int) number of bases in which range there should be no INDEL" << std::endl <<
          "   -o        prefix of output file" <<std::endl;
    InputParser inputParser (argc, argv);
    if(inputParser.cmdOptionExists("-h") ||inputParser.cmdOptionExists("--help")  ){
        std::cerr << usage.str();
        return 1;
    } else if( inputParser.cmdOptionExists("-v") && inputParser.cmdOptionExists("-o") ) {
        std::string variants = inputParser.getCmdOption("-v");
        std::string outputPrefix = inputParser.getCmdOption("-o");
        int rangeLength;
        if( inputParser.cmdOptionExists("-l") ){
            rangeLength = std::stoi( inputParser.getCmdOption("-l") );
        }else{
            rangeLength = 3;
        }
        countNumberSNPAndIndel(variants, outputPrefix, rangeLength);
        return 0;
    }else {
        std::cerr << usage.str();
        return 1;
    }
}

int myGenerateRandomSdi( int argc, char** argv  ){
    std::stringstream usage;
    usage <<  "Usage:    * generateRandomSdi -v variants" << std::endl <<
          "   -h        produce help message" << std::endl <<
          "   -v        variant calling result in vcf/sdi format" << std::endl <<
          "   -r        reference genome in fasta format" << std::endl <<
          "   -o        prefix of output file" <<std::endl;
    InputParser inputParser (argc, argv);
    if(inputParser.cmdOptionExists("-h") ||inputParser.cmdOptionExists("--help")  ){
        std::cerr << usage.str();
        return 1;
    } else if( inputParser.cmdOptionExists("-v") && inputParser.cmdOptionExists("-o") && inputParser.cmdOptionExists("-r") ) {
        std::string variants = inputParser.getCmdOption("-v");
        std::string outputPrefix = inputParser.getCmdOption("-o");
        std::string reference = inputParser.getCmdOption("-r");
        generateRandomSdi(variants, reference, outputPrefix);
        return 0;
    }else {
        std::cerr << usage.str();
        return 1;
    }
}