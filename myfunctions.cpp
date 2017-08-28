/*
 * =====================================================================================
 *
 *       Filename:  myfunctions.cpp
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  05/31/2017 09:51:42
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Baoxing Song (songbx.me), song@mpipz.mpg.de
 *   Organization:  MPIPZ
 *
 * =====================================================================================
 */
#include <cstdlib>
/*************************************************************************




 ************************************************************************/

#include<iostream>
#include<map>
#include "myfunctions.h"
#include <fstream>
#include <sstream>
#include <regex>
void hereOutPutLiftOrOrthologousResult(std::map<std::string, std::vector<Gene> >& genes, std::string& outputGffFile );


int getPseudoGenomeSequence(std::map<std::string, Fasta>& referenceSequences, std::map<std::string, std::vector<Variant> >& variantsMap,
                            std::map<std::string, Fasta>& targetSequences) {
    std::set<std::string> chromosomes;
    for (std::map<std::string, Fasta>::iterator iterFasta = referenceSequences.begin();
         iterFasta != referenceSequences.end(); iterFasta++) {
        chromosomes.insert(iterFasta->first);
    }
    return getPseudoGenomeSequence(referenceSequences, variantsMap, targetSequences, chromosomes);
}
int getPseudoGenomeSequence(std::map<std::string, Fasta>& referenceSequences, std::map<std::string, std::vector<Variant> >& variantsMap,
                            std::map<std::string, Fasta>& targetSequences, std::string& chromosome){
    std::set<std::string> chromosomes;
    chromosomes.insert(chromosome);
    return getPseudoGenomeSequence(referenceSequences, variantsMap, targetSequences, chromosomes);
}
int getPseudoGenomeSequence(std::map<std::string, Fasta>& referenceSequences, std::map<std::string, std::vector<Variant> >& variantsMap,
                            std::map<std::string, Fasta>& targetSequences, std::set<std::string>& chromosomes){
    for( std::map<std::string, Fasta>::iterator iterFasta=referenceSequences.begin(); iterFasta!=referenceSequences.end(); iterFasta++ ) {
        std::string chromosome = iterFasta->first;
        if( chromosomes.find(chromosome)!=chromosomes.end() ){
            std::string referenceSequence = iterFasta->second.getSequence();
            std::stringstream sequencestream;
            int currentPosition = 1;
            for (std::vector<Variant>::iterator iterSdi = variantsMap[chromosome].begin(); iterSdi != variantsMap[chromosome].end(); iterSdi++) {
                int position = (*iterSdi).getPosition();
                std::string ref = (*iterSdi).getReference();
                std::string alter = (*iterSdi).getAlternative();
                int changingLength = (*iterSdi).getChanginglength();
                if (ref != "-") {
                    if (referenceSequence.substr(position - 1, ref.size()).compare(ref) !=0) {
                        std::cerr << "the record does not confirm with reference sequence at: " << chromosome << "\t"
                                  << position << std::endl;
                        exit(1);
                    }
                }
                if (currentPosition < position) {
                    sequencestream << referenceSequence.substr(currentPosition - 1, position - currentPosition);
                } else if( currentPosition == position+1 && changingLength>0 && ref.compare("-")==0 && (*(iterSdi-1)).getChanginglength()==0 && (*(iterSdi-1)).getReference().size()==1 ){
                    currentPosition = position+1;
                    sequencestream.seekp(-1,sequencestream.cur);
                    sequencestream << alter << (*(iterSdi-1)).getAlternative();
                    continue;
                }else if (currentPosition > position) {
                    std::cerr << "the sdi file is not well sorted, it should be sorted with coordinate " <<std::endl
                              << chromosome << ": currentPosition:"<< currentPosition << " position:" << position << std::endl;
                    exit(1);
                }

                if (changingLength == 0) {
                    currentPosition = position + ref.size();
                    sequencestream << alter;
                }
                if (changingLength > 0) {
                    if ( ref.compare("-")==00 ) {
                        sequencestream << alter;
                        currentPosition = position;
                    } else {
                        sequencestream << alter;
                        currentPosition = position + ref.size();
                    }
                }
                if (changingLength < 0) {
                    if (alter.compare("-")==0 ) {
                        currentPosition = position - changingLength;
                    } else {
                        sequencestream << alter;
                        currentPosition = position + ref.size();
                    }
                }
            }
            int totalSize = referenceSequence.size();
            sequencestream << referenceSequence.substr(currentPosition - 1, totalSize - currentPosition + 1);
            std::string targetSequence = sequencestream.str();
            Fasta fasta(chromosome, targetSequence);
            targetSequences[chromosome]=fasta;
        }
    }
    return 0;
}
int getPseudoGenomeSequence(std::map<std::string, Fasta>& referenceSequences, std::string& sdiFile,
                            std::map<std::string, Fasta>& targetSequences, std::string& vcfFix){
    std::map<std::string, std::vector<Variant> > sdiMaps;
    readSdiFile(sdiFile, sdiMaps, vcfFix);
    return getPseudoGenomeSequence(referenceSequences, sdiMaps, targetSequences);
}
int getPseudoGenomeSequence(std::string& referenceGenomeFastaFile,
                            std::string& sdiFile, std::map<std::string, Fasta>& targetSequences, std::string& vcfFix){
    std::map<std::string, Fasta> referenceSequences;
    readFastaFile(referenceGenomeFastaFile, referenceSequences);
    return getPseudoGenomeSequence(referenceSequences, sdiFile, targetSequences, vcfFix);
}
int getPseudoGenomeSequence(std::string& referenceGenomeFastaFile,
                            std::string& sdiFile, std::string& outputFile, std::string& vcfFix){
    std::map<std::string, Fasta> targetSequences;
    int resultCode = getPseudoGenomeSequence(referenceGenomeFastaFile, sdiFile, targetSequences, vcfFix);
    std::ofstream ofile;
    ofile.open(outputFile);
    for( std::map<std::string, Fasta>::iterator i=targetSequences.begin() ; i!= targetSequences.end(); i ++ ){
        std::string name = i->first;
        std::string sequence = i->second.getSequence();
        int linewidth=60;
        writeFasta(ofile, name, sequence,  linewidth);
    }
    ofile.close();
    return resultCode;
}
void writeFasta(std::ostream& out, std::string& seqname, std::string& sequence, int& linewidth) {
    out << ">" << seqname << std::endl;
    std::size_t pos = 0;
    while (pos < sequence.length()) {
        out << sequence.substr(pos, linewidth) << std::endl;
        pos += linewidth;
    }
}

int myCoordinateLiftOver( std::string& sdiFile, std::string& chromosome, int& position, std::string& vcfFix ){
    std::map<std::string, std::vector<Variant> > sdiMaps;
    readSdiFile(sdiFile, sdiMaps, vcfFix);
    return getChangedFromBasement( chromosome, position, sdiMaps );
}

void myGffCoordinateLiftOver( std::string& sdiFile, std::string& gffFile, std::string& outputFile, std::string& vcfFix ){
    std::map<std::string, std::vector<Variant> > sdiMaps;
    readSdiFile(sdiFile, sdiMaps, vcfFix);

    std::ofstream ofile;
    ofile.open(outputFile);
    std::ifstream infile(gffFile);

    //std::regex reg("^([\\s\\S]+)$");
    std::regex reg("^([\\s\\S]*?)\t([\\s\\S]*?)\t([\\s\\S]*?)\t(\\S*?)\t(\\S*?)\t([\\s\\S]+)$");
    std::string line;
    while (std::getline(infile, line)){
        std::smatch match;
        regex_search(line, match, reg);
        //if( match.empty() || line[0]=='#'){
        if( match.empty() ){
            ofile << line << std::endl;
            //std::cout << "not match: " << line << std::endl;
        }else{
            int start = std::stoi(match[4]);
            int end = std::stoi(match[5]);
            int startL = getChangedFromBasement( match[1], start, sdiMaps );
            int endL = getChangedFromBasement( match[1], end, sdiMaps );
            ofile << match[1] << "\t" << match[2] << "\t" << match[3] <<
                  "\t" << startL << "\t" << endL <<"\t" << match[6] << std::endl;
        }

//        std::vector<std::string> splits = split(line, '\\s');
//        if( splits.size() <9 || line[0]=='#' ){
//            ofile << line << std::endl;
//            std::cout << line << std::endl;
//        }else{
//            int start = std::stoi(splits[3]);
//            int end = std::stoi(splits[4]);
//            int startL = getChangedFromBasement( splits[0], start, sdiMaps );
//            int endL = getChangedFromBasement( splits[0], end, sdiMaps );
//            ofile << splits[0] << "\t" << splits[1] << "\t" << splits[2] <<
//                  "\t" << startL << "\t" << endL <<"\t" << splits[5] << "\t" <<splits[6] << "\t"<<
//                splits[7] << "\t" << splits[8] << std::endl;
//        }
    }
    ofile.close();
}
void getSequences(std::string& gffFile, std::string& genome, std::string& outputProteinSequences,
                  std::string& outputCdsSequences, std::string& outputGenomeSequences){
    std::string cdsParentRegex="([\\s\\S]*)Parent=([\\s\\S]*?)[;,][\\s\\S]*$";
    getSequences(gffFile, genome, outputProteinSequences,
                 outputCdsSequences, outputGenomeSequences, cdsParentRegex);
}
void getSequences(std::string& gffFile, std::string& genomeFile, std::string& outputProteinSequences,
                  std::string& outputCdsSequences, std::string& outputGenomeSequences, std::string& regex){
    std::map<std::string, std::vector<Transcript> > transcriptHashSet;
    readGffFile (gffFile, transcriptHashSet, regex);
    std::map<std::string, Fasta> genome;
    readFastaFile(genomeFile, genome);
    std::ofstream oPfile;
    oPfile.open(outputProteinSequences);
    std::ofstream oCfile;
    oCfile.open(outputCdsSequences);
    std::ofstream oGfile;
    oGfile.open(outputGenomeSequences);
    NucleotideCodeSubstitutionMatrix nucleotideCodeSubstitutionMatrix;
    for( std::map<std::string, std::vector<Transcript> >::iterator it1=transcriptHashSet.begin();
            it1!=transcriptHashSet.end(); ++it1 ){
        if( genome.find(it1->first) != genome.end() ){
            for ( std::vector<Transcript>::iterator it2=it1->second.begin(); it2!=it1->second.end(); ++it2 ) {
                (*it2).updateInfor(genome);
                checkOrfState( (*it2), genome, nucleotideCodeSubstitutionMatrix);

                std::string cdsSequence = (*it2).getCdsSequence();

                oPfile << ">" << (*it2).getName() << std::endl;
                oCfile << ">" << (*it2).getName() << std::endl;
                oGfile << ">" << (*it2).getName() << std::endl;

                oCfile << (*it2).getCdsSequence() << std::endl;
                oGfile << (*it2).getGeneomeSequence() << std::endl;

                std::string proteinSequence = nA2AA(cdsSequence, nucleotideCodeSubstitutionMatrix);
                oPfile << proteinSequence << std::endl;
                if( (*it2).getIfOrfShift() ){
                    std::cerr << (*it2).getName() << (*it2).getMetaInformation() << std::endl;
                }
            }
        }
    }
    oPfile.close();
    oCfile.close();
    oGfile.close();
}

void myReAnnotationLiftoverSingleLine( std::string& referenceGenomeFile, std::string& inputGffFile, 
        std::string& variantsFile, std::string& outputGffFile, std::string& regex, int &maxThread, std::string& regexG, int & lengthThread, std::string & vcfFix ){
    std::map<std::string, std::vector<Gene> > genes;
    reAnnotationSingleLine( referenceGenomeFile, inputGffFile, variantsFile, genes, maxThread, regex, regexG, lengthThread, vcfFix );
    hereOutPutLiftOrOrthologousResult(genes, outputGffFile);
}

void myReAnnotationLiftoverAndOrthologous( std::string& referenceGenomeFile, std::string& inputGffFile,
                                       std::string& variantsFile, std::string& outputGffFile, std::string& regex, int &maxThread, std::string& regexG, int& lengthThread, std::string & vcfFix  ){

    std::map<std::string, std::vector<Gene> > genes;
    reAnnotationAndExonerate( referenceGenomeFile, inputGffFile, variantsFile, genes, maxThread, regex, regexG, outputGffFile, lengthThread, vcfFix );
//    std::cout << "234" << std::endl;
    hereOutPutLiftOrOrthologousResult(genes, outputGffFile );
}

void myReAnnotationAndExonerateAndNovo( std::string& referenceGenomeFile, std::string& inputGffFile,std::string& novoGffFilePath,
                                         std::string& variantsFile, std::string& outputGffFile, std::string& regex, int &maxThread,
                                         std::string& regexG,std::string& novoRegex, std::string& novoRegexG , int& lengthThread, std::string & vcfFix){

    std::map<std::string, std::vector<Gene> > genes;
    reAnnotationAndExonerateAndNovo( referenceGenomeFile, inputGffFile, novoGffFilePath, variantsFile, genes, maxThread, regex, regexG, novoRegex, novoRegexG, outputGffFile, lengthThread, vcfFix);
    hereOutPutLiftOrOrthologousResult(genes, outputGffFile );
}



void hereOutPutLiftOrOrthologousResult(std::map<std::string, std::vector<Gene> >& genes, std::string& outputGffFile ){
    NucleotideCodeSubstitutionMatrix nucleotideCodeSubstitutionMatrix;
    std::ofstream ofile;
    ofile.open(outputGffFile);

    for( std::map<std::string, std::vector<Gene> >::iterator it=genes.begin(); it!=genes.end(); ++it ){
        std::sort(it->second.begin(), it->second.end(), compare_gene());
        for( std::vector<Gene>::iterator it2=it->second.begin(); it2!=it->second.end(); ++it2 ){
            std::string st = "+";

            if( NEGATIVE == (*it2).getStrand() ){
                st="-";
            }
            std::string geneResource = (*it2).getSource();
            if( geneResource.length()<1 ){
                geneResource="LIFTOVER";
            }
            ofile << it->first << "\t"+geneResource+"\tgene\t" << (*it2).getStart() << "\t" << (*it2).getEnd()
                  <<"\t.\t"<< st <<"\t.\tID="<< (*it2).getName()<< std::endl;
            for( std::vector<Transcript>::iterator it3 = (*it2).getTranscriptVector().begin();
                 it3!=(*it2).getTranscriptVector().end(); ++it3  ){

                std::string transcriptResource = (*it3).getSource();
                if( transcriptResource.length()<1 ){
                    transcriptResource="LIFTOVER";
                }

                ofile << it->first << "\t"+transcriptResource+"\tmRNA\t" << (*it3).getStart() << "\t" <<
                      (*it3).getEnd() << "\t.\t"<< st <<"\t.\tID="<< (*it3).getName() << ";Parent=" << (*it2).getName() << std::endl;
                //int cdsId = 1;
                for( std::vector<Cds>::iterator it4=(*it3).getCdsVector().begin();
                     it4!=(*it3).getCdsVector().end(); ++it4 ){
                    ofile << it->first << "\t"+transcriptResource+"\tCDS\t" << (*it4).getStart() << "\t" <<
                          (*it4).getEnd() << "\t.\t" <<st << "\t.\tParent=" << (*it3).getName()<< std::endl;
                   // ++cdsId;
                }
                ofile << "#metainformation: " << (*it3).getMetaInformation() << std::endl;
                ofile << "#genome sequence: " << (*it3).getGeneomeSequence() << std::endl;
                ofile << "#CDS sequence: " << (*it3).getCdsSequence() << std::endl;
                std::string cdsSequence = (*it3).getCdsSequence();
                //std::string proteinSequence = nA2AA(cdsSequence, nucleotideCodeSubstitutionMatrix );
                //ofile << "#protein sequence: " << proteinSequence << std::endl;
            }
        }
    }
    ofile.close();
}


std::ostream &printSdi(std::ostream& out, const Variant& variant){
    out << variant.getChromosome() << "\t" << variant.getPosition() << "\t" <<
        variant.getChanginglength() << "\t" << variant.getReference() << "\t" << variant.getAlternative() <<std::endl;
    return out;
}

bool caseInsensitiveStringCompare(const std::string& str1, const std::string& str2) {
    if (str1.size() != str2.size()) {
        return false;
    }
    for (std::string::const_iterator c1 = str1.begin(), c2 = str2.begin(); c1 != str1.end(); ++c1, ++c2) {
        if (tolower(*c1) != tolower(*c2)) {
            return false;
        }
    }
    return true;
}

void countNumberOfTwoneighborSNP( std::string& sdiFile, std::string & outputPrefix, int & rangeLength, std::string& vcfFix){

    NucleotideCodeSubstitutionMatrix nucleotideCodeSubstitutionMatrix;
    std::set<std::string>& legalNasString = nucleotideCodeSubstitutionMatrix.getLegalNasString();

    std::ofstream ofilealldouble;
    ofilealldouble.open(outputPrefix+".alldouble.sdi");

    std::ofstream ofileaacc;
    ofileaacc.open(outputPrefix+".aacc.sdi");

    std::ofstream ofileaacg;
    ofileaacg.open(outputPrefix+".aacg.sdi");

    std::ofstream ofileaact;
    ofileaact.open(outputPrefix+".aact.sdi");

    std::ofstream ofileaagc;
    ofileaagc.open(outputPrefix+".aagc.sdi");

    std::ofstream ofileaagg;
    ofileaagg.open(outputPrefix+".aagg.sdi");

    std::ofstream ofileaagt;
    ofileaagt.open(outputPrefix+".aagt.sdi");

    std::ofstream ofileaatc;
    ofileaatc.open(outputPrefix+".aatc.sdi");

    std::ofstream ofileaatg;
    ofileaatg.open(outputPrefix+".aatg.sdi");

    std::ofstream ofileaatt;
    ofileaatt.open(outputPrefix+".aatt.sdi");

    std::ofstream ofileacca;
    ofileacca.open(outputPrefix+".acca.sdi");

    std::ofstream ofileaccg;
    ofileaccg.open(outputPrefix+".accg.sdi");

    std::ofstream ofileacct;
    ofileacct.open(outputPrefix+".acct.sdi");

    std::ofstream ofileacga;
    ofileacga.open(outputPrefix+".acga.sdi");

    std::ofstream ofileacgg;
    ofileacgg.open(outputPrefix+".acgg.sdi");

    std::ofstream ofileacgt;
    ofileacgt.open(outputPrefix+".acgt.sdi");

    std::ofstream ofileacta;
    ofileacta.open(outputPrefix+".acta.sdi");

    std::ofstream ofileactg;
    ofileactg.open(outputPrefix+".actg.sdi");

    std::ofstream ofileagca;
    ofileagca.open(outputPrefix+".agca.sdi");

    std::ofstream ofileagcc;
    ofileagcc.open(outputPrefix+".agcc.sdi");

    std::ofstream ofileagct;
    ofileagct.open(outputPrefix+".agct.sdi");

    std::ofstream ofileagga;
    ofileagga.open(outputPrefix+".agga.sdi");

    std::ofstream ofileaggc;
    ofileaggc.open(outputPrefix+".aggc.sdi");

    std::ofstream ofileagta;
    ofileagta.open(outputPrefix+".agta.sdi");

    std::ofstream ofileagtc;
    ofileagtc.open(outputPrefix+".agtc.sdi");

    std::ofstream ofileatca;
    ofileatca.open(outputPrefix+".atca.sdi");

    std::ofstream ofileatcc;
    ofileatcc.open(outputPrefix+".atcc.sdi");

    std::ofstream ofileatcg;
    ofileatcg.open(outputPrefix+".atcg.sdi");

    std::ofstream ofileatga;
    ofileatga.open(outputPrefix+".atga.sdi");

    std::ofstream ofileatgc;
    ofileatgc.open(outputPrefix+".atgc.sdi");

    std::ofstream ofileatta;
    ofileatta.open(outputPrefix+".atta.sdi");

    std::ofstream ofilecagc;
    ofilecagc.open(outputPrefix+".cagc.sdi");

    std::ofstream ofilecagg;
    ofilecagg.open(outputPrefix+".cagg.sdi");

    std::ofstream ofilecatc;
    ofilecatc.open(outputPrefix+".catc.sdi");

    std::ofstream ofilecatg;
    ofilecatg.open(outputPrefix+".catg.sdi");

    std::ofstream ofileccga;
    ofileccga.open(outputPrefix+".ccga.sdi");

    std::ofstream ofileccgg;
    ofileccgg.open(outputPrefix+".ccgg.sdi");

    std::ofstream ofileccta;
    ofileccta.open(outputPrefix+".ccta.sdi");

    std::ofstream ofilecgga;
    ofilecgga.open(outputPrefix+".cgga.sdi");

    std::ofstream ofilecggc;
    ofilecggc.open(outputPrefix+".cggc.sdi");

    std::ofstream ofilecgta;
    ofilecgta.open(outputPrefix+".cgta.sdi");

    std::ofstream ofilegatc;
    ofilegatc.open(outputPrefix+".gatc.sdi");

    std::ofstream ofilegcta;
    ofilegcta.open(outputPrefix+".gcta.sdi");

    std::ofstream ofileunclassified;
    ofileunclassified.open(outputPrefix+".unclassified.sdi");


    std::ofstream ofileRaatt;
    ofileRaatt.open(outputPrefix+".Raatt.sdi");
    std::ofstream ofileRacgt;
    ofileRacgt.open(outputPrefix+".Racgt.sdi");
    std::ofstream ofileRatat;
    ofileRatat.open(outputPrefix+".Ratat.sdi");
    std::ofstream ofileRagct;
    ofileRagct.open(outputPrefix+".Ragct.sdi");
    std::ofstream ofileRcatg;
    ofileRcatg.open(outputPrefix+".Rcatg.sdi");
    std::ofstream ofileRccgg;
    ofileRccgg.open(outputPrefix+".Rccgg.sdi");
    std::ofstream ofileRcgcg;
    ofileRcgcg.open(outputPrefix+".Rcgcg.sdi");
    std::ofstream ofileRgatc;
    ofileRgatc.open(outputPrefix+".Rgatc.sdi");
    std::ofstream ofileRgcgc;
    ofileRgcgc.open(outputPrefix+".Rgcgc.sdi");
    std::ofstream ofileRtata;
    ofileRtata.open(outputPrefix+".Rtata.sdi");

    std::map<std::string, std::vector<Variant> > sdiMaps;
    //std::cout << "495" <<std::endl;
    readSdiFile(sdiFile, sdiMaps, vcfFix);
    //std::cout << "497" <<std::endl;
    for( std::map<std::string, std::vector<Variant> >::iterator it=sdiMaps.begin();
            it!=sdiMaps.end(); ++it){
        size_t thisChromosomeRecordsNumber = sdiMaps[it->first].size();
        size_t i=1;
        std::vector<size_t> wantedIds;
      //  std::cout << "501" <<std::endl;
        if( sdiMaps[it->first][i].getPosition() == sdiMaps[it->first][i-1].getPosition()+1
                && sdiMaps[it->first][i].getChanginglength()==0 &&
                sdiMaps[it->first][i-1].getChanginglength() == 0) {
            if( legalNasString.find(sdiMaps[it->first][i].getReference())!= legalNasString.end() &&
                legalNasString.find(sdiMaps[it->first][i-1].getReference())!= legalNasString.end()
                && legalNasString.find(sdiMaps[it->first][i].getAlternative())!= legalNasString.end() &&
                    legalNasString.find(sdiMaps[it->first][i-1].getAlternative())!= legalNasString.end() ){
                bool addIt = true;
                for( int j = 0; j <= rangeLength; ++j){
                    if (sdiMaps[it->first][i + 1].overlap(sdiMaps[it->first][i].getPosition() + j)) { //
                        addIt = false;
                    }
                }
                if( addIt ){
                    wantedIds.push_back(i);
                }
            }
        }
        //std::cout << "519" <<std::endl;
        for( i=2; i < thisChromosomeRecordsNumber-1; ++i ){
            if( sdiMaps[it->first][i].getPosition() == sdiMaps[it->first][i-1].getPosition()+1
                && sdiMaps[it->first][i].getChanginglength()==0 &&
                sdiMaps[it->first][i-1].getChanginglength() == 0){
                if( legalNasString.find(sdiMaps[it->first][i].getReference())!= legalNasString.end() &&
                    legalNasString.find(sdiMaps[it->first][i-1].getReference())!= legalNasString.end()
                    && legalNasString.find(sdiMaps[it->first][i].getAlternative())!= legalNasString.end() &&
                    legalNasString.find(sdiMaps[it->first][i-1].getAlternative())!= legalNasString.end()) {
                    bool addIt = true;
                    for( int j = 0; j <= rangeLength; ++j){
                        if (sdiMaps[it->first][i + 1].overlap(sdiMaps[it->first][i].getPosition() + j)) { //
                            addIt = false;
                        }
                        if (sdiMaps[it->first][i - 2].overlap(sdiMaps[it->first][i - 1].getPosition() - j)) {
                            addIt = false;
                        }
                    }
                    if( addIt ){
                        wantedIds.push_back(i);
                    }
                }
            }
        }
        //std::cout << "543" <<std::endl;
        i = thisChromosomeRecordsNumber-1;
        if( sdiMaps[it->first][i].getPosition() == sdiMaps[it->first][i-1].getPosition()+1
            && sdiMaps[it->first][i].getChanginglength()==0 &&
            sdiMaps[it->first][i-1].getChanginglength() == 0) {
            if( legalNasString.find(sdiMaps[it->first][i].getReference())!= legalNasString.end() &&
                legalNasString.find(sdiMaps[it->first][i-1].getReference())!= legalNasString.end()
                && legalNasString.find(sdiMaps[it->first][i].getAlternative())!= legalNasString.end() &&
                legalNasString.find(sdiMaps[it->first][i-1].getAlternative())!= legalNasString.end()){
                bool addIt = true;
                for( int j = 0; j <= rangeLength; ++j){
                    if (sdiMaps[it->first][i - 2].overlap(sdiMaps[it->first][i - 1].getPosition() - j)) {
                        addIt = false;
                    }
                }
                if( addIt ){
                    wantedIds.push_back(i);
                }
            }
        }
        for( std::vector<size_t>::iterator itj=wantedIds.begin();
                itj!=wantedIds.end(); ++itj){
            size_t j = (*itj);
            printSdi(ofilealldouble, sdiMaps[it->first][j - 1]);
            printSdi(ofilealldouble, sdiMaps[it->first][j]);
             //  care about both the reference and alternative
            if( ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "A") &&
                          caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "A") &&
                          caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "C") &&
                          caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "C") ) ||

                        ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "C") &&
                          caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "C") &&
                          caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "A") &&
                          caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "A") ) ||

                        ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "T") &&
                          caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "T") &&
                          caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "G") &&
                          caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "G") ) ||

                        ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "G") &&
                          caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "G") &&
                          caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "T") &&
                          caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "T") )
            ){
                printSdi(ofileaacc, sdiMaps[it->first][j - 1]);
                printSdi(ofileaacc, sdiMaps[it->first][j]);
            }else if( ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "C") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "G") ) ||

                      ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "C") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "G") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "A") ) ||

                      ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "T") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "T") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "C") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "G") ) ||

                      ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "C") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "G") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "T") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "T") )
                    ){
                printSdi(ofileaacg, sdiMaps[it->first][j - 1]);
                printSdi(ofileaacg, sdiMaps[it->first][j]);
            }else if( ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "C") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "T") ) ||

                      ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "C") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "T") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "A") ) ||

                      ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "T") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "T") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "G") ) ||

                      ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "G") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "T") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "T") )
                    ){
                printSdi(ofileaact, sdiMaps[it->first][j - 1]);
                printSdi(ofileaact, sdiMaps[it->first][j]);
            }else if( ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "G") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "C") ) ||

                      ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "G") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "C") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "A") ) ||

                      ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "T") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "T") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "G") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "C") ) ||

                      ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "G") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "C") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "T") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "T") )
                    ){
                printSdi(ofileaagc, sdiMaps[it->first][j - 1]);
                printSdi(ofileaagc, sdiMaps[it->first][j]);
            }else if( ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "G") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "G") ) ||

                      ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "G") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "G") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "A") ) ||

                      ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "T") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "T") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "C") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "C") ) ||

                      ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "C") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "C") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "T") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "T") )
                    ){
                printSdi(ofileaagg, sdiMaps[it->first][j - 1]);
                printSdi(ofileaagg, sdiMaps[it->first][j]);
            }else if( ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "G") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "T") ) ||

                      ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "G") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "T") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "A") ) ||

                      ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "T") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "T") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "C") ) ||

                      ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "C") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "T") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "T") )
                    ){
                printSdi(ofileaagt, sdiMaps[it->first][j - 1]);
                printSdi(ofileaagt, sdiMaps[it->first][j]);
            }else if( ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "T") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "C") ) ||

                      ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "T") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "C") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "A") ) ||

                      ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "T") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "T") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "G") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "A") ) ||

                      ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "G") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "T") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "T") )
                    ){
                printSdi(ofileaatc, sdiMaps[it->first][j - 1]);
                printSdi(ofileaatc, sdiMaps[it->first][j]);
            }else if( ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "T") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "G") ) ||

                      ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "T") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "G") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "A") ) ||

                      ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "T") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "T") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "C") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "A") ) ||

                      ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "C") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "T") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "T") )
                    ){
                printSdi(ofileaatg, sdiMaps[it->first][j - 1]);
                printSdi(ofileaatg, sdiMaps[it->first][j]);
            }else if( ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "T") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "T") ) ||

                      ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "T") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "T") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "A") )
                    ){
                printSdi(ofileaatt, sdiMaps[it->first][j - 1]);
                printSdi(ofileaatt, sdiMaps[it->first][j]);
            }else if( ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "C") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "C") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "A") ) ||

                      ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "C") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "C") ) ||

                      ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "G") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "T") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "T") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "G") ) ||

                      ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "T") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "G") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "G") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "T") )
                    ){
                printSdi(ofileacca, sdiMaps[it->first][j - 1]);
                printSdi(ofileacca, sdiMaps[it->first][j]);
            }else if( ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "C") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "C") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "G") ) ||

                      ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "C") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "G") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "C") ) ||

                      ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "G") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "T") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "C") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "G") ) ||

                      ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "C") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "G") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "G") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "T") )
                    ){
                printSdi(ofileaccg, sdiMaps[it->first][j - 1]);
                printSdi(ofileaccg, sdiMaps[it->first][j]);
            }else if( ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "C") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "C") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "T") ) ||

                      ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "C") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "T") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "C") ) ||

                      ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "G") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "T") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "G") ) ||

                      ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "G") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "G") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "T") )
                    ){
                printSdi(ofileacct, sdiMaps[it->first][j - 1]);
                printSdi(ofileacct, sdiMaps[it->first][j]);
            }else if( ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "C") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "G") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "A") ) ||

                      ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "G") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "C") ) ||

                      ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "G") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "T") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "T") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "C") ) ||

                      ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "T") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "C") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "G") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "T") )
                    ){
                printSdi(ofileacga, sdiMaps[it->first][j - 1]);
                printSdi(ofileacga, sdiMaps[it->first][j]);
            }else if( ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "C") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "G") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "G") ) ||

                      ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "G") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "G") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "C") ) ||

                      ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "G") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "T") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "C") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "C") ) ||

                      ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "C") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "C") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "G") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "T") )
                    ){
                printSdi(ofileacgg, sdiMaps[it->first][j - 1]);
                printSdi(ofileacgg, sdiMaps[it->first][j]);
            }else if( ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "C") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "G") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "T") ) ||

                      ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "G") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "T") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "C") )
                    ){
                printSdi(ofileacgt, sdiMaps[it->first][j - 1]);
                printSdi(ofileacgt, sdiMaps[it->first][j]);
            }else if( ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "C") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "T") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "A") ) ||

                      ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "T") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "C") ) ||

                      ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "G") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "T") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "T") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "A") ) ||

                      ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "T") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "G") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "T") )
                    ){
                printSdi(ofileacta, sdiMaps[it->first][j - 1]);
                printSdi(ofileacta, sdiMaps[it->first][j]);
            }else if( ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "C") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "T") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "G") ) ||

                      ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "T") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "G") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "C") ) ||

                      ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "G") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "T") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "C") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "A") ) ||

                      ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "C") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "G") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "T") )
                    ){
                printSdi(ofileactg, sdiMaps[it->first][j - 1]);
                printSdi(ofileactg, sdiMaps[it->first][j]);
            }else if( ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "A") &&
                       caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "G") &&
                       caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "C") &&
                       caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "A") ) ||

                     ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "C") &&
                       caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "A") &&
                       caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "A") &&
                       caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "G") ) ||

                     ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "C") &&
                       caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "T") &&
                       caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "T") &&
                       caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "G") ) ||

                     ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "T") &&
                       caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "G") &&
                       caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "C") &&
                       caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "T") )
                    ){
                printSdi(ofileagca, sdiMaps[it->first][j - 1]);
                printSdi(ofileagca, sdiMaps[it->first][j]);
            }else if( ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "G") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "C") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "C") ) ||

                      ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "C") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "C") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "G") ) ||

                      ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "C") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "T") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "G") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "G") ) ||

                      ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "G") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "G") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "C") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "T") )
                    ){
                printSdi(ofileagcc, sdiMaps[it->first][j - 1]);
                printSdi(ofileagcc, sdiMaps[it->first][j]);
            }else if( ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "G") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "C") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "T") ) ||

                      ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "C") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "T") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "G") )
                    ){
                printSdi(ofileagct, sdiMaps[it->first][j - 1]);
                printSdi(ofileagct, sdiMaps[it->first][j]);
            }else if( ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "G") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "G") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "A") ) ||

                      ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "G") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "G") ) ||

                      ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "C") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "T") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "T") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "C") ) ||

                      ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "T") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "C") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "C") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "T") )
                    ){
                printSdi(ofileagga, sdiMaps[it->first][j - 1]);
                printSdi(ofileagga, sdiMaps[it->first][j]);
            }else if( ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "G") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "G") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "C") ) ||

                      ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "G") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "C") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "G") ) ||

                      ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "C") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "T") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "G") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "C") ) ||

                      ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "G") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "C") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "C") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "T") )
                    ){
                printSdi(ofileaggc, sdiMaps[it->first][j - 1]);
                printSdi(ofileaggc, sdiMaps[it->first][j]);
            }else if( ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "G") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "T") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "A") ) ||

                      ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "T") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "G") ) ||

                      ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "C") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "T") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "T") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "A") ) ||

                      ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "T") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "C") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "T") )
                    ){
                printSdi(ofileagta, sdiMaps[it->first][j - 1]);
                printSdi(ofileagta, sdiMaps[it->first][j]);
            }else if( ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "A") &&
                       caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "G") &&
                       caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "T") &&
                       caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "C") ) ||

                     ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "T") &&
                       caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "C") &&
                       caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "A") &&
                       caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "G") ) ||

                     ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "C") &&
                       caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "T") &&
                       caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "G") &&
                       caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "A") ) ||

                     ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "G") &&
                       caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "A") &&
                       caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "C") &&
                       caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "T") )
                    ){
                printSdi(ofileagtc, sdiMaps[it->first][j - 1]);
                printSdi(ofileagtc, sdiMaps[it->first][j]);
            }else if( ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "T") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "C") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "A") ) ||

                      ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "C") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "T") ) ||

                      ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "T") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "T") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "G") ) ||

                      ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "T") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "G") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "T") )
                    ){
                printSdi(ofileatca, sdiMaps[it->first][j - 1]);
                printSdi(ofileatca, sdiMaps[it->first][j]);
            }else if( ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "T") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "C") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "C") ) ||

                      ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "C") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "C") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "T") ) ||

                      ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "T") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "G") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "G") ) ||

                      ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "G") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "G") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "T") )
                    ){
                printSdi(ofileatcc, sdiMaps[it->first][j - 1]);
                printSdi(ofileatcc, sdiMaps[it->first][j]);
            }else if( ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "T") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "C") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "G") ) ||

                      ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "C") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "G") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "T") )
                    ){
                printSdi(ofileatcg, sdiMaps[it->first][j - 1]);
                printSdi(ofileatcg, sdiMaps[it->first][j]);
            }else if( ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "T") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "G") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "A") ) ||

                      ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "G") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "T") ) ||

                      ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "T") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "T") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "C") ) ||

                      ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "T") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "C") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "T") )
                    ){
                printSdi(ofileatga, sdiMaps[it->first][j - 1]);
                printSdi(ofileatga, sdiMaps[it->first][j]);
            }else if( ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "T") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "G") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "C") ) ||

                      ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "G") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "C") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "T") )
                    ){
                printSdi(ofileatgc, sdiMaps[it->first][j - 1]);
                printSdi(ofileatgc, sdiMaps[it->first][j]);
            }else if( ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "T") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "T") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "A") ) ||

                      ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "T") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "T") )
                    ){
                printSdi(ofileatta, sdiMaps[it->first][j - 1]);
                printSdi(ofileatta, sdiMaps[it->first][j]);
            }else if( ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "C") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "G") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "C") ) ||

                      ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "G") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "C") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "C") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "A") ) ||

                      ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "T") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "G") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "G") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "C") ) ||

                      ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "G") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "C") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "T") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "G") )
                    ){
                printSdi(ofilecagc, sdiMaps[it->first][j - 1]);
                printSdi(ofilecagc, sdiMaps[it->first][j]);
            }else if( ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "C") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "G") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "G") ) ||

                      ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "G") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "G") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "C") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "A") ) ||

                      ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "T") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "G") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "C") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "C") ) ||

                      ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "C") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "C") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "T") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "G") )
                    ){
                printSdi(ofilecagg, sdiMaps[it->first][j - 1]);
                printSdi(ofilecagg, sdiMaps[it->first][j]);
            }else if( ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "C") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "T") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "C") ) ||

                      ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "T") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "C") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "C") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "A") ) ||

                      ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "T") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "G") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "G") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "A") ) ||

                      ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "G") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "T") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "G") )
                    ){
                printSdi(ofilecatc, sdiMaps[it->first][j - 1]);
                printSdi(ofilecatc, sdiMaps[it->first][j]);
            } else if( ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "C") &&
                         caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "A") &&
                         caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "T") &&
                         caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "G") ) ||

                       ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "T") &&
                         caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "G") &&
                         caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "C") &&
                         caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "A") )
                    ){
                printSdi(ofilecatg, sdiMaps[it->first][j - 1]);
                printSdi(ofilecatg, sdiMaps[it->first][j]);
            } else if ( ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "C") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "C") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "G") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "A") ) ||

                      ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "G") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "C") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "C") ) ||

                      ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "G") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "G") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "T") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "C") ) ||

                      ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "T") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "C") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "G") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "G") )
                    ){
                printSdi(ofileccga, sdiMaps[it->first][j - 1]);
                printSdi(ofileccga, sdiMaps[it->first][j]);
            }else if( ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "C") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "C") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "G") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "G") ) ||

                      ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "G") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "G") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "C") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "C") )
                    ){
                printSdi(ofileccgg, sdiMaps[it->first][j - 1]);
                printSdi(ofileccgg, sdiMaps[it->first][j]);
            }else if( ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "C") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "C") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "T") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "A") ) ||

                      ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "T") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "C") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "C") ) ||

                      ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "G") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "G") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "T") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "A") ) ||

                      ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "T") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "G") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "G") )
                    ){
                printSdi(ofileccta, sdiMaps[it->first][j - 1]);
                printSdi(ofileccta, sdiMaps[it->first][j]);
            }else if( ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "C") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "G") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "G") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "A") ) ||

                      ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "G") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "C") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "G") ) ||

                      ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "C") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "G") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "T") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "C") ) ||

                      ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "T") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "C") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "C") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "G") )
                    ){
                printSdi(ofilecgga, sdiMaps[it->first][j - 1]);
                printSdi(ofilecgga, sdiMaps[it->first][j]);
            }else if( ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "C") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "G") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "G") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "C") ) ||

                      ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "G") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "C") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "C") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "G") )
                    ){
                printSdi(ofilecggc, sdiMaps[it->first][j - 1]);
                printSdi(ofilecggc, sdiMaps[it->first][j]);
            }else if( ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "C") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "G") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "T") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "A") ) ||

                      ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "T") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "C") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "G") )
                    ){
                printSdi(ofilecgta, sdiMaps[it->first][j - 1]);
                printSdi(ofilecgta, sdiMaps[it->first][j]);
            }else if( ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "G") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "T") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "C") ) ||

                      ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "T") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "C") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "G") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "A") )
                    ){
                printSdi(ofilegatc, sdiMaps[it->first][j - 1]);
                printSdi(ofilegatc, sdiMaps[it->first][j]);
            }else if( ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "G") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "C") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "T") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "A") ) ||

                      ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "T") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getAlternative(), "G") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getAlternative(), "C") )
                    ){
                printSdi(ofilegcta, sdiMaps[it->first][j - 1]);
                printSdi(ofilegcta, sdiMaps[it->first][j]);
            }else{
                printSdi(ofileunclassified, sdiMaps[it->first][j - 1]);
                printSdi(ofileunclassified, sdiMaps[it->first][j]);
            }

            //  only care about what is the reference
            if( ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "A") &&
                  caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "A")  ) ||

                ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "T") &&
                  caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "T") ) ){
                printSdi(ofileRaatt, sdiMaps[it->first][j - 1]);
                printSdi(ofileRaatt, sdiMaps[it->first][j]);
            } else if ( ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "A") &&
                          caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "C")  ) ||

                        ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "G") &&
                          caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "T") ) ){
                printSdi(ofileRacgt, sdiMaps[it->first][j - 1]);
                printSdi(ofileRacgt, sdiMaps[it->first][j]);
            } else if ( ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "A") &&
                        caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "T")  ) ){
                printSdi(ofileRatat, sdiMaps[it->first][j - 1]);
                printSdi(ofileRatat, sdiMaps[it->first][j]);
            } else if ( ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "A") &&
                          caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "G")  ) ||

                        ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "C") &&
                          caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "T") ) ){
                printSdi(ofileRagct, sdiMaps[it->first][j - 1]);
                printSdi(ofileRagct, sdiMaps[it->first][j]);
            } else if ( ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "C") &&
                         caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "A")  ) ||

                       ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "T") &&
                         caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "G") ) ){
                printSdi(ofileRcatg, sdiMaps[it->first][j - 1]);
                printSdi(ofileRcatg, sdiMaps[it->first][j]);
            } else if ( ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "C") &&
                         caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "C")  ) ||

                       ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "G") &&
                         caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "G") ) ){
                printSdi(ofileRccgg, sdiMaps[it->first][j - 1]);
                printSdi(ofileRccgg, sdiMaps[it->first][j]);
            } else if ( ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "C") &&
                          caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "G")  )  ){
                printSdi(ofileRcgcg, sdiMaps[it->first][j - 1]);
                printSdi(ofileRcgcg, sdiMaps[it->first][j]);
            } else if ( ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "G") &&
                          caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "A")  ) ||

                        ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "T") &&
                          caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "C") ) ){
                printSdi(ofileRgatc, sdiMaps[it->first][j - 1]);
                printSdi(ofileRgatc, sdiMaps[it->first][j]);
            } else if ( ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "G") &&
                          caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "C")  ) ){
                printSdi(ofileRgcgc, sdiMaps[it->first][j - 1]);
                printSdi(ofileRgcgc, sdiMaps[it->first][j]);
            }else if ( ( caseInsensitiveStringCompare(sdiMaps[it->first][j-1].getReference(), "T") &&
                         caseInsensitiveStringCompare(sdiMaps[it->first][j].getReference(), "A")  ) ){
                printSdi(ofileRtata, sdiMaps[it->first][j - 1]);
                printSdi(ofileRtata, sdiMaps[it->first][j]);
            }
        }
    }
    ofilealldouble.close();


    ofileaacc.close();
    ofileaacg.close();
    ofileaact.close();
    ofileaagc.close();
    ofileaagg.close();
    ofileaagt.close();
    ofileaatc.close();
    ofileaatg.close();
    ofileaatt.close();
    ofileacca.close();
    ofileaccg.close();
    ofileacct.close();
    ofileacga.close();
    ofileacgg.close();
    ofileacgt.close();
    ofileacta.close();
    ofileactg.close();
    ofileagca.close();
    ofileagcc.close();
    ofileagct.close();
    ofileagga.close();
    ofileaggc.close();
    ofileagta.close();
    ofileagtc.close();
    ofileatca.close();
    ofileatcc.close();
    ofileatcg.close();
    ofileatga.close();
    ofileatgc.close();
    ofileatta.close();
    ofilecagc.close();
    ofilecagg.close();
    ofilecatc.close();
    ofilecatg.close();
    ofileccga.close();
    ofileccgg.close();
    ofileccta.close();
    ofilecgga.close();
    ofilecggc.close();
    ofilecgta.close();
    ofilegatc.close();
    ofilegcta.close();
    ofileunclassified.close();


    ofileRaatt.close();
    ofileRacgt.close();
    ofileRatat.close();
    ofileRagct.close();
    ofileRcatg.close();
    ofileRccgg.close();
    ofileRcgcg.close();
    ofileRgatc.close();
    ofileRgcgc.close();
    ofileRtata.close();

}


void countNumberSNPAndIndel( std::string& sdiFile, std::string & outputPrefix, int & rangeLength, std::string& vcfFix) {

    NucleotideCodeSubstitutionMatrix nucleotideCodeSubstitutionMatrix;
    std::set<std::string> &legalNasString = nucleotideCodeSubstitutionMatrix.getLegalNasString();
    std::map<std::string, std::vector<Variant> > sdiMaps;
    readSdiFile(sdiFile, sdiMaps, vcfFix);

    std::ofstream ofileallsnps;
    ofileallsnps.open(outputPrefix+".allsnps.sdi");
    std::ofstream ofileallindels;
    ofileallindels.open(outputPrefix+".allindels.sdi");

    std::ofstream ofilelegalsnps;
    ofilelegalsnps.open(outputPrefix+".legalsnps.sdi");

    std::ofstream ofilelegalsnpsNoNearByIndel;
    ofilelegalsnpsNoNearByIndel.open(outputPrefix+".legalsnpsNoNearByIndel.sdi");

    std::ofstream ofilelegalsnpsIsolate;
    ofilelegalsnpsIsolate.open(outputPrefix+".legalsnpsIsolate.sdi");

    std::ofstream ofileillegalsnpsNoNearByIndel;
    ofileillegalsnpsNoNearByIndel.open(outputPrefix+".illegalsnpsNoNearByIndel.sdi");

    std::ofstream ofileillegalsnpsIsolate;
    ofileillegalsnpsIsolate.open(outputPrefix+".illegalsnpsIsolate.sdi");

    std::ofstream ofileillegalsnps;
    ofileillegalsnps.open(outputPrefix+".illegalsnps.sdi");

    for( std::map<std::string, std::vector<Variant> >::iterator it=sdiMaps.begin();
         it!=sdiMaps.end(); ++it) {
        size_t thisChromosomeRecordsNumber = sdiMaps[it->first].size();

        for( size_t i=0; i < thisChromosomeRecordsNumber; ++i ){
            if(  sdiMaps[it->first][i].getChanginglength()==0 && sdiMaps[it->first][i].getReference().size()==1 ){
                if( legalNasString.find(sdiMaps[it->first][i].getReference())!= legalNasString.end()
                    && legalNasString.find(sdiMaps[it->first][i].getAlternative())!= legalNasString.end() ) {
                    printSdi(ofilelegalsnps, sdiMaps[it->first][i]);

                    //ofilelegalsnpsNoNearByIndel begin
                    bool addIt = true;
                    for( int j = 0; j <= rangeLength; ++j){
                        if ( i>0 && sdiMaps[it->first][i-1].getChanginglength()!=0 && sdiMaps[it->first][i-1].overlap(sdiMaps[it->first][i].getPosition() - j)) { //
                            addIt = false;
                        }
                        if (i<(thisChromosomeRecordsNumber-1) && sdiMaps[it->first][i+1].getChanginglength()!=0 && sdiMaps[it->first][i +1].overlap(sdiMaps[it->first][i].getPosition() + j)) {
                            addIt = false;
                        }
                    }
                    if( addIt ){
                        printSdi(ofilelegalsnpsNoNearByIndel, sdiMaps[it->first][i]);
                    }//ofilelegalsnpsNoNearByIndel end


                    //ofilelegalsnpsIsolate begin
                    addIt = true;
                    for( int j = 0; j <= rangeLength; ++j){
                        if ( i>0  && sdiMaps[it->first][i-1].overlap(sdiMaps[it->first][i].getPosition() - j)) { //
                            addIt = false;
                        }
                        if (i<(thisChromosomeRecordsNumber-1) && sdiMaps[it->first][i +1].overlap(sdiMaps[it->first][i].getPosition() + j)) {
                            addIt = false;
                        }
                    }
                    if( addIt ){
                        printSdi(ofilelegalsnpsIsolate, sdiMaps[it->first][i]);
                    }
                    //ofilelegalsnpsIsolate end

                }else{
                    printSdi(ofileillegalsnps, sdiMaps[it->first][i]);

                    //ofileillegalsnpsNoNearByIndel begin
                    bool addIt = true;
                    for( int j = 0; j <= rangeLength; ++j){
                        if ( i>0 && sdiMaps[it->first][i-1].getChanginglength()!=0 && sdiMaps[it->first][i-1].overlap(sdiMaps[it->first][i].getPosition() - j)) { //
                            addIt = false;
                        }
                        if (i<(thisChromosomeRecordsNumber-1) && sdiMaps[it->first][i+1].getChanginglength()!=0 && sdiMaps[it->first][i +1].overlap(sdiMaps[it->first][i].getPosition() + j)) {
                            addIt = false;
                        }
                    }
                    if( addIt ){
                        printSdi(ofileillegalsnpsNoNearByIndel, sdiMaps[it->first][i]);
                    }//ofileillegalsnpsNoNearByIndel end


                    //ofileillegalsnpsIsolate begin
                    addIt = true;
                    for( int j = 0; j <= rangeLength; ++j){
                        if ( i>0  && sdiMaps[it->first][i-1].overlap(sdiMaps[it->first][i].getPosition() - j)) { //
                            addIt = false;
                        }
                        if (i<(thisChromosomeRecordsNumber-1) && sdiMaps[it->first][i +1].overlap(sdiMaps[it->first][i].getPosition() + j)) {
                            addIt = false;
                        }
                    }
                    if( addIt ){
                        printSdi(ofileillegalsnpsIsolate, sdiMaps[it->first][i]);
                    }
                    //ofileillegalsnpsIsolate end

                }
                printSdi(ofileallsnps, sdiMaps[it->first][i]);


            }else{
                printSdi(ofileallindels, sdiMaps[it->first][i]);
            }
        }
    }
    ofileallsnps.close();
    ofileallindels.close();
    ofilelegalsnps.close();
    ofileillegalsnps.close();
    ofilelegalsnpsIsolate.close();
    ofilelegalsnpsNoNearByIndel.close();
    ofileillegalsnpsNoNearByIndel.close();
    ofileillegalsnpsIsolate.close();
}


int ranDomItDeletion( Variant & obvervedVariant, std::map<std::string, std::vector<Variant> >& randomSdiMaps, std::string& chr, int& chrSize ){
    bool notInserted = true; // the record is not inserted
    int timesTried = 0;
    while( notInserted ){

        int position = (rand() % chrSize) + 1;
        std::string ref = obvervedVariant.getReference();
        std::string alt = obvervedVariant.getAlternative();
        Variant variant(chr, position, ref, alt); // the (*it2) variant records but assigned with a random position
        if( variant.getChanginglength() <0 && variant.getPosition()-variant.getChanginglength() > chrSize+1  ){
            // if this is deletion and deletion is beyond the chromosome, it does not make sense, so avoid such records
        }else{
            bool hasOverLap = false;
            for( std::vector<Variant>::iterator it3=randomSdiMaps[chr].begin(); it3!=randomSdiMaps[chr].end(); ++it3 ){
                // if the new record conflicts with any of the previous records, re-do position random
                int start1 = it3->getPosition();
                int end1 = it3->getPosition() + abs(it3->getChanginglength());
                int start2 = variant.getPosition();
                int end2 = variant.getPosition() + abs(variant.getChanginglength());
                if(  end1 < start2 || end2 < start1  ){
                    continue; // for such case, there is definitly no confliction
                }
                if(
                    ( (it3->getPosition()<=variant.getPosition() && variant.getPosition()<=it3->getPosition()-it3->getChanginglength()) ||
                      ( variant.getPosition()<=it3->getPosition() && it3->getPosition()<=variant.getPosition()-variant.getChanginglength() ) ) ){
                    // deletion VS deletion
                    //two deletion could not neighbor with each other, else they should be marged as one single records
                    hasOverLap = true;
                    break;
                }
            }
            if( ! hasOverLap  ){
                randomSdiMaps[chr].push_back(variant); // insert the variants with random position into random records data structure
                //printSdi(std::cout, variant);
                notInserted = false;
                return 0; // good
            }
        }
        timesTried++;
        if( timesTried>100000 ){
            //after 100000 fails, give up and redo everything
            return 1;
        }
    }
    return 0; // good
}

int ranDomItSnp( Variant & obvervedVariant, std::map<std::string, std::vector<Variant> >& randomSdiMaps, std::string& chr, int& chrSize ){
    bool notInserted = true; // the record is not inserted
    int timesTried = 0;
    while( notInserted ){

        int position = (rand() % chrSize) + 1;
        std::string ref = obvervedVariant.getReference();
        std::string alt = obvervedVariant.getAlternative();
        Variant variant(chr, position, ref, alt); // the (*it2) variant records but assigned with a random position
        if( variant.getChanginglength() <0 && variant.getPosition()-variant.getChanginglength() > chrSize+1  ){
            // if this is deletion and deletion is beyond the chromosome, it does not make sense, so avoid such records
        }else{
            bool hasOverLap = false;
            for( std::vector<Variant>::iterator it3=randomSdiMaps[chr].begin(); it3!=randomSdiMaps[chr].end(); ++it3 ){
                // if the new record conflicts with any of the previous records, re-do position random
                int start1 = it3->getPosition();
                int end1 = it3->getPosition() + abs(it3->getChanginglength());
                int start2 = variant.getPosition();
                int end2 = variant.getPosition() + abs(variant.getChanginglength());
                if(  end1 < start2 || end2 < start1  ){
                    continue; // for such case, there is definitly no confliction
                }
                if( it3->getChanginglength()==0 && variant.getPosition()==it3->getPosition() ){
                    //SNP VS SNP
                    //two SNP could not share the same position
                    hasOverLap = true;
                    break;
                }else if( it3->getChanginglength()<0  && it3->getPosition()<=variant.getPosition() && variant.getPosition() <it3->getPosition()-variant.getChanginglength() ){
                    //deletion VS SNP
                    hasOverLap = true;
                    break;
                }// there is no confliction between insertion and SNP
            }
            if( ! hasOverLap  ){
                randomSdiMaps[chr].push_back(variant); // insert the variants with random position into random records data structure
//                printSdi(std::cout, variant);
                notInserted = false;
                return 0; // good
            }
        }
        timesTried++;
        if( timesTried>100000 ){
            //after 100000 fails, give up and redo everything
            return 1;
        }
    }
    return 0; // good
}


int ranDomItInsertion( Variant & obvervedVariant, std::map<std::string, std::vector<Variant> >& randomSdiMaps, std::string& chr, int& chrSize ){
    bool notInserted = true; // the record is not inserted
    int timesTried = 0;
    while( notInserted ){
        int position = (rand() % chrSize) + 1;
        std::string ref = obvervedVariant.getReference();
        std::string alt = obvervedVariant.getAlternative();
        Variant variant(chr, position, ref, alt); // the (*it2) variant records but assigned with a random position
        if( variant.getChanginglength() <0 && variant.getPosition()-variant.getChanginglength() > chrSize+1  ){
            // if this is deletion and deletion is beyond the chromosome, it does not make sense, so avoid such records
        }else{
            bool hasOverLap = false;
            for( std::vector<Variant>::iterator it3=randomSdiMaps[chr].begin(); it3!=randomSdiMaps[chr].end(); ++it3 ){
                // if the new record conflicts with any of the previous records, re-do position random
                int start1 = it3->getPosition();
                int end1 = it3->getPosition() + abs(it3->getChanginglength());
                int start2 = variant.getPosition();
                int end2 = variant.getPosition() + abs(variant.getChanginglength());
                if(  end1 < start2 || end2 < start1  ){
                    continue; // for such case, there is definitly no confliction
                }

                if( it3->getChanginglength()>0 && variant.getChanginglength()>0 && it3->getPosition()==position ){
                    //insertion VS insertion
                    //two insertion could not at the same position, else they should be marged as one single records
                    hasOverLap = true;
                    break;
                }else if( it3->getChanginglength()<0 && variant.getPosition()>it3->getPosition() &&  variant.getPosition()< it3->getPosition()-it3->getChanginglength() ){
                    //deletion VS insertion 1
                    //could no insert at deleted sequence
                    hasOverLap = true;
                    break;
                }// there is no confliction between insertion and SNP
            }
            if( ! hasOverLap  ){
                randomSdiMaps[chr].push_back(variant); // insert the variants with random position into random records data structure
//                printSdi(std::cout, variant);
                notInserted = false;
                return 0; // good
            }
        }
        timesTried++;
        if( timesTried>100000 ){
            //after 100000 fails, give up and redo everything
            return 1;
        }
    }
    return 0; // good
}

int ranDomIt( Variant & obvervedVariant, std::map<std::string, std::vector<Variant> >& randomSdiMaps, std::string& chr, int& chrSize ){
    bool notInserted = true; // the record is not inserted
    int timesTried = 0;
    while( notInserted ){

        int position = (rand() % chrSize) + 1;
        std::string ref = obvervedVariant.getReference();
        std::string alt = obvervedVariant.getAlternative();
        Variant variant(chr, position, ref, alt); // the (*it2) variant records but assigned with a random position
        if( variant.getChanginglength() <0 && variant.getPosition()-variant.getChanginglength() > chrSize+1  ){
            // if this is deletion and deletion is beyond the chromosome, it does not make sense, so avoid such records
        }else{
            bool hasOverLap = false;
            for( std::vector<Variant>::iterator it3=randomSdiMaps[chr].begin(); it3!=randomSdiMaps[chr].end(); ++it3 ){
                // if the new record conflicts with any of the previous records, re-do position random
                int start1 = it3->getPosition();
                int end1 = it3->getPosition() + abs(it3->getChanginglength());
                int start2 = variant.getPosition();
                int end2 = variant.getPosition() + abs(variant.getChanginglength());
                if(  end1 < start2 || end2 < start1  ){
                    continue; // for such case, there is definitly no confliction
                }
                if( it3->getChanginglength()<0 && variant.getChanginglength()<0 &&
                    ( (it3->getPosition()<=variant.getPosition() && variant.getPosition()<=it3->getPosition()-it3->getChanginglength()) ||
                      ( variant.getPosition()<=it3->getPosition() && it3->getPosition()<=variant.getPosition()-variant.getChanginglength() ) ) ){
                    // deletion VS deletion
                    //two deletion could not neighbor with each other, else they should be marged as one single records
                    hasOverLap = true;
                    break;
                }else if( it3->getChanginglength()>0 && variant.getChanginglength()>0 && it3->getPosition()==position ){
                    //insertion VS insertion
                    //two insertion could not at the same position, else they should be marged as one single records
                    hasOverLap = true;
                    break;
                }else if( it3->getChanginglength()<0 && variant.getChanginglength()>0 && variant.getPosition()>it3->getPosition() &&  variant.getPosition()< it3->getPosition()-it3->getChanginglength() ){
                    //deletion VS insertion 1
                    //could no insert at deleted sequence
                    hasOverLap = true;
                    break;
                }else if( variant.getChanginglength()<0 && it3->getChanginglength()>0 && it3->getPosition()>variant.getPosition() &&  it3->getPosition()< variant.getPosition()-variant.getChanginglength() ){
                    //deletion VS insertion 2
                    hasOverLap = true;
                    break;
                }else if( variant.getChanginglength()==0 && it3->getChanginglength()==0 && variant.getPosition()==it3->getPosition() ){
                    //SNP VS SNP
                    //two SNP could not share the same position
                    hasOverLap = true;
                    break;
                }else if( variant.getChanginglength()<0 && it3->getChanginglength()==0 && variant.getPosition()<=it3->getPosition() && it3->getPosition() <variant.getPosition()-variant.getChanginglength() ){
                    //deletion VS SNP 1
                    //could no SNP at deleted sequence
                    hasOverLap = true;
                    break;
                }else if( it3->getChanginglength()<0 && variant.getChanginglength()==0 && it3->getPosition()<=variant.getPosition() && variant.getPosition() <it3->getPosition()-variant.getChanginglength() ){
                    //deletion VS SNP 2
                    hasOverLap = true;
                    break;
                }// there is no confliction between insertion and SNP
            }
            if( ! hasOverLap  ){
                randomSdiMaps[chr].push_back(variant); // insert the variants with random position into random records data structure
                //printSdi(std::cout, variant);
                notInserted = false;
                return 0; // good
            }
        }
        timesTried++;
        if( timesTried>100000 ){
            //after 100000 fails, give up and redo everything
            return 1;
        }
    }
    return 0; // good
}

void generateRandomSdi( std::string& sdiFile, std::string& referenceGenomeFile, std::string & outputPrefix, std::string& vcfFix ){

    std::map<std::string, std::vector<Variant> > sdiMaps;
    readSdiFile(sdiFile, sdiMaps, vcfFix);

    std::map<std::string, Fasta> genome;
    readFastaFile(referenceGenomeFile, genome);

    std::map<std::string, int> chrSizeMap;
    std::map<std::string, std::vector<Variant> > randomSdiMaps;
    for( std::map<std::string, Fasta>::iterator it=genome.begin(); it!=genome.end(); ++it ){
        chrSizeMap[it->first]=it->second.getSequence().size();
        randomSdiMaps[it->first]=std::vector<Variant>();
    }
    std::cout << "begin to random" << std::endl;
    for( std::map<std::string, std::vector<Variant> >::iterator it=sdiMaps.begin(); it!=sdiMaps.end(); ++it) {
        std::string chr = it->first;
        int chrSize = chrSizeMap[chr];
        std::cout << "begin to random " << chr << std::endl;

        reDoThisChrLable:
            int deletionNumber = 0;

            for( std::vector<Variant>::iterator it2=it->second.begin(); it2!=it->second.end(); ++it2 ){ // the observed variant records
                if( it2->getChanginglength() < -10000 ){
                    int results = ranDomItDeletion( (*it2), randomSdiMaps, chr, chrSize );
                    if( 1==results ){
                        randomSdiMaps[chr].clear();
                        std::cout << "begin to re-random " << chr << std::endl;
                        goto reDoThisChrLable;
                    }
                    ++deletionNumber;
                    if( deletionNumber % 5000 ==0 ) {
                        std::cout << deletionNumber << " deletion finished" << std::endl;
                    }
                }
            }

            for( std::vector<Variant>::iterator it2=it->second.begin(); it2!=it->second.end(); ++it2 ){ // the observed variant records
                if( it2->getChanginglength() < -1000 &&it2->getChanginglength() >= -10000  ){
                    int results = ranDomItDeletion( (*it2), randomSdiMaps, chr, chrSize );
                    if( 1==results ){
                        randomSdiMaps[chr].clear();
                        std::cout << "begin to re-random " << chr << std::endl;
                        goto reDoThisChrLable;
                    }
                    ++deletionNumber;
                    if( deletionNumber % 5000 ==0 ) {
                        std::cout << deletionNumber << " deletion finished" << std::endl;
                    }
                }
            }
            for( std::vector<Variant>::iterator it2=it->second.begin(); it2!=it->second.end(); ++it2 ){ // the observed variant records
                if( it2->getChanginglength() < 0 && it2->getChanginglength() >= -1000 ){
                    int results = ranDomItDeletion( (*it2), randomSdiMaps, chr, chrSize );
                    if( 1==results ){
                        randomSdiMaps[chr].clear();
                        std::cout << "begin to re-random " << chr << std::endl;
                        goto reDoThisChrLable;
                    }
                    ++deletionNumber;
                    if( deletionNumber % 5000 ==0 ) {
                        std::cout << deletionNumber << " deletion finished" << std::endl;
                    }
                }
            }

            int snpNumber = 0;
            for( std::vector<Variant>::iterator it2=it->second.begin(); it2!=it->second.end(); ++it2 ){ // the observed variant records
                if( it2->getChanginglength() == 0 ){
                    int results = ranDomItSnp( (*it2), randomSdiMaps, chr, chrSize );
                    if( 1==results ){
                        randomSdiMaps[chr].clear();
                        std::cout << "begin to re-random " << chr << std::endl;
                        goto reDoThisChrLable;
                    }
                    ++snpNumber;
                    if( snpNumber % 5000 ==0 ) {
                        std::cout << snpNumber << " SNP finished" << std::endl;
                    }
                }
            }

            int insertNumber = 0;
            for( std::vector<Variant>::iterator it2=it->second.begin(); it2!=it->second.end(); ++it2 ){ // the observed variant records
                if( it2->getChanginglength() > 0 ) {
                    int results = ranDomItInsertion((*it2), randomSdiMaps, chr, chrSize);
                    if (1 == results) {
                        randomSdiMaps[chr].clear();
                        std::cout << "begin to re-random " << chr << std::endl;
                        goto reDoThisChrLable;
                    }
                    ++insertNumber;
                    if (insertNumber % 5000 == 0) {
                        std::cout << insertNumber << " insert finished" << std::endl;
                    }
                }
            }
    }

    std::ofstream ofileRandomSdi;
    ofileRandomSdi.open(outputPrefix+".random.sdi");
    for(std::map<std::string, std::vector<Variant> >::iterator it=randomSdiMaps.begin(); it!=randomSdiMaps.end(); ++it){
        std::sort(it->second.begin(), it->second.end(), compare_sdi_record());
        for( std::vector<Variant>::iterator it2=it->second.begin(); it2!=it->second.end(); ++it2 ){
            printSdi(ofileRandomSdi, (*it2));
        }
    }
    ofileRandomSdi.close();
}

//
//**
// *
// * AA TT
// * AC GT
// * AT AT
// * AG CT
// *
// * CA TG
// * CC GG
// * CG CG
// * CT AG
// *
// * GA TC
// * GC GC
// * GG CC
// * GT AC
// *
// * TA TA
// * TC GA
// * TG CA
// * TT AA
// *
// * THEN WHAT IS LEFT
// *
// * AA/TT
// * AC/GT
// * AT/AT
// * AG/CT
// * CA/TG
// * CC/GG
// * CG/CG
// * GA/TC
// * GC/GC
// * TA/TA
// * /

/*
double snp combinations:
AA - CC  CC - AA  TT - GG  GG - TT
AA - CG  CG - AA  TT - CG  CG - TT
AA - CT  CT - AA  TT - AG  AG - TT
AA - GC  GC - AA  TT - GC  GC - TT
AA - GG  GG - AA  TT - CC  CC - TT
AA - GT  GT - AA  TT - AC  AC - TT
AA - TC  TC - AA  TT - GA  GA - TT
AA - TG  TG - AA  TT - CA  CA - TT
AA - TT  TT - AA

AC - CA  CA - AC  TG - GT  GT - TG
AC - CG  CG - AC  CG - GT  GT - CG
AC - CT  CT - AC  AG - GT  GT - AG
AC - GA  GA - AC  TC - GT  GT - TC
AC - GG  GG - AC  CC - GT  GT - CC
AC - GT  GT - AC
AC - TA  TA - AC  TA - GT  GT - TA
AC - TG  TG - AC  CA - GT  GT - CA

AG - CA  CA - AG  CT - TG  TG - CT
AG - CC  CC - AG  CT - GG  GG - CT
AG - CT  CT - AG
AG - GA  GA - AG  CT - TC  TC - CT
AG - GC  GC - AG  CT - GC  GC - CT
AG - TA  TA - AG  CT - TA  TA - CT
AG - TC  TC - AG  CT - GA  GA - CT

AT - CA  CA - AT  AT - TG  TG - AT
AT - CC  CC - AT  AT - GG  GG - AT
AT - CG  CG - AT
AT - GA  GA - AT  AT - TC  TC - AT
AT - GC  GC - AT
AT - TA  TA - AT

CA - GC  GC - CA  TG - GC  GC - TG
CA - GG  GG - CA  TG - CC  CC - TG
CA - TC  TC - CA  TG - GA  GA - TG
CA - TG  TG - CA

CC - GA  GA - CC  GG - TC  TC - GG
CC - GG  GG - CC
CC - TA  TA - CC  GG - TA  TA - GG

CG - GA  GA - CG  CG - TC  TC - CG
CG - GC  GC - CG
CG - TA  TA - CG

GA - TC  TC - GA

GC - TA  TA - GC


 */