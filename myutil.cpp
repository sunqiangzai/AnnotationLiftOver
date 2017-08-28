// =====================================================================================
// 
//       Filename:  myutil.cpp
// 
//    Description:  
// 
//        Version:  1.0
//        Created:  04/12/2017 04:53:23 PM
//       Revision:  none
//       Compiler:  g++
// 
//         Author:  Baoxing Song (songbx.me), song@mpipz.mpg.de
//        Company:  MPIPZ
// 
// =====================================================================================
/*************************************************************************


************************************************************************/

#include <string>
#include <vector>
#include <iostream>
#include <sstream>
#include <fstream>
#include <regex>
#include "myutil.h"
#include <sys/stat.h>
#include <unistd.h>

std::vector<std::string> &split(const std::string &s, char delim,std::vector<std::string> &elems) {
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        if (item.length() > 0) {
            elems.push_back(item);
        }
    }
    return elems;
}


std::vector<std::string> split(const std::string &s, char delim) {
    std::vector<std::string> elems;
    split(s, delim, elems);
    return elems;
}


void songToUpCase( std::string& str  ){
    transform(str.begin(), str.end(), str.begin(),::toupper);
}
//*  input a string, the transformed would also be returned */
void songToLowCase( std::string& str   ){
    transform(str.begin(), str.end(), str.begin(), ::tolower);
}
/* this implemention is very slow, maybe should find a better way to redo it */
std::string& songStrRemoveBlank( std::string& str  ){
    std::string pattern="\\s";
    std::string pattern2="";
    return songStrReplaceAll( str, pattern, pattern2 );
}

std::string& songStrReplaceAll( std::string& str, std::string& pattern, std::string& pattern2  ){
    std::regex vowel_re(pattern);
    str=std::regex_replace(str, vowel_re, pattern2);
    return str;
}

std::string getReverseComplementary(std::string& sequence) {
    std::stringstream reversecomplementary;
    for (int i = sequence.length() - 1; i >= 0; i--) {
        char c = sequence[i];
        if ('A' == c) {
            c = 'T';
        } else if ('T' == c) {
            c = 'A';
        } else if ('U' == c) {
            c = 'A';
        } else if ('C' == c) {
            c = 'G';
        } else if ('G' == c) {
            c = 'C';
        } else if ('R' == c) {
            c = 'Y';
        } else if ('Y' == c) {
            c = 'R';
        } else if ('K' == c) {
            c = 'M';
        } else if ('M' == c) {
            c = 'K';
        } else if ('B' == c) {
            c = 'V';
        } else if ('V' == c) {
            c = 'B';
        } else if ('D' == c) {
            c = 'H';
        } else if ('H' == c) {
            c = 'D';
        }
        reversecomplementary<< c;
    }
    return reversecomplementary.str();
//    std::string rc = reversecomplementary.str();
//    return rc;
}
void readFastaFile( const std::string& filePath, std::map<std::string, Fasta>& sequences ){
    std::ifstream infile(filePath);
    if( ! infile.good()){
        std::cerr << "error in opening fasta file " << filePath << std::endl;
        exit (1);
    }
    std::regex reg("^>(\\S+)");
    std::string name="";
    std::stringstream sequencestream;
    std::string line="";
    while (std::getline(infile, line)){
        std::smatch match;
        regex_search(line, match, reg);
        if( match.empty()   ){
            sequencestream << line;
        }else{
            if( name.size()>0 ){
                std::string sequence = sequencestream.str();
                Fasta f (name, sequence);
                sequences[name]=f;
            }
            name=match[1];
            sequencestream.str(std::string());
        }
    }
    if( name.size()>0 ){
        std::string sequence = sequencestream.str();
        Fasta f (name, sequence);
        if( f.getSequence().size()>0 ){
            sequences[name]=f;
        }
    }
}

std::string getSubsequence(std::map<std::string, Fasta>& sequences, std::string seqName, int start, int end){
    return getSubsequence(sequences, seqName, start, end, POSITIVE);
}
std::string getSubsequence(std::map<std::string, Fasta>& sequences, std::string seqName, int start, int end, STRAND strand){
    if( sequences.find(seqName)!=sequences.end() ){
        //std::cout << "114 " << start << " " << end << std::endl;
        std::string seq = sequences[seqName].getSubsequence(start, end);
        if( strand == POSITIVE ){
            return seq;
        }else{
//            std::cout << "118" << std::endl;
            std::string revSeq = getReverseComplementary(seq);
  //          std::cout << "120" << std::endl;
            return revSeq;
        }
    }
    return "";
}
void readSdiFile(const std::string& filePath, std::map<std::string, std::vector<Variant> >& variantsMap){
    readSdiFile(filePath, variantsMap, "" );
}
void readSdiFile (const std::string& filePath, std::map<std::string, std::vector<Variant> >& variantsMap, const std::string& chromosome){
    std::ifstream infile(filePath);
    if( ! infile.good()){
        std::cerr << "error in opening variants file " << filePath << std::endl;
        exit (1);
    }
//    std::regex reg("^(\\S+)\\s+(\\S+)\\s+(\\S+)\\s+(\\S+)\\s+(\\S+)");
    std::string line="";
//    std::cout << "143" << std::endl;
    while (std::getline(infile, line)){
        if( line.compare(0, 1, "#")==0 ){
            //std::cout << line << std::endl;
            continue;
        }

        if( chromosome.length()>0) {
            if( line.compare(0, chromosome.size(), chromosome)!=0 ){
                continue;
            }
        }

//        std::cout << "good " << line << std::endl;
//        std:: smatch match;
//        std::cout << line << std::endl;
//        regex_search(line, match, reg);
//        std::cout << line << std::endl;
        std::vector<std::string> splits = split(line, '\t');
        if( splits.size() >=5 ){
  //          std::cout << line << std::endl;
            std::string chromosome = splits[0];
            int position = std::stoi(splits[1]);
            std::string reference = splits[3];
            std::string alternative = splits[4];
            Variant variant(chromosome, position, reference, alternative);
            if( variantsMap.find(chromosome) == variantsMap.end() ){
                variantsMap[chromosome]=std::vector<Variant>();
            }
            variantsMap[chromosome].push_back(variant);
//
//            println(std::cout, variant);
        }
//        if( match.empty() ){
//            //std::cout << "empty " << line << std::endl;
//        }else{
//            std::cout << line << std::endl;
//            std::string chromosome = match[1];
//            int position = std::stoi(match[2]);
//            std::string reference = match[4];
//            std::string alternative = match[5];
//            Variant variant(chromosome, position, reference, alternative);
//            if( variantsMap.find(chromosome) == variantsMap.end() ){
//                variantsMap[chromosome]=std::vector<Variant>();
//            }
//            variantsMap[chromosome].push_back(variant);
//        }
    }
    //std::cout << "164" << std::endl;
    for(std::map<std::string, std::vector<Variant> >::iterator it=variantsMap.begin(); it!=variantsMap.end(); ++it){
        int lastTotalChanged=0;
        for( std::vector<Variant>::iterator it2=(*it).second.begin(); it2!=(*it).second.end(); it2++ ){
            (*it2).setLastTotalChanged(lastTotalChanged);
            lastTotalChanged +=(*it2).getChanginglength();
        }
    }
}

int getChangedFromBasement(std::string chromosomeName, int basement, std::map<std::string, std::vector<Variant> >& variantsMap){
    if( variantsMap.find(chromosomeName)!=variantsMap.end() && variantsMap[chromosomeName].size()>0 ){
        int start = 0 ;
        int end = variantsMap[chromosomeName].size()-1;
        int lastStart = start;
//        std::cout << "166" << std::endl;
        if(variantsMap[chromosomeName][start].getPosition() >= basement){
            return basement;//allChanged;
        }
  //      std::cout << "169" << std::endl;
        if(variantsMap[chromosomeName][end].getPosition() <= basement){
            return basement+variantsMap[chromosomeName][end].getChanginglength()+variantsMap[chromosomeName][end].getLastTotalChanged();//allChanged;
        }
    //    std::cout << "173" << std::endl;
        while(!((variantsMap[chromosomeName][start].getPosition() < basement) && (variantsMap[chromosomeName][start+1].getPosition() >= basement))){
            if((variantsMap[chromosomeName][start].getPosition() < basement)){
                lastStart = start;
                if(1 == (end - start)){
                    start = end;
                }else{
                    start = (start+end)/2;
                }
            }else{
                end = start;
                start = lastStart;
            }
        }
        if((variantsMap[chromosomeName][start].getChanginglength() < 0) &&  basement<=variantsMap[chromosomeName][start].getPosition()-variantsMap[chromosomeName][start].getChanginglength()){
            return variantsMap[chromosomeName][start].getPosition() + variantsMap[chromosomeName][start].getLastTotalChanged();//lastChangedPoint;
        }
        return basement+variantsMap[chromosomeName][start+1].getLastTotalChanged();//allChanged;
    }else{
        return basement;
    }
}

void readGffFile (const std::string& filePath, std::map<std::string, std::vector<Transcript> >& transcriptHashSet){
    std::string cdsParentRegex="([\\s\\S]*)Parent=([\\s\\S]*?)[;,][\\s\\S]*$";
    readGffFile (filePath, transcriptHashSet, cdsParentRegex);
}

void readGffFile (const std::string& filePath, std::map<std::string, std::vector<Transcript> >& transcriptHashSet, std::string& cdsParentRegex){
    std::map<std::string, Transcript> transcriptHashMap;
    std::ifstream infile(filePath);
    if( ! infile.good()){
        std::cerr << "error in opening GFF/GTF file " << filePath << std::endl;
        exit (1);
    }
    std::regex reg("^(\\S*)\t([\\s\\S]*)\tCDS\t(\\S*)\t(\\S*)\t(\\S*)\t(\\S*)\t(\\S*)\t"+cdsParentRegex);
    std::string line="";
    while (std::getline(infile, line)){
        std::smatch match;
        regex_search(line, match, reg);

        if( match.empty() || line[0]=='#'){
        }else{
            int start = stoi(match[3]);
            int end = stoi(match[4]);
            if(start>end){
                int temp=start;
                start = end;
                end = temp;
            }
            std::string information = match[9];
            if(transcriptHashMap.find(information) != transcriptHashMap.end() ){
            }else{
                std::string chromosomeName = match[1];
                STRAND strand;
                if( match[6].compare("-") == 0){
                    strand = NEGATIVE;
                }else{
                    strand = POSITIVE;
                }
                Transcript transcript1 (information, chromosomeName, strand);
                transcriptHashMap[information] = transcript1;
            }
            Cds cds(start, end);
            //cds.setTranscript(transcriptHashMap[information]);
            transcriptHashMap[information].addCds(cds);
        }
    }

    for (std::map<std::string, Transcript>::iterator it=transcriptHashMap.begin(); it!=transcriptHashMap.end(); ++it){
        if(transcriptHashSet.find(it->second.getChromeSomeName()) == transcriptHashSet.end() ){
            transcriptHashSet[it->second.getChromeSomeName()]=std::vector<Transcript>();
        }
        it->second.updateInfor();
        transcriptHashSet[it->second.getChromeSomeName()].push_back(it->second);
    }
}


void annotationLiftOver(std::map<std::string, std::vector<Transcript> >& refTranscriptHashSet,
                        std::map<std::string, Transcript >& targetTranscriptHashMap,
                        std::map<std::string, std::vector<Variant> >& variantsMap,
                        std::map<std::string, Fasta>& targetGenome, std::map<std::string, Fasta>& referenceGenome,
                        std::string chromosome){
    std::set<std::string> chromosomes;
    chromosomes.insert(chromosome);
    annotationLiftOver( refTranscriptHashSet,
                        targetTranscriptHashMap,
                        variantsMap,
                        targetGenome, referenceGenome,
                        chromosomes);
}

void annotationLiftOver(std::map<std::string, std::vector<Transcript> >& refTranscriptHashSet,
                        std::map<std::string, Transcript >& targetTranscriptHashMap,
                        std::map<std::string, std::vector<Variant> >& variantsMap,
                        std::map<std::string, Fasta>& targetGenome, std::map<std::string, Fasta>& referenceGenome){
    std::set<std::string> chromosomes;
    for(std::map<std::string,std::vector<Transcript> >::iterator it1=refTranscriptHashSet.begin(); it1!=refTranscriptHashSet.end(); it1++){
        std::string chromosome = it1->first;
        chromosomes.insert(chromosome);
    }
    annotationLiftOver(refTranscriptHashSet, targetTranscriptHashMap, variantsMap, targetGenome, referenceGenome, chromosomes);
}

void annotationLiftOver(std::map<std::string, std::vector<Transcript> >& refTranscriptHashSet,
                        std::map<std::string, Transcript >& targetTranscriptHashMap,
                        std::map<std::string, std::vector<Variant> >& variantsMap,
                        std::map<std::string, Fasta>& targetGenome, std::map<std::string, Fasta>& referenceGenome,
                        std::set<std::string>& chromosomes){
    NucleotideCodeSubstitutionMatrix nucleotideCodeSubstitutionMatrix;
    for(std::map<std::string,std::vector<Transcript> >::iterator it1=refTranscriptHashSet.begin(); it1!=refTranscriptHashSet.end(); it1++){
        std::string chromosome = it1->first;
        if( chromosomes.find(chromosome)!=chromosomes.end() ){
            for( std::vector<Transcript>::iterator it2=it1->second.begin(); it2!=it1->second.end(); it2++){
                Transcript referenceTranscript = (*it2);
                Transcript targetTranscript((*it2).getName(), chromosome, (*it2).getStrand());
                for( std::vector<Cds>::iterator it3=referenceTranscript.getCdsVector().begin(); it3!=referenceTranscript.getCdsVector().end();it3++ ){
                    int liftStart = getChangedFromBasement(chromosome, (*it3).getStart(), variantsMap);
                    int liftEnd = getChangedFromBasement(chromosome, (*it3).getEnd(), variantsMap);
                    Cds cds(liftStart, liftEnd);
                    targetTranscript.addCds(cds);
                }
                targetTranscript.setSource("LIFTOVER");
                targetTranscript.updateInfor(targetGenome);
                checkOrfState(targetTranscript, targetGenome, nucleotideCodeSubstitutionMatrix);
                targetTranscriptHashMap[targetTranscript.getName()]=targetTranscript;
            }
        }
    }
}

void checkOrfState( Transcript& targetTranscript, Transcript &referenceTranscript,
                    std::map<std::string, Fasta>& targetGenome,
                    std::map<std::string, Fasta>& referenceGenome,
                    NucleotideCodeSubstitutionMatrix& nucleotideCodeSubstitutionMatrix){
    std::stringstream metaInformation;
    bool orfShift = false;
    if( ifSpliceSitesOk(targetTranscript, referenceTranscript, targetGenome, referenceGenome, nucleotideCodeSubstitutionMatrix) ){
        metaInformation << "_spliceSitesConserved";
    }else{
        metaInformation << "_spliceSitesDestroyed";
        orfShift = true;
    }
    std::string cdsSequenceString = targetTranscript.getCdsSequence();
    if (cdsSequenceString.length() < 3) {
        metaInformation << "_exonLengthLessThan3";
        orfShift = true;
    } else {
        metaInformation << "_exonLengthMoreThan3";
        if (ifLengthDivisibleByThree(cdsSequenceString)) {
            metaInformation << "_exonLengthIsDivisibleBy3";
        } else {
            metaInformation << "_exonLengthIsNotMultipleOf3";
            orfShift = true;
        }
        if (ifNewStopCodon(cdsSequenceString, nucleotideCodeSubstitutionMatrix)) {
            metaInformation << "_prematureStopCodon";
            orfShift = true;
        } else {
            metaInformation << "_noPrematureStopCodon";
        }
        if (ifEndWithStopCodon(cdsSequenceString, nucleotideCodeSubstitutionMatrix)) {
            metaInformation << "_endWithStopCodon";
        } else {
            metaInformation << "_notEndWithStopCodon";
            orfShift = true;
        }
        if (ifStartWithStartCodon(cdsSequenceString, nucleotideCodeSubstitutionMatrix)) {
            metaInformation << "_startWithStartCodon";
        } else {
            metaInformation << "_notStartWithStartCodon";
            orfShift = true;
        }
    }
    if( !orfShift ){
        metaInformation << "_ConservedFunction";
    }
    metaInformation << "_reference" << referenceTranscript.getStart() << "-" << referenceTranscript.getEnd();
    metaInformation << "_local" << targetTranscript.getStart() << "-" << targetTranscript.getEnd();
    if( POSITIVE ==  targetTranscript.getStrand()){
        metaInformation << "_positive";
    }else{
        metaInformation << "_negative";
    }
    std::string tempMetaInformation = metaInformation.str();
    targetTranscript.setMetaInformation(tempMetaInformation);
    targetTranscript.setIfOrfShift(orfShift);
}

void checkOrfState( Transcript& targetTranscript,
                    std::map<std::string, Fasta>& targetGenome,
                    NucleotideCodeSubstitutionMatrix& nucleotideCodeSubstitutionMatrix){
    std::stringstream metaInformation;
    bool orfShift = false;
    if( ifSpliceSitesOk(targetTranscript,targetGenome, nucleotideCodeSubstitutionMatrix) ){
        metaInformation << "_spliceSitesConserved";
    }else{
        metaInformation << "_spliceSitesDestroyed";
        orfShift = true;
    }
    std::string cdsSequenceString = targetTranscript.getCdsSequence();
    if (cdsSequenceString.length() < 3) {
        metaInformation << "_exonLengthLessThan3";
        orfShift = true;
    } else {
        metaInformation << "_exonLengthMoreThan3";
        if (ifLengthDivisibleByThree(cdsSequenceString)) {
            metaInformation << "_exonLengthIsDivisibleBy3";
        } else {
            metaInformation << "_exonLengthIsNotMultipleOf3";
            orfShift = true;
        }
        if (ifNewStopCodon(cdsSequenceString, nucleotideCodeSubstitutionMatrix)) {
            metaInformation << "_prematureStopCodon";
            orfShift = true;
        } else {
            metaInformation << "_noPrematureStopCodon";
        }
        if (ifEndWithStopCodon(cdsSequenceString, nucleotideCodeSubstitutionMatrix)) {
            metaInformation << "_endWithStopCodon";
        } else {
            metaInformation << "_notEndWithStopCodon";
            orfShift = true;
        }
        if (ifStartWithStartCodon(cdsSequenceString, nucleotideCodeSubstitutionMatrix)) {
            metaInformation << "_startWithStartCodon";
        } else {
            metaInformation << "_notStartWithStartCodon";
            orfShift = true;
        }
    }
    if( !orfShift ){
        metaInformation << "_ConservedFunction";
    }
    metaInformation << "_local" << targetTranscript.getStart() << "-" << targetTranscript.getEnd();
    if( POSITIVE ==  targetTranscript.getStrand()){
        metaInformation << "_positive";
    }else{
        metaInformation << "_negative";
    }
    std::string tempMetaInformation = metaInformation.str();
    targetTranscript.setMetaInformation(tempMetaInformation);
    targetTranscript.setIfOrfShift(orfShift);
}

bool ifSpliceSitesOk(Transcript& targetTranscript, Transcript &referenceTranscript,
                     std::map<std::string, Fasta>& targetGenome,
                     std::map<std::string, Fasta>& referenceGenome, NucleotideCodeSubstitutionMatrix& nucleotideCodeSubstitutionMatrix) {
    bool ifSelectTaur10 = true;
    for (size_t i = 1; i < referenceTranscript.getCdsVector().size(); i++) {
        int le;
        int ts;
        std::string s1;
        std::string s2;
        int let;
        int tst;
        std::string s1t;
        std::string s2t;
        if (referenceTranscript.getStrand() == POSITIVE) {
            le = referenceTranscript.getCdsVector()[i-1].getEnd();
            ts = referenceTranscript.getCdsVector()[i].getStart();
            s1 = getSubsequence(referenceGenome, referenceTranscript.getChromeSomeName(), le + 1, le + 2, POSITIVE);
            s2 = getSubsequence(referenceGenome,referenceTranscript.getChromeSomeName(), ts - 2, ts - 1, POSITIVE);

            let = targetTranscript.getCdsVector()[i-1].getEnd();
            tst = targetTranscript.getCdsVector()[i].getStart();

            s1t = getSubsequence(targetGenome, targetTranscript.getChromeSomeName(), let + 1, let + 2, POSITIVE);
            s2t = getSubsequence(targetGenome, targetTranscript.getChromeSomeName(), tst - 2, tst - 1, POSITIVE);
        } else {
            le = referenceTranscript.getCdsVector()[i-1].getEnd();
            ts = referenceTranscript.getCdsVector()[i].getStart();
            s1 = getSubsequence(referenceGenome, referenceTranscript.getChromeSomeName(), ts - 2, ts - 1, NEGATIVE);
            s2 = getSubsequence(referenceGenome, referenceTranscript.getChromeSomeName(), le + 1, le + 2, NEGATIVE);

            let = targetTranscript.getCdsVector()[i - 1].getEnd();
            tst = targetTranscript.getCdsVector()[i].getStart();
            s1t = getSubsequence(targetGenome, targetTranscript.getChromeSomeName(), tst - 2, tst - 1, NEGATIVE);
            s2t = getSubsequence(targetGenome, targetTranscript.getChromeSomeName(), let + 1, let + 2, NEGATIVE);
        }

        if( ! checkSpliceSites(s1, s2, s1t, s2t, nucleotideCodeSubstitutionMatrix)  ){
            return false;
        }
    }
    return ifSelectTaur10;
}

bool ifSpliceSitesOk(Transcript& targetTranscript,
                     std::map<std::string, Fasta>& targetGenome,
                      NucleotideCodeSubstitutionMatrix& nucleotideCodeSubstitutionMatrix) {
    bool ifSelectTaur10 = true;
    for (size_t i = 1; i < targetTranscript.getCdsVector().size(); i++) {
        int let;
        int tst;
        std::string s1t;
        std::string s2t;
        if (targetTranscript.getStrand() == POSITIVE) {
            let = targetTranscript.getCdsVector()[i-1].getEnd();
            tst = targetTranscript.getCdsVector()[i].getStart();
            s1t = getSubsequence(targetGenome, targetTranscript.getChromeSomeName(), let + 1, let + 2, POSITIVE);
            s2t = getSubsequence(targetGenome, targetTranscript.getChromeSomeName(), tst - 2, tst - 1, POSITIVE);
        } else {
            let = targetTranscript.getCdsVector()[i - 1].getEnd();
            tst = targetTranscript.getCdsVector()[i].getStart();
            s1t = getSubsequence(targetGenome, targetTranscript.getChromeSomeName(), tst - 2, tst - 1, NEGATIVE);
            s2t = getSubsequence(targetGenome, targetTranscript.getChromeSomeName(), let + 1, let + 2, NEGATIVE);
            //std::cerr << let << "\t" << tst << "\t" << s1t << "\t" << s2t << std::endl;
        }

        if( ! checkSpliceSites( s1t, s2t, nucleotideCodeSubstitutionMatrix)  ){
            return false;
        }
    }
    return ifSelectTaur10;
}


bool checkSpliceSites( std::string& s1t, std::string& s2t, NucleotideCodeSubstitutionMatrix& nucleotideCodeSubstitutionMatrix){
    s1t = gtIUPACcodesTranslation(s1t);
    s2t = agIUPACcodesTranslation(s2t);
    std::vector<std::string> dornors;
    std::vector<std::string> acceptors;
//https://metasystems.riken.jp/conjoing/faqs
//    GT-AG, GC-AG, CT-AG, GG-AG, GA-AG, CA-AG, GT-CG, GT-GG,
// TG-AG, GT-TG, AT-AG, TC-AG, GT-AC, CC-AG, TT-AG,
// GT-AA, AA-AG, GT-AT, GT-GC, AC-AG, AG-AG, CT-AT, GT-CC,
// GA-CT, GT-CT, GT-GA, CC-GC, GG-GT, GT-GT, GT-TC, CT-TT, and GT-TT.
    dornors.push_back("GT");    acceptors.push_back("AG");
    dornors.push_back("GC");    acceptors.push_back("AG");
    dornors.push_back("CT");    acceptors.push_back("AG");
    dornors.push_back("GG");    acceptors.push_back("AG");
    dornors.push_back("GA");    acceptors.push_back("AG");
    dornors.push_back("CA");    acceptors.push_back("AG");
    dornors.push_back("GT");    acceptors.push_back("CG");
    dornors.push_back("GT");    acceptors.push_back("GG");
    dornors.push_back("TG");    acceptors.push_back("AG");
    dornors.push_back("GT");    acceptors.push_back("TG");
    dornors.push_back("AT");    acceptors.push_back("AG");
    dornors.push_back("TC");    acceptors.push_back("AG");
    dornors.push_back("GT");    acceptors.push_back("AC");
    dornors.push_back("CC");    acceptors.push_back("AG");
    dornors.push_back("TT");    acceptors.push_back("AG");
    dornors.push_back("GT");    acceptors.push_back("AA");
    dornors.push_back("AA");    acceptors.push_back("AG");
    dornors.push_back("GT");    acceptors.push_back("AT");
    dornors.push_back("GT");    acceptors.push_back("GC");
    dornors.push_back("AC");    acceptors.push_back("AG");
    dornors.push_back("AG");    acceptors.push_back("AG");
    dornors.push_back("CT");    acceptors.push_back("AT");
    dornors.push_back("GT");    acceptors.push_back("CC");
    dornors.push_back("GA");    acceptors.push_back("CT");
    dornors.push_back("GT");    acceptors.push_back("CT");
    dornors.push_back("GT");    acceptors.push_back("GA");
    dornors.push_back("CC");    acceptors.push_back("GC");
    dornors.push_back("GG");    acceptors.push_back("GT");
    dornors.push_back("GT");    acceptors.push_back("GT");
    dornors.push_back("GT");    acceptors.push_back("TC");
    dornors.push_back("CT");    acceptors.push_back("TT");
    dornors.push_back("GT");    acceptors.push_back("TT");

    //https://www.ncbi.nlm.nih.gov/pmc/articles/PMC84117/
    dornors.push_back("AT");    acceptors.push_back("AC");
    //https://www.ncbi.nlm.nih.gov/pmc/articles/PMC113136/
    dornors.push_back("CT");    acceptors.push_back("GC");

    //https://www.ncbi.nlm.nih.gov/pmc/articles/PMC84117/
    dornors.push_back("AT");    acceptors.push_back("AA");

//tair10
    // dornors.push_back("GT");    acceptors.push_back("TA");
    // dornors.push_back("TA");    acceptors.push_back("AG");
    // dornors.push_back("TA");    acceptors.push_back("AT");
    // dornors.push_back("AT");    acceptors.push_back("AT");
    // dornors.push_back("GC");    acceptors.push_back("CT");
    // dornors.push_back("CG");    acceptors.push_back("AG");
    // dornors.push_back("TG");    acceptors.push_back("AT");


    for( size_t size =0; size<dornors.size(); ++size ){
        std::string dornor=dornors[size];
        std::string acceptor=acceptors[size];
        std::set<std::string> allDornors;
        std::set<std::string> allAcceptors;
        nucleotideCodeSubstitutionMatrix.getAllPossibleWithIupac(dornor, allDornors);
        nucleotideCodeSubstitutionMatrix.getAllPossibleWithIupac(acceptor, allAcceptors);
        if( allDornors.find(s1t)!=allDornors.end() && allAcceptors.find(s2t)!=allAcceptors.end() ){
            return true;
        }
    }
//    std::cerr << s1t << " " << s2t <<std::endl;
    return false;
}

bool checkSpliceSites(std::string& s1, std::string& s2, std::string& s1t, std::string& s2t, NucleotideCodeSubstitutionMatrix& nucleotideCodeSubstitutionMatrix){
    s2 = agIUPACcodesTranslation(s2);
    s1 = gtIUPACcodesTranslation(s1);
    s2t = agIUPACcodesTranslation(s2t);
    s1t = gtIUPACcodesTranslation(s1t);

    if ( checkSpliceSites(s1t, s2t, nucleotideCodeSubstitutionMatrix) ) {
        return true;
    } else {
        std::set<std::string> allS1s;
        nucleotideCodeSubstitutionMatrix.getAllPossibleWithIupac(s1, allS1s);
        if (allS1s.find(s1t) == allS1s.end()) {
            std::cout << s1t << "   " << s2t << std::endl;
            return false;
        }
        std::set<std::string> allS2s;
        nucleotideCodeSubstitutionMatrix.getAllPossibleWithIupac(s2, allS2s);
        if (allS2s.find(s2t) == allS2s.end() ) {
            //std::cout << s1t << "  " << s2t << std::endl;
            return false;
        }
    }
    return true;
}

std::string agIUPACcodesTranslation(std::string& ag) {
    std::string agString = ag;
    if ( agString.length()==2 && ('A' == ag[0] || 'R' == ag[0] || 'W' == ag[0]
        || 'M' == ag[0] || 'D' == ag[0]
        || 'H' == ag[0] || 'V' == ag[0] || 'N' == ag[0])
        && ('G' == ag[1] || 'K' == ag[1]
        || 'R' == ag[1] || 'S' == ag[1]
        || 'B' == ag[1] || 'D' == ag[1]
        || 'V' == ag[1] || 'N' == ag[1])) {
        agString = "AG";
    }
    return agString;
}

std::string gtIUPACcodesTranslation(std::string& gt) {
    std::string gtString = gt;
    if ( gtString.length()==2 && ('G' == gtString[0] || 'K' == gtString[0]
        || 'R' == gtString[0] || 'S' == gtString[0]
        || 'B' == gtString[0] || 'D' == gtString[0]
        || 'V' == gtString[0] || 'N' == gtString[0])
        && ('T' == gtString[1] || 'K' == gtString[1]
        || 'Y' == gtString[1]
        || 'W' == gtString[1]
        || 'U' == gtString[1]
        || 'B' == gtString[1]
        || 'D' == gtString[1]
        || 'H' == gtString[1] || 'N' == gtString[1])) {
        gtString = "GT";
    }
    return gtString;
}
bool ifLengthDivisibleByThree(std::string& sequence){
    return 0 == sequence.length() % 3;
}

bool ifNewStopCodon(std::string& cdsSequence, NucleotideCodeSubstitutionMatrix& nucleotideCodeSubstitutionMatrix) {
    for (size_t j = 0; j < cdsSequence.length() - 3; j += 3) {
        std::string threeNaInFrame = cdsSequence.substr(j,3);
        if (nucleotideCodeSubstitutionMatrix.getMustStopCodons().find(threeNaInFrame) !=
                nucleotideCodeSubstitutionMatrix.getMustStopCodons().end()) {// new stop codon
            //std::cout << "443 " << cdsSequence << " " << threeNaInFrame << " " << j << std::endl;
            return true;
        }
    }
    return false;
}

bool ifEndWithStopCodon(std::string& cdsSequence, NucleotideCodeSubstitutionMatrix& nucleotideCodeSubstitutionMatrix) {
    std::string threeNaInFrame = cdsSequence.substr(cdsSequence.length() - 3, 3);
    //std::cout << "472 " << threeNaInFrame << std::endl << cdsSequence << std::endl;
    return nucleotideCodeSubstitutionMatrix.getPossibleStopCodons().find(threeNaInFrame)
           != nucleotideCodeSubstitutionMatrix.getPossibleStopCodons().end();
}

bool ifStartWithStartCodon(std::string& cdsSequence, NucleotideCodeSubstitutionMatrix& nucleotideCodeSubstitutionMatrix) {
    std::string threeNaInFrame = cdsSequence.substr(0, 3);
    //std::cout << "479 " << threeNaInFrame << std::endl << cdsSequence << std::endl;
    return nucleotideCodeSubstitutionMatrix.getPossibleStartCodons().find(threeNaInFrame)
        != nucleotideCodeSubstitutionMatrix.getPossibleStartCodons().end();
}
// orf checking related function end

//translation related begin
// this function is private here
void na2aaLocal( std::string& seq, std::stringstream& ssSb, NucleotideCodeSubstitutionMatrix& nucleotideCodeSubstitutionMatrix){
    seq=songStrRemoveBlank(seq);
    songToUpCase(seq);
    std::string pattern="U";
    std::string pattern2="T";
    seq = songStrReplaceAll(seq, pattern, pattern2);

    for (size_t jj = 0; jj <= seq.length() - 3; jj += 3) {
        std::string codeSeq = seq.substr(jj, 3);
        BEGINMIDDLEEND position=MIDDLE;
        if( jj==0  ){
            position=BEGIN;
        }else if( jj==seq.length() - 3  ){
            position=END;
        }
        ssSb << nucleotideCodeSubstitutionMatrix.getGeneticCode( codeSeq, position );
    }
}

std::string nA2AA(std::string& seq){
    NucleotideCodeSubstitutionMatrix nucleotideCodeSubstitutionMatrix;
    std::string pattern="-";
    std::string pattern2="";
    seq=songStrReplaceAll(seq, pattern, pattern2);
    std::stringstream ssSb;
    na2aaLocal( seq, ssSb, nucleotideCodeSubstitutionMatrix);
    return ssSb.str();
}

std::string nA2AANoDeletedIndel(std::string& seq){
    NucleotideCodeSubstitutionMatrix nucleotideCodeSubstitutionMatrix;
    std::stringstream ssSb;
    na2aaLocal( seq, ssSb, nucleotideCodeSubstitutionMatrix);
    return ssSb.str();
}

std::string nA2AA(std::string& seq, NucleotideCodeSubstitutionMatrix& nucleotideCodeSubstitutionMatrix){
    std::string pattern="-";
    std::string pattern2="";
    seq=songStrReplaceAll(seq, pattern, pattern2);
    std::stringstream ssSb;
    na2aaLocal( seq, ssSb, nucleotideCodeSubstitutionMatrix);
    return ssSb.str();
}

std::string nA2AANoDeletedIndel(std::string& seq, NucleotideCodeSubstitutionMatrix& nucleotideCodeSubstitutionMatrix){
    std::stringstream ssSb;
    na2aaLocal( seq, ssSb, nucleotideCodeSubstitutionMatrix);
    return ssSb.str();
}
//translation related end

bool if_file_exists (const std::string& name) {
    struct stat buffer;
    return (stat (name.c_str(), &buffer) == 0);
}

long GetFileSize(std::string& filename) {
    struct stat stat_buf;
    int rc = stat(filename.c_str(), &stat_buf);
    return rc == 0 ? stat_buf.st_size : -1;
}
