// =====================================================================================
//
//       Filename:  reAnnotationAndMsa.cpp
//
//    Description:
//
//        Version:  1.0
//        Created:  04/12/2017 02:13:56 PM
//       Revision:  none
//       Compiler:  g++
//
//         Author:  Baoxing Song (songbx.me), song@mpipz.mpg.de
//        Company:  MPIPZ
//
// =====================================================================================
 
#include <stdlib.h>
#include "reAnnotationAndMsa.h"
#include "alignNeedlemanForTranscript.h"
#include <iostream>
#include <pthread.h>
#include "model.h"
#include <unistd.h>
#include <fstream>
#include <sstream>
#include "parameters.h"
#include <regex>
// parallele related things begin
pthread_mutex_t mutex;
struct myPara{
    std::string chromosomeName;
    std::string accessionId;
    Transcript* transcript;
    MyThreadCount* myThreadCount;
    std::map<std::string, Fasta>* targetGenome;
    myPara(std::string _accessionId, Transcript& _transcript, MyThreadCount& _myThreadCount, std::map<std::string,
            Fasta>& _targetGenome, std::string _chromosomeName){
        accessionId=_accessionId;
        transcript=&_transcript;
        myThreadCount=&_myThreadCount;
        targetGenome=&_targetGenome;
        chromosomeName=_chromosomeName;
    }
};

MyThreadCount myThreadCount;
pthread_attr_t attr;

void * transcriptRealignment( void *arg );
// parallele related things end


std::map<std::string, std::map<std::string, Transcript> > accessionsTargetTranscriptHashSet; // accessionsId transcriptName, transcript
NucleotideCodeSubstitutionMatrix nucleotideCodeSubstitutionMatrix;
std::map<std::string, Fasta> referenceGenome;
std::map<std::string, std::vector<Transcript> > referenceTranscriptHashSet;
std::string referenceID= "Col";
void reAnnotationSingleLine( std::string& referenceGenomeFilePath, std::string referenceGffFilePath, std::string sdiFile,
                    std::map<std::string, std::vector<Gene> >& genes, int maxThread, std::string& regex, std::string& regexG){
    pthread_attr_init( &attr );
    pthread_attr_setscope( &attr, PTHREAD_SCOPE_SYSTEM );
    std::cout << "62" << std::endl;
    readFastaFile(referenceGenomeFilePath, referenceGenome);
    std::cout << "64 " << referenceGffFilePath << " " << regex << std::endl;
    readGffFile (referenceGffFilePath, referenceTranscriptHashSet, regex);
    std::cout << "66" << std::endl;
    std::map<std::string, std::vector<Variant> > variantsMaps;
    readSdiFile (sdiFile,  variantsMaps);
    std::cout << "69" << std::endl;
    std::map<std::string, Fasta> targetGenome;
    std::cout << "71" << std::endl;
    getPseudoGenomeSequence(referenceGenome, variantsMaps, targetGenome);
    std::cout << "73" << std::endl;
    std::map<std::string, std::vector<Transcript> > targetTranscriptsTHashSet;
    annotationLiftOver( referenceTranscriptHashSet, targetTranscriptsTHashSet, variantsMaps, targetGenome, referenceGenome);
    std::cout << "76" << std::endl;
    std::string accessionId = "temp";
    std::cout << "78" << std::endl;
    for( std::map<std::string, Fasta>::iterator it1=referenceGenome.begin(); it1!=referenceGenome.end(); it1++ ){
        if( referenceTranscriptHashSet.find(it1->first)!= referenceTranscriptHashSet.end() ){
            std::string chromosomeName = it1->first;
            if( targetGenome.find(chromosomeName)!=targetGenome.end() ){
            }else{
                std::cerr << "reAnnotationAndMsa line37, there are bugs" << std::endl;
            }
            for( std::vector<Transcript>::iterator it3=referenceTranscriptHashSet[chromosomeName].begin();
                    it3!=referenceTranscriptHashSet[chromosomeName].end();++it3 ){
                (*it3).updateInfor(referenceGenome);
                accessionsTargetTranscriptHashSet[referenceID][(*it3).getName()] = (*it3);
            }
            for( std::vector<Transcript>::iterator it3=targetTranscriptsTHashSet[chromosomeName].begin(); it3!=targetTranscriptsTHashSet[chromosomeName].end(); it3++) {
                accessionsTargetTranscriptHashSet[accessionId][(*it3).getName()] = (*it3);
                if ((*it3).getIfOrfShift()) {

                    bool isThisThreadUnrun = true;
                    while (isThisThreadUnrun) {
                        if (myThreadCount.getCount() < maxThread) {
                            //std::cout << 62 << std::endl;
                            pthread_t thread;
                            struct myPara *para = new struct myPara(accessionId, (*it3), myThreadCount, targetGenome,
                                                                    chromosomeName);
                            //std::cout << "111 " << accessionId << " " << chromosomeName << std::endl;
                            pthread_mutex_lock(&mutex);
                            std::cout << "98 " << para->accessionId << " " << para->chromosomeName << std::endl;
                            myThreadCount.plusOne();
                            pthread_mutex_unlock(&mutex);
                            pthread_create(&thread, &attr, transcriptRealignment, (void *) para);
                            //pthread_join(&thread, NULL);
                            //std::cout << "70" << std::endl;
                            isThisThreadUnrun = false;
                            break;
                        } else {
                            //std::cout << "77" << std::endl;
                            usleep(10);
                            //std::cout << "79" << std::endl;
                        }
                    }
                }
            }
            while(myThreadCount.hasNext()){// wait for all the thread
                usleep(50);
            }
        }
    }
    pthread_mutex_destroy(&mutex);

    std::map< std::string, std::map<std::string, Gene > > tempGenes;
    for( std::map<std::string, Transcript>::iterator it=accessionsTargetTranscriptHashSet[accessionId].begin();
            it!=accessionsTargetTranscriptHashSet[accessionId].end(); ++it){
        std::string transcriptId = it->first;
        std::string chromosome = it->second.getChromeSomeName();
        if(tempGenes.find(chromosome) == tempGenes.end() ){
            tempGenes[chromosome]=std::map<std::string, Gene >();
        }
        std::cerr << "129" << std::endl;
        std::regex reg(regexG);
        std::smatch match;
        regex_search(transcriptId, match, reg);
        std::cerr << "133" << std::endl;
        if( match.empty()  ){

        }else{
            std::string geneName = match[1];
            if( tempGenes[chromosome].find(geneName) == tempGenes[chromosome].end() ){
                Gene gene(geneName, it->second.getStrand());
                tempGenes[chromosome].insert( std::pair<std::string, Gene>(geneName, gene) );
            }
            tempGenes[chromosome][geneName].checkOverLapAndAddTranscript(it->second);
        }
    }
    for( std::map< std::string, std::map<std::string, Gene > >::iterator it=tempGenes.begin();
        it!=tempGenes.end(); ++it ){
        if( genes.find(it->first) == genes.end()  ){
            genes[it->first]=std::vector<Gene>();
        }
        for( std::map<std::string, Gene >::iterator it2=it->second.begin();
            it2!=it->second.end(); ++it2){
            genes[it->first].push_back(it2->second);
        }
    }
}

void reAnnotationAndMsa( std::string& referenceGenomeFilePath,
    std::string referenceGffFilePath, std::map<std::string, std::string> sdiFilePaths, int maxThread){

    pthread_attr_init( &attr );
    pthread_attr_setscope( &attr, PTHREAD_SCOPE_SYSTEM );

    readFastaFile(referenceGenomeFilePath, referenceGenome);
    readGffFile (referenceGffFilePath, referenceTranscriptHashSet);
    for( std::map<std::string, Fasta>::iterator it1=referenceGenome.begin(); it1!=referenceGenome.end(); it1++ ){
        if( referenceTranscriptHashSet.find(it1->first)!= referenceTranscriptHashSet.end() ){
            // check chromosome by chromosome, for RAM saving purpose. Do the keys of following maps are name of lines
            std::map<std::string, Fasta> accessionsTargetSequences;// the Fasta of current chromosome for each accession

            std::map<std::string, std::vector<Variant> > accessionsVariantsMaps; //accessions variants
            std::string chromosomeName = it1->first;
            for( std::map<std::string, std::string>::iterator it2=sdiFilePaths.begin(); it2!=sdiFilePaths.end(); it2++ ){
                std::string accessionId = it2->first;
                std::map<std::string, std::vector<Variant> > variantsTMaps;
                readSdiFile (it2->second,  variantsTMaps, chromosomeName);// sdiFile, variantsTmaps chromosomeName

                if( variantsTMaps.find(chromosomeName)!=variantsTMaps.end() ){
                    accessionsVariantsMaps[accessionId] = variantsTMaps[chromosomeName];
                }else{
                    accessionsVariantsMaps[accessionId] = std::vector<Variant>();
                }

                //std::cout << "25" << std::endl;
                std::map<std::string, Fasta> targetGenome;
                getPseudoGenomeSequence(referenceGenome, variantsTMaps, targetGenome, it1->first);
                if( targetGenome.find(chromosomeName)!=targetGenome.end() ){
                    accessionsTargetSequences[chromosomeName]=targetGenome.begin()->second;
                }else{
                    std::cerr << "reAnnotationAndMsa line37, there are bugs" << std::endl;
                }
                //std::cout << "39 " << targetGenome.begin()->second.getSequence().length() << std::endl;

                std::map<std::string, std::vector<Transcript> > targetTranscriptsTHashSet;
                annotationLiftOver( referenceTranscriptHashSet, targetTranscriptsTHashSet, variantsTMaps, targetGenome, referenceGenome, chromosomeName);
                //std::cout << "43 " << std::endl;
                accessionsTargetTranscriptHashSet[it2->first]=std::map<std::string, Transcript>();


                for( std::vector<Transcript>::iterator it3=targetTranscriptsTHashSet[chromosomeName].begin(); it3!=targetTranscriptsTHashSet[chromosomeName].end(); it3++){
                    accessionsTargetTranscriptHashSet[accessionId][(*it3).getName()]=(*it3);

                    if( (*it3).getIfOrfShift() && accessionId.compare(referenceID)!=0 ){

                        bool isThisThreadUnrun=true;
                        while(isThisThreadUnrun){
                            if(myThreadCount.getCount() < maxThread){
                                //std::cout << 62 << std::endl;
                                pthread_t thread;
                                struct myPara * para= new struct myPara(accessionId, (*it3), myThreadCount, targetGenome, chromosomeName);
                                //std::cout << "111 " << accessionId << " " << chromosomeName << std::endl;
                                pthread_mutex_lock(&mutex);
                                //std::cout << "113 " << para->accessionId << " " << para->chromosomeName << std::endl;
                                myThreadCount.plusOne();
                                pthread_mutex_unlock(&mutex);
                                pthread_create(&thread, &attr, transcriptRealignment, (void *)para);
                                //pthread_join(&thread, NULL);
                                //std::cout << "70" << std::endl;
                                isThisThreadUnrun=false;
                                break;
                            }else{
                                //std::cout << "77" << std::endl;
                                usleep(10);
                                //std::cout << "79" << std::endl;
                            }
                        }
                    }else{
                        (*it3).cleanSaveRam();
                    }
                }
                while(myThreadCount.hasNext()){// wait for all the thread
                    usleep(50);
                }
            }
        }
    }
    pthread_mutex_destroy(&mutex);
}

void * transcriptRealignment( void *arg ){
    struct myPara *pstru = (struct myPara *)arg;
    Transcript* it3 =  pstru->transcript;
    std::string accessionId=pstru->accessionId;
    std::map<std::string, Fasta>& targetGenome = (*pstru->targetGenome);
    std::string chromosomeName=pstru->chromosomeName;
    pthread_mutex_lock(&mutex);
    std::cout << "146 " << accessionId << " " << chromosomeName << std::endl;
    pthread_mutex_unlock(&mutex);


   // pthread_mutex_lock(&mutex);
    //std::cout << accessionId << " " << (*it3).getName() << (*it3).getMetaInformation() << std::endl;
    //pthread_mutex_unlock(&mutex);


    std::string refGenomeSequence;
    pthread_mutex_lock(&mutex);
    refGenomeSequence = accessionsTargetTranscriptHashSet[referenceID][(*it3).getName()].getGeneomeSequence();
    pthread_mutex_unlock(&mutex);
   
   
    if( (*it3).getStrand() == POSITIVE ){
        size_t startCodonPosition=1;
        size_t stopCodonPosition=refGenomeSequence.length()-2;
        std::vector<SpliceSitePosition> splitSitePositions;
        if( (*it3).getCdsVector().size()>1 ){
            for( size_t i=1; i<(*it3).getCdsVector().size(); ++i ){
                SpliceSitePosition spliceSitePosition(
                        accessionsTargetTranscriptHashSet[referenceID][(*it3).getName()].getCdsVector()[i-1].getEnd()-accessionsTargetTranscriptHashSet[referenceID][(*it3).getName()].getStart()+2,
                        accessionsTargetTranscriptHashSet[referenceID][(*it3).getName()].getCdsVector()[i].getStart()-accessionsTargetTranscriptHashSet[referenceID][(*it3).getName()].getStart());
                splitSitePositions.push_back(spliceSitePosition);
            }
        }
        pthread_mutex_lock(&mutex);
       // std::cout << "167" << std::endl;
        std::string dna_b = getSubsequence(targetGenome, chromosomeName, (*it3).getStart()-refGenomeSequence.length(), (*it3).getEnd()+refGenomeSequence.length(), POSITIVE);
        //std::cout << "172" << std::endl;
        pthread_mutex_unlock(&mutex);
        
        NeedlemanWunschForTranscript nw (refGenomeSequence, dna_b, startCodonPosition, stopCodonPosition, splitSitePositions);
        
        pthread_mutex_lock(&mutex);
        nw.print_results();
        pthread_mutex_unlock(&mutex);

        Transcript targetTranscript((*it3).getName(), chromosomeName, (*it3).getStrand());
        int targetPosition=0;
        int referencePosition=0;
        //prepare for the new Cds After re-alignment begin
        std::vector<int> targetCdsStarts;
        std::vector<int> targetCdsEnds;
        for( size_t tp = 0; tp<nw.getAlignment_a().length(); ++tp ){
            if( nw.getAlignment_b()[tp] != '-' ){
                ++targetPosition;
            }
            if( nw.getAlignment_a()[tp] != '-' ){
                ++referencePosition;
                for( std::vector<Cds>::iterator it4=accessionsTargetTranscriptHashSet[referenceID][(*it3).getName()].getCdsVector().begin();
                     it4!=accessionsTargetTranscriptHashSet[referenceID][(*it3).getName()].getCdsVector().end(); it4++){
                    if( referencePosition+accessionsTargetTranscriptHashSet[referenceID][(*it3).getName()].getStart()-1 == (*it4).getStart() ){
                        if( nw.getAlignment_b()[tp] == '-' ){
                            targetCdsStarts.push_back((*it3).getStart()-refGenomeSequence.length() + targetPosition); //(*i3) is the target transcript, refGenomeSequence.length() is the extend length
                        }else{
                            targetCdsStarts.push_back((*it3).getStart()-refGenomeSequence.length() + targetPosition -1);
                        }
                    }
                    if( referencePosition+accessionsTargetTranscriptHashSet[referenceID][(*it3).getName()].getStart()-1 == (*it4).getEnd() ){
                        if( nw.getAlignment_b()[tp] == '-' ){
                            targetCdsEnds.push_back((*it3).getStart()-refGenomeSequence.length() + targetPosition);
                        }else{
                            targetCdsEnds.push_back((*it3).getStart()-refGenomeSequence.length() + targetPosition -1);
                        }
                    }
                }
            }
        }
        //prepare for the new Cds After re-alignment end
        for(size_t i5=0; i5<targetCdsStarts.size(); i5++ ){
            Cds cds(targetCdsStarts[i5], targetCdsEnds[i5]);
            targetTranscript.addCds(cds);
        }
        targetTranscript.updateInfor(targetGenome);
        checkOrfState( targetTranscript, accessionsTargetTranscriptHashSet[referenceID][(*it3).getName()],
                       targetGenome, referenceGenome, nucleotideCodeSubstitutionMatrix);
        pthread_mutex_lock(&mutex);
        accessionsTargetTranscriptHashSet[accessionId][targetTranscript.getName()]=targetTranscript; // update the map data structure with modified transcript
        std::cout << accessionId << " " << targetTranscript.getName() << " " << targetTranscript.getMetaInformation() << std::endl <<
                  targetTranscript.getGeneomeSequence() << std::endl <<
                  targetTranscript.getCdsSequence() << std::endl;
        pthread_mutex_unlock(&mutex);
        //targetTranscript.cleanSaveRam();// empty the transcriprt sequence and genome sequence, for RAM saving purpose
        
    }else{
        // pthread_mutex_lock(&mutex);
        // std::cout << accessionId << " " << (*it3).getName() << (*it3).getMetaInformation() << std::endl;
        // pthread_mutex_unlock(&mutex);
        size_t startCodonPosition=1;
        size_t stopCodonPosition=refGenomeSequence.length()-2;
        std::vector<SpliceSitePosition> spliceSitePositions;
        if( (*it3).getCdsVector().size()>1 ){
            for( size_t i=(*it3).getCdsVector().size()-1; i>0; --i ){
                SpliceSitePosition spliceSitePosition(
                        accessionsTargetTranscriptHashSet[referenceID][(*it3).getName()].getEnd()-accessionsTargetTranscriptHashSet[referenceID][(*it3).getName()].getCdsVector()[i].getStart()+2,
                        accessionsTargetTranscriptHashSet[referenceID][(*it3).getName()].getEnd()-accessionsTargetTranscriptHashSet[referenceID][(*it3).getName()].getCdsVector()[i-1].getEnd()-1);
                spliceSitePositions.push_back(spliceSitePosition);
                //std::cout << "135 " << spliteSitePosition.getDonorSplitSitePosition() << " " << spliteSitePosition.getacceptorSplitSitePosition() << std::endl;
            }
        }
        std::string dna_b = getSubsequence(targetGenome, chromosomeName, (*it3).getStart()-refGenomeSequence.length(), (*it3).getEnd()+refGenomeSequence.length(), NEGATIVE);
        NeedlemanWunschForTranscript nw (refGenomeSequence, dna_b, startCodonPosition, stopCodonPosition, spliceSitePositions);
        pthread_mutex_lock(&mutex);
        nw.print_results();
        pthread_mutex_unlock(&mutex);


        Transcript targetTranscript((*it3).getName(), chromosomeName, (*it3).getStrand());
        int targetPosition=0;
        int referencePosition=0;
        //prepare for the new Cds After re-alignment begin
        std::vector<int> targetCdsStarts;
        std::vector<int> targetCdsEnds;
        for( size_t tp = 0; tp<nw.getAlignment_a().length(); ++tp ){
            if( nw.getAlignment_b()[tp] != '-' ){
                ++targetPosition;
            }
            if( nw.getAlignment_a()[tp] != '-' ){
                ++referencePosition;
                for( std::vector<Cds>::iterator it4=accessionsTargetTranscriptHashSet[referenceID][(*it3).getName()].getCdsVector().begin();
                     it4!=accessionsTargetTranscriptHashSet[referenceID][(*it3).getName()].getCdsVector().end(); it4++){
                    if( referencePosition+accessionsTargetTranscriptHashSet[referenceID][(*it3).getName()].getStart()-1 == (*it4).getStart() ){
                        if( nw.getAlignment_b()[tp] == '-' ){
                            targetCdsStarts.push_back((*it3).getEnd() + refGenomeSequence.length() - targetPosition);
                        }else{
                            targetCdsStarts.push_back((*it3).getEnd() + refGenomeSequence.length() - targetPosition +1);
                        }
                    }
                    if( referencePosition+accessionsTargetTranscriptHashSet[referenceID][(*it3).getName()].getStart()-1 == (*it4).getEnd() ){
                        if( nw.getAlignment_b()[tp] == '-' ){
                            targetCdsEnds.push_back((*it3).getEnd() + refGenomeSequence.length() - targetPosition );
                        }else{
                            targetCdsEnds.push_back((*it3).getEnd() + refGenomeSequence.length() - targetPosition +1);
                        }
                    }
                }
            }
        }
        //prepare for the new Cds After re-alignment end
        for(size_t i5=0; i5<targetCdsStarts.size(); i5++ ){
            Cds cds(targetCdsStarts[i5], targetCdsEnds[i5]);
            targetTranscript.addCds(cds);
        }
        targetTranscript.updateInfor(targetGenome);
        checkOrfState( targetTranscript, accessionsTargetTranscriptHashSet[referenceID][(*it3).getName()],
                       targetGenome, referenceGenome, nucleotideCodeSubstitutionMatrix);
        pthread_mutex_lock(&mutex);
        accessionsTargetTranscriptHashSet[accessionId][targetTranscript.getName()]=targetTranscript; // update the map data structure with modified transcript
        std::cout << accessionId << " " << targetTranscript.getName() << " " << targetTranscript.getMetaInformation() << std::endl <<
                  targetTranscript.getGeneomeSequence() << std::endl <<
                  targetTranscript.getCdsSequence() << std::endl;
        pthread_mutex_unlock(&mutex);
        //targetTranscript.cleanSaveRam();// empty the transcriprt sequence and genome sequence, for RAM saving purpose
    }
    
    //(*it3).cleanSaveRam();
    delete(pstru);
    pthread_mutex_lock(&mutex);
    (*pstru->myThreadCount).countDown();
    pthread_mutex_unlock(&mutex);
    return NULL;
}


void runExonerate(std::string& cdsSequence, std::string& refSequence){
    std::cout << cdsSequence << std::endl << refSequence << std::endl;
    std::string tempFolder = createdTempFloder();
    std::string tempFile = tempFolder + "/" + generateUUID();
    std::string cdsFile = tempFile + "cds.fasta";
    std::string refFile = tempFile + "ref.fasta";

    std::ofstream ofile;
    ofile.open(cdsFile);
    ofile << ">cds" << std::endl << cdsSequence << std::endl;
    ofile.close();

    ofile.open(refFile);
    ofile << ">ref" << std::endl << refSequence << std::endl;
    ofile.close();

    std::string command = "exonerate --maxintron 30000 --model est2genome -i -10 --score 10 --bestn 1 --minintron 10 "+tempFile + "cds.fasta " + tempFile + "ref.fasta --showtargetgff true >" +tempFile;
    std::cout << command << std::endl;
    system(&command[0]);
    readExonerateEstResult(tempFile);
    std::string cleanFileCommand = "rm " + cdsFile + "; rm " + refFile + "; rm " + tempFile;
    //system(&cleanFileCommand[0]);

//    system("exonerate --bestn 1 --maxintron 30000 --intronpenalty -10 --model protein2genome --percent 10 --score 10 --minintron 10 ../../col_0/protein/$file ../genomeSequence/$file >$file");
    return ;
}

void readExonerateEstResult( std::string& fileLocation ){
    std::ifstream infile(fileLocation);
    std::regex reg("\\d\\s*:\\s*([\\w\\.\\s\\>-]+)\\s*:\\s*\\d");
    std::stringstream sequencestream;
    std::string line="";
    int lineNumber = 0;
    while ( std::getline(infile, line) ){
        std::smatch match;
        regex_search(line, match, reg);
        if( match.empty()   ){

        }else{
            lineNumber++;
            if( lineNumber % 2 == 0 ){
                sequencestream << match[1];
            }
        }
    }
    infile.close();
    std::stringstream metaInformation;
    std::string cdsSequenceString = sequencestream.str();
    std::regex vowel_re("\\s");
    cdsSequenceString=std::regex_replace(cdsSequenceString, vowel_re, "");
    std::regex vowel_re2("[a-z]+\\.+[a-z]+");
    cdsSequenceString=std::regex_replace(cdsSequenceString, vowel_re2, "");

    bool orfShift = false;
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
        std::ifstream infile(fileLocation);
        std::regex reg("^(\\S*)\t(.*)\texon\t(\\S*)\t(\\S*)\t(\\S*)\t\\+");
        std::string line="";
        while (std::getline(infile, line)){
            std::smatch match;
            regex_search(line, match, reg);

            if( match.empty() ){
            }else{
                int start = stoi(match[3]);
                int end = stoi(match[4]);
                if(start>end){
                    int temp=start;
                    start = end;
                    end = temp;
                }
                Cds cds(start, end);
                std::cout << start << " " << end << std::endl;
            }
        }
        infile.close();
    }
    std::cout << metaInformation.str() << std::endl;
    std::cout << cdsSequenceString << std::endl;
    return;
}



void readAugustusGff( std::string& fileLocation ){
    std::map<std::string, std::vector<Transcript> > transcriptHashSet;
    std::string cdsRegex = "ID=(.*?);Parent=(.*?)$";
    readGffFile (fileLocation, transcriptHashSet, cdsRegex);
    for( std::map<std::string, std::vector<Transcript> >::iterator it=transcriptHashSet.begin(); it!=transcriptHashSet.end(); ++it){
        std::cout << it->first << std::endl;
    }
}
