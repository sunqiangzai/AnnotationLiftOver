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
// parallel related things begin

struct myPara{
    Transcript* tartgetTranscript;
    Transcript* referenceTranscript;
    MyThreadCount* myThreadCount;
    std::map<std::string, Fasta>* targetGenome;
    std::map<std::string, Fasta>* referenceGenome;
    NucleotideCodeSubstitutionMatrix* nucleotideCodeSubstitutionMatrix;
    std::string chromosomeName;
    std::map<std::string, Transcript>* targetTranscriptsHashMap;
    myPara(Transcript& _tartgetTranscript, Transcript& _referenceTranscript,  MyThreadCount& _myThreadCount, std::map<std::string,
            Fasta>& _targetGenome, std::map<std::string, Fasta>& _referenceGenome,
           NucleotideCodeSubstitutionMatrix& _nucleotideCodeSubstitutionMatrix, std::string _chromosomeName,
           std::map<std::string, Transcript>& _targetTranscriptsHashMap){
        tartgetTranscript=&_tartgetTranscript;
        referenceTranscript=& _referenceTranscript;
        myThreadCount=&_myThreadCount;
        targetGenome=&_targetGenome;
        referenceGenome=&_referenceGenome;
        nucleotideCodeSubstitutionMatrix=&_nucleotideCodeSubstitutionMatrix;
        chromosomeName=_chromosomeName;
        targetTranscriptsHashMap=&_targetTranscriptsHashMap;
    }
};

struct myPara2{
    std::string cdsSequence;
    std::string targetSequence;
    NucleotideCodeSubstitutionMatrix* nucleotideCodeSubstitutionMatrix;
    int startTarget;
    int endTarget;
    STRAND strand;
    MyThreadCount* myThreadCount;
    std::map<std::string, Transcript>* homologousGoodTranscripts;
    std::string transcriptName;
    myPara2(std::string _cdsSequence, std::string _targetSequence, NucleotideCodeSubstitutionMatrix& _nucleotideCodeSubstitutionMatrix,
            int _startTarget, int _endTarget, STRAND _strand, MyThreadCount& _myThreadCount, std::map<std::string,
            Transcript>* _homologousGoodTranscripts, std::string _transcriptName ){
        cdsSequence=_cdsSequence;
        targetSequence=_targetSequence;
        nucleotideCodeSubstitutionMatrix=&_nucleotideCodeSubstitutionMatrix;
        startTarget=_startTarget;
        endTarget=_endTarget;
        strand=_strand;
        myThreadCount=&_myThreadCount;
        homologousGoodTranscripts=_homologousGoodTranscripts;
        transcriptName=_transcriptName;
    }
};


pthread_attr_t attr;
pthread_mutex_t mutex;

void * transcriptRealignment( void *arg );
void * runExonerateEst( void *arg );
void * runExonerateProtein( void *arg );

// parallel related things end


void TranscriptsTogenes(std::string& regexG, std::map<std::string, std::vector<Gene> >& genes,
                        std::map<std::string, Transcript>& targetTranscriptsHashMap){
    std::map< std::string, std::map<std::string, Gene > > tempGenes;
    for( std::map<std::string, Transcript>::iterator it=targetTranscriptsHashMap.begin();
         it!=targetTranscriptsHashMap.end(); ++it){
        std::string transcriptId = it->first;
        std::string chromosome = it->second.getChromeSomeName();
        if(tempGenes.find(chromosome) == tempGenes.end() ){
            tempGenes[chromosome]=std::map<std::string, Gene >();
        }
        std::regex reg(regexG);
        std::smatch match;
        regex_search(transcriptId, match, reg);
        if( match.empty() ){

        }else{
            std::string geneName = match[1];
            if( tempGenes[chromosome].find(geneName) == tempGenes[chromosome].end() ){
                Gene gene(geneName, it->second.getStrand());
                tempGenes[chromosome].insert( std::pair<std::string, Gene>(geneName, gene) );
            }
            tempGenes[chromosome][geneName].addTranscript(it->second);
            tempGenes[chromosome][geneName].setSource(it->second.getSource());
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
    std::cerr << "125" << std::endl;
}


void reAnnotationSingleLine( std::string& referenceGenomeFilePath, std::string referenceGffFilePath, std::string sdiFile,
                    std::map<std::string, std::vector<Gene> >& genes, int maxThread, std::string& regex, std::string& regexG){
    pthread_attr_init( &attr );
    pthread_attr_setscope( &attr, PTHREAD_SCOPE_SYSTEM );
    NucleotideCodeSubstitutionMatrix nucleotideCodeSubstitutionMatrix;
    std::map<std::string, Fasta> referenceGenome;
    std::map<std::string, std::vector<Transcript> > referenceTranscriptHashSet;
    readFastaFile(referenceGenomeFilePath, referenceGenome);
    readGffFile (referenceGffFilePath, referenceTranscriptHashSet, regex);
    std::map<std::string, std::vector<Variant> > variantsMaps;
    readSdiFile (sdiFile,  variantsMaps);
    std::map<std::string, Fasta> targetGenome;
    getPseudoGenomeSequence(referenceGenome, variantsMaps, targetGenome);
    std::map<std::string, Transcript> targetTranscriptsHashMap;// transcriptName Transcript
    annotationLiftOver( referenceTranscriptHashSet, targetTranscriptsHashMap, variantsMaps, targetGenome, referenceGenome);
    std::map<std::string, Transcript > referenceTranscriptHashMap; // transcriptName Transcript
    reAnnotationSingleLine( referenceTranscriptHashSet, referenceGenome, targetTranscriptsHashMap,
                            targetGenome, nucleotideCodeSubstitutionMatrix, maxThread);
    TranscriptsTogenes(regexG, genes, targetTranscriptsHashMap);

    pthread_mutex_destroy(&mutex);
}


void reAnnotationSingleLine( std::map<std::string, std::vector<Transcript> >& referenceTranscriptHashSet,
                             std::map<std::string, Fasta>& referenceGenome, std::map<std::string, Transcript>& targetTranscriptsHashMap,
                             std::map<std::string, Fasta>& targetGenome, NucleotideCodeSubstitutionMatrix& nucleotideCodeSubstitutionMatrix,
                             int maxThread){
    MyThreadCount myThreadCount;
    //std::cerr << "131" << std::endl;
    for( std::map<std::string, Fasta>::iterator it1=referenceGenome.begin(); it1!=referenceGenome.end(); ++it1 ){
        if( referenceTranscriptHashSet.find(it1->first)!= referenceTranscriptHashSet.end() ){
            std::string chromosomeName = it1->first;
            for( std::vector<Transcript>::iterator it2=referenceTranscriptHashSet[chromosomeName].begin(); it2!=referenceTranscriptHashSet[chromosomeName].end(); ++it2) {
                (*it2).updateInfor(referenceGenome);
                Transcript * it3 = & targetTranscriptsHashMap[it2->getName()];
                if( (*it3).getIfOrfShift() ){
                    bool isThisThreadUnrun = true;
                    while (isThisThreadUnrun) {
                        if (myThreadCount.getCount() < maxThread) {
                            pthread_t thread;
                            struct myPara *para = new struct myPara( (*it3), (*it2), myThreadCount, targetGenome, referenceGenome, nucleotideCodeSubstitutionMatrix, chromosomeName, targetTranscriptsHashMap);
                            pthread_mutex_lock(&mutex);
                            myThreadCount.plusOne();
                            pthread_mutex_unlock(&mutex);
                            pthread_create(&thread, &attr, transcriptRealignment, (void *) para);
                            isThisThreadUnrun = false;
                            break;
                        } else {
                            usleep(10);
                        }
                    }
                }
            }
        }
    }
    while(myThreadCount.hasNext()){// wait for all the thread
        usleep(50);
    }
}


void reAnnotationAndMsa( std::string& referenceGenomeFilePath, std::string referenceGffFilePath, std::map<std::string, std::string> sdiFilePaths, int maxThread, std::string& regex){

    pthread_attr_init( &attr );
    pthread_attr_setscope( &attr, PTHREAD_SCOPE_SYSTEM );

    std::map<std::string, Fasta> referenceGenome;
    std::map<std::string, std::vector<Transcript> > referenceTranscriptHashSet;

    readFastaFile(referenceGenomeFilePath, referenceGenome);
    readGffFile (referenceGffFilePath, referenceTranscriptHashSet, regex);

    NucleotideCodeSubstitutionMatrix nucleotideCodeSubstitutionMatrix;

    for( std::map<std::string, Fasta>::iterator it1=referenceGenome.begin(); it1!=referenceGenome.end(); it1++ ){
        if( referenceTranscriptHashSet.find(it1->first)!= referenceTranscriptHashSet.end() ){
            // check chromosome by chromosome, for RAM saving purpose. Do the keys of following maps are name of lines
            std::map<std::string, Fasta> accessionsTargetSequences;// the Fasta of current chromosome for each accession

            std::map<std::string, std::vector<Variant> > accessionsVariantsMaps; //accessions variants
            std::string chromosomeName = it1->first;
            std::map<std::string, std::map<std::string, Transcript> > accessionsTargetTranscriptHashSet; //accession, transcriptName transcript
            for( std::map<std::string, std::string>::iterator it2=sdiFilePaths.begin(); it2!=sdiFilePaths.end(); it2++ ){
                std::string accessionId = it2->first;
                std::map<std::string, std::vector<Variant> > variantsTMaps;
                readSdiFile (it2->second,  variantsTMaps, chromosomeName);// sdiFile, variantsTmaps chromosomeName

                if( variantsTMaps.find(chromosomeName)!=variantsTMaps.end() ){
                    accessionsVariantsMaps[accessionId] = variantsTMaps[chromosomeName];
                }else{
                    accessionsVariantsMaps[accessionId] = std::vector<Variant>();
                }

                std::map<std::string, Fasta> targetGenome;
                getPseudoGenomeSequence(referenceGenome, variantsTMaps, targetGenome, it1->first);
                if( targetGenome.find(chromosomeName)!=targetGenome.end() ){
                    accessionsTargetSequences[chromosomeName]=targetGenome.begin()->second;
                }else{
                    std::cerr << "reAnnotationAndMsa line 213, there are bugs" << std::endl;
                }
                //std::cout << "39 " << targetGenome.begin()->second.getSequence().length() << std::endl;
                std::map<std::string, Transcript> targetTranscriptsHashMap;// transcriptName Transcript
                annotationLiftOver( referenceTranscriptHashSet, targetTranscriptsHashMap, variantsTMaps, targetGenome, referenceGenome, chromosomeName);
                MyThreadCount myThreadCount;
                for( std::vector<Transcript>::iterator it4=referenceTranscriptHashSet[chromosomeName].begin(); it4!=referenceTranscriptHashSet[chromosomeName].end(); ++it4) {
                    Transcript * it3 = & targetTranscriptsHashMap[it4->getName()];

                    if( (*it3).getIfOrfShift()  ){

                        bool isThisThreadUnrun=true;
                        while(isThisThreadUnrun){
                            if(myThreadCount.getCount() < maxThread){
                                //std::cout << 62 << std::endl;
                                pthread_t thread;
                                struct myPara *para = new struct myPara( (*it3), (*it4), myThreadCount, targetGenome,
                                                                         referenceGenome, nucleotideCodeSubstitutionMatrix, 
                                                                         chromosomeName, targetTranscriptsHashMap);
                                pthread_mutex_lock(&mutex);
                                myThreadCount.plusOne();
                                pthread_mutex_unlock(&mutex);
                                pthread_create(&thread, &attr, transcriptRealignment, (void *)para);
                                isThisThreadUnrun=false;
                                break;
                            }else{
                                usleep(10);
                            }
                        }
                    }
                }
                while(myThreadCount.hasNext()){// wait for all the thread
                    usleep(50);
                }
                for( std::map<std::string, Transcript>::iterator it5 = targetTranscriptsHashMap.begin();
                        it5!=targetTranscriptsHashMap.end(); ++it5){
                    it5->second.cleanSaveRam();
                }
                accessionsTargetTranscriptHashSet[accessionId]=targetTranscriptsHashMap;
            }
        }
    }
    pthread_mutex_destroy(&mutex);
}

void * transcriptRealignment( void *arg ){
    struct myPara *pstru = (struct myPara *)arg;
    Transcript* it3 =  pstru->tartgetTranscript;
    Transcript* itr =  pstru->referenceTranscript;
    NucleotideCodeSubstitutionMatrix& nucleotideCodeSubstitutionMatrix = (*pstru->nucleotideCodeSubstitutionMatrix);
    std::map<std::string, Fasta>& targetGenome = (*pstru->targetGenome);
    std::map<std::string, Fasta>& referenceGenome = (*pstru->referenceGenome);
    std::string chromosomeName=pstru->chromosomeName;
    std::map<std::string, Transcript>& targetTranscriptsHashMap = (*pstru->targetTranscriptsHashMap);
    std::string refGenomeSequence;
    pthread_mutex_lock(&mutex);
    refGenomeSequence = (*itr).getGeneomeSequence();
    pthread_mutex_unlock(&mutex);
    //std::cout << "281" << std::endl;
    if( (*it3).getStrand() == POSITIVE ){
        size_t startCodonPosition=1;
        size_t stopCodonPosition=refGenomeSequence.length()-2;
        std::vector<SpliceSitePosition> splitSitePositions;
        if( (*it3).getCdsVector().size()>1 ){
            for( size_t i=1; i<(*it3).getCdsVector().size(); ++i ){
                SpliceSitePosition spliceSitePosition(
                        (*itr).getCdsVector()[i-1].getEnd()-(*itr).getStart()+2,
                        (*itr).getCdsVector()[i].getStart()-(*itr).getStart());
                splitSitePositions.push_back(spliceSitePosition);
            }
        }
        pthread_mutex_lock(&mutex);
        std::string dna_b = getSubsequence(targetGenome, chromosomeName, (*it3).getStart()-refGenomeSequence.length(), (*it3).getEnd()+refGenomeSequence.length(), POSITIVE);
        pthread_mutex_unlock(&mutex);
        NeedlemanWunschForTranscript nw (refGenomeSequence, dna_b, startCodonPosition, stopCodonPosition, splitSitePositions);
        pthread_mutex_lock(&mutex);
        nw.print_results();
        pthread_mutex_unlock(&mutex);

        Transcript targetTranscript((*it3).getName(), chromosomeName, (*it3).getStrand());
        targetTranscript.setSource("REALIGNMENT");
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
                for( std::vector<Cds>::iterator it4=(*itr).getCdsVector().begin();
                     it4!=(*itr).getCdsVector().end(); ++it4){
                    if( referencePosition+(*itr).getStart()-1 == (*it4).getStart() ){
                        if( nw.getAlignment_b()[tp] == '-' ){
                            targetCdsStarts.push_back((*it3).getStart()-refGenomeSequence.length() + targetPosition); //(*i3) is the target transcript, refGenomeSequence.length() is the extend length
                        }else{
                            targetCdsStarts.push_back((*it3).getStart()-refGenomeSequence.length() + targetPosition -1);
                        }
                    }
                    if( referencePosition+(*itr).getStart()-1 == (*it4).getEnd() ){
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
        checkOrfState( targetTranscript, (*itr),
                       targetGenome, referenceGenome, nucleotideCodeSubstitutionMatrix);
        pthread_mutex_lock(&mutex);
        targetTranscriptsHashMap[targetTranscript.getName()]=targetTranscript; // update the map data structure with modified transcript
//        std::cout << accessionId << " " << targetTranscript.getName() << " " << targetTranscript.getMetaInformation() << std::endl <<
//                  targetTranscript.getGeneomeSequence() << std::endl <<
//                  targetTranscript.getCdsSequence() << std::endl;
        pthread_mutex_unlock(&mutex);
        //targetTranscript.cleanSaveRam();// empty the transcriprt sequence and genome sequence, for RAM saving purpose
        
    }else{
        size_t startCodonPosition=1;
        size_t stopCodonPosition=refGenomeSequence.length()-2;
        std::vector<SpliceSitePosition> spliceSitePositions;
        if( (*it3).getCdsVector().size()>1 ){
            for( size_t i=(*it3).getCdsVector().size()-1; i>0; --i ){
                SpliceSitePosition spliceSitePosition(
                                           (*itr).getEnd()-(*itr).getCdsVector()[i].getStart()+2,
                                           (*itr).getEnd()-(*itr).getCdsVector()[i-1].getEnd()-1);
                spliceSitePositions.push_back(spliceSitePosition);
            }
        }
        pthread_mutex_lock(&mutex);
        std::string dna_b = getSubsequence(targetGenome, chromosomeName, (*it3).getStart()-refGenomeSequence.length(), (*it3).getEnd()+refGenomeSequence.length(), NEGATIVE);
        pthread_mutex_unlock(&mutex);
        NeedlemanWunschForTranscript nw (refGenomeSequence, dna_b, startCodonPosition, stopCodonPosition, spliceSitePositions);
        pthread_mutex_lock(&mutex);
        nw.print_results();
        pthread_mutex_unlock(&mutex);


        Transcript targetTranscript((*it3).getName(), chromosomeName, (*it3).getStrand());
        targetTranscript.setSource("REALIGNMENT");
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
                for( std::vector<Cds>::iterator it4=(*itr).getCdsVector().begin();
                     it4!=(*itr).getCdsVector().end(); it4++){
                    if( referencePosition+(*itr).getStart()-1 == (*it4).getStart() ){
                        if( nw.getAlignment_b()[tp] == '-' ){
                            targetCdsStarts.push_back((*it3).getEnd() + refGenomeSequence.length() - targetPosition);
                        }else{
                            targetCdsStarts.push_back((*it3).getEnd() + refGenomeSequence.length() - targetPosition +1);
                        }
                    }
                    if( referencePosition+(*itr).getStart()-1 == (*it4).getEnd() ){
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
        checkOrfState( targetTranscript, (*itr),
                       targetGenome, referenceGenome, nucleotideCodeSubstitutionMatrix);
        pthread_mutex_lock(&mutex);
        targetTranscriptsHashMap[targetTranscript.getName()]=targetTranscript; // update the map data structure with modified transcript
//        std::cout << accessionId << " " << targetTranscript.getName() << " " << targetTranscript.getMetaInformation() << std::endl <<
//                  targetTranscript.getGeneomeSequence() << std::endl <<
//                  targetTranscript.getCdsSequence() << std::endl;
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
//
//int main(){
//    std::string referenceGenomeFilePath = "/Users/song/Dropbox/geneStructureBasedRealignment/testData/Col.fa";
//    std::string referenceGffFilePath = "/Users/song/Dropbox/geneStructureBasedRealignment/Chr1_TAIR10_GFF3_genes.gff";
//    std::string sdiFile = "/Users/song/Dropbox/geneStructureBasedRealignment/testData/PA10000.sdi";
//    std::map<std::string, std::vector<Gene> > genes;
//    int maxThread=3;
//    std::string regex="(.*)Parent=(.*?)[;,].*$";
//    std::string regexG="(.+?)\\.";
//    reAnnotationAndExonerate(referenceGenomeFilePath, referenceGffFilePath, sdiFile, genes,
//                             maxThread, regex, regexG);
//    return 0;
//}

void reAnnotationAndExonerate( std::string& referenceGenomeFilePath, std::string referenceGffFilePath, std::string sdiFile,
                             std::map<std::string, std::vector<Gene> >& genes, int maxThread, std::string& regex, std::string& regexG){
    NucleotideCodeSubstitutionMatrix nucleotideCodeSubstitutionMatrix;
    std::map<std::string, Fasta> referenceGenome;
    std::map<std::string, std::vector<Transcript> > referenceTranscriptHashSet;
    readFastaFile(referenceGenomeFilePath, referenceGenome);
    readGffFile (referenceGffFilePath, referenceTranscriptHashSet, regex);

    std::map<std::string, std::vector<Variant> > variantsMaps;
    readSdiFile (sdiFile,  variantsMaps);
    std::map<std::string, Fasta> targetGenome;
    getPseudoGenomeSequence(referenceGenome, variantsMaps, targetGenome);
    std::map<std::string, Transcript> targetTranscriptsHashMap;// transcriptName Transcript
    annotationLiftOver( referenceTranscriptHashSet, targetTranscriptsHashMap, variantsMaps, targetGenome, referenceGenome);

    std::cerr << "annotationLiftOver done" << std::endl;

    std::map<std::string, Transcript > referenceTranscriptHashMap; // transcriptName Transcript

    reAnnotationSingleLine( referenceTranscriptHashSet, referenceGenome, targetTranscriptsHashMap,
                            targetGenome, nucleotideCodeSubstitutionMatrix, maxThread);

    std::cerr << "realignment done" << std::endl;
    std::map<std::string, Transcript> homologousGoodTranscripts;

    std::map<std::string, Transcript> *homologousGoodTranscriptsP = & homologousGoodTranscripts;

    pthread_attr_init( &attr );
    pthread_attr_setscope( &attr, PTHREAD_SCOPE_SYSTEM );
    MyThreadCount myThreadCount;

    for( std::map<std::string, std::vector<Transcript> >::iterator itchromosome=referenceTranscriptHashSet.begin();
            itchromosome!=referenceTranscriptHashSet.end(); ++itchromosome){
        for( std::vector<Transcript>::iterator itTranscript=referenceTranscriptHashSet[itchromosome->first].begin();
             itTranscript!=referenceTranscriptHashSet[itchromosome->first].end(); ++itTranscript){
            if( targetTranscriptsHashMap[(*itTranscript).getName()].getIfOrfShift() ){
                std::string refGenomeSequence=(*itTranscript).getGeneomeSequence();
                std::string cdsSequence=(*itTranscript).getCdsSequence();
                //std::cerr << "469 " << (*itTranscript).getName() << std::endl;
                int startTarget = targetTranscriptsHashMap[(*itTranscript).getName()].getStart()-refGenomeSequence.length();
                int endTarget = targetTranscriptsHashMap[(*itTranscript).getName()].getEnd()+refGenomeSequence.length();
                std::string targetSequence = getSubsequence(targetGenome, (*itTranscript).getChromeSomeName(),
                                                            startTarget, endTarget, (*itTranscript).getStrand());
                std::string transcriptName = (*itTranscript).getName();
                std::string tchromeSomeName = (*itTranscript).getChromeSomeName();
                STRAND tstrand = (*itTranscript).getStrand();
                Transcript targetTranscript(transcriptName, tchromeSomeName, tstrand );
                homologousGoodTranscripts[transcriptName]=(targetTranscript);
                //int lastIndex = homologousGoodTranscripts.size()-1;
                //Transcript& targetTranscript2 = homologousGoodTranscripts[lastIndex];
                bool isThisThreadUnrun = true;
                while (isThisThreadUnrun) {
                    if (myThreadCount.getCount() < maxThread) {
                        pthread_t thread;
                        struct myPara2 *para2 = new struct myPara2(cdsSequence, targetSequence, nucleotideCodeSubstitutionMatrix,
                        	startTarget, endTarget, (*itTranscript).getStrand(), myThreadCount, homologousGoodTranscriptsP, transcriptName );
                        pthread_mutex_lock(&mutex);
                        myThreadCount.plusOne();
                        pthread_mutex_unlock(&mutex);
                        pthread_create(&thread, &attr, runExonerateEst, (void *) para2);
                        isThisThreadUnrun = false;
                        break;
                    } else {
                        usleep(10);
                    }
                }
            }
        }
    }

    while(myThreadCount.hasNext()){// wait for all the thread
        usleep(50);
    }
    for( std::map<std::string, Transcript>::iterator it = homologousGoodTranscripts.begin();
         it!=homologousGoodTranscripts.end(); ++it){
        if( ! it->second.getIfOrfShift() ){
            targetTranscriptsHashMap[it->second.getName()] = it->second;
        }
    }
    std::cerr << "homologous CDS based annotation done" << std::endl;

    std::map<std::string, Transcript> homologousGoodTranscripts2;
    std::map<std::string, Transcript> *homologousGoodTranscriptsP2 = & homologousGoodTranscripts2;
    for( std::map<std::string, std::vector<Transcript> >::iterator itchromosome=referenceTranscriptHashSet.begin();
         itchromosome!=referenceTranscriptHashSet.end(); ++itchromosome){
        for( std::vector<Transcript>::iterator itTranscript=referenceTranscriptHashSet[itchromosome->first].begin();
             itTranscript!=referenceTranscriptHashSet[itchromosome->first].end(); ++itTranscript){
            if( targetTranscriptsHashMap[(*itTranscript).getName()].getIfOrfShift() ){
                std::string refGenomeSequence=(*itTranscript).getGeneomeSequence();
                std::string cdsSequence=(*itTranscript).getCdsSequence();

                int startTarget = targetTranscriptsHashMap[(*itTranscript).getName()].getStart()-refGenomeSequence.length();
                int endTarget = targetTranscriptsHashMap[(*itTranscript).getName()].getEnd()+refGenomeSequence.length();
                std::string targetSequence = getSubsequence(targetGenome, (*itTranscript).getChromeSomeName(),
                                                            startTarget, endTarget, (*itTranscript).getStrand());

                Transcript targetTranscript((*itTranscript).getName(), (*itTranscript).getChromeSomeName(),
                                            (*itTranscript).getStrand());
                homologousGoodTranscripts2[(*itTranscript).getName()]=targetTranscript;
                bool isThisThreadUnrun = true;
                while (isThisThreadUnrun) {
                    if (myThreadCount.getCount() < maxThread) {
                        pthread_t thread;
                        std::string protenSequene = nA2AA( cdsSequence, nucleotideCodeSubstitutionMatrix);
                        //int lastIndex = homologousGoodTranscripts2.size()-1;
                        struct myPara2 *para2 = new struct myPara2(protenSequene, targetSequence, nucleotideCodeSubstitutionMatrix,
                                                                   startTarget, endTarget, (*itTranscript).getStrand(), myThreadCount,
                                                                   homologousGoodTranscriptsP2, (*itTranscript).getName() );
                        pthread_mutex_lock(&mutex);
                        myThreadCount.plusOne();
                        pthread_mutex_unlock(&mutex);
                        pthread_create(&thread, &attr, runExonerateProtein, (void *) para2);
                        isThisThreadUnrun = false;
                        break;
                    } else {
                        usleep(10);
                    }
                }
            }else{
//                std::cerr << "552 " <<  std::endl;
            }
        }
    }
    while(myThreadCount.hasNext()){// wait for all the thread
        usleep(50);
    }
    for( std::map<std::string, Transcript>::iterator it = homologousGoodTranscripts2.begin();
         it!=homologousGoodTranscripts2.end(); ++it){
        if( ! it->second.getIfOrfShift() ){
            targetTranscriptsHashMap[it->second.getName()] = it->second;
        }
    }
    std::cerr << "homologous protein based annotation done" << std::endl;
    for( std::map<std::string, Transcript>::iterator it=targetTranscriptsHashMap.begin();
            it!=targetTranscriptsHashMap.end(); ++it){
        std::cerr << it->second.getName() <<  std::endl;
    }
    TranscriptsTogenes(regexG, genes, targetTranscriptsHashMap);
//    for( std::vector<Gene>::iterator it = genes.begin(); it!=genes.end(); ++it){
//        std::cerr << it->second.getName() <<  std::endl;
//    }
    for( std::map<std::string, Transcript>::iterator it = targetTranscriptsHashMap.begin(); it!=targetTranscriptsHashMap.end();++it ){
        it->second.updateInfor(targetGenome);
    }

    pthread_mutex_destroy(&mutex);
    std::cerr << "transcript structure to gene structure done" << std::endl;
}


void * runExonerateEst( void *arg ) {
    struct myPara2 *pstru = (struct myPara2 *)arg;

    std::string cdsSequence=(*pstru).cdsSequence;
    std::string targetSequence=(*pstru).targetSequence;
    NucleotideCodeSubstitutionMatrix* nucleotideCodeSubstitutionMatrix = pstru->nucleotideCodeSubstitutionMatrix;
    std::map<std::string, Transcript>* homologousGoodTranscripts = pstru ->homologousGoodTranscripts;
    Transcript t("Chr0333", "test", POSITIVE);
    (*homologousGoodTranscripts)["test"]=t;
    std::string transcriptName = pstru->transcriptName;

//    std::cerr << "600003 " << (*homologousGoodTranscripts)[transcriptName].getName() << std::endl;
    int startTarget=(*pstru).startTarget;
    int endTarget=(*pstru).endTarget;
    STRAND strand=(*pstru).strand;
//    std::cerr << "570" << std::endl;
//    std::cerr << "6000088888 " << (*homologousGoodTranscripts)[transcriptName].getName() << std::endl;
//    std::cout << cdsSequence << std::endl << targetSequence << std::endl;
    if( cdsSequence.length()>0 &&  targetSequence.length()>0){
//        std::cerr << "611111 " << (*homologousGoodTranscripts)[transcriptName].getName() << std::endl;
        pthread_mutex_lock(&mutex);
        std::string tempFolder = createdTempFloder();
        pthread_mutex_unlock(&mutex);
        std::string tempFile = tempFolder + "/" + generateUUID();
        std::string cdsFile = tempFile + "cds.fasta";
        std::string targetFile = tempFile + "target.fasta";

        std::ofstream ofile;
        ofile.open(cdsFile);
        ofile << ">cds" << std::endl << cdsSequence << std::endl;
        ofile.close();

        ofile.open(targetFile);
        ofile << ">target" << std::endl << targetSequence << std::endl;
        ofile.close();
        //std::cerr << "622228 " << targetTranscript->getName() << std::endl;


        std::string command =
                "exonerate --maxintron 30000 --model est2genome -i -10 --score 10 --bestn 1 --minintron 10 " + cdsFile +
                " " + targetFile + " --showtargetgff true >" + tempFile;
//        std::cerr << command << std::endl;

        system(&command[0]);
        readExonerateEstResult(tempFile, (*nucleotideCodeSubstitutionMatrix), startTarget, endTarget, strand, homologousGoodTranscripts, transcriptName);
        std::string cleanFileCommand = "rm " + cdsFile + "; rm " + targetFile + "; rm " + tempFile;
//        std::cerr << cleanFileCommand << std::endl;
        //system(&cleanFileCommand[0]);
    }

    pthread_mutex_lock(&mutex);
//    std::cerr << "line 532 " << std::endl;
//    std::cerr << "6444443838 " << (*homologousGoodTranscripts)[transcriptName].getName() << std::endl;
    pstru->myThreadCount->countDown();
//    std::cerr << "line 599 " << pstru->myThreadCount->getCount() << std::endl;
//
//    std::cerr << "line 636 " << pstru->myThreadCount->getCount() << std::endl;
    pthread_mutex_unlock(&mutex);
    return NULL;
}

void readExonerateEstResult( std::string& fileLocation, NucleotideCodeSubstitutionMatrix& nucleotideCodeSubstitutionMatrix,
                             int startTarget, int endTarget, STRAND strand, std::map<std::string, Transcript>* homologousGoodTranscripts, std::string& transcriptName ){

    if( (!if_file_exists(fileLocation)) || GetFileSize(fileLocation)<1 ){
        return;
    }
    std::ifstream infile(fileLocation);
//    std::cerr << "648 " << (*homologousGoodTranscripts)[transcriptName].getName() << std::endl;
    std::regex reg("\\d\\s*:\\s*([\\w\\.\\s\\>-]+)\\s*:\\s*\\d");
//    std::cerr << "651 " << (*homologousGoodTranscripts)[transcriptName].getName() << std::endl;
    std::stringstream sequencestream;
    std::string line;
    int lineNumber = 0;
//    std::cerr << "652 " << (*homologousGoodTranscripts)[transcriptName].getName() << std::endl;
//    return;
    while ( std::getline(infile, line) ){
        std::smatch match;
        regex_search(line, match, reg);
        if( match.empty() ){

        }else{
            lineNumber++;
            if( lineNumber % 2 == 0 ){
                sequencestream << match[1];
            }
        }
    }
    infile.close();
//    std::cerr << "765 " << (*homologousGoodTranscripts)[transcriptName].getName() << std::endl;

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
    //std::cerr << "703 " << (*homologousGoodTranscripts)[transcriptName].getName() << std::endl;
    //return;
    if( !orfShift && POSITIVE ==strand ){
        metaInformation << "_ConservedFunction";
        std::string mi = metaInformation.str();

        pthread_mutex_lock(&mutex);
        (*homologousGoodTranscripts)[transcriptName].setIfOrfShift(false);
        (*homologousGoodTranscripts)[transcriptName].setSource("CDSALIGNMENT");
        (*homologousGoodTranscripts)[transcriptName].setMetaInformation(mi);
        pthread_mutex_unlock(&mutex);

        std::ifstream infile2(fileLocation);
        std::regex reg2("^(\\S*)\t(.*)\texon\t(\\S*)\t(\\S*)\t(\\S*)\t\\+");
        std::string line2="";
        while (std::getline(infile2, line2)){
            std::smatch match;
            regex_search(line2, match, reg2);

            if( match.empty() ){
            }else{
                int start = stoi(match[3]);
                int end = stoi(match[4]);
                if(start>end){
                    int temp=start;
                    start = end;
                    end = temp;
                }
                start = start + startTarget -1;
                end = end + startTarget - 1;
                Cds cds(start, end);
                pthread_mutex_lock(&mutex);
                (*homologousGoodTranscripts)[transcriptName].addCds(cds);
                pthread_mutex_unlock(&mutex);
                std::cerr << start << " " << end << std::endl;
            }
        }
        infile2.close();
    }else if( !orfShift ){
        metaInformation << "_ConservedFunction";
        std::string mi = metaInformation.str();

        pthread_mutex_lock(&mutex);
        (*homologousGoodTranscripts)[transcriptName].setIfOrfShift(false);
        (*homologousGoodTranscripts)[transcriptName].setSource("CDSALIGNMENT");
        (*homologousGoodTranscripts)[transcriptName].setMetaInformation(mi);
        pthread_mutex_unlock(&mutex);

        std::ifstream infile2(fileLocation);
        std::regex reg2("^(\\S*)\t(.*)\texon\t(\\S*)\t(\\S*)\t(\\S*)\t\\-");
        std::string line2="";
        while (std::getline(infile2, line2)){
            std::smatch match;
            regex_search(line2, match, reg2);

            if( match.empty() ){
            }else{
                int start = stoi(match[3]);
                int end = stoi(match[4]);
                if(start>end){
                    int temp=start;
                    start = end;
                    end = temp;
                }
                start = start + startTarget -1;
                end = end + startTarget - 1;
                Cds cds(start, end);

                pthread_mutex_lock(&mutex);
                (*homologousGoodTranscripts)[transcriptName].addCds(cds);
                pthread_mutex_unlock(&mutex);
                std::cout << start << " " << end << std::endl;
            }
        }
        infile2.close();

    }else{
        std::string mi = metaInformation.str();
        pthread_mutex_lock(&mutex);
        (*homologousGoodTranscripts)[transcriptName].setIfOrfShift(true);
        (*homologousGoodTranscripts)[transcriptName].setMetaInformation(mi);
        pthread_mutex_unlock(&mutex);
    }
//    pthread_mutex_lock(&mutex);
//    std::cerr << metaInformation.str() << std::endl;
////    std::cerr << cdsSequenceString << " 729" << std::endl;
//    pthread_mutex_unlock(&mutex);
    return;
}




void* runExonerateProtein(void *arg ){
    struct myPara2 *pstru = (struct myPara2 *)arg;
    std::string protein=(*pstru).cdsSequence;
    std::string targetSequence=(*pstru).targetSequence;
    NucleotideCodeSubstitutionMatrix* nucleotideCodeSubstitutionMatrix = pstru->nucleotideCodeSubstitutionMatrix;
    std::map<std::string, Transcript>* homologousGoodTranscripts = pstru ->homologousGoodTranscripts;
    std::string transcriptName = pstru->transcriptName;
    int startTarget=(*pstru).startTarget;
    int endTarget=(*pstru).endTarget;
    STRAND strand=(*pstru).strand;

    if(protein.length()>0){
        std::string tempFolder = createdTempFloder();
        std::string tempFile = tempFolder + "/" + generateUUID();
        std::string proteinFile = tempFile + "protein.fasta";
        std::string targetFile = tempFile + "target.fasta";

        std::ofstream ofile;
        ofile.open(proteinFile);
        ofile << ">cds" << std::endl << protein << std::endl;
        ofile.close();

        ofile.open(targetFile);
        ofile << ">target" << std::endl << targetSequence << std::endl;
        ofile.close();
        std::string command = "exonerate --bestn 1 --maxintron 30000 --intronpenalty -10 --model protein2genome --percent 10 --score 10 --minintron 10 "+proteinFile + " " + targetFile + " --showtargetgff true >" +tempFile;
        std::cout << command << std::endl;
        system(&command[0]);
        readExonerateProteinResult(tempFile, (*nucleotideCodeSubstitutionMatrix), startTarget, endTarget, strand, homologousGoodTranscripts, transcriptName) ;
        std::string cleanFileCommand = "rm " + proteinFile + "; rm " + targetFile + "; rm " + tempFile;
        //system(&cleanFileCommand[0]);
    }
    pthread_mutex_lock(&mutex);
    pstru->myThreadCount->countDown();
    pthread_mutex_unlock(&mutex);
    return NULL;
}



void readExonerateProteinResult( std::string& fileLocation, NucleotideCodeSubstitutionMatrix& nucleotideCodeSubstitutionMatrix,
                             int startTarget, int endTarget, STRAND strand, std::map<std::string, Transcript>* homologousGoodTranscripts,
                                 std::string& transcriptName){
    if( (!if_file_exists(fileLocation)) || GetFileSize(fileLocation)<1 ){
        return;
    }
    std::ifstream infile(fileLocation);
    std::regex reg("^\\s*\\d+\\s*:\\s*([\\w\\.\\s\\>\\-\\{\\}\\*]+)\\s*:\\s*\\d+\\s*$");
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
    std::regex vowel_re("\\{");
    cdsSequenceString=std::regex_replace(cdsSequenceString, vowel_re, "");
    std::regex vowel_re2("\\}");
    cdsSequenceString=std::regex_replace(cdsSequenceString, vowel_re2, "");
    std::regex vowel_re3("\\s");
    cdsSequenceString=std::regex_replace(cdsSequenceString, vowel_re3, "");
    std::regex vowel_re4("[a-z]+\\.+[a-z]+");
    cdsSequenceString=std::regex_replace(cdsSequenceString, vowel_re4, "");
    std::regex vowel_re5("\\*");
    cdsSequenceString=std::regex_replace(cdsSequenceString, vowel_re5, "");

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
    if( !orfShift && POSITIVE ==strand ){

        metaInformation << "_ConservedFunction";
        std::string mi = metaInformation.str();
        pthread_mutex_lock(&mutex);
        (*homologousGoodTranscripts)[transcriptName].setSource("PROTEINALIGNMENT");
        (*homologousGoodTranscripts)[transcriptName].setIfOrfShift(false);
        (*homologousGoodTranscripts)[transcriptName].setMetaInformation(mi);
        pthread_mutex_unlock(&mutex);

        std::ifstream infile(fileLocation);
        std::regex reg("^(\\S*)\t(.*)\tcds\t(\\S*)\t(\\S*)\t(\\S*)\t\\+");
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
                start = start + startTarget -1;
                end = end + startTarget - 1;
                Cds cds(start, end);
                pthread_mutex_lock(&mutex);
                (*homologousGoodTranscripts)[transcriptName].addCds(cds);
                pthread_mutex_unlock(&mutex);
                std::cout << start << " " << end << std::endl;
            }
        }
        infile.close();
    }else if( !orfShift ){
        metaInformation << "_ConservedFunction";
        std::string mi = metaInformation.str();

        pthread_mutex_lock(&mutex);
        (*homologousGoodTranscripts)[transcriptName].setSource("PROTEINALIGNMENT");
        (*homologousGoodTranscripts)[transcriptName].setIfOrfShift(false);
        (*homologousGoodTranscripts)[transcriptName].setMetaInformation(mi);
        pthread_mutex_unlock(&mutex);

        std::ifstream infile(fileLocation);
        std::regex reg("^(\\S*)\t(.*)\tcds\t(\\S*)\t(\\S*)\t(\\S*)\t\\-");
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
                start = start + startTarget -1;
                end = end + startTarget - 1;
                Cds cds(start, end);
                pthread_mutex_lock(&mutex);
                (*homologousGoodTranscripts)[transcriptName].addCds(cds);
                pthread_mutex_unlock(&mutex);
                std::cout << start << " " << end << std::endl;
            }
        }
        infile.close();

    }else{
        std::string mi = metaInformation.str();
        pthread_mutex_lock(&mutex);
        (*homologousGoodTranscripts)[transcriptName].setIfOrfShift(true);
        (*homologousGoodTranscripts)[transcriptName].setMetaInformation(mi);
        pthread_mutex_unlock(&mutex);
    }
//    pthread_mutex_lock(&mutex);
//    std::cout << metaInformation.str() << std::endl;
//    std::cout << cdsSequenceString << std::endl;
//    pthread_mutex_unlock(&mutex);
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
