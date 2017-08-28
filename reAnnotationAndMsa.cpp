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
#include <thread>
#include <unistd.h>
#include <fstream>
#include <sstream>
#include "parameters.h"
#include <regex>
#include <mutex>
#include <chrono>
std::mutex gmutex;


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
                STRAND strand = it->second.getStrand();
                Gene gene(geneName, strand);
                tempGenes[chromosome].insert( std::pair<std::string, Gene>(geneName, gene) );
            }
            tempGenes[chromosome][geneName].addTranscript(it->second);
            std::string source = it->second.getSource();
            tempGenes[chromosome][geneName].setSource(source);
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
//    std::cerr << "125" << std::endl;
}


void reAnnotationSingleLine( std::string& referenceGenomeFilePath, std::string referenceGffFilePath, std::string sdiFile,
                    std::map<std::string, std::vector<Gene> >& genes, int maxThread, std::string& regex, std::string& regexG, int & lengthThread){

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
//    std::map<std::string, Transcript > referenceTranscriptHashMap; // transcriptName Transcript
    reAnnotationSingleLine( referenceTranscriptHashSet, referenceGenome, targetTranscriptsHashMap,
                            targetGenome, nucleotideCodeSubstitutionMatrix, maxThread, lengthThread);
    TranscriptsTogenes(regexG, genes, targetTranscriptsHashMap);

}


void reAnnotationSingleLine( std::map<std::string, std::vector<Transcript> >& referenceTranscriptHashSet,
                             std::map<std::string, Fasta>& referenceGenome, std::map<std::string, Transcript>& targetTranscriptsHashMap,
                             std::map<std::string, Fasta>& targetGenome, NucleotideCodeSubstitutionMatrix& nucleotideCodeSubstitutionMatrix,
                             int& maxThread, int & lengthThread){
    std::atomic_int number_of_runing_threads(0);
    for( std::map<std::string, Fasta>::iterator it1=referenceGenome.begin(); it1!=referenceGenome.end(); ++it1 ){
        if( referenceTranscriptHashSet.find(it1->first)!= referenceTranscriptHashSet.end() ){
            std::string chromosomeName = it1->first;
            for( std::vector<Transcript>::iterator it2=referenceTranscriptHashSet[chromosomeName].begin(); it2!=referenceTranscriptHashSet[chromosomeName].end(); ++it2) {
                (*it2).updateInfor(referenceGenome);
                Transcript * it3 = & targetTranscriptsHashMap[it2->getName()];
                if( (*it3).getIfOrfShift() ){
                    bool isThisThreadUnrun = true;
                    while (isThisThreadUnrun) {
                        if (number_of_runing_threads < maxThread) {
                            std::thread t(transcriptRealignment, std::ref(*it3), std::ref(*it2), std::ref(nucleotideCodeSubstitutionMatrix),  std::ref(targetGenome),
                                          std::ref(referenceGenome),
                                          chromosomeName, std::ref(targetTranscriptsHashMap), std::ref(number_of_runing_threads), std::ref(lengthThread));
                            ++number_of_runing_threads;
                            t.detach();
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
    while( number_of_runing_threads >0 ){// wait for all the thread
        usleep(50);
    }
}


void reAnnotationAndMsa( std::string& referenceGenomeFilePath, std::string referenceGffFilePath, std::map<std::string, std::string> sdiFilePaths, int& maxThread, std::string& regex, int & lengthThread ){

    std::atomic_int number_of_runing_threads(0);

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
                std::string chromosomeName = it1->first;
                getPseudoGenomeSequence(referenceGenome, variantsTMaps, targetGenome, chromosomeName);
                if( targetGenome.find(chromosomeName)!=targetGenome.end() ){
                    accessionsTargetSequences[chromosomeName]=targetGenome.begin()->second;
                }else{
                    std::cerr << "reAnnotationAndMsa line 213, there are bugs" << std::endl;
                }
                //std::cout << "39 " << targetGenome.begin()->second.getSequence().length() << std::endl;
                std::map<std::string, Transcript> targetTranscriptsHashMap;// transcriptName Transcript
                annotationLiftOver( referenceTranscriptHashSet, targetTranscriptsHashMap, variantsTMaps, targetGenome, referenceGenome, chromosomeName);
                for( std::vector<Transcript>::iterator it4=referenceTranscriptHashSet[chromosomeName].begin(); it4!=referenceTranscriptHashSet[chromosomeName].end(); ++it4) {
                    Transcript * it3 = & targetTranscriptsHashMap[it4->getName()];

                    if( (*it3).getIfOrfShift()  ){

                        bool isThisThreadUnrun=true;
                        while(isThisThreadUnrun){
                            if(number_of_runing_threads < maxThread){
                                std::thread t(transcriptRealignment, std::ref(*it3), std::ref(*it4), std::ref(nucleotideCodeSubstitutionMatrix),  std::ref(targetGenome),
                                              std::ref(referenceGenome),
                                              chromosomeName, std::ref(targetTranscriptsHashMap), std::ref(number_of_runing_threads), std::ref(lengthThread));
                                t.detach();
                                ++number_of_runing_threads;
                                isThisThreadUnrun=false;
                                break;
                            }else{
                                usleep(5);
                            }
                        }
                    }
                }
                while(number_of_runing_threads > 0){// wait for all the thread
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
}

void transcriptRealignment( Transcript& tartgetTranscript, Transcript& referenceTranscript,
                            NucleotideCodeSubstitutionMatrix& nucleotideCodeSubstitutionMatrix,
        std::map<std::string, Fasta>& targetGenome,
        std::map<std::string, Fasta>& referenceGenome, std::string chromosomeName, std::map<std::string,
        Transcript>& targetTranscriptsHashMap, std::atomic_int & number_of_runing_threads, int& lengthThread  ){

    std::string refGenomeSequence = referenceTranscript.getGeneomeSequence();
    std::string dna_b = getSubsequence(targetGenome, chromosomeName, tartgetTranscript.getStart()-refGenomeSequence.length(), tartgetTranscript.getEnd()+refGenomeSequence.length(), tartgetTranscript.getStrand());
    if( dna_b.length() <= lengthThread &&  refGenomeSequence.length()<=lengthThread ){

        if( tartgetTranscript.getStrand() == POSITIVE ){
            int startCodonPosition=1;
            int stopCodonPosition=refGenomeSequence.length()-2;
            std::vector<SpliceSitePosition> splitSitePositions;
            if( referenceTranscript.getCdsVector().size()>1 ){
                for( size_t i=1; i<referenceTranscript.getCdsVector().size(); ++i ){
                    SpliceSitePosition spliceSitePosition(
                            referenceTranscript.getCdsVector()[i-1].getEnd()-referenceTranscript.getStart()+2,
                            referenceTranscript.getCdsVector()[i].getStart()-referenceTranscript.getStart());
                    splitSitePositions.push_back(spliceSitePosition);
                }
            }

            NeedlemanWunschForTranscript nw (refGenomeSequence, dna_b, startCodonPosition, stopCodonPosition, splitSitePositions);

            Transcript targetTranscript(tartgetTranscript.getName(), chromosomeName, tartgetTranscript.getStrand());
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
                    for( std::vector<Cds>::iterator it4=referenceTranscript.getCdsVector().begin();
                         it4!=referenceTranscript.getCdsVector().end(); ++it4){
                        if( referencePosition+referenceTranscript.getStart()-1 == (*it4).getStart() ){
                            if( nw.getAlignment_b()[tp] == '-' ){
                                targetCdsStarts.push_back(tartgetTranscript.getStart()-refGenomeSequence.length() + targetPosition-1); //(*i3) is the target transcript, refGenomeSequence.length() is the extend length
                            }else{
                                targetCdsStarts.push_back(tartgetTranscript.getStart()-refGenomeSequence.length() + targetPosition -1);
                            }
                        }
                        if( referencePosition+referenceTranscript.getStart()-1 == (*it4).getEnd() ){
                            if( nw.getAlignment_b()[tp] == '-' ){
                                targetCdsEnds.push_back(tartgetTranscript.getStart()-refGenomeSequence.length() + targetPosition-1);
                            }else{
                                targetCdsEnds.push_back(tartgetTranscript.getStart()-refGenomeSequence.length() + targetPosition -1);
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
//            targetTranscript.updateInfor(targetGenome);
//            checkOrfState( targetTranscript, referenceTranscript,
//                           targetGenome, referenceGenome, nucleotideCodeSubstitutionMatrix);
            checkOrfState( targetTranscript, targetGenome, nucleotideCodeSubstitutionMatrix);

            if( !targetTranscript.getIfOrfShift() ){
                gmutex.lock();
                targetTranscriptsHashMap[targetTranscript.getName()]=targetTranscript; // update the map data structure with modified transcript
                gmutex.unlock();
           }//else{
    //            gmutex.lock();
    //            targetTranscriptsHashMap[targetTranscript.getName()].setSource("LIFTOVER");
    //            gmutex.unlock();
    //        }
        }else{
            int startCodonPosition=1;
            int stopCodonPosition=refGenomeSequence.length()-2;
            std::vector<SpliceSitePosition> spliceSitePositions;
            if( referenceTranscript.getCdsVector().size()>1 ){
                for( size_t i=referenceTranscript.getCdsVector().size()-1; i>0; --i ){
                    SpliceSitePosition spliceSitePosition(
                                                referenceTranscript.getEnd()-referenceTranscript.getCdsVector()[i].getStart()+2,
                                                referenceTranscript.getEnd()-referenceTranscript.getCdsVector()[i-1].getEnd()-1);
                    spliceSitePositions.push_back(spliceSitePosition);
                }
            }
            NeedlemanWunschForTranscript nw (refGenomeSequence, dna_b, startCodonPosition, stopCodonPosition, spliceSitePositions);
            gmutex.lock();
            nw.print_results();
            gmutex.unlock();
            Transcript targetTranscript(tartgetTranscript.getName(), chromosomeName, tartgetTranscript.getStrand());
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
                    for( std::vector<Cds>::iterator it4=referenceTranscript.getCdsVector().begin();
                         it4!=referenceTranscript.getCdsVector().end(); it4++){
                        if( referencePosition+referenceTranscript.getStart()-1 == (*it4).getStart() ){
                            if( nw.getAlignment_b()[tp] == '-' ){
                                targetCdsStarts.push_back(tartgetTranscript.getEnd() + refGenomeSequence.length() - targetPosition+1);
                            }else{
                                targetCdsStarts.push_back(tartgetTranscript.getEnd() + refGenomeSequence.length() - targetPosition +1);
                            }
                        }
                        if( referencePosition+referenceTranscript.getStart()-1 == (*it4).getEnd() ){
                            if( nw.getAlignment_b()[tp] == '-' ){
                                targetCdsEnds.push_back(tartgetTranscript.getEnd() + refGenomeSequence.length() - targetPosition +1);
                            }else{
                                targetCdsEnds.push_back(tartgetTranscript.getEnd() + refGenomeSequence.length() - targetPosition +1);
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
            checkOrfState( targetTranscript, targetGenome, nucleotideCodeSubstitutionMatrix);
//            checkOrfState( targetTranscript, referenceTranscript,
//                           targetGenome, referenceGenome, nucleotideCodeSubstitutionMatrix);

            if( !targetTranscript.getIfOrfShift() ){
                gmutex.lock();
                targetTranscriptsHashMap[targetTranscript.getName()]=targetTranscript; // update the map data structure with modified transcript
                gmutex.unlock();
            }
            //targetTranscript.cleanSaveRam();// empty the transcriprt sequence and genome sequence, for RAM saving purpose
        }
    }
    --number_of_runing_threads;
    return;
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
                             std::map<std::string, std::vector<Gene> >& genes, int maxThread, std::string& regex, std::string& regexG, std::string & outputGffFile, int & lengthThread){

    NucleotideCodeSubstitutionMatrix nucleotideCodeSubstitutionMatrix;
    std::map<std::string, Fasta> referenceGenome;
    std::map<std::string, std::vector<Transcript> > referenceTranscriptHashSet;
    readFastaFile(referenceGenomeFilePath, referenceGenome);
    readGffFile (referenceGffFilePath, referenceTranscriptHashSet, regex);

    std::set<std::string> toRemoveChromosomes;
    for( std::map<std::string, Fasta>::iterator it=referenceGenome.begin();
            it!=referenceGenome.end(); ++it){
        if( referenceTranscriptHashSet.find(it->first)==referenceTranscriptHashSet.end() ){
            toRemoveChromosomes.insert(it->first);
        }
    }
    for( std::map<std::string, std::vector<Transcript> >::iterator it=referenceTranscriptHashSet.begin();
            it!=referenceTranscriptHashSet.end(); ++it){
        if( referenceGenome.find(it->first)==referenceGenome.end() ){
            toRemoveChromosomes.insert(it->first);
        }
    }
    for( std::set<std::string>::iterator it=toRemoveChromosomes.begin();
            it!=toRemoveChromosomes.end(); ++it){
        if( referenceGenome.find(*it) != referenceGenome.end() ){
            referenceGenome.erase(*it);
        }
        if( referenceTranscriptHashSet.find(*it) != referenceTranscriptHashSet.end() ){
            referenceTranscriptHashSet.erase(*it);
        }
    }

    std::map<std::string, std::vector<Variant> > variantsMaps;
    readSdiFile (sdiFile, variantsMaps);
    std::map<std::string, Fasta> targetGenome;
    getPseudoGenomeSequence(referenceGenome, variantsMaps, targetGenome);
    std::map<std::string, Transcript> targetTranscriptsHashMap;// transcriptName Transcript
    annotationLiftOver( referenceTranscriptHashSet, targetTranscriptsHashMap, variantsMaps, targetGenome, referenceGenome);
    //std::cerr << "annotationLiftOver done. maxThread:" << maxThread << std::endl;

    std::vector<std::thread> threads;

    std::string tempFolder = createdTempFloder(); // create tempFolder
    std::atomic_int number_of_runing_threads(0);
    for(std::map<std::string, std::vector<Transcript> >::iterator it=referenceTranscriptHashSet.begin();
        it!=referenceTranscriptHashSet.end(); ++it){
        for( std::vector<Transcript>::iterator it2=it->second.begin(); it2!=it->second.end();++it2 ){
            (*it2).updateInfor(referenceGenome);
            Transcript * it3 = & targetTranscriptsHashMap[it2->getName()];
            if( (*it3).getIfOrfShift() ){
                bool isThisThreadUnrun = true;
                while (isThisThreadUnrun) {
                    if (number_of_runing_threads < maxThread) {
                        std::thread t(transcriptRealignmentAndExonerate, std::ref(*it3), std::ref(*it2), std::ref(nucleotideCodeSubstitutionMatrix),  std::ref(targetGenome),
                                      std::ref(referenceGenome),
                                      it->first, std::ref(targetTranscriptsHashMap), std::ref(number_of_runing_threads), std::ref(outputGffFile), std::ref(lengthThread) );
                        ++number_of_runing_threads;
                        t.detach();
                        isThisThreadUnrun = false;
//                        break;
                    } else {
                        usleep(20000);
                    }
                }
            }
        }

    }
    while( number_of_runing_threads >0 ){// wait for all the thread
        usleep(100000);
    }
    //std::cerr << "realignment done" << std::endl;
    TranscriptsTogenes(regexG, genes, targetTranscriptsHashMap);
    //std::cerr << "transcript structure to gene structure done" << std::endl;
}

void reAnnotationAndExonerateAndNovo( std::string& referenceGenomeFilePath, std::string referenceGffFilePath, std::string novoGffFilePath, std::string sdiFile,
                               std::map<std::string, std::vector<Gene> >& genes, int maxThread, std::string& regex, std::string& regexG, std::string novoRegex, std::string& novoRegexG, std::string & outputGffFile, int & lengthThread){

    NucleotideCodeSubstitutionMatrix nucleotideCodeSubstitutionMatrix;
    std::map<std::string, Fasta> referenceGenome;
    std::map<std::string, std::vector<Transcript> > referenceTranscriptHashSet;
    readFastaFile(referenceGenomeFilePath, referenceGenome);
    readGffFile (referenceGffFilePath, referenceTranscriptHashSet, regex);

    std::set<std::string> toRemoveChromosomes;
    for( std::map<std::string, Fasta>::iterator it=referenceGenome.begin();
         it!=referenceGenome.end(); ++it){
        if( referenceTranscriptHashSet.find(it->first)==referenceTranscriptHashSet.end() ){
            toRemoveChromosomes.insert(it->first);
        }
    }
    for( std::map<std::string, std::vector<Transcript> >::iterator it=referenceTranscriptHashSet.begin();
         it!=referenceTranscriptHashSet.end(); ++it){
        if( referenceGenome.find(it->first)==referenceGenome.end() ){
            toRemoveChromosomes.insert(it->first);
        }
    }
    for( std::set<std::string>::iterator it=toRemoveChromosomes.begin();
         it!=toRemoveChromosomes.end(); ++it){
        if( referenceGenome.find(*it) != referenceGenome.end() ){
            referenceGenome.erase(*it);
        }
        if( referenceTranscriptHashSet.find(*it) != referenceTranscriptHashSet.end() ){
            referenceTranscriptHashSet.erase(*it);
        }
    }

    std::map<std::string, std::vector<Variant> > variantsMaps;
    readSdiFile (sdiFile, variantsMaps);
    std::map<std::string, Fasta> targetGenome;
    getPseudoGenomeSequence(referenceGenome, variantsMaps, targetGenome);
    std::map<std::string, Transcript> targetTranscriptsHashMap;// transcriptName Transcript
    annotationLiftOver( referenceTranscriptHashSet, targetTranscriptsHashMap, variantsMaps, targetGenome, referenceGenome);
    std::cerr << "annotationLiftOver done. maxThread:" << maxThread << std::endl;

    std::vector<std::thread> threads;

    std::string prefixUuid = outputGffFile;
    std::regex reg("[\\s\\S]*[(\\)(\\/)]([\\s\\S]*)");
    std::smatch match;
    regex_search(outputGffFile, match, reg);
    if( match.empty()){

    }else{
        prefixUuid=match[3];
    }
    std::string tempFolder = createdTempFloder(); // create tempFolder
    std::atomic_int number_of_runing_threads(0);
    for(std::map<std::string, std::vector<Transcript> >::iterator it=referenceTranscriptHashSet.begin();
        it!=referenceTranscriptHashSet.end(); ++it){
        for( std::vector<Transcript>::iterator it2=it->second.begin(); it2!=it->second.end();++it2 ){
            Transcript * it3 = & targetTranscriptsHashMap[it2->getName()];
            if( (*it3).getIfOrfShift() ){
                (*it2).updateInfor(referenceGenome);
                bool isThisThreadUnrun = true;
                while (isThisThreadUnrun) {
                    if (number_of_runing_threads < maxThread) {
                        std::thread t(transcriptRealignmentAndExonerate, std::ref(*it3), std::ref(*it2), std::ref(nucleotideCodeSubstitutionMatrix),  std::ref(targetGenome),
                                      std::ref(referenceGenome),
                                      it->first, std::ref(targetTranscriptsHashMap), std::ref(number_of_runing_threads), std::ref(prefixUuid), std::ref(lengthThread) );
                        ++number_of_runing_threads;
                        t.detach();
                        isThisThreadUnrun = false;
                        break;
                    } else {
                        usleep(20000);
                    }
                }
            }
        }
    }
    while( number_of_runing_threads >0 ){// wait for all the thread
        usleep(100000);
    }
    std::cerr << "realignment done" << std::endl;

    //replace old orf-shift annotation with novo annotation begin
    std::map<std::string, std::vector<Transcript> > novoTranscriptHashSet;
    readGffFile (novoGffFilePath, novoTranscriptHashSet, novoRegex);

    std::map<std::string, Transcript> adaptedNovoTranscriptAnnotation1;
    std::map<std::string, Transcript> adaptedNovoTranscriptAnnotation2;
    for( std::map<std::string, std::vector<Transcript> >::iterator it=novoTranscriptHashSet.begin(); it!=novoTranscriptHashSet.end(); ++it){
        if( referenceGenome.find(it->first) != referenceGenome.end() ){
            for(std::vector<Transcript>::iterator it2=novoTranscriptHashSet[it->first].begin();
                it2!=novoTranscriptHashSet[it->first].end(); ++it2){
                bool ifHasOverLap=false;
                for( std::map<std::string, Transcript>::iterator it3=targetTranscriptsHashMap.begin();
                        it3!=targetTranscriptsHashMap.end(); ++it3){
                    if( it2->ifOverLap(it3->second) ){
                        ifHasOverLap=true;
                        if( it3->second.getIfOrfShift() ){
                            it2->updateInfor(targetGenome);
                            checkOrfState( (*it2), targetGenome, nucleotideCodeSubstitutionMatrix);
                            if( it2->getIfOrfShift() ){
                                std::cout << "There are ORF in de novo gene structure not complete\n" << it2->getName() << " " << it2->getMetaInformation()  << std::endl;
                            }else{
                                it2->setSource("DENOVOOVERLAP");
                                it2->setName(it3->second.getName());
                                adaptedNovoTranscriptAnnotation1[it2->getName()]=(*it2);
                            }
                        }
                    }
                }
                if( !ifHasOverLap ){
                    it2->updateInfor(targetGenome);
                    checkOrfState( (*it2), targetGenome, nucleotideCodeSubstitutionMatrix);
                    if( it2->getIfOrfShift() ){
                        std::cout << "There are ORF in de novo gene structure not complete\n" << it2->getName() << " " << it2->getMetaInformation()  << std::endl;
                    }else {
                        it2->setSource("DENOVO");
                        adaptedNovoTranscriptAnnotation2[it2->getName()] = (*it2);
                    }
                }
            }
        }
    }
    //replace old orf-shift annotation with novo annotation end

    std::cerr << "de novo annotation included" << std::endl;

    // remove those ORF-shift annotation begin
    std::vector<std::string> transcriptToRemove;
    for( std::map<std::string, Transcript>::iterator it=targetTranscriptsHashMap.begin();
         it!=targetTranscriptsHashMap.end(); ++it){
        if( it->second.getIfOrfShift() ){
            transcriptToRemove.push_back(it->first);
        }
    }

    for( std::vector<std::string>::iterator it=transcriptToRemove.begin();
            it!=transcriptToRemove.end(); ++it){
        if( targetTranscriptsHashMap.find(*it) != targetTranscriptsHashMap.end() ){
            targetTranscriptsHashMap.erase(*it);
        }
    }

    for(  std::map<std::string, Transcript>::iterator it=adaptedNovoTranscriptAnnotation1.begin();
         it!=adaptedNovoTranscriptAnnotation1.end(); ++it){
        targetTranscriptsHashMap[it->first]=it->second;
    }
    // remove those ORF-shift annotation end
    std::cerr << "ORF shift transcripts removed" << std::endl;

    TranscriptsTogenes(regexG, genes, targetTranscriptsHashMap);

    std::map<std::string, std::vector<Gene>> genes2;
    TranscriptsTogenes(novoRegexG, genes2, adaptedNovoTranscriptAnnotation2);

    for( std::map<std::string, std::vector<Gene>>::iterator it=genes2.begin();
            it!=genes2.end(); ++it){
        if( genes.find(it->first) == genes.end() ){
            genes[it->first]=std::vector<Gene>();
        }
        for( std::vector<Gene>::iterator it2=genes2[it->first].begin();
             it2!=genes2[it->first].end(); ++it2 ){
            genes[it->first].push_back(*it2);
        }
    }
    std::cerr << "transcript structure to gene structure done" << std::endl;
}

void transcriptRealignmentAndExonerate( Transcript& tartgetTranscript, Transcript& referenceTranscript,
                            NucleotideCodeSubstitutionMatrix& nucleotideCodeSubstitutionMatrix,
                            std::map<std::string, Fasta>& targetGenome,
                            std::map<std::string, Fasta>& referenceGenome, std::string chromosomeName, std::map<std::string,
        Transcript>& targetTranscriptsHashMap, std::atomic_int & number_of_runing_threads, std::string & prefixUuid, int & lengthThread ){
    std::string refGenomeSequence = referenceTranscript.getGeneomeSequence();
    Transcript targetTranscript(tartgetTranscript.getName(), chromosomeName, tartgetTranscript.getStrand());
    int startTarget = tartgetTranscript.getStart()-refGenomeSequence.length();
    int endTarget = tartgetTranscript.getEnd()+refGenomeSequence.length();

    if( endTarget > targetGenome[chromosomeName].getSequence().length() ){
        endTarget = targetGenome[chromosomeName].getSequence().length();
    }

    int startCodonPosition=1;
    int stopCodonPosition=refGenomeSequence.length()-2;
    std::vector<SpliceSitePosition> spliceSitePositions;

    targetTranscript.setSource("REALIGNMENT");
    int targetPosition=0;
    int referencePosition=0;
    //prepare for the new Cds After re-alignment begin
    std::vector<int> targetCdsStarts;
    std::vector<int> targetCdsEnds;
    std::string dna_b = getSubsequence(targetGenome, chromosomeName, startTarget, endTarget, referenceTranscript.getStrand());
    if( dna_b.length()<=lengthThread && refGenomeSequence.length()<=lengthThread ){ //  if the sequence is too long, don't try to align it
        if( referenceTranscript.getStrand() == POSITIVE ){
            if( referenceTranscript.getCdsVector().size()>1 ){
                for( size_t i=1; i<referenceTranscript.getCdsVector().size(); ++i ){
                    SpliceSitePosition spliceSitePosition(
                            referenceTranscript.getCdsVector()[i-1].getEnd()-referenceTranscript.getStart()+2,
                            referenceTranscript.getCdsVector()[i].getStart()-referenceTranscript.getStart());
                    spliceSitePositions.push_back(spliceSitePosition);
                }
            }

            NeedlemanWunschForTranscript nw (refGenomeSequence, dna_b, startCodonPosition, stopCodonPosition, spliceSitePositions);

            for( size_t tp = 0; tp<nw.getAlignment_a().length(); ++tp ){
                if( nw.getAlignment_b()[tp] != '-' ){
                    ++targetPosition;
                }
                if( nw.getAlignment_a()[tp] != '-' ){
                    ++referencePosition;
                    for( std::vector<Cds>::iterator it4=referenceTranscript.getCdsVector().begin();
                         it4!=referenceTranscript.getCdsVector().end(); ++it4){
                        if( referencePosition+referenceTranscript.getStart()-1 == (*it4).getStart() ){
                            //if( nw.getAlignment_b()[tp] == '-' ){
                                targetCdsStarts.push_back(startTarget + targetPosition-1); //(*i3) is the target transcript, refGenomeSequence.length() is the extend length
//                            }else{
//                                targetCdsStarts.push_back(startTarget + targetPosition -1);
//                            }
                        }
                        if( referencePosition+referenceTranscript.getStart()-1 == (*it4).getEnd() ){
                            //if( nw.getAlignment_b()[tp] == '-' ){
                                targetCdsEnds.push_back(startTarget + targetPosition-1);
//                            }else{
//                                targetCdsEnds.push_back(startTarget + targetPosition -1);
//                            }
                        }
                    }
                }
            }
        }else{
            if( referenceTranscript.getCdsVector().size()>1 ){
                for( size_t i=referenceTranscript.getCdsVector().size()-1; i>0; --i ){
                    SpliceSitePosition spliceSitePosition(
                            referenceTranscript.getEnd()-referenceTranscript.getCdsVector()[i].getStart()+2,
                            referenceTranscript.getEnd()-referenceTranscript.getCdsVector()[i-1].getEnd()-1);
                    spliceSitePositions.push_back(spliceSitePosition);
                }
            }
            NeedlemanWunschForTranscript nw (refGenomeSequence, dna_b, startCodonPosition, stopCodonPosition, spliceSitePositions);
//            gmutex.lock();
//            nw.print_results();
//            gmutex.unlock();

            for( size_t tp = 0; tp<nw.getAlignment_a().length(); ++tp ){
                if( nw.getAlignment_b()[tp] != '-' ){
                    ++targetPosition;
                }
                if( nw.getAlignment_a()[tp] != '-' ){
                    ++referencePosition;
                    for( std::vector<Cds>::iterator it4=referenceTranscript.getCdsVector().begin();
                         it4!=referenceTranscript.getCdsVector().end(); it4++){
                        //if( referencePosition+referenceTranscript.getStart()-1 == (*it4).getStart() ){
                        if( - referencePosition+referenceTranscript.getEnd()+1 == (*it4).getStart() ){
                            //if( nw.getAlignment_b()[tp] == '-' ){
                                targetCdsStarts.push_back(endTarget - targetPosition +1);
//                            }else{
//                                targetCdsStarts.push_back(endTarget - targetPosition +1);
//                            }
                        }
                        //if( referencePosition+referenceTranscript.getStart()-1 == (*it4).getEnd() ){
                        if( - referencePosition+referenceTranscript.getEnd()+1 == (*it4).getEnd() ){
                            //if( nw.getAlignment_b()[tp] == '-' ){
                                targetCdsEnds.push_back(endTarget - targetPosition +1);
//                            }else{
//                                targetCdsEnds.push_back(endTarget - targetPosition +1);
//                            }
                        }
                    }
                }
            }
        }
        for(size_t i5=0; i5<targetCdsStarts.size(); i5++ ){
            Cds cds(targetCdsStarts[i5], targetCdsEnds[i5]);
            //std::cout << "CDS " << targetCdsStarts[i5] << " " << targetCdsEnds[i5] << std::endl;
            targetTranscript.addCds(cds);
        }
        targetTranscript.updateInfor(targetGenome);
//        checkOrfState( targetTranscript, referenceTranscript,
//                       targetGenome, referenceGenome, nucleotideCodeSubstitutionMatrix);
        checkOrfState( targetTranscript, targetGenome,  nucleotideCodeSubstitutionMatrix);

        if( targetTranscript.getIfOrfShift() ){
            std::string cdsSequence=referenceTranscript.getCdsSequence();
            if(cdsSequence.length() >0 && cdsSequence.length() < lengthThread*4){
                std::string transcriptName = referenceTranscript.getName();
                std::string targetSequence = getSubsequence(targetGenome, referenceTranscript.getChromeSomeName(),
                                                            startTarget, endTarget, POSITIVE);

                runExonerateEst(transcriptName, cdsSequence, targetSequence,
                                nucleotideCodeSubstitutionMatrix, targetTranscriptsHashMap, startTarget, endTarget,
                                referenceTranscript.getStrand(), chromosomeName, prefixUuid, targetGenome);
            }
        }else{
            gmutex.lock();
            targetTranscriptsHashMap[targetTranscript.getName()]=targetTranscript; // update the map data structure with modified transcript
            gmutex.unlock();
        }
    }else{
        std::string cdsSequence=referenceTranscript.getCdsSequence();
        if(cdsSequence.length() >0 && cdsSequence.length() < lengthThread*4 ){
            std::string transcriptName = referenceTranscript.getName();
            std::string targetSequence = getSubsequence(targetGenome, referenceTranscript.getChromeSomeName(),
                                                        startTarget, endTarget, POSITIVE);

            runExonerateEst(transcriptName, cdsSequence, targetSequence,
                            nucleotideCodeSubstitutionMatrix, targetTranscriptsHashMap, startTarget, endTarget,
                            referenceTranscript.getStrand(), chromosomeName, prefixUuid, targetGenome);
        }
    }
    --number_of_runing_threads;
    return;
}

void runExonerateEst(std::string& transcriptName, std::string& cdsSequence, std::string& targetSequence,
                     NucleotideCodeSubstitutionMatrix& nucleotideCodeSubstitutionMatrix,
                     std::map<std::string, Transcript>& targetTranscriptsHashMap, int& startTarget, int& endTarget,
                     STRAND strand, std::string& tchromeSomeName, std::string& prefixUuid , std::map<std::string, Fasta>& targetGenome ){
    try{
        std::string tempFolder = getTempFolder();
        gmutex.lock();
        std::string tempFile = tempFolder + "/" + generateUUID( prefixUuid );
        gmutex.unlock();
        std::string cdsFile = tempFile + "cds.fasta";
        std::string targetFile = tempFile + "target.fasta";

        std::ofstream ofile2;
        ofile2.open(cdsFile);
        ofile2 << ">cds" << std::endl << cdsSequence << std::endl;
        ofile2.close();

        std::ofstream ofile;
        ofile.open(targetFile);
        ofile << ">target" << std::endl << targetSequence << std::endl;
        ofile.close();

        std::string command =
                "exonerate --maxintron 30000 --model est2genome -i -10 --score 10 --bestn 1 --minintron 10 " +
                cdsFile +
                " " + targetFile + " --showtargetgff true >" + tempFile;
        system(&command[0]);
        readExonerateEstResult(transcriptName, cdsSequence, targetSequence,
                               nucleotideCodeSubstitutionMatrix,
                               targetTranscriptsHashMap, startTarget, endTarget,
                               strand, tchromeSomeName, tempFile, targetGenome);
        std::string cleanFileCommand = "rm " + cdsFile;
        system(&cleanFileCommand[0]);
        cleanFileCommand = "rm " + targetFile;
        system(&cleanFileCommand[0]);
        cleanFileCommand = "rm " + tempFile;
        system(&cleanFileCommand[0]);
    }catch (...){
        return;
    }
    return;
}


void readExonerateEstResult( std::string& transcriptName, std::string& cdsSequence, std::string& targetSequence,
                             NucleotideCodeSubstitutionMatrix& nucleotideCodeSubstitutionMatrix,
                             std::map<std::string, Transcript>& targetTranscriptsHashMap, int& startTarget, int& endTarget,
                             STRAND& strand, std::string& tchromeSomeName, std::string& fileLocation, std::map<std::string, Fasta>& targetGenome ){

    if( (!if_file_exists(fileLocation)) || GetFileSize(fileLocation)<1 ){
        return;
    }

    Transcript targetTranscript(transcriptName, tchromeSomeName, strand );
    int cdsNumber = 0;
    if( POSITIVE ==strand ){
        std::ifstream infile2(fileLocation);
        std::regex reg2("^(target)\\t(\\S*)\\texon\\t(\\S*)\\t(\\S*)\\t(\\S*)\\t\\+");
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
                targetTranscript.addCds(cds);
                ++cdsNumber;
            }
        }
        infile2.close();
    }else{
        std::ifstream infile2(fileLocation);
        std::regex reg2("^(target)\\t(\\S*)\\texon\\t(\\S*)\\t(\\S*)\\t(\\S*)\\t\\-");
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
                targetTranscript.addCds(cds);
                ++cdsNumber;
            }
        }
        infile2.close();
    }
    if( cdsNumber > 0 && targetTranscript.getCdsVector().size() > 0 ){
//        std::cout << transcriptName << " exonerate est updateInfor begin" << std::endl;

        targetTranscript.updateInfor(targetGenome);

//        std::cout << transcriptName << " exonerate est checkOrfState begin" << std::endl;
//        checkOrfState( targetTranscript, referenceTranscript,
//                       targetGenome, referenceGenome, nucleotideCodeSubstitutionMatrix);
        checkOrfState( targetTranscript, targetGenome, nucleotideCodeSubstitutionMatrix);
//        std::cout << transcriptName << " exonerate est checkOrfState finished" << std::endl;

        if( targetTranscript.getIfOrfShift() ){

//            if( targetTranscript.getMetaInformation().find("spliceSitesDestroyed")!=std::string::npos ){
//                std::cout << "879 there is something wrong with the splice sites of cds based alignment " << targetTranscript.getMetaInformation() << std::endl;
//            }else{
//                std::cout << "881 " << targetTranscript.getMetaInformation() << std::endl;
//            }
            std::string protenSequene = nA2AA(cdsSequence, nucleotideCodeSubstitutionMatrix);
            runExonerateProtein(transcriptName, protenSequene, targetSequence,
                                nucleotideCodeSubstitutionMatrix,
                                targetTranscriptsHashMap, startTarget, endTarget,
                                strand, tchromeSomeName, fileLocation, targetGenome );
        }else{
            targetTranscript.setSource("CDSALIGNMENT");
            gmutex.lock();
            targetTranscriptsHashMap[transcriptName]=targetTranscript;
            gmutex.unlock();
        }
    }else{
        std::string protenSequene = nA2AA(cdsSequence, nucleotideCodeSubstitutionMatrix);
        runExonerateProtein(transcriptName, protenSequene, targetSequence,
                            nucleotideCodeSubstitutionMatrix,
                            targetTranscriptsHashMap, startTarget, endTarget,
                            strand, tchromeSomeName, fileLocation, targetGenome );
    }
    return;
}

void runExonerateProtein(std::string& transcriptName, std::string& protenSequene, std::string& targetSequence,
                             NucleotideCodeSubstitutionMatrix& nucleotideCodeSubstitutionMatrix,
                             std::map<std::string, Transcript>& targetTranscriptsHashMap, int& startTarget, int& endTarget,
                        STRAND& strand, std::string& tchromeSomeName, std::string& fileLocation , std::map<std::string, Fasta>& targetGenome
                         ) {
    std::string proteinFile = fileLocation + "protein.fasta";
    std::string targetFile = fileLocation + "target.fasta";
    std::string tempFile=fileLocation + "protein.aln";

    std::ofstream ofile;
    ofile.open(proteinFile);
    ofile << ">protein" << std::endl << protenSequene << std::endl;
    ofile.close();

    std::string command = "exonerate --bestn 1 --maxintron 30000 --intronpenalty -10 --model protein2genome --percent 10 --score 10 --minintron 10 "+proteinFile + " " + targetFile + " --showtargetgff true >" +tempFile;
    system(&command[0]);
    readExonerateProteinResult(tempFile, nucleotideCodeSubstitutionMatrix, startTarget, endTarget, strand, targetTranscriptsHashMap, transcriptName, tchromeSomeName, targetGenome) ;
    std::string cleanFileCommand = "rm " + proteinFile;
    system(&cleanFileCommand[0]);
    cleanFileCommand = "rm " + tempFile;
    system(&cleanFileCommand[0]);
    return;
}

void readExonerateProteinResult( std::string& fileLocation, NucleotideCodeSubstitutionMatrix& nucleotideCodeSubstitutionMatrix,
                                 int& startTarget, int& endTarget, STRAND& strand, std::map<std::string, Transcript>& targetTranscriptsHashMap,
                                 std::string& transcriptName, std::string& tchromeSomeName, std::map<std::string, Fasta>& targetGenome){
    if( (!if_file_exists(fileLocation)) || GetFileSize(fileLocation)<1 ){
        return;
    }

    Transcript targetTranscript(transcriptName, tchromeSomeName, strand );
    targetTranscript.setSource("PROTEINALIGNMENT");
    int cdsNumber = 0;
    if( POSITIVE ==strand ){
        std::ifstream infile(fileLocation);
        std::regex reg("^(target)\t([\\s\\S]*)\tcds\t(\\S*)\t(\\S*)\t(\\S*)\t\\+");
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
                targetTranscript.addCds(cds);
                ++cdsNumber;
                //std::cout << start << " " << end << std::endl;
            }
        }
        infile.close();
    }else {
        std::ifstream infile(fileLocation);
        std::regex reg("^(target)\t([\\s\\S]*)\tcds\t(\\S*)\t(\\S*)\t(\\S*)\t\\-");
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
                targetTranscript.addCds(cds);
                ++cdsNumber;
            }
        }
        infile.close();
    }
    if( cdsNumber > 0 && targetTranscript.getCdsVector().size() > 0 ) {
        targetTranscript.updateInfor(targetGenome);
        checkOrfState(targetTranscript, targetGenome, nucleotideCodeSubstitutionMatrix);
        if (!targetTranscript.getIfOrfShift()) {
            gmutex.lock();
            targetTranscriptsHashMap[transcriptName] = targetTranscript;
            gmutex.unlock();
        }
//        else{
//            if( targetTranscript.getMetaInformation().find("spliceSitesDestroyed")!=std::string::npos ){
//                std::cout << "991 there is something wrong with the splice sites of cds based alignment " << targetTranscript.getMetaInformation() << std::endl;
//            }else{
//                std::cout << "993 " << targetTranscript.getMetaInformation() << std::endl;
//            }
//        }
    }
    return;
}

void readAugustusGff( std::string& fileLocation ){
    std::map<std::string, std::vector<Transcript> > transcriptHashSet;
    std::string cdsRegex = "ID=([\\s\\S]*?);Parent=([\\s\\S]*?)$";
    readGffFile (fileLocation, transcriptHashSet, cdsRegex);
    for( std::map<std::string, std::vector<Transcript> >::iterator it=transcriptHashSet.begin(); it!=transcriptHashSet.end(); ++it){
        std::cout << it->first << std::endl;
    }
}

















void reAnnotationAndExonerate2( std::string& referenceGenomeFilePath, std::string referenceGffFilePath, std::string sdiFile,
                               std::map<std::string, std::vector<Gene> >& genes, int maxThread, std::string& regex, std::string& regexG, std::string& outputGffFile, int & lengthThread){
    NucleotideCodeSubstitutionMatrix nucleotideCodeSubstitutionMatrix;
    std::map<std::string, Fasta> referenceGenome;
    std::map<std::string, std::vector<Transcript> > referenceTranscriptHashSet;
    readFastaFile(referenceGenomeFilePath, referenceGenome);
    readGffFile (referenceGffFilePath, referenceTranscriptHashSet, regex);

    std::map<std::string, std::vector<Variant> > variantsMaps;
    readSdiFile (sdiFile, variantsMaps);
    std::map<std::string, Fasta> targetGenome;
    getPseudoGenomeSequence(referenceGenome, variantsMaps, targetGenome);
    std::map<std::string, Transcript> targetTranscriptsHashMap;// transcriptName Transcript
    annotationLiftOver( referenceTranscriptHashSet, targetTranscriptsHashMap, variantsMaps, targetGenome, referenceGenome);

    std::cerr << "annotationLiftOver done" << std::endl;

    std::map<std::string, Transcript > referenceTranscriptHashMap; // transcriptName Transcript

    for(std::map<std::string, std::vector<Transcript> >::iterator it=referenceTranscriptHashSet.begin();
        it!=referenceTranscriptHashSet.end(); ++it){
        for( std::vector<Transcript>::iterator it2=it->second.begin(); it2!=it->second.end();++it2 ){
            it2->updateInfor(referenceGenome);
        }
    }
    reAnnotationSingleLine( referenceTranscriptHashSet, referenceGenome, targetTranscriptsHashMap,
                            targetGenome, nucleotideCodeSubstitutionMatrix, maxThread, lengthThread);

    std::cerr << "realignment done" << std::endl;
    std::map<std::string, Transcript> homologousGoodTranscripts;
    std::string tempFolder = createdTempFloder(); // create tempFolder
    std::atomic_int number_of_runing_threads(0);
    //time_t rawtime;
    //std::cerr << time(&rawtime) << " start time" << std::endl;

    for( std::map<std::string, std::vector<Transcript> >::iterator itchromosome=referenceTranscriptHashSet.begin();
         itchromosome!=referenceTranscriptHashSet.end(); ++itchromosome){
        if( referenceGenome.find(itchromosome->first) != referenceGenome.end() && targetGenome.find(itchromosome->first) != targetGenome.end()){
            for( std::vector<Transcript>::iterator itTranscript=referenceTranscriptHashSet[itchromosome->first].begin();
                 itTranscript!=referenceTranscriptHashSet[itchromosome->first].end(); ++itTranscript){
                if( targetTranscriptsHashMap[(*itTranscript).getName()].getIfOrfShift() ){
                    std::string refGenomeSequence=(*itTranscript).getGeneomeSequence();
                    std::string cdsSequence=(*itTranscript).getCdsSequence();
                    int startTarget = targetTranscriptsHashMap[(*itTranscript).getName()].getStart()-refGenomeSequence.length();
                    int endTarget = targetTranscriptsHashMap[(*itTranscript).getName()].getEnd()+refGenomeSequence.length();
                    std::string targetSequence = getSubsequence(targetGenome, (*itTranscript).getChromeSomeName(),
                                                                startTarget, endTarget, POSITIVE);

                    if (cdsSequence.length() > 0 && targetSequence.length() > 0) {
                        std::string tempFile = tempFolder + "/" + generateUUID(outputGffFile);
                        std::string cdsFile = tempFile + "cds.fasta";
                        std::string targetFile = tempFile + "target.fasta";

                        std::ofstream ofile;
                        ofile.open(cdsFile);
                        ofile << ">cds" << std::endl << cdsSequence << std::endl;
                        ofile.close();

                        std::ofstream ofile2;
                        ofile2.open(targetFile);
                        ofile2 << ">target" << std::endl << targetSequence << std::endl;
                        ofile.close();
                        std::string command =
                                "exonerate --maxintron 30000 --model est2genome -i -10 --score 10 --bestn 1 --minintron 10 " +
                                cdsFile +
                                " " + targetFile + " --showtargetgff true >" + tempFile;

                        bool isThisThreadUnrun = true;
                        while (isThisThreadUnrun) {
                            if (number_of_runing_threads < maxThread) {
                                std::thread t(runSystemCommand, command, std::ref(number_of_runing_threads));
                                t.detach();
                                ++number_of_runing_threads;
                                isThisThreadUnrun = false;
                                break;
                            } else {
                                usleep(10);
                            }
                        }
                    }
                }
            }
        }else{
            referenceTranscriptHashSet.erase(itchromosome->first);
        }
    }
    while(number_of_runing_threads > 0 ){// wait for all the thread
        usleep(50);
    }
    std::chrono::milliseconds ms2 = std::chrono::duration_cast< std::chrono::milliseconds >(
            std::chrono::system_clock::now().time_since_epoch()
    );

    //std::cerr << time(&rawtime) << " end time" << std::endl;

//    for( std::map<std::string, Transcript>::iterator it = homologousGoodTranscripts.begin();
//         it!=homologousGoodTranscripts.end(); ++it){
//        if( ! it->second.getIfOrfShift() ){
//            it->second.updateInfor(targetGenome);
//            targetTranscriptsHashMap[it->second.getName()] = it->second;
//        }
//    }
    std::cerr << "homologous CDS based annotation done" << std::endl;
/*
    std::map<std::string, Transcript> homologousGoodTranscripts2;

    for( std::map<std::string, std::vector<Transcript> >::iterator itchromosome=referenceTranscriptHashSet.begin();
         itchromosome!=referenceTranscriptHashSet.end(); ++itchromosome){
        if( referenceGenome.find(itchromosome->first) != referenceGenome.end()
            && targetGenome.find(itchromosome->first) != targetGenome.end()) {
            for (std::vector<Transcript>::iterator itTranscript = referenceTranscriptHashSet[itchromosome->first].begin();
                 itTranscript != referenceTranscriptHashSet[itchromosome->first].end(); ++itTranscript) {
                if (targetTranscriptsHashMap[(*itTranscript).getName()].getIfOrfShift()) {
                    std::string refGenomeSequence = (*itTranscript).getGeneomeSequence();
                    std::string cdsSequence = (*itTranscript).getCdsSequence();

                    int startTarget =
                            targetTranscriptsHashMap[(*itTranscript).getName()].getStart() - refGenomeSequence.length();
                    int endTarget =
                            targetTranscriptsHashMap[(*itTranscript).getName()].getEnd() + refGenomeSequence.length();
                    std::string targetSequence = getSubsequence(targetGenome, (*itTranscript).getChromeSomeName(),
                                                                startTarget, endTarget, POSITIVE);

                    std::string protenSequene = nA2AA(cdsSequence, nucleotideCodeSubstitutionMatrix);
                    if( protenSequene.length()>0 && targetSequence.length() > 0 ) {

                        std::string tempFile = tempFolder + "/" + generateUUID();
                        std::string proteinFile = tempFile + "protein.fasta";
                        std::string targetFile = tempFile + "target.fasta";

                        std::ofstream ofile;
                        ofile.open(proteinFile);
                        ofile << ">cds" << std::endl << protenSequene << std::endl;
                        ofile.close();

                        std::ofstream ofile2;
                        ofile2.open(targetFile);
                        ofile2 << ">target" << std::endl << targetSequence << std::endl;
                        ofile2.close();
                        std::string command = "exonerate --bestn 1 --maxintron 30000 --intronpenalty -10 --model protein2genome --percent 10 --score 10 --minintron 10 "+proteinFile + " " + targetFile + " --showtargetgff true >" +tempFile;
                        bool isThisThreadUnrun = true;
                        while (isThisThreadUnrun) {
                            if (number_of_runing_threads < maxThread) {
                                std::thread t(runExonerateEst, command, std::ref(number_of_runing_threads));
                                t.detach();
                                ++number_of_runing_threads;
                                isThisThreadUnrun = false;
                                break;
                            } else {
                                usleep(10);
                            }
                        }
                    }
                }
            }
        }else{
            referenceTranscriptHashSet.erase(itchromosome->first);
            if( targetTranscriptsHashMap.find(itchromosome->first) != targetTranscriptsHashMap.end() ){
                targetTranscriptsHashMap.erase(itchromosome->first);
            }
        }
    }
    while( number_of_runing_threads > 0 ){// wait for all the thread
        usleep(50);
    }
    for( std::map<std::string, Transcript>::iterator it = homologousGoodTranscripts2.begin();
         it!=homologousGoodTranscripts2.end(); ++it){
        if( ! it->second.getIfOrfShift() ){
            it->second.updateInfor(targetGenome);
            targetTranscriptsHashMap[it->second.getName()] = it->second;
        }
    }*/
    std::cerr << "homologous protein based annotation done" << std::endl;
    TranscriptsTogenes(regexG, genes, targetTranscriptsHashMap);
    std::cerr << "transcript structure to gene structure done" << std::endl;
}

void runSystemCommand( std::string command, std::atomic_int & number_of_runing_threads ){
    std::cerr << command << " " << number_of_runing_threads << std::endl;
    system(&command[0]);
    --number_of_runing_threads;
}



void readExonerateEstResult2( std::string& transcriptName, std::string& cdsSequence, std::string& targetSequence,
                              NucleotideCodeSubstitutionMatrix& nucleotideCodeSubstitutionMatrix,
                              std::map<std::string, Transcript>& targetTranscriptsHashMap, int& startTarget, int& endTarget,
                              STRAND& strand, std::string& tchromeSomeName, std::string& fileLocation, std::map<std::string, Fasta>& targetGenome,
                              Transcript& referenceTranscript, std::map<std::string, Fasta>& referenceGenome){

    if( (!if_file_exists(fileLocation)) || GetFileSize(fileLocation)<1 ){
        return;
    }

    std::ifstream infile(fileLocation);
    std::regex reg("\\d\\s*:\\s*([\\w\\.\\s\\>-]+)\\s*:\\s*\\d");
    std::stringstream sequencestream;
    std::string line;
    int lineNumber = 0;
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

    std::stringstream metaInformation;
    std::string cdsSequenceString = sequencestream.str();
    std::regex vowel_re("\\s");
    cdsSequenceString=std::regex_replace(cdsSequenceString, vowel_re, "");
    std::regex vowel_re2("[a-z]+\\.+[a-z]+");
    cdsSequenceString=std::regex_replace(cdsSequenceString, vowel_re2, "");
    std::regex vowel_re3("-");
    cdsSequenceString=std::regex_replace(cdsSequenceString, vowel_re3, "");
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

    Transcript targetTranscript(transcriptName, tchromeSomeName, strand );

    if( !orfShift && POSITIVE ==strand ){
        metaInformation << "_ConservedFunction";
        std::string mi = metaInformation.str();

        targetTranscript.setIfOrfShift(false);
        targetTranscript.setSource("CDSALIGNMENT");
        targetTranscript.setMetaInformation(mi);

        std::ifstream infile2(fileLocation);

        std::regex reg2("^(\\S*)\\t([\\s\\S]*)\\texon\\t(\\S*)\\t(\\S*)\\t(\\S*)\\t\\+");
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
                targetTranscript.addCds(cds);
            }
        }
        infile2.close();
        targetTranscript.updateInfor(targetGenome);
        gmutex.lock();
        targetTranscriptsHashMap[transcriptName]=targetTranscript;
        gmutex.unlock();
    }else if( !orfShift ){
        metaInformation << "_ConservedFunction";
        std::string mi = metaInformation.str();
        targetTranscript.setIfOrfShift(false);
        targetTranscript.setSource("CDSALIGNMENT");
        targetTranscript.setMetaInformation(mi);

        std::ifstream infile2(fileLocation);
        std::regex reg2("^(\\S*)\\t([\\s\\S]*)\\texon\\t(\\S*)\\t(\\S*)\\t(\\S*)\\t\\-");
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
                targetTranscript.addCds(cds);
            }
        }
        infile2.close();
        targetTranscript.updateInfor(targetGenome);
        gmutex.lock();
        targetTranscriptsHashMap[transcriptName]=targetTranscript;
        gmutex.unlock();
    }else{
        std::string protenSequene = nA2AA(cdsSequence, nucleotideCodeSubstitutionMatrix);
        runExonerateProtein(transcriptName, protenSequene, targetSequence,
                            nucleotideCodeSubstitutionMatrix,
                            targetTranscriptsHashMap, startTarget, endTarget,
                            strand, tchromeSomeName, fileLocation, targetGenome );
    }
    return;
}



void readExonerateProteinResult2( std::string& fileLocation, NucleotideCodeSubstitutionMatrix& nucleotideCodeSubstitutionMatrix,
                                  int& startTarget, int& endTarget, STRAND& strand, std::map<std::string, Transcript>& targetTranscriptsHashMap,
                                  std::string& transcriptName, std::string& tchromeSomeName, std::map<std::string, Fasta>& targetGenome){
    if( (!if_file_exists(fileLocation)) || GetFileSize(fileLocation)<1 ){
        return;
    }

    std::ifstream infile(fileLocation);
    std::regex reg("^\\s*\\d+\\s*:\\s*([\\w\\.\\s\\>-\\{\\}\\*]+)\\s*:\\s*\\d+\\s*$");
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
    std::cout << cdsSequenceString << std::endl;
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
    std::regex vowel_re6("-");
    cdsSequenceString=std::regex_replace(cdsSequenceString, vowel_re6, "");
    std::cout << cdsSequenceString << std::endl;
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

    Transcript targetTranscript(transcriptName, tchromeSomeName, strand );

    if( !orfShift && POSITIVE ==strand ){
        metaInformation << "_ConservedFunction";
        std::string mi = metaInformation.str();
        targetTranscript.setSource("PROTEINALIGNMENT");
        targetTranscript.setIfOrfShift(false);
        targetTranscript.setMetaInformation(mi);

        std::ifstream infile(fileLocation);
        std::regex reg("^(\\S*)\t([\\s\\S]*)\tcds\t(\\S*)\t(\\S*)\t(\\S*)\t\\+");
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
                targetTranscript.addCds(cds);
                //std::cout << start << " " << end << std::endl;
            }
        }
        infile.close();
        targetTranscript.updateInfor(targetGenome);
        gmutex.lock();
        targetTranscriptsHashMap[transcriptName]=targetTranscript;
        gmutex.unlock();
    }else if( !orfShift ){
        metaInformation << "_ConservedFunction";
        std::string mi = metaInformation.str();

        targetTranscript.setSource("PROTEINALIGNMENT");
        targetTranscript.setIfOrfShift(false);
        targetTranscript.setMetaInformation(mi);

        std::ifstream infile(fileLocation);
        std::regex reg("^(\\S*)\t([\\s\\S]*)\tcds\t(\\S*)\t(\\S*)\t(\\S*)\t\\-");
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
                targetTranscript.addCds(cds);
            }
        }
        infile.close();
        targetTranscript.updateInfor(targetGenome);
        gmutex.lock();
        targetTranscriptsHashMap[transcriptName]=targetTranscript;
        gmutex.unlock();
    }
    return;
}




//for MSA begin
void songPopulationMsa( std::string referenceGenomeFilePath, std::string referenceGffFilePath, std::string accessionIdList,
    std::string targetSdiPath, std::string targetGffPath, int windownsSize, int overlapSize, std::string& regex){

    std::map<std::string, Fasta> referenceGenome;
    std::map<std::string, std::vector<Transcript> > referenceTranscriptHashSet;
    readFastaFile(referenceGenomeFilePath, referenceGenome);
    readGffFile (referenceGffFilePath, referenceTranscriptHashSet, regex);

    std::ifstream infile(accessionIdList);
    std::string line="";
    while (std::getline(infile, line)){

    }
    return;
}


