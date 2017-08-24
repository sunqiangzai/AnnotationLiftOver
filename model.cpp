// =====================================================================================
// 
//       Filename:  model.cpp
// 
//    Description:  
// 
//        Version:  1.0
//        Created:  04/12/2017 02:45:02 PM
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
#include <iostream>
#include "Cds.h"
#include "model.h"
#include <sstream>
#include "myutil.h"

//Fasta class begin
Fasta::Fasta(std::string& _name, std::string& _sequence){
    this->name=_name;
    this->sequence=_sequence;
}
Fasta::Fasta(){

}
void Fasta::setName(std::string& _name){
    name=_name;
}
void Fasta::setSequence(std::string& _sequence){
    sequence=_sequence;
}
Fasta::Fasta( const Fasta& f ){
    this->name=f.name;
    this->sequence=f.sequence;
}
std::string Fasta::getSubsequence( size_t& _start,  size_t& _end){
    size_t start = _start;
    size_t end = _end;
    if( start > end ){
        size_t temp = start;
        start=end;
        end=temp;
    }
    if( start < 1 ){
        start = 1;
    }
    if( end>sequence.length() ){
        end = sequence.length();
    }
    return sequence.substr(start-1, end-start+1);
}
std::string Fasta::getSubsequence( int& _start,  int& _end){
    size_t start = _start;
    size_t end = _end;
    if( start > end ){
        size_t temp = start;
        start=end;
        end=temp;
    }
    if( start < 1 ){
        start = 1;
    }
    if( end>sequence.length() ){
        end = sequence.length();
    }
    return sequence.substr(start-1, end-start+1);
}
std::ostream& print(std::ostream& out, const Fasta& f){
    out << ">" << f.getName() << std::endl << f.getSequence() << std::endl;
    return out;
}
//Fasta class end

//variant record class begin
//Variant::Variant(std::string _chromosome, int _position , int _changingLength,
//    std::string _reference, std::string _alternative): chromosome(_chromosome),
//    position(_position), changingLength(_changingLength),reference(_reference),
//    alternative(_alternative){
//}
Variant::Variant(std::string& _chromosome, int& _position, std::string& _reference,
        std::string& _alternative){
    int alternativeSize = _alternative.size();
    if( _alternative.compare("-") ==0  ){
        alternativeSize = 0;
    }else{
        std::size_t found = _alternative.find("-");
        if (found!=std::string::npos){
            std::cerr << "the record contains strange - at: " << _chromosome << "\t" << _position << std::endl;
            exit(1);
        }
    }
    int referenceSize = _reference.size();
    if( _reference.compare("-") ==0  ){
        referenceSize=0;
    }else {
        std::size_t found = _reference.find("-");
        if (found!=std::string::npos){
            std::cerr << "the record contains strange - at: " << _chromosome << "\t" << _position << std::endl;
            exit(1);
        }
    }
    changingLength = alternativeSize - referenceSize;
    chromosome=_chromosome;
    position=_position;
    reference=_reference;
    alternative=_alternative;
//    Variant(_chromosome, _position, changingLength, _reference, _alternative);
}

Variant::Variant (const Variant& variant){
    this->chromosome=variant.chromosome;
    this->position = variant.position;
    this->changingLength = variant.changingLength;
    this->reference = variant.reference;
    this->alternative = variant.alternative;

}
bool Variant::overlap(int _position){
    if( this->changingLength>=0 ){
        if( this->position == _position ){
            return true;
        }
    }else{
        if( (this->position-this->changingLength)>=_position && this->position<=_position  ){
            return true;
        }
    }
    return false;
}
//
//bool Variant::overlap (const Variant& _variant){
//    if( _variant.changingLength >=0  ){
//        if( this->changingLength>=0 ){
//            if( this->getPosition() == _variant.getPosition() ){
//                return true;
//            }else{
//                return false;
//            }
//        }else{
//            if( _variant.changingLength ){
//
//            }else{
//
//            }
//        }
//    }
//}

std::ostream &print(std::ostream& out, const Variant& variant){
    out << variant.getChromosome() << "\t" << variant.getPosition() << "\t" <<
        variant.getChanginglength() << "\t" << variant.getReference() << "\t" << variant.getAlternative();
    return out;
}
std::ostream &println(std::ostream& out, const Variant& variant){
    out << variant.getChromosome() << "\t" << variant.getPosition() << "\t" <<
        variant.getChanginglength() << "\t" << variant.getReference() << "\t" << variant.getAlternative() << std::endl;
    return out;
}
//variant record class end



//gff related class begin

//trabscript begin
//Transcript::Transcript(std::string name, std::set<Cds> cdsHashSet, STRAND strand, std::string chromeSomeName):
//    _name(name),_cdsHashSet(cdsHashSet), _strand(strand), _chromeSomeName(chromeSomeName){
//}
// Transcript::Transcript(){

Transcript::Transcript(std::string name, std::string& chromeSomeName, STRAND strand){
    this->_name=name;
    this->_chromeSomeName=chromeSomeName;
    this->_strand=strand;
    this->_start =  std::numeric_limits<int>::max();
    this->_end=0;
    _cdsHashSet=std::set<Cds>();
}
Transcript::Transcript(){
    this->_strand=POSITIVE;
    this->_start =  std::numeric_limits<int>::max();
    this->_end=0;
    _cdsHashSet=std::set<Cds>();
}
std::string Transcript::getName(){
    return _name;
}
void Transcript::setName(std::string name){
    this->_name=name;
}
std::set<Cds>& Transcript::getCdsHashSet(){
    return _cdsHashSet;
}
std::vector<Cds>& Transcript::getCdsVector(){
    return _cdsVector;
}
STRAND Transcript::getStrand(){
    return _strand;
}
std::string Transcript::getChromeSomeName(){
    return _chromeSomeName;
}
void Transcript::addCds(Cds& cds){
    this->_cdsHashSet.insert(cds);
    this->_cdsVector.push_back(cds);
}
int Transcript::getStart(){
    return _start;
}
int Transcript::getEnd(){
    return _end;
}
void Transcript::updateInfor() {
    this->_start =  std::numeric_limits<int>::max();
    this->_end=0;
    for( std::vector<Cds>::iterator it=this->_cdsVector.begin(); it!=this->_cdsVector.end(); it++ ){
        if( (*it).getStart() < this->_start) {
            this->_start = (*it).getStart();
        }
        if( (*it).getEnd() > this->_end ){
            this->_end = (*it).getEnd();
        }
    }
    std::sort(this->_cdsVector.begin(), this->_cdsVector.end(), [](Cds a, Cds b) {
        return a < b;
    });
}

void Transcript::updateInfor(std::map<std::string, Fasta>& genome) {
    this->updateInfor();
    this->_geneomeSequence = getSubsequence(genome, _chromeSomeName, _start, _end, _strand);
    std::string reserve;
    reserve.reserve(this->_geneomeSequence.size());
    std::stringstream cdsss(reserve);
    if( POSITIVE == _strand ){
        for( size_t i=0; i<_cdsVector.size(); i++ ){
            //std::cout << "193 " << this->_name << " " << _chromeSomeName << " " << _cdsVector[i].getStart() << " " << _cdsVector[i].getEnd() << " " << _strand << std::endl;
            cdsss << getSubsequence(genome, _chromeSomeName, _cdsVector[i].getStart(), _cdsVector[i].getEnd(), _strand);
        }
    }else{
        for( size_t i=_cdsVector.size(); i>0; i-- ){
            //std::cout << "199 " << i-1 << " " << this->_name << " " << _chromeSomeName << " " << _cdsVector[i-1].getStart() << " " << _cdsVector[i-1].getEnd() << " " << _strand << std::endl;
            cdsss << getSubsequence(genome, _chromeSomeName, _cdsVector[i-1].getStart(), _cdsVector[i-1].getEnd(), _strand);
        }
    }
    this->_cdsSequence = cdsss.str();
}
std::string Transcript::getGeneomeSequence(){
    return this->_geneomeSequence;
}
std::string Transcript::getCdsSequence(){
    return this->_cdsSequence;
}
void Transcript::setMetaInformation(std::string& metaInformation){
    this->_metaInformation=metaInformation;
}
std::string& Transcript::getMetaInformation(){
    return this->_metaInformation;
}
bool& Transcript::getIfOrfShift(){
    return _ifOrfShift;
}
void Transcript::cleanCdsVector(){
    this->_cdsVector.clear();
}
bool Transcript::ifOverLap( Transcript& transcript ){
    if( this->_strand == transcript._strand ){
        if( this->_start <= transcript._start &&
                transcript._start <= this->_end){
            return true;
        }
        if( this->_start <= transcript._end &&
            transcript._end <= this->_end){
            return true;
        }
        if( transcript._start <= this->_start &&
                this->_start <= transcript._end){
            return true;
        }
        if( transcript._start <= this->_end &&
            this->_end <= transcript._end){
            return true;
        }
    }
    return false;
}
void Transcript::setIfOrfShift(bool ifOrfShift){
    this->_ifOrfShift=ifOrfShift;
}
std::string Transcript::getSource(){
    return this->_source;
}
void Transcript::setSource(std::string source){
    this->_source=source;
}
//transcript end

//gene begin
Gene::Gene(std::string& name, STRAND& strand){
    this->_start = std::numeric_limits<int>::max();
    this->_end=0;
    this->_name=name;
    this->_strand=strand;
}
Gene::Gene() {
    this->_start = std::numeric_limits<int>::max();
    this->_end=0;
}
std::string Gene::getName(){
    return this->_name;
}
STRAND Gene::getStrand(){
    return this->_strand;
}
int Gene::getStart(){
    return this->_start;
}
int Gene::getEnd(){
    return this->_end;
}
std::string Gene::getChromeSomeName(){
    return this->_chromeSomeName;
}
std::vector<Transcript>& Gene::getTranscriptVector(){
    return this->_transcriptVector;
}
bool Gene::checkOverLapAndAddTranscript(Transcript& transcript){
    if( this->_strand == transcript.getStrand() ){
        if( (this->_start < transcript.getStart() && transcript.getStart() > this->_end )
            && (this->_start < transcript.getEnd() && transcript.getEnd() > this->_end ) ){
            this->addTranscript(transcript);
            return true;
        }
    }
    return false;
}

void Gene::addTranscript(Transcript& transcript){
    this->_transcriptVector.push_back(transcript);
    this->updateStartEnd(transcript);
}
void Gene::updateStartEnd( Transcript& transcript  ){
    if( _start > transcript.getStart()  ){
       _start = transcript.getStart();
    } 
    if(_end < transcript.getEnd()){
        _end=transcript.getEnd();
    }
}
std::string Gene::getSource(){
    return this->_source;
}
void Gene::setSource(std::string& source){
    this->_source=source;
}
//gene end
//gff related class end
