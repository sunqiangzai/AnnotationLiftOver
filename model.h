// =====================================================================================
// 
//       Filename:  model.h
// 
//    Description:  Here is the common model used for bioinformatics data
//    analysis, and operation function for those models
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
/*************************************************************************
The operation function for the specific model should be just below the class.
The functions crossing different classes should not defind in the service.h

For the specific class the default constructor should be defined clearly or forbiden

 ************************************************************************/

#ifndef _MODEL_H
#define _MODEL_H


#include <string>
#include <set>
#include "Cds.h"
#include <vector>
#include <map>
#include <limits>


enum STRAND
{
    POSITIVE, NEGATIVE
};

//Fasta class begin
class Fasta{
    private:
        std::string name;
        std::string sequence;
    public:
//        Fasta()=delete; //this is a new technique in C++11. It allows forbidding a function without depending on any other trick. It also clearly express your intention in the code.
//        Fasta()=default; //automatically generate a default constructor
        Fasta();
//        Fasta()= delete;
        Fasta(std::string& _name, std::string& _sequence);//When you declare any other constructor, the compiler will not generate the default constructor for you
        std::string getName() const {return name;}
        std::string getSequence() const {return sequence;}
        void setName(std::string& _name);
        void setSequence(std::string& _sequence);
        Fasta(const Fasta& f);//copy constructure
        std::string getSubsequence( size_t& start,  size_t& end);
        std::string getSubsequence( int& start,  int& end);
// virtual ~fasta() = 0;//TODO: to be improved
};
//this is a nomember function for class fasta, inputs are constant reference for
//memory saving, return the osream back, so than print (cout, f) << endl; could
//be used
std::ostream &print(std::ostream&, const Fasta&); 
//Fasta class end

// variant record class begin
class Variant{
    private:
        std::string chromosome;
        int position;
        int changingLength;
        std::string reference;
        std::string alternative;
        int lastTotalChanged;
    public:
        Variant()=delete;
//        Variant(std::string _chromosome, int _position, int _changingLength, std::string _reference, std::string _alternative);
        Variant(std::string& _chromosome, int& _position, std::string& _reference, std::string& _alternative);
        std::string getChromosome() const{return chromosome; }
        int getPosition() const {return position;}
        void setLastTotalChanged(int _lastTotalChanged){this->lastTotalChanged=_lastTotalChanged;}
        int getLastTotalChanged( ){return lastTotalChanged;}
        int getChanginglength() const {return changingLength;}
        std::string getReference() const {return reference;}
        std::string getAlternative() const {return alternative;}
        Variant(const Variant& variant);
        //bool overlap(const Variant& variant);
        bool overlap(int position);
};
std::ostream &print(std::ostream&, const Variant&);

struct compare_sdi_record
{
    inline bool operator() ( Variant& variant1, Variant& variant2 )
    {
        std::string chr1 = variant1.getChromosome();
        std::string chr2 = variant2.getChromosome();
        if( chr1.compare(chr2) ==0 ){
            if( variant1.getPosition() == variant2.getPosition() ){
                return 0-(variant1.getChanginglength() - variant2.getChanginglength());
            }else{
                return variant1.getPosition() < variant2.getPosition();
            }
        }else{
            return (chr1 < chr1);
        }
    }
};
//variant record class end

//gff related class begin
class Transcript{
    private:
        std::string _name;
        std::set<Cds> _cdsHashSet;
        std::vector<Cds> _cdsVector;
        std::string _chromeSomeName;
        STRAND _strand;
        int _start;
        int _end;
        std::string _metaInformation; //for ORF lift over
        std::string _geneomeSequence;
        std::string _cdsSequence;
        bool _ifOrfShift;
        std::string _source;
    public:
//        Transcript(std::string name, std::set<Cds> cdsHashSet, STRAND strand, string chromeSomeName);
        Transcript();
        Transcript(std::string name, std::string& chromeSomeName, STRAND strand);
        std::string getName();
        void setName(std::string name);
        std::set<Cds>& getCdsHashSet();
        std::vector<Cds>& getCdsVector();
        STRAND getStrand();
        std::string getChromeSomeName();
        void addCds(Cds& cds);
        int getStart();
        int getEnd();
        void updateInfor();
        void updateInfor(std::map<std::string, Fasta>& genome);
        std::string getGeneomeSequence();
        std::string getCdsSequence();
        void setMetaInformation(std::string& metaInformation);
        std::string& getMetaInformation();
        bool& getIfOrfShift();
        bool ifOverLap(Transcript& transcript);
        void setIfOrfShift(bool ifOrfShift);

        void cleanCdsVector();

        void cleanSaveRam(){
            _geneomeSequence.clear();
            _cdsSequence.clear();
        }
        std::string getSource();
        void setSource(std::string source);
        bool operator==( const Transcript& transcript ) const {
            if( _chromeSomeName != transcript._chromeSomeName ){
                return false;
            }else{
                if( _cdsHashSet.size() == transcript._cdsHashSet.size() ){
                    for (std::set<Cds>::iterator it = _cdsHashSet.begin(); it != _cdsHashSet.end(); ++it){
                        Cds i = *it;
                        if( transcript._cdsHashSet.find(i) == transcript._cdsHashSet.end()){
                            return false;
                        }
                    }
                    return true;
                }
                return false;
            }
        }
        bool operator!=( const Transcript& transcript ) const {
            if( _chromeSomeName != transcript._chromeSomeName ){
                return true;
            }else{
                if( _cdsHashSet.size() == transcript._cdsHashSet.size() ){
                    for (std::set<Cds>::iterator it = _cdsHashSet.begin(); it != _cdsHashSet.end(); ++it){
                        Cds i = *it;
                        if( transcript._cdsHashSet.find(i) == transcript._cdsHashSet.end()){
                            return true;
                        }
                    }
                    return false;
                }
                return true;
            }
        }

        bool operator<( const Transcript& transcript ) const {
            bool ifEqual = true;
            if( _chromeSomeName != transcript._chromeSomeName ){
                ifEqual =  false;
            }else{
                if( _cdsHashSet.size() == transcript._cdsHashSet.size() ){
                    for (std::set<Cds>::iterator it = _cdsHashSet.begin(); it != _cdsHashSet.end(); ++it){
                        Cds i = *it;
                        if( transcript._cdsHashSet.find(i) == transcript._cdsHashSet.end()){
                            ifEqual =  false;
                        }
                    }
                    ifEqual =  true;
                }
                ifEqual =  false;
            }
            if( ifEqual ){
                return false;
            }
            if( _chromeSomeName != transcript._chromeSomeName ){
                return _chromeSomeName < transcript._chromeSomeName;
            }else{
                if( _start != transcript._start ){
                    if( _start < transcript._start ){
                        return true;
                    }else{
                        return false;
                    }
                }else if ( _end != transcript._end ){
                    if(_end < transcript._end) {
                        return true;
                    }else{
                        return false;
                    }
                } else{
                    return true; //not perfecrt, but it should work for std::set<Transcripts>
                }
            }

        }
        bool operator>( const Transcript& transcript ) const{
            bool ifEqual = true;
            if( _chromeSomeName != transcript._chromeSomeName ){
                ifEqual =  false;
            }else{
                if( _cdsHashSet.size() == transcript._cdsHashSet.size() ){
                    for (std::set<Cds>::iterator it = _cdsHashSet.begin(); it != _cdsHashSet.end(); ++it){
                        Cds i = *it;
                        if( transcript._cdsHashSet.find(i) == transcript._cdsHashSet.end()){
                            ifEqual =  false;
                        }
                    }
                    ifEqual =  true;
                }
                ifEqual =  false;
            }
            if( ifEqual ){
                return false;
            }
            if( _chromeSomeName != transcript._chromeSomeName ){
                return _chromeSomeName > transcript._chromeSomeName;
            }else{
                if( _start != transcript._start ){
                    if( _start > transcript._start ){
                        return true;
                    }else{
                        return false;
                    }
                }else if ( _end != transcript._end ){
                    if(_end > transcript._end) {
                        return true;
                    }else{
                        return false;
                    }
                } else{
                    return false; //not perfecrt, but it should work for std::set<Transcripts>
                }
            }
        }
};

//gene begin
class Gene{
private:
    std::string _name;
    STRAND _strand;
    int _start;
    int _end;
    std::string _chromeSomeName;
    std::vector<Transcript> _transcriptVector;
    void updateStartEnd( Transcript& transcript);
    std::string _source;
public:
    Gene(std::string& name, STRAND& strand);
    Gene();
    std::string getName();
    STRAND getStrand();
    int getStart();
    int  getEnd();
    std::string getChromeSomeName();
    std::vector<Transcript>& getTranscriptVector();
    bool checkOverLapAndAddTranscript(Transcript& transcript);
    void addTranscript(Transcript& transcript);
    std::string getSource();
    void setSource(std::string& source);
};
struct compare_gene
{
    inline bool operator() ( Gene& gene1, Gene& gene2 )
    {
        std::string chr1 = gene1.getChromeSomeName();
        std::string chr2 = gene2.getChromeSomeName();
        if( gene1.getChromeSomeName().compare(gene2.getChromeSomeName()) ==0 ){
            return gene1.getStart() < gene2.getStart();
        }else{
            return (chr1 < chr1);
        }
    }
};
//gene end
//gff related class end

#endif
