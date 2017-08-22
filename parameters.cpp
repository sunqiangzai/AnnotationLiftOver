//
// Created by baoxing on 6/12/17.
//

#include "parameters.h"

int _targetExtendLengthP = 100;
int _alignmentExonMatchP = 6;
int _alignmentExonMismatchP = -1;
int _alignmentExonOpenGapP = -4;
int _alignmentExonExtendGapP = -2;
int _alignmentIntronMatchP = 1;
int _alignmentIntronMismatchP = -1;
int _alignmentIntronOpenGapP = -2;
int _alignmentIntronExtendGapP = -1;
int _alignmentStartStopCodonMatchP = 10;
int _alignmentStartStopCodonMismatchP = -10;
int _alignmentStartStopCodonOpenGapP = -10;
int _alignmentStartStopCodonExtendGapP = -10;
int _alignmentSpliceSitesMatchP = 0;
int _alignmentSpliceSitesMismatchP = 0;
// There are a lot of splice donor and acceptor pairs. It not easy to check the pairs during alignment.
// Here we do not give penalty for splice sites, but give penalty to gap.
// If the exon regions could be well aligned, their should be not problem for splice sites.
// Then we could check the splice sites pairs for OFS state analysis.
int _alignmentSpliceSitesOpenGapP = -10;
int _alignmentSpliceSitesExtendGapP = -10;
int _myMsaMatchP = 5;
int _myMsaMisMatchP = -3;
int _myMsaGapP = -1;
int _myMsaOpenGapP = -3;
std::string _tempFolder = "./temp";



void setTargetExtendLength(int targetExtendLengthP){
    _targetExtendLengthP=targetExtendLengthP;
}
int getTargetExtendLength(){
    return _targetExtendLengthP;
}

void setAlignmentExonMatchP(int alignmentExonMatchP){
    _alignmentExonMatchP=alignmentExonMatchP;
}
int getAlignmentExonMatchP(){
    return _alignmentExonMatchP;
}

void setAlignmentExonMismatchP(int alignmentExonMismatchP){
    _alignmentExonMismatchP=alignmentExonMismatchP;
}
int getAlignmentExonMismatchP(){
    return _alignmentExonMismatchP;
}

void setAlignmentExonOpenGapintP (int alignmentExonOpenGapP){
    _alignmentExonOpenGapP=alignmentExonOpenGapP;
}
int getAlignmentExonOpenGapP(){
    return _alignmentExonOpenGapP;
}

void setAlignmentExonExtendGapintP (int alignmentExonExtendGapP){
    _alignmentExonExtendGapP=alignmentExonExtendGapP;
}
int getAlignmentExonExtendGapP(){
    return _alignmentExonExtendGapP;
}

void setAlignmentIntronMatchP(int alignmentIntronMatchP){
    _alignmentIntronMatchP=alignmentIntronMatchP;
}
int getAlignmentIntronMatchP(){
    return _alignmentIntronMatchP;
}

void setAlignmentIntronMismatchP(int alignmentIntronMismatchP){
    _alignmentIntronMismatchP=alignmentIntronMismatchP;
}
int getAlignmentIntronMismatchP(){
    return _alignmentIntronMismatchP;
}

void setAlignmentIntronOpenGapintP (int alignmentIntronOpenGapP){
    _alignmentIntronOpenGapP=alignmentIntronOpenGapP;
}
int getAlignmentIntronOpenGapP(){
    return _alignmentIntronOpenGapP;
}

void setAlignmentIntronExtendGapintP (int alignmentIntronExtendGapP){
    _alignmentIntronExtendGapP=alignmentIntronExtendGapP;
}
int getAlignmentIntronExtendGapP(){
    return _alignmentIntronExtendGapP;
}

void setAlignmentStartStopCodonMatchP(int alignmentStartStopCodonMatchP){
    _alignmentStartStopCodonMatchP=alignmentStartStopCodonMatchP;
}
int getAlignmentStartStopCodonMatchP(){
    return _alignmentStartStopCodonMatchP;
}

void setAlignmentStartStopCodonMismatchP(int alignmentStartStopCodonMismatchP){
    _alignmentStartStopCodonMismatchP=alignmentStartStopCodonMismatchP;
}
int getAlignmentStartStopCodonMismatchP(){
    return _alignmentStartStopCodonMismatchP;
}

void setAlignmentStartStopCodonOpenGapintP (int alignmentStartStopCodonOpenGapP){
    _alignmentStartStopCodonOpenGapP=alignmentStartStopCodonOpenGapP;
}
int getAlignmentStartStopCodonOpenGapP(){
    return _alignmentStartStopCodonOpenGapP;
}

void setAlignmentStartStopCodonExtendGapintP (int alignmentStartStopCodonExtendGapP){
    _alignmentStartStopCodonExtendGapP=alignmentStartStopCodonExtendGapP;
}
int getAlignmentStartStopCodonExtendGapP(){
    return _alignmentStartStopCodonExtendGapP;
}

void setAlignmentSpliceSitesMatchP(int alignmentSpliceSitesMatchP){
    _alignmentSpliceSitesMatchP=alignmentSpliceSitesMatchP;
}
int getAlignmentSpliceSitesMatchP(){
    return _alignmentSpliceSitesMatchP;
}

// There are a lot of splice donor and acceptor pairs. It not easy to check the pairs during alignment.
// Here we do not give penalty for splice sites, but give penalty to gap.
// If the exon regions could be well aligned, their should be not problem for splice sites.
// Then we could check the splice sites pairs for OFS state analysis.
void setAlignmentSpliceSitesMismatchP(int alignmentSpliceSitesMismatchP){
    _alignmentSpliceSitesMismatchP=alignmentSpliceSitesMismatchP;
}
int getAlignmentSpliceSitesMismatchP(){
    return _alignmentSpliceSitesMismatchP;
}

void setAlignmentSpliceSitesOpenGapintP (int alignmentSpliceSitesOpenGapP){
    _alignmentSpliceSitesOpenGapP=alignmentSpliceSitesOpenGapP;
}
int getAlignmentSpliceSitesOpenGapP(){
    return _alignmentSpliceSitesOpenGapP;
}

void setAlignmentSpliceSitesExtendGapintP (int alignmentSpliceSitesExtendGapP){
    _alignmentSpliceSitesExtendGapP=alignmentSpliceSitesExtendGapP;
}
int getAlignmentSpliceSitesExtendGapP(){
    return _alignmentSpliceSitesExtendGapP;
}


void setMyMsaMatchP (int myMsaMatchP){
    _myMsaMatchP=myMsaMatchP;
}
int getMyMsaMatchP(){
    return _myMsaMatchP;
}

void setMyMsaMisMatchP (int myMsaMisMatchP){
    _myMsaMisMatchP=myMsaMisMatchP;
}
int getMyMsaMisMatchP(){
    return _myMsaMisMatchP;
}

void setMyMsaGapP (int myMsaGapP){
    _myMsaGapP=myMsaGapP;
}
int getMyMsaGapP(){
    return _myMsaGapP;
}


void setMyMsaOpenGapP (int myMsaOpenGapP){
    _myMsaOpenGapP=myMsaOpenGapP;
}
int getMyMsaOpenGapP(){
    return _myMsaOpenGapP;
}


void setTempFolder(std::string& tempFolder ){
    _tempFolder=tempFolder;
}
std::string getTempFolder( ){
    return _tempFolder;
}

std::string createdTempFloder(){
    std::string tempFolder = getTempFolder();
    struct stat info;
    bool ifTempFolderCreated = false;
    while( !ifTempFolderCreated ){
        if( stat( &tempFolder[0], &info ) != 0 ){
            std::string command = "mkdir " + getTempFolder();
            system(&command[0]);
        } else if( info.st_mode & S_IFDIR ){  // S_ISDIR() doesn't exist on my windows
            ifTempFolderCreated=true;
        } else {
            std::cout << "could not access temp folder: " +  tempFolder << std::endl;
            exit(2);
        }
    }
    return _tempFolder;
}
std::string getTempFloder(){
    return _tempFolder;
}

std::string generateUUID(){
    return sole::uuid0().str();
}

std::string generateUUID(std::string& prefix){
    return prefix + "." + sole::uuid0().str();
}