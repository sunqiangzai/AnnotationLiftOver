//
// Created by baoxing on 6/12/17.
//

#ifndef VARIANTSANNOTATIONSOFTWARE_PARAMETERS_H
#define VARIANTSANNOTATIONSOFTWARE_PARAMETERS_H

#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <iostream>
#include "sole.h"

void setTargetExtendLength(int targetExtendLengthP);
int getTargetExtendLength();

void setAlignmentExonMatchP(int alignmentExonMatchP);
int getAlignmentExonMatchP();

void setAlignmentExonMismatchP(int alignmentExonMismatchP);
int getAlignmentExonMismatchP();

void setAlignmentExonOpenGapintP (int alignmentExonOpenGapP);
int getAlignmentExonOpenGapP();

void setAlignmentExonExtendGapintP (int alignmentExonExtendGapP);
int getAlignmentExonExtendGapP();

void setAlignmentIntronMatchP(int alignmentIntronMatchP);
int getAlignmentIntronMatchP();

void setAlignmentIntronMismatchP(int alignmentIntronMismatchP);
int getAlignmentIntronMismatchP();

void setAlignmentIntronOpenGapintP (int alignmentIntronOpenGapP);
int getAlignmentIntronOpenGapP();

void setAlignmentIntronExtendGapintP (int alignmentIntronExtendGapP);
int getAlignmentIntronExtendGapP();

void setAlignmentStartStopCodonMatchP(int alignmentStartStopCodonMatchP);
int getAlignmentStartStopCodonMatchP();

void setAlignmentStartStopCodonMismatchP(int alignmentStartStopCodonMismatchP);
int getAlignmentStartStopCodonMismatchP();

void setAlignmentStartStopCodonOpenGapintP (int alignmentStartStopCodonOpenGapP);
int getAlignmentStartStopCodonOpenGapP();

void setAlignmentStartStopCodonExtendGapintP (int alignmentStartStopCodonExtendGapP);
int getAlignmentStartStopCodonExtendGapP();

void setAlignmentSpliceSitesMatchP(int alignmentSpliceSitesMatchP);
int getAlignmentSpliceSitesMatchP();

void setAlignmentSpliceSitesMismatchP(int alignmentSpliceSitesMismatchP);
int getAlignmentSpliceSitesMismatchP();

void setAlignmentSpliceSitesOpenGapintP (int alignmentSpliceSitesOpenGapP);
int getAlignmentSpliceSitesOpenGapP();

void setAlignmentSpliceSitesExtendGapintP (int alignmentSpliceSitesExtendGapP);
int getAlignmentSpliceSitesExtendGapP();

void setMyMsaMatchP (int myMsaMatchP);
int getMyMsaMatchP();

void setMyMsaMisMatchP (int myMsaMisMatchP);
int getMyMsaMisMatchP();

void setMyMsaGapP (int myMsaGapP);
int getMyMsaGapP();

void setMyMsaOpenGapP (int myMsaOpenGapP);
int getMyMsaOpenGapP();

void setTempFolder(std::string& tempFolder );
std::string getTempFolder( );

std::string createdTempFloder();

std::string generateUUID();
std::string getTempFloder();
std::string generateUUID(std::string& prefix);
#endif