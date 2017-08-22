/*
 * =====================================================================================
 *
 *       Filename:  controlLayer.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  06/25/2017 09:38:17
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (), 
 *   Organization:  
 *
 * =====================================================================================
 */
/*************************************************************************




 ************************************************************************/

#ifndef _CONTROLLAYER_H
#define _CONTROLLAYER_H
int getGenomeSequence(int argc, char** argv);
int coordinateLiftOver( int argc, char** argv );
int gffCoordinateLiftOver( int argc, char** argv );
int getSequences(int argc, char** argv);
int annotationLiftOver( int argc, char** argv );
int annotationLiftOverAndOrth( int argc, char** argv);
int reAnnotationAndExonerateAndNovo( int argc, char** argv );
int myCountNumberOfTwoneighborSNP( int argc, char** argv );
int mycountNumberSNPAndIndel( int argc, char** argv );
int myGenerateRandomSdi( int argc, char** argv );
#endif
