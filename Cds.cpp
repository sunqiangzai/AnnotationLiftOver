/*
 * =====================================================================================
 *
 *       Filename:  Cds.cpp
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  06/05/2017 15:47:39
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Baoxing Song (songbx.me), song@mpipz.mpg.de
 *   Organization:  MPIPZ
 *
 * =====================================================================================
 */
#include <stdlib.h>
/*************************************************************************




 ************************************************************************/

#include <iostream>
#include "Cds.h"
//cds begin
Cds::Cds(int start, int end){
    if( start < end ){
        _start=start;
        _end=end;
    }else{
        _end=start;
        _start=end;
    }
}
int Cds::getStart(){
    return _start;
}
int Cds::getEnd(){
    return _end;
}

// void Cds::setTranscript(Transcript& transcript){
//     this->_transcript = transcript;
// }
// Transcript* Cds::getTranscript(){
//     return _transcript;
// }
//cds end
