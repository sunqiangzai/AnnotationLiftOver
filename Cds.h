/*
 * =====================================================================================
 *
 *       Filename:  Cds.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  06/05/2017 15:47:30
 *       Revision:  none
 *       Compiler:  gcc
 *
 *        Author:  Baoxing Song (songbx.me), song@mpipz.mpg.de
 *   Organization:  MPIPZ
 *
 * =====================================================================================
 */
/*************************************************************************




 ************************************************************************/

#ifndef _CDS_H
#define _CDS_H
#include <set>
class Cds{
    private:
//        Transcript _transcript;
        int _start;
        int _end;
    public:
        Cds(int start, int end);
        int getStart();
        int getEnd();
        void temp();
		Cds()=delete;
       	//void setTranscript(Transcript& transcript);
       	//Transcript* getTranscript();
		
	bool operator<( const Cds& cds ) const{
	    if( _start < cds._start){
	        return true;
	    }else if( _start == cds._start && _end<cds._end ){
	        return true;
	    }
	    return false;
	}
	bool operator>(const Cds& cds )const {
	    if( _start > cds._start) {
	        return true;
	    }else if( _start == cds._start && _end > cds._end ){
	        return true;
	    }
	    return false;
	}
	bool operator==(const Cds& cds ) const{
	    if( _start==cds._start && _end==cds._start ){
	    	return true;
	    }
	    return false;
	}
	bool operator!=(const Cds& cds ) const{
	    if( _start==cds._start && _end==cds._start ){
	    	return false;
	    }
	    return true;
	}
};
#endif
