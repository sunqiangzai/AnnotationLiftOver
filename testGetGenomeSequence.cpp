// =====================================================================================
// 
//       Filename:  main.cpp
// 
//    Description:  
// 
//        Version:  1.0
//        Created:  04/12/2017 02:59:09 PM
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
#include "myfunctions.h"
#include "iostream"
int main(int argc, char** argv){
    if( argc == 4 ){
        std::cout << argc << " " << argv[0] << std::endl;
        getPseudoGenomeSequence(argv[1], argv[2], argv[3]);
    }else{
        std::cout << argc << " referenceGenomeFile sdiFile outPutFile" << std::endl;
    }
    return 0;
}

