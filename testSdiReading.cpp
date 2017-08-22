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
#include <iostream>
#include "model.h"
#include <vector>
#include <sstream>
#include <fstream>
#include <regex>
#include <algorithm>
#include <map>
#include "myutil.h"
using namespace std;
int main(){
    map<string, vector<Variant> > sdiMaps;
    string inputFile = "./testData/PA10000.sdi";
    readSdiFile(inputFile, sdiMaps);
    map<string, vector<Variant> >::iterator  iter;
    for(iter = sdiMaps.begin(); iter != sdiMaps.end(); iter++){
        cout << iter ->first << endl;
        for(std::vector<Variant>::iterator vit = iter->second.begin(); vit != iter->second.end(); ++vit) {
            cout << "\t" << (*vit).getChromosome() << "\t" << (*vit).getPosition()
               << "\t" << (*vit).getChanginglength() << "\t"
              <<(*vit).getReference() << "\t" << (*vit).getAlternative() << endl;
//            print (cout, vit) << endl;
        }
    }
    return 0;
}

