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
#include "myutil.h"
using namespace std;
int main(){
    map<string, Fasta> sequences;
    Fasta f("name", "ATCG");
    cout << ">" << f.getName() << endl;
    cout << f.getSequence() << endl;
    std::string infile="/biodata/dep_tsiantis/grp_gan/song/rdINDELallHere/inputData/fullSdiFile/col_0.fa";
    readFastaFile(infile, sequences);
//    for( vector<Fasta>::size_type i=0 ; i< sequences.size(); i ++ ){
//        cout << sequences[i+0].getName() << " good" << endl;
//    }
    return 0;
}

