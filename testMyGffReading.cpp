/*
 * =====================================================================================
 *
 *       Filename:  testMyGffReading.cpp
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  06/05/2017 16:05:11
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (), 
 *   Organization:  
 *
 * =====================================================================================
 */
#include <stdlib.h>
/*************************************************************************




 ************************************************************************/

#include <iostream>
#include "model.h"
#include "myutil.h"
#include "Cds.h"
#include <map>
#include <set>
using namespace std;

int main()
{
    std::map<std::string, std::vector<Transcript> > variantsMap;
    readGffFile ("../../TAIR10_GFF3_genes.gff", variantsMap);
    std::cout << variantsMap.size() << std::endl;
    for( std::map<std::string, std::vector<Transcript> >::iterator it=variantsMap.begin(); it!=variantsMap.end(); it++  ){
    	std::cout << it->first << std::endl;
    	for( std::vector<Transcript>::iterator it2=it->second.begin(); it2!=it->second.end();it2++){
    		std::cout << (*it2).getName() << endl;
    	}
    }
    return 0;
}
