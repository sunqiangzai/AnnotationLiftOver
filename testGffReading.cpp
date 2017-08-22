/*
 * =====================================================================================
 *
 *       Filename:  testGffReading.cpp
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  06/05/2017 00:39:34
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

#include<iostream>
using namespace std;

#include <iostream>
#include <seqan/sequence.h>
#include <seqan/stream.h>
#include <seqan/store.h>

using namespace seqan;

int main()
{
    CharString fileName = getAbsolutePath("/Users/song/TAIR10_GFF3_genes.gff");
    GffFileIn file(toCString(fileName));

    FragmentStore<> store;
    readRecords(store, file);
    // Create iterator
    Iterator<FragmentStore<>, AnnotationTree<> >::Type it;
    it = begin(store, AnnotationTree<>());
    // Iterate to the first annotation of type "exon"
    while (!atEnd(it) && getType(it) != "exon")
        goNext(it);
    // Output:
    std::cout << "  type: " << getType(it) << std::endl;
    std::cout << "  begin position: " << getAnnotation(it).beginPos << std::endl;
    std::cout << "  end position: " << getAnnotation(it).endPos << std::endl;
    std::cout << "  id: " << value(it) << std::endl;
    std::cout << "  parent id: " << getAnnotation(it).parentId << std::endl;
    std::cout << "  parent name: " << getParentName(it) << std::endl;
    return 0;
}