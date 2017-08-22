//
// Created by baoxing on 6/12/17.
//
#include "reAnnotationAndMsa.h"

int main(  ){
    std::string referenceGenomeFilePath="testData/Col.fa";
    std::string referenceGffFilePath="./TAIR10_GFF3_genes.gff";
    std::map<std::string, std::string> sdiFilePaths;
    sdiFilePaths["PA10001"]="./testData/PA10001.sdi";
    sdiFilePaths["PA10000"]="./testData/PA10000.sdi";
    sdiFilePaths["Col"]="./testData/Col.sdi";
    reAnnotationAndMsa( referenceGenomeFilePath, referenceGffFilePath, sdiFilePaths);
    return 0;
}
