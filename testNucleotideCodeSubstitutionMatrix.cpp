//
// Created by baoxing on 6/13/17.
//

#include "nucleotideCodeSubstitutionMatrix.h"
#include <iostream>
int main (){
    NucleotideCodeSubstitutionMatrix nucleotideCodeSubstitutionMatrix;

    if(nucleotideCodeSubstitutionMatrix.getPossibleStopCodons().find("TAA") == nucleotideCodeSubstitutionMatrix.getPossibleStopCodons().end()  ){
        std::cout << "10 " << nucleotideCodeSubstitutionMatrix.getPossibleStopCodons().size() << std::endl;
    }else{
        std::cout << "12" << std::endl;
    }
    if(nucleotideCodeSubstitutionMatrix.getPossibleStartCodons().find("ATG") == nucleotideCodeSubstitutionMatrix.getPossibleStartCodons().end()  ){
        std::cout << "16 " << nucleotideCodeSubstitutionMatrix.getPossibleStartCodons().size() << std::endl;
    }else{
        std::cout << "18" << std::endl;
    }
    return 0;
}