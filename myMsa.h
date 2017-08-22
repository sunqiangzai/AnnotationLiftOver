/*
 * =====================================================================================
 *
 *       Filename:  myMsa.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  06/14/2017 09:49:02
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

#ifndef _MYMSA_H
#define _MYMSA_H
#include "model.h"
#include "parameters.h"
class MyMsa{
    private:
        std::vector<std::string> _names;
        std::vector<std::string> _seqs;
        std::vector<int> _starts;
        std::vector<int> _ends;
        std::vector<std::vector<Variant> > _variants;
        std::string chr;
    public:
        std::vector<std::string>& getNames();
        std::vector<std::string>& getSeqs();
        std::vector<int>& getStarts();
        std::vector<int>& getEnds();
        std::vector<std::vector<Variant> >& getVariants();
        std::string& getChr();

        void setNames( std::vector<std::string>& names);
        void setSeqs( std::vector<std::string>& seqs);
        void setStarts( std::vector<int>& starts);
        void setEnds( std::vector<int>& ends);
        void setVariants(std::vector<std::vector<Variant> >& variants);
        void setChr(std::string & chr);
};

std::vector<std::string>& myMsaAlignment(MyMsa& myMsa);
std::vector<std::string>& mafftMsaAlignment(MyMsa& myMsa);
#endif
