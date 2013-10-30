/* 
 * File:   alignment_t.h
 * Author: dieter
 *
 * Created on February 28, 2013, 10:59 AM
 */

#pragma once

#include <string>
#include <assert.h>
#include <sstream>
#include <bitset>
#include <memory>

#include "sparseSA.h"
#include "utils.h"
#include "options.h"
#include "mate_t.h"


using namespace std;

struct Alignment_t {
    
public:
    typedef std::shared_ptr<Alignment_t> Ptr;
    typedef vector<char> CigarChars;
    typedef vector<int>  CigarLengths;
    
    string cigar; //TODO: remove this fields, only used when printed
    string NMtag; //TODO: remove this fields, only used when printed
    CigarChars   cigarChars; //Change these to fixed-length values
    CigarLengths cigarLengths; //Change these to fixed-length values (make sure to increase them when necessary): make own string and vector classes
    string rname; //leave out, only for printing
    long globPos; // position in index (concat of all reference sequences)
    long pos; // position in reference sequence
    bitset<11> flag;
    int mapq;
    int editDist;
    int alignmentScore;
    int refLength; //length in reference sequence spanned by this alignment
    //paired-end fields
    vector<Mate_t> mateInfo;
    bool fwStrand; // forward or reverse strand 
    
    Alignment_t();
    Alignment_t(const Alignment_t & o);
    ~Alignment_t() {}
    
    //functions

    bool paired();
    bool concordant() const;
    void setLocalPos(const SparseSA& sa);
    int pairedCount();
    void createCigar(bool newVersion);
    void setFieldsFromCigar(const Dp_scores & scores);
    void addMate(const Alignment_t::Ptr o, bool concordant, bool upstream);
    
    // Add different parts to the cigar string
    void addSymbol(const char symbol, const int length);    
    void addMatch(const int length);
    void addMisMatch(const int length);
    void addDeletion(const int length);
    void addInsertion(const int length);
    void addSkippedRegion(const int length);
    
    void deleteLastSymbol();
    
    // Add the cigar string from alignment at the end in the same order
    void addAlignmentForward(Alignment_t & alignment);
    // Add the cigar string from alignment at the end in the reversed order
    void addAlignmentReversed(Alignment_t & alignment);
    
    // Add the string of alignment characters to the alignment
    void addStringForward(string & s);
    // Add the string of alignment characters to the alignment in reverse order
    void addStringReversed(string & s);
  
    // Check if the cigarChars are correct, e.g. combine the same adjacent characters
    void postProcess();
    
    void print(const string &ref, const string &fwQuery,const string &rcQuery, const size_t wraplength);
};

