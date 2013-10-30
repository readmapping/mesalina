/* 
 * File:   read_t.h
 * Author: dieter
 *
 * Created on February 28, 2013, 11:08 AM
 */

#pragma once

#include <string>
#include <assert.h>
#include <iostream>
#include <sstream>

#include "sparseSA.h"
#include "dp.h"
#include "alignment_t.h"


using namespace std;

struct Read_t {
    string qname;//TODO should be reference
    vector< Alignment_t::Ptr > alignments;
    string sequence;//TODO should be reference
    string rcSequence; // Reverse Complement
    string qual;//TODO should be reference
    string rQual;
    int    pairedAlignmentCount;
    
    Read_t();
    Read_t( string name, string &seq, string &qw, bool nucleotidesOnly );
    ~Read_t();
    
    void init( bool nucleotidesOnly );
    void postprocess( const Dp_scores & scores, const SparseSA& sa );
    string emptyAlingment( bool paired, bool mateFailed, bool upstream );
    int alignmentCount();
    void removeLastAlignment();
    string printUnpairedAlignment(int i);
    string printUnpairedAlignmentSplice(int i);
    string printPairedAlignments(int i);
    
    void addAlignment(Alignment_t::Ptr alignment);
    void showAlignment(const size_t i, const size_t wraplength = 50) const;
    
};
