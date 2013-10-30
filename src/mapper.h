/* 
 * File:   mapper.h
 * Author: mvyvermn
 *
 * Created on 2 augustus 2012, 16:18
 */

#ifndef MAPPER_H
#define	MAPPER_H

#include <string>
#include <assert.h>
#include <iostream>
#include <sstream>
#include <bitset>

#include "sparseSA.h"
#include "utils.h"
#include "options.h"
#include "alignment_t.h"
#include "read_t.h"

using namespace std;

// interval in match results + bases covering the result
struct Lis_t {
  Lis_t(): begin(0), end(0), len(0), fw(1), alignment(NULL), extended(0) {}
  Lis_t(vector<Match_t> * matches, int b, int e, int l, bool fw_): matches(matches), begin(b), end(e), len(l), fw(fw_), alignment(NULL), extended(0) {}
  int begin; // position in reference sequence
  int end; // position in query
  int len; // length of match
  bool fw;
  bool extended;
  vector<Match_t> * matches;
  Alignment_t::Ptr alignment;
};



extern void calculateSeeds(const SparseSA& sa, const string& P, int min_len, int alignmentCount, vector<Match_t>& matches, bool tryHarder);

extern void unpairedMatch(const SparseSA& sa, DynProg& dp_, Read_t & read,const align_opt & alnOptions, bool print);
extern void inexactMatch(const SparseSA& sa, DynProg& dp_, Read_t& read, const align_opt & alnOptions, bool fwStrand, bool print);
//TODO: calculate global position in above function

//PAIRED END FUNCTIONS
//extern bool isConcordant(const alignment_t& mate1, const alignment_t& mate2, const paired_opt& options);
extern void pairedMatch(const SparseSA& sa, DynProg& dp_, Read_t & mate1, Read_t & mate2, const align_opt & alnOptions, const paired_opt & pairedOpt, bool print);

#endif	/* MAPPER_H */

