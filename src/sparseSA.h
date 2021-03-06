/*
 * Copyright 2012, Michael Vyverman <michael.vyverman@ugent.be>
 *
 * This file is part of ALFALFA.
 * 
 * ALFALFA is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * ALFALFA is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with ALFALFA.  If not, see <http://www.gnu.org/licenses/>.
 * 
 * Part of the code originates from "A Practical Algorithm for Finding 
 * Maximal Exact Matches in Large Sequence Data Sets Using Sparse Suffix Arrays"
 * By Khan et al. Copyright (c) 2009 
 * You should have received a copy of the copyright notice with this code.
 */

#ifndef sparseSA_h__
#define sparseSA_h__

#include <vector>
#include <string>
#include <algorithm>
#include <limits>
#include <iostream>

using namespace std;

// Stores the LCP array in an unsigned char (0-255).  Values larger
// than or equal to 255 are stored in a sorted array.
// Simulates a vector<int> LCP;
struct Vec_uchar {
  struct Item_t{
    Item_t(){}
    Item_t(size_t i, int v) { idx = i; val = v; }
    size_t idx; int val;
    bool operator < (Item_t t) const { return idx < t.idx;  }
  };
  vector<unsigned char> vec;  // LCP values from 0-65534
  vector<Item_t> M;
  void resize(size_t N) { vec.resize(N); }
  // Vector X[i] notation to get LCP values.
  int operator[] (size_t idx) const{
    if(vec[idx] == numeric_limits<unsigned char>::max())
      return lower_bound(M.begin(), M.end(), Item_t(idx,0))->val;
    else
      return vec[idx];
  }
  // Actually set LCP values, distingushes large and small LCP
  // values.
  void set(size_t idx, int v) {
    if(v >= numeric_limits<unsigned char>::max()) {
      vec.at(idx) = numeric_limits<unsigned char>::max();
      M.push_back(Item_t(idx, v));
    }
    else { vec.at(idx) = (unsigned char)v; }
  }
  // Once all the values are set, call init. This will assure the
  // values >= 255 are sorted by index for fast retrieval.
  void init() { sort(M.begin(), M.end()); cerr << "M.size()=" << M.size() << endl; }
};

// Match find by findMEM.
struct Match_t {
  Match_t(): ref(0), query(0), len(0) {}
  Match_t(long r, long q, long l): ref(r), query(q), len(l) {}
  long ref; // position in reference sequence
  long query; // position in query
  long len; // length of match
};

// depth : [start...end]
struct Interval_t {
  Interval_t(): start(1), end(0), depth(-1) {}
  Interval_t(long s, long e, long d): start(s), end(e), depth(d) {}
  void reset(long e) { start = 0; end = e; depth = 0; }
  long depth, start, end;
  long size() { return end - start + 1; }
};

struct SparseSA {
  vector<string> &descr; // Descriptions of concatenated sequences.
  vector<long> &startpos; // Lengths of concatenated sequences.
  long maxdescrlen; // Maximum length of the sequence description, used for formatting.
  bool _4column; // Use 4 column output format.

  long N; //!< Length of the sequence.
  long logN; // ceil(log(N))
  long NKm1; // N/K - 1
  string &S; //!< Reference to sequence data.
  vector<unsigned int> SA;  // Suffix array.
  vector<int> ISA;  // Inverse suffix array.
  Vec_uchar LCP; // Simulates a vector<int> LCP.
  vector<int> CHILD; //child table

  long K; // suffix sampling, K = 1 every suffix, K = 2 every other suffix, K = 3, every 3rd suffix
  bool hasChild;
  bool hasSufLink;

  // Maps a hit in the concatenated sequence set to a position in that sequence.
  void from_set(long hit, long &seq, long &seqpos) const {
    // Use binary search to locate index of sequence and position
    // within sequence.
    vector<long>::iterator it = upper_bound(startpos.begin(), startpos.end(), hit);
    seq = distance(startpos.begin(), it) - 1;
    it--;
    seqpos = hit - *it;
  }

  // Get the index of the start and the end of the chromosome in the concatenated sequence,
  void getChromBounds(const long hit, long& chrStart, long& chrEnd) const {
      long seq, seqpos;
      this->from_set(hit, seq, seqpos);
        
      chrStart = this->startpos[seq];
      if(this->startpos.size() > 1 && seq < this->startpos.size() - 1)
          chrEnd = this->startpos[seq+1] - 1; //Minus one for the separation character
      else
          chrEnd =  S.length() - 1; //Minus one for the separation character
  }
  
  // Constructor builds sparse suffix array.
  SparseSA(string &S_, vector<string> &descr_, vector<long> &startpos_, bool __4column, long K_);

  // Modified Kasai et all for LCP computation.
  void computeLCP();
  //Modified Abouelhoda et all for CHILD Computation.
  void computeChild();

  // Radix sort required to construct transformed text for sparse SA construction.
  void radixStep(int *t_new, int *SA, long &bucketNr, long *BucketBegin, long l, long r, long h);

  // Prints match to cout.
  void print_match(const Match_t m) const;
  void print_match(const Match_t m, vector<Match_t> &buf) const; // buffered version
  void print_match(const string meta, vector<Match_t> &buf, bool rc) const; // buffered version

  // Binary search for left boundry of interval.
  inline long bsearch_left(char c, long i, long s, long e) const;
  // Binary search for right boundry of interval.
  inline long bsearch_right(char c, long i, long s, long e) const;

  // Simple suffix array search.
  inline bool search(const string &P, long &start, long &end) const;

  // Simple top down traversal of a suffix array.
  inline bool top_down(char c, long i, long &start, long &end) const;
  inline bool top_down_faster(char c, long i, long &start, long &end) const;
  inline bool top_down_child(char c, Interval_t &cur) const;

  // Traverse pattern P starting from a given prefix and interval
  // until mismatch or min_len characters reached.
  inline void traverse(const string &P, long prefix, Interval_t &cur, int min_len) const;
  inline void traverse_faster(const string &P,const long prefix, Interval_t &cur, int min_len) const;

  // Simulate a suffix link.
  inline bool suffixlink(Interval_t &m) const;

  // Expand ISA/LCP interval. Used to simulate suffix links.
  inline bool expand_link(Interval_t &link) const {
    long thresh = 2 * link.depth * logN, exp = 0; // Threshold link expansion.
    long start = link.start;
    long end = link.end;
    while(LCP[start] >= link.depth) {
      exp++;
      if(exp >= thresh) return false;
      start--;
    }
    while(end < NKm1 && LCP[end+1] >= link.depth) {
      exp++;
      if(exp >= thresh) return false;
      end++;
    }
    link.start = start; link.end = end;
    return true;
  }

  // Given a position i in S, finds a left maximal match of minimum
  // length within K steps.
  inline void find_Lmaximal(const string &P, long prefix, long i, long len, vector<Match_t> &matches, int min_len, bool print) const;

  // Given an interval where the given prefix is matched up to a
  // mismatch, find all MEMs up to a minimum match depth.
  void collectMEMs(const string &P, long prefix, const Interval_t mli,
  Interval_t xmi, vector<Match_t> &matches, int min_len, bool print) const;

  void collectSMAMs(const string &P, long prefix, const Interval_t mli,
  Interval_t xmi, vector<Match_t> &matches, int min_len, int maxCount, bool print) const;

  // Find all MEMs given a prefix pattern offset k.
  void findMEM(long k, const string &P, vector<Match_t> &matches, int min_len, bool print) const;

  // Find all MEMs given a prefix pattern offset k.
  void findSMAM(long k, const string &P, vector<Match_t> &matches, int min_len, int maxCount, bool print) const;

  // NOTE: min_len must be > 1
  void findMAM(const string &P, vector<Match_t> &matches, int min_len, bool print) const;
  inline bool is_leftmaximal(const string &P, long p1, long p2) const;

  // Maximal Almost-Unique Match (MAM). Match is unique in the indexed
  // sequence S. as computed by MUMmer version 2 by Salzberg
  // et. al. Note this is a "one-sided" query. It "streams" the query
  // P throught he index.  Consequently, repeats can occur in the
  // pattern P.
  void MAM(const string &P, vector<Match_t> &matches, int min_len, bool print) const{
    if(K != 1) return;  // Only valid for full suffix array.
    findMAM(P, matches, min_len, print);
  }

  // Find Maximal Exact Matches (MEMs)
  void MEM(const string &P, vector<Match_t> &matches, int min_len, bool print, int num_threads = 1) const;

  // Find Maximal Exact Matches (MEMs)
  void SMAM(const string &P, vector<Match_t> &matches, int min_len, int maxCount, bool print, int num_threads = 1) const;

  // Maximal Unique Match (MUM)
  void MUM(const string &P, vector<Match_t> &unique, int min_len, bool print) const;
  
  //save index to files
  void save(const string &prefix);
  
  //load index from file
  bool load(const string &prefix);
  
  //construct
  void construct();

};


#endif // __sparseSA_hpp__

