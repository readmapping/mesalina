/* 
 * File:   splicing.h
 * Author: mvyvermn
 *
 * Created on 20 november 2012, 18:11
 */

#ifndef SPLICING_H
#define	SPLICING_H

#include "mapper.h"
#include "types.h"
#include "dp.h"
#include <getopt.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>


using namespace std;

struct SpliceOptions_t{//commentary + sort + constructor
    SpliceOptions_t(){ initOptions(); }
    void initOptions(){
        query_threads   = 1;
        nucleotidesOnly = false;
        minMemLength    = 20;
        maxBranching    = 20;
        minCoverage     = 40;
        queryfile       = ref_fasta = outputName = "";
        indexLocation   = "";
        maxSpliceLength = 400000;
        printAln        = false;
    }
    int query_threads;
    //I/O options
    string queryfile;
    string ref_fasta;
    string indexLocation;// index prefix
    string outputName;
    //Sequence options
    bool nucleotidesOnly;//read only nucleotides (remove N's etc. from query)
    //alignment options
    int minMemLength;// minimum MEM length
    int maxBranching;//maximum number of MEMs for a single query position
    int minCoverage;//minimal query coverage for a cluster
    int maxSpliceLength;//max splice gap
    bool printAln; // print the alignment of each read
};

//i indexLocation
//x reference 
//U unpaired read sequence
//C minCoverage
//L minimum length
//n nucleotides only
//q query threads
//k maxbranching
//h help
//X maxSpliceLength
//A printAln
static const char * short_optionsS = "i:x:U:C:L:nq:k:hX:S:A";

static struct option long_optionsS[] = {
    {(char*)"threads",          required_argument, 0,            'q'},
    {(char*)"seedminlength",    required_argument, 0,            'L'},
    {(char*)"maxBranching",     required_argument, 0,            'k'},
    {(char*)"min-query-coverage",      required_argument, 0,     'C'},
    {(char*)"wildcards",        no_argument,       0,            'n'},
    {(char*)"help",             no_argument,       0,            'h'},
    {(char*)"index",            required_argument, 0,            'i'},
    {(char*)"max-splice",       required_argument, 0,            'X'},
    {(char*)"printAln",         no_argument,       0,            'A'},
    {(char*)0, 0, 0, 0} // terminator
};

static void usageSplice(const string prog) {
  cerr << "Usage: " << prog << " splice [options]" << endl;
  cerr << "the options for are ordered by functionality: " << endl;
  cerr << endl;
  cerr << "I/O OPTIONS " << endl;
  cerr << "-x (string)                reference sequence in multi-fasta" << endl;
  cerr << "-i/--index (string)        prefix or path of the index to load" << endl;
  cerr << "-U                         query file with unpaired reads (fasta or fastq)" << endl;
  cerr << "-S                         output file name (will be sam) [referenceName.sam]" << endl;
  cerr << endl;
  cerr << "PERFORMANCE OPTIONS " << endl;
  cerr << "-p/--threads (int)         number of threads [1]" << endl;
  cerr << endl;
  cerr << "ALIGNMENT OPTIONS " << endl;
  cerr << "-L/--seedminlength (int)   minimum length of the seeds used [20]" << endl;
  cerr << "-k/--maxBranching (int)    max number of seeds per query position [20]" << endl;
  cerr << "-C/--min-query-coverage (int)  minimum percent of bases of read the seeds have to cover [25]" << endl;
  cerr << "-X/--max-splice            maximum allowed size of an intron [400000]" << endl;
  cerr << "-n/--wildcards             treat Ns as wildcard characters" << endl;
  cerr << endl;
  cerr << "MISC OPTIONS " << endl;
  cerr << "-A/--printAln              print the alignment" <<endl;
  cerr << "-h/--help                  print this statement" << endl;
  exit(1);
}

static void processSpliceParameters(int argc, char* argv[], SpliceOptions_t& opt, const string program){
    cerr << "COMMAND: " << argv[1] << endl;
    cerr << "parsing options: ..." << endl;
    int option_index = 0;
    int c;
    while ((c = getopt_long(
    argc-1, argv+1,
    short_optionsS, long_optionsS, &option_index)) != -1) {//memType cannot be chosen
        switch (c) {
            case 'x': opt.ref_fasta = optarg; break;
            case 'U': opt.queryfile = optarg; break;
            case 'S': opt.outputName = optarg; break;
            case 'k': opt.maxBranching = atoi(optarg); break;
            case 'L': opt.minMemLength = atoi(optarg); break;
            case 'n': opt.nucleotidesOnly = 1; break;
            case 'q': opt.query_threads = atoi(optarg); break;
            case 'C': opt.minCoverage = atoi(optarg); break;
            case 'X': opt.maxSpliceLength = atoi(optarg); break;
            case 'h': usageSplice(program); break;
            case 'i': opt.indexLocation = optarg; break;
            case 'A': opt.printAln = true; break;
            case -1: /* Done with options. */
            break;
            case 0: if (long_options[option_index].flag != 0) break;
            default: usageSplice(program);
            throw 1;
        }
    }
    //add the reference query and output files
    if(opt.ref_fasta.empty() || opt.indexLocation.empty() || opt.queryfile.empty()){
        fprintf(stderr, "ALFALFA splice requires input reference file, index location and query file \n");
        exit(1);
    }
}

extern void spliceMap(const SparseSA& sa, Read_t& read, const SpliceOptions_t & spliceOptions, bool fwStrand);


// Fill the vector with the correct intron types
void getIntronTypes(IntronTypes &intronTypes, const string& ref,
                    const Boundaries &offset, const bool & fwStrand, const bool left);

// Structure to keep track of the possible maximum values in sandwich dp
struct Scores {

public:
    int scoreGTAG;
    int scoreGCAG;
    int scoreATAC;
    int scoreNON;
    
    int colGTAG;
    int colGCAG;
    int colATAC;
    int colNON;
    
    void resetScores() {
        this->scoreGTAG = MIN_INF;
        this->scoreGCAG = MIN_INF;
        this->scoreATAC = MIN_INF;
        this->scoreNON  = MIN_INF;
        
        this->colGTAG = 0;
        this->colGCAG = 0;
        this->colATAC = 0;
        this->colNON  = 0;
    }
    
    Scores(){resetScores();}
};


extern void processChains(const Chains& chains, const SparseSA& sa, 
            Read_t& read, bool fwStrand,  const SpliceOptions_t & spliceOptions);

int intronScore(const string& dinucLeft, const string& dinucRight, bool forward);

void sandwichDP(DynProg& dpLeft, DynProg& dpRight,        // The two dynamic programing objects
                Boundaries& offsetL, Boundaries& offsetR, // The boundary settings for both sides
                Alignment_t::Ptr alignment,  // The two output objects for both sides
                const string& ref, const Read_t &read, bool fwStrand, const Dp_type& type);

// Function used to extend the alignment to the left, used for the first seed
void extendLeft(DynProg& dp, Boundaries& offset,
                Alignment_t::Ptr alignment,
                const long chrStart, const int peelback, 
                const Match_t & match,
                const string& ref, const Read_t &read, bool fwStrand, const Dp_type& type);

// Function used to extend the alignment to the right, used for the last seed
void extendRight(DynProg& dp, Boundaries& offset,
                Alignment_t::Ptr alignment,
                const long chrEnd, const int peelback, 
                const Match_t & match,
                const string& ref, const Read_t &read, bool fwStrand, const Dp_type& type);

void getBestSplicePositionA(int& bestRowL, int& bestColL, int& bestRowR, int& bestColR, // Used to store the best splicesite positions
                           const DynProg& dpLeft, const DynProg& dpRight,               // The dp objects for both ends of the gap
                           const Boundaries& offsetL, const Boundaries& offsetR,        // The offsets of the alligned strings
                           const string& ref,                                           // The reference string, needed to determine the introntype
                           const bool& fwStrand);

void getBestSplicePositionB(int& bestRowL, int& bestColL, int& bestRowR, int& bestColR, // Used to store the best splicesite positions
                           const DynProg& dpLeft, const DynProg& dpRight,               // The dp objects for both ends of the gap
                           const Boundaries& offsetL, const Boundaries& offsetR,        // The offsets of the alligned strings
                           const string& ref,                                           // The reference string, needed to determine the introntype
                           const bool& fwStrand);

void getChainingA(Chains &chains, vector<Match_t> &matches, const long querylength, const SpliceOptions_t & spliceOptions);

#endif	/* SPLICING_H */

