/**
 *
 *
 */



#include "splicing.h"
#include "utils.h"


//#define DEBUG
#ifdef DEBUG
#define debug(x) x
#else
#define debug(x) ((void)0)
#endif

//#define DEBUG_DP
#ifdef DEBUG_DP
#define debug_dp(x) x
#else
#define debug_dp(x) ((void)0)
#endif

// for left and right extention of the alignment
//#define DEBUG_DP2
#ifdef DEBUG_DP2
#define debug_dp2(x) x
#else
#define debug_dp2(x) ((void)0)
#endif

//#define DEBUG_CHAIN
#ifdef DEBUG_CHAIN
#define debug_ch(x) x
#else
#define debug_ch(x) ((void)0)
#endif

//#define DEBUG_SAM
#ifdef DEBUG_SAM
#define debug_sam(x) x
#else
#define debug_sam(x) ((void)0)
#endif


//#define DEBUG_SEEDS
#ifdef DEBUG_SEEDS
#define debug_seeds(x) x
#else
#define debug_seeds(x) ((void)0)
#endif

using namespace std;

void spliceMap(const SparseSA& sa, Read_t& read, const SpliceOptions_t & spliceOptions, bool fwStrand) {
    string* P = &read.sequence;
    long Plength = (long) P->length();
    if (!fwStrand)
        P = &read.rcSequence;
    int min_len = spliceOptions.minMemLength;
    vector<Match_t> matches; // will contain all the seeds
    //calc seeds and sort according to genome position
    calculateSeeds(sa, *P, min_len, spliceOptions.maxBranching, matches, false);

#ifdef DEBUG_SEEDS
    cerr << endl << "Seeds before any processing: " << endl;
    cerr << "Number of seeds: " << matches.size() << endl;
    for (vector<Match_t>::iterator match = matches.begin(); match != matches.end(); match++)
        cerr << match->ref << "\t" << match->query << "\t" << match->len << endl;
    cerr << endl;
#endif

    
    ////////////////////////////////////////////////////////////////////////////
    //  SEED CHAINING
    ////////////////////////////////////////////////////////////////////////////

   
    Chains chains; // Each cluster will contain a selected set of seeds, which will be extended to the query
    getChainingA(chains, matches, Plength, spliceOptions);
    
    string strand = fwStrand ? "forward" : "reverse";
#ifdef DEBUG_SEEDS
    cerr << "clusters found for query " << read.qname << " on the " << strand << " strand." << endl;
    cerr << "number of clusters: " << chains.size() << endl;
    for (int idx = 0; idx < chains.size(); idx++) {
        cerr << "cluster " << (idx + 1) << endl;
        for (int idx2 = 0; idx2 < chains[idx].size(); idx2++) {
            cerr << chains[idx][idx2].ref << "\t" << chains[idx][idx2].query << "\t" << chains[idx][idx2].len << endl;
        }
        cerr << endl;
    }
#endif    
    ////////////////////////////////////////////////////////////////////////////
    // PROCESS THE CHAINS
    ////////////////////////////////////////////////////////////////////////////

    processChains(chains, sa, read, fwStrand, spliceOptions);

}

#define MIN_INTRON_LENGTH 9

void processChains(const Chains& chains, const SparseSA& sa, Read_t& read, bool fwStrand,  const SpliceOptions_t & spliceOptions) {

//    cerr << "queryLength: " << read.sequence.size() << endl;

    long chrStart, chrEnd;

   
    // TODO: should be initialized somewhere else
    int peelback = 1; // This number of nucleotides will be used to peel back around the gap
//    int min_len = spliceOptions.minMemLength;
//    
//    if (peelback < 2*min_len)
//        peelback = min_len/2;
    
    long queryGap, refGap; //The length of the gaps between the seeds
    Dp_scores scores(3, -3, -10, -3); //Values are taken from GMAP
    const Dp_type type;

    DynProg dpLeft(255, true, scores);
    DynProg dpRight(255, true, scores);
    Boundaries offsetL;
    Boundaries offsetR;
    Dp_output outputL;
    Dp_output outputR;

    bool bad = false;
    
    // Loop over the different possible chains
    for (Chains::const_iterator chain = chains.cbegin(); chain != chains.cend(); chain++) {

        bad = false;
//        Alignment_t * alignment = new Alignment_t();
        Alignment_t::Ptr alignment ( new Alignment_t() );
        if(!fwStrand) {
            alignment->flag.set(4,true);
            alignment->fwStrand = false;
        }

        Chain::const_iterator seed = chain->cbegin();

        // Get the boundaries of the current chromosome in the concatenated sequence
        sa.getChromBounds(seed->ref, chrStart, chrEnd);

        // There is only one seed in this chain
        if (chain->size() == 1) {
            extendLeft(dpLeft, offsetL, alignment,
                   chrStart, peelback,
                    *seed,
                    sa.S, read, fwStrand, type);
            
            extendRight(dpRight, offsetR,
                alignment,
                chrEnd, peelback,
                *seed,
                sa.S, read, fwStrand, type);
        }

        else {
        ////////////////////////////////////////////////////////////////////////
        // Process the first seed
        ////////////////////////////////////////////////////////////////////////

        extendLeft(dpLeft, offsetL, alignment,
                   chrStart, peelback,
                   *seed,
                   sa.S, read, fwStrand, type);
 
        ////////////////////////////////////////////////////////////////////////
        // Process the sequential seeds
        ////////////////////////////////////////////////////////////////////////

        for (seed; seed != chain->end() - 1 ; seed++) {
            // Calculate the gap sizes between this and the next seed
            queryGap = (seed + 1)->query - seed->query - seed->len;
            refGap   = (seed + 1)->ref - seed->ref - seed->len;

            debug(cerr << endl << "querygap: " << queryGap << " refgap: " << refGap << endl);

            // TODO: something better !!!
            if (queryGap - refGap > MIN_INTRON_LENGTH){
                bad = true;
                break;
            }
            // Check for the intron, using sandwich dp
            else if (refGap - queryGap > MIN_INTRON_LENGTH) {
                debug_dp(cerr << "Sandwich DP" << endl);
                // For now use the same length of the genomic segment as the length of the query part
                offsetL.setValues( seed->ref + seed->len - peelback,             /*refBegin*/
//                                   seed->ref + seed->len + 1.2*queryGap + peelback,  /*refEnd*/
                                   seed->ref + seed->len + queryGap + peelback,  /*refEnd*/
                                   seed->query + seed->len - peelback ,           /*queryBegin*/
                                   (seed + 1)->query + peelback                  /*queryEnd*/);

                // use the same query segment for the right part
                offsetR.setValues( (seed + 1)->ref - queryGap - peelback ,         /*refBegin*/
                                   (seed + 1)->ref + peelback,                    /*refEnd*/
                                   seed->query + seed->len - peelback ,            /*queryBegin*/
                                   (seed + 1)->query + peelback                   /*queryEnd*/);

                // The alignment is updated inside this function
                sandwichDP(dpLeft, dpRight, offsetL, offsetR, alignment, sa.S, read, fwStrand, type);

                // add the seed as match
                alignment->addMatch((seed+1)->len - 2*peelback);
                
            }
            else if (refGap == 1 && queryGap == 1) {
                debug_dp(cerr << "1 mismatch" << endl);
                // Add the peelbackmatch from the previous seed
                alignment->addMatch(peelback);
                // Add this one mismatch
                alignment->addMisMatch(1);
                alignment->addMatch((seed+1)->len - peelback);
            } else {
                debug_dp(cerr << "Normal DP" << endl);
                // Normal dynamic programming
                int bandsize = 10;
                bandsize = queryGap/4 > bandsize ? queryGap/2 : bandsize ;
                              
                long refstrLB = seed->ref + seed->len - peelback-1L;
                long refstrRB = (seed + 1)-> ref + peelback;
                long queryLB  = seed->query + seed->len - peelback-1L;
                long queryRB  = (seed+1)-> query + peelback;
                Boundaries offset (refstrLB + 1L, refstrRB - 1L, queryLB + 1, queryRB - 1);
                Dp_type types;
                Dp_output output;
                dpLeft.dpBandStatic( sa.S, fwStrand ? read.sequence : read.rcSequence, offset, types, ERRORSTRING, output, bandsize, false);

                

#ifdef DEBUG_DP
                cerr << "alignment of normal DP:" <<endl;
                for(size_t i = 0 ; i < output.cigarChars.size(); i++)
                    cerr << output.cigarChars[i] << ", " << output.cigarLengths[i] <<endl;
                cerr << "end alignment" <<endl;
#endif
                
                alignment->cigarChars.insert(alignment->cigarChars.end(), output.cigarChars.begin(), output.cigarChars.end());
                alignment->cigarLengths.insert(alignment->cigarLengths.end(), output.cigarLengths.begin(), output.cigarLengths.end());
                alignment->addMatch((seed+1)->len - 2*peelback);
                alignment->alignmentScore += dpLeft.scores.match*(seed->len -2*peelback) + output.dpScore;
//                curEditDist += output.editDist;
                output.clear();
            }
           //}
//            for(size_t i = 0; i < alignment->cigarChars.size(); i++)
//                cout << alignment->cigarChars[i] << ", " << alignment->cigarLengths[i]<<endl;
//            cout <<endl;
        }

        ////////////////////////////////////////////////////////////////////////
        // Process the last seed from this chain
        ////////////////////////////////////////////////////////////////////////

        extendRight(dpRight, offsetR,
                alignment,
                chrEnd, peelback,
                *seed,
                sa.S, read, fwStrand, type);
        } // chain contains more than one seed

        
        // TODO: Add some checks whether this is a good alignment
        if( !bad ) {

            alignment->postProcess();
            alignment->setFieldsFromCigar(scores);
            read.addAlignment(alignment);
            
            if(spliceOptions.printAln) {
                cout << endl << "Alignment for query " << read.qname << ":"<< endl;
                alignment->print(sa.S, read.sequence, read.rcSequence, 60);
            }
                  
#ifdef DEBUG
            for(size_t i = 0; i < alignment->cigarChars.size() ; i++)
                cerr << alignment->cigarChars.at(i) << ", " << alignment->cigarLengths.at(i) <<endl;
#endif

            read.postprocess(scores, sa);
            read.printUnpairedAlignment(0);

        }

    } // Loop over different chains
}




// Function used to extend the alignment to the left, used for the first seed
void extendLeft(DynProg& dp, Boundaries& offset,
                Alignment_t::Ptr alignment,
                const long chrStart, const int peelback,
                const Match_t & seed,
                const string& ref, const Read_t &read, bool fwStrand, const Dp_type& type){

    long editDist = peelback;

    long refstrLB = seed.ref;
    long refstrRB = refstrLB + seed.len -1L;
    long queryLB  = seed.query;
    long queryRB  = queryLB + seed.len-1L;
    bool clipping = false;

    long curEditDist = 0;
    Dp_output output;

    if(seed.query > 0 && seed.ref > chrStart) {
        debug_dp2(cerr << "Extending left 1" <<endl);
        
        // Try to extend the alignment to the left, without looking for introns
        long alignmentBoundLeft = max(refstrLB-queryLB-editDist,chrStart);
        offset.setValues(alignmentBoundLeft, refstrLB-1, 0, queryLB-1);
        const string * query = fwStrand ? &read.sequence : &read.rcSequence;

        Dp_type types;
        types.freeRefB = true;
        types.freeQueryB = clipping;
        dp.banded    = true;
        dp.bandLeft  = editDist-curEditDist;
        dp.bandRight = editDist-curEditDist;
        dp.bandSize  = 0;
        dp.updateMatrix(type);
        dp.dpBandStatic( ref, fwStrand ? read.sequence : read.rcSequence, offset, types, ERRORSTRING, output, editDist-curEditDist, false);


//        if(output.cigarChars[output.cigarChars[0]] == 'I')
//                queryLB += (long)output.cigarLengths[0];
//        if(clipping && offset.queryB> 0){
//            alignment->cigarChars.push_back('S');
//            alignment->cigarLengths.push_back(offset.queryB);
//        }

        alignment->cigarChars.insert(alignment->cigarChars.end(),output.cigarChars.begin(),output.cigarChars.end());
        alignment->cigarLengths.insert(alignment->cigarLengths.end(),output.cigarLengths.begin(),output.cigarLengths.end());
        alignment->globPos = offset.refB+1L;
    //            alignment->alignmentScore += output.dpScore;
    ////            curEditDist += output.editDist;
        output.clear();

        // Add the rest of the seed as match
        if (seed.len - peelback > 0)
            alignment->addMatch(seed.len - peelback); // one peelback at the end of the seed
    }

    else if (seed.query > 0) { // The seed starts at the beginning of the chromosome
        debug_dp2(cerr << "Extending left 2" << endl);
        // Add the first nucleotides of the query as insertions
        alignment->addInsertion(seed.query);
        //Add the the rest of the seed as match
        if ( seed.len - peelback - seed.query > 0 )
            alignment->addMatch(seed.len - peelback - seed.query);
        alignment->globPos = seed.ref + 1L;
    }
    else { // We can add the seed as a full match
        debug_dp2(cerr << "Extending left 3" <<endl);
        if ( seed.len - peelback > 0 )
            alignment->addMatch(seed.len - peelback);
        alignment->globPos = seed.ref + 1L;
    }
}

// Function used to extend the alignment to the right, used for the last seed
void extendRight(DynProg& dp, Boundaries& offset,
                Alignment_t::Ptr alignment,
                const long chrEnd, const int peelback,
                const Match_t & seed,
                const string& ref, const Read_t &read, bool fwStrand, const Dp_type& type){

    
if(seed.query + seed.len < read.sequence.length() && seed.ref + seed.len < chrEnd ) {
    debug_dp2(cerr << "rightext: 1" <<endl);
    // Add the end of the seed as a match
    alignment->addMatch(peelback);

    // Try to extend the alignment to the right, without looking for introns

    const string * query = fwStrand ? &read.sequence : &read.rcSequence;
    Dp_output output;
    long refstrRB = seed.ref + seed.len - 1L;
    long queryRB = seed.query + seed.len - 1L;
    long curEditDist = 0;
    long editDist = peelback;
    int bandsize = 10;
    long refAlignRB = refstrRB + ((long)read.sequence.length() - queryRB) + 1L + editDist - curEditDist;
    bool clipping = false;
        if (refAlignRB >= chrEnd)
            refAlignRB = chrEnd - 1L;
        offset.setValues(refstrRB + 1L, refAlignRB, queryRB + 1, query->length() - 1);
        Dp_type types;
        types.freeRefE = true;
        types.freeQueryE = clipping;


        dp.dpBandStatic( ref, *query, offset, types, ERRORSTRING, output, bandsize, false);
        if(offset.queryE > queryRB){
            int addToLength = offset.queryE-queryRB;
            if(output.cigarChars[output.cigarChars.size()-1] == 'I')
                addToLength -= output.cigarLengths[output.cigarLengths.size()-1];
            alignment->refLength += addToLength;
        }
        alignment->cigarChars.insert(alignment->cigarChars.end(), output.cigarChars.begin(), output.cigarChars.end());
        alignment->cigarLengths.insert(alignment->cigarLengths.end(), output.cigarLengths.begin(), output.cigarLengths.end());
        alignment->alignmentScore += output.dpScore;
        if (clipping && offset.queryE < query->length() - 1) {
            alignment->cigarChars.push_back('S');
            alignment->cigarLengths.push_back(query->length() - 1 - offset.queryE);
        }
        curEditDist += output.editDist;
        output.clear();

}
else if(seed.query + seed.len < read.sequence.length()) {
    debug_dp2(cerr<< "rightext: 2" <<endl);
    //Add the length of the seed as full match to the alignment
    alignment->addMatch(peelback);
    // Add the last nucleotides of the query as insertions
    if ( read.sequence.length() - seed.query - seed.len > 0 )
        alignment->addInsertion(read.sequence.length() - seed.query - seed.len);
}
else {
    debug_dp2(cerr << "rightext: 3" <<endl);
    // Add the seed as a full match
    alignment->addMatch(peelback);

}
}

void sandwichDP(DynProg& dpLeft, DynProg& dpRight,
        Boundaries& offsetL, Boundaries& offsetR,
        Alignment_t::Ptr alignment,
        const string& ref, const Read_t &read, bool fwStrand, const Dp_type& type) {
    debug_dp(cerr << endl << "Begin Sandwich DP\n=======================" << endl);

    size_t qlength  = offsetL.queryE - offsetL.queryB;
    size_t rlengthL = offsetL.refE - offsetL.refB;
    size_t rlengthR = offsetR.refE - offsetR.refB;

    const string * const query = &(fwStrand ? read.sequence : read.rcSequence);

   
    debug_dp(cerr << "Left\n");
    debug_dp(cerr << "qB: " <<offsetL.queryB << "\tqE: " << offsetL.queryE<<endl);
    debug_dp(cerr << "rB: " <<offsetL.refB << "\trE: " << offsetL.refE<<endl);
    debug_dp(cerr << "Right\n");
    debug_dp(cerr << "qB: " <<offsetR.queryB << "\tqE: " << offsetR.queryE<<endl);
    debug_dp(cerr << "rB: " <<offsetR.refB << "\trE: " << offsetR.refE<<endl);

    ////////////////////////////////////////////////////////////////////////////
    // Compute the dp matrices for the left and right end of the gap
    ////////////////////////////////////////////////////////////////////////////

    dpLeft.banded = true;
    dpLeft.bandLeft = qlength / 2;
    dpLeft.bandRight = qlength / 2;
    dpLeft.L1 = rlengthL;
    dpLeft.L2 = qlength;
    dpLeft.bandSize = qlength;

    dpLeft.updateMatrix(type);
    dpLeft.dpFillStatic(ref, *query, true, offsetL, type, true);
    //debug_dp(dpLeft.print_matrices(ref, *query, offsetL, false, true));


    dpRight.banded = true;
    dpRight.bandLeft = qlength / 2;
    dpRight.bandRight = qlength / 2;
    dpRight.L1 = rlengthR;
    dpRight.L2 = qlength;
    dpRight.bandSize = qlength;

    dpRight.updateMatrix(type);
    dpRight.dpFillStatic(ref, *query, true, offsetR, type, false);
    //debug_dp(dpRight.print_matrices(ref, *query, offsetR, false, false));


//    cerr << "dp.. .L2: " << dpLeft.L2 << endl;

    
    debug_dp(cerr << "qlength: " << qlength << "\ndp.. .L2: " << dpLeft.L2 << endl);
    debug_dp(cerr << "dpleft.L1: " << dpLeft.L1 << "\ndpRight.L1: " << dpRight.L1 << endl);

    ////////////////////////////////////////////////////////////////////////////
    // Get the best splicesite positions
    ////////////////////////////////////////////////////////////////////////////
    int bestRowL = 0;
    int bestColL = 0;
    int bestRowR = 0;
    int bestColR = 0;

    // These vectors will contain the introntype of the nucleotides following
    // this position.
    // The introntypes for the right will be ordered like the dp matrix
//    IntronTypes intronTypesLeft (rlengthL - 2);
//    IntronTypes intronTypesRight (rlengthR - 2);
//
//    getIntronTypes(intronTypesLeft, ref, offsetL, fwStrand, true);
//    getIntronTypes(intronTypesRight, ref, offsetR, fwStrand, false);

    getBestSplicePositionA(bestRowL, bestColL, bestRowR, bestColR,
                          dpLeft, dpRight,
                          offsetL, offsetR,
                          ref,
                          fwStrand);
//    cerr << "After:" <<endl;
//cerr << "bestRowL: "<< bestRowL << " bestColL: " << bestColL <<endl;
//cerr << "bestRowR: "<< bestRowR << " bestColR: " << bestColR <<endl;
    debug_dp(cerr << "bestRowL: "<< bestRowL << " bestColL: " << bestColL <<endl);
    debug_dp(cerr << "bestRowR: "<< bestRowR << " bestColR: " << bestColR <<endl);

    ////////////////////////////////////////////////////////////////////////////
    // Traceback to get the alignment
    ////////////////////////////////////////////////////////////////////////////

    Dp_output outputL;
    Dp_output outputR;

    stringstream ssL;
    int i = bestRowL;
    int j = bestColL;
    dpLeft.dpTraceBackStatic(i, j, type, outputL,
            ssL, offsetL, true, ref, *query, true);

    stringstream ssR;
    i = bestRowR;
    j = bestColR;
    dpRight.dpTraceBackStatic(i, j, type, outputR,
            ssR, offsetR, true, ref, *query, false);
    

#ifdef DEBUG_DP
    cerr <<endl;
    cerr << query->substr(offsetL.queryB, bestRowL) <<"       "<< query->substr(offsetR.queryE-bestRowR, bestRowR) <<endl;
    cerr <<ssL.str()<< "       "  << ssR.str() <<endl;
    cerr <<ref.substr(offsetL.refB, bestColL+2)<<"..." << ref.substr(offsetR.refE-bestColR-2,  bestColR+2) <<endl<<endl;
#endif
    

    ////////////////////////////////////////////////////////////////////////////
    // Output
    ////////////////////////////////////////////////////////////////////////////

    // Add the end of the left exon
    string s = ssL.str();
    alignment->addStringReversed(s);

    // Add the intron
    alignment->addSkippedRegion(offsetR.refE - bestColR - (offsetL.refB + bestColL));

    // Add the beginning of the right exon
    // The stringstream is alrdeady in the correct order !
    s = ssR.str();
    alignment->addStringForward( s );


#ifdef DEBUG
    cerr << "=== Splicesite ===" << endl;
    for(size_t i = 0 ; i < outputL.cigarChars.size(); i++)
        cerr << outputL.cigarChars[i] << outputL.cigarLengths[i];
    cerr << endl;
    if (fwStrand) {
        cerr << "splicesite in query: " << bestRowL + offsetL.queryB << " on forward strand" << endl;
        cerr << "intron in ref: " << offsetL.refB + bestColL + 1 << " .. " << offsetR.refE - bestColR << endl;
    } else {
        cerr << "splicesite in query: " << bestRowL + offsetL.queryB << " on reverse strand" << endl;
        cerr << "equivalent position in forward query: " << read.sequence.length() - bestRowL - offsetL.queryB << endl;
        cerr << "intron in ref: - " << offsetR.refE - bestColR + 1 + 1 << " .. " << offsetL.refB + bestColL - 1 << endl;
    }
    cerr << "intronlength:  " << offsetR.refE - bestColR - (offsetL.refB + bestColL) << endl;
    cerr << "type:          " << ref.substr(offsetL.refB + bestColL, 2) << "-" << ref.substr(offsetR.refE - bestColR - 2, 2) << endl;
    cerr << "==================" << endl;
#endif
    debug_dp(cerr << endl << "end Sandwich DP\n=======================" << endl << endl);
}

// TODO
void getIntronTypes(IntronTypes &intronTypes, const string& ref,
                    const Boundaries &offset, const bool & fwStrand, const bool left)
{
 /*
    Introntypes:
        Forward    Reverse
        GT-AG       CT-AC
        GC-AG       CT-GC
        AT-AC       GT-AT
          -           -
  */
    for(size_t i = 0; i < intronTypes.size(); i++ ){
        if(left && fwStrand)
        {
            if(ref[offset.refB + i + 1] == 'G' && ref[offset.refB + i + 2] == 'T')
                intronTypes[i] = INTRON_GTAG;
            else if(ref[offset.refB + i + 1] == 'G' && ref[offset.refB + i + 2] == 'C')
                intronTypes[i] = INTRON_GCAG;
            else if(ref[offset.refB + i + 1] == 'A' && ref[offset.refB + i + 2] == 'T')
                intronTypes[i] = INTRON_ATAC;
            else
                intronTypes[i] = INTRON_NON;
        }
        else if (left && !fwStrand)
        {
            // !!! problem with CT !!!
            if(ref[offset.refB + i + 1] == 'C' && ref[offset.refB + i + 2] == 'T')
                intronTypes[i] = INTRON_GTAG;
            else if(ref[offset.refB + i + 1] == 'C' && ref[offset.refB + i + 2] == 'T')
                intronTypes[i] = INTRON_GCAG;
            else if(ref[offset.refB + i + 1] == 'G' && ref[offset.refB + i + 2] == 'T')
                intronTypes[i] = INTRON_ATAC;
            else
                intronTypes[i] = INTRON_NON;
        }
        else if (!left && fwStrand)
        {
            // !!! problem with CT !!!
            if(ref[offset.refE - i - 1] == 'G' && ref[offset.refE - i - 2] == 'A')
                intronTypes[i] = INTRON_GTAG;
            else if(ref[offset.refE - i - 1] == 'G' && ref[offset.refE - i - 2] == 'A')
                intronTypes[i] = INTRON_GCAG;
            else if(ref[offset.refE - i - 1] == 'C' && ref[offset.refE - i - 2] == 'A')
                intronTypes[i] = INTRON_ATAC;
            else
                intronTypes[i] = INTRON_NON;
        }
        else if (!left && !fwStrand)
        {
            // !!! problem with CT !!!
            if(ref[offset.refE - i - 1] == 'C' && ref[offset.refE - i - 2] == 'A')
                intronTypes[i] = INTRON_GTAG;
            else if(ref[offset.refE - i - 1] == 'C' && ref[offset.refE - i - 2] == 'G')
                intronTypes[i] = INTRON_GCAG;
            else if(ref[offset.refE - i - 1] == 'T' && ref[offset.refE - i - 2] == 'A')
                intronTypes[i] = INTRON_ATAC;
            else
                intronTypes[i] = INTRON_NON;
        }
    }

}

#define INTRON_CANONICAL +30;
#define INTRON_SEMICANONICAL +10;
#define INTRON_NONINTRON 0;

// Can be done without the need to create 2 new strings for the dinucleotides

int intronScore(const string& dinucLeft, const string& dinucRight, bool forward) {
    if (forward && dinucLeft == "gt" && dinucRight == "ag") {
        return INTRON_CANONICAL;
    } else if (!forward && dinucLeft == "ct" && dinucRight == "ac") {
        return INTRON_CANONICAL;
    } else if (forward && dinucLeft == "gc" && dinucRight == "ag") {
        return INTRON_SEMICANONICAL;
    } else if (!forward && dinucLeft == "ct" && dinucRight == "gc") {
        return INTRON_SEMICANONICAL;
    } else if (forward && dinucLeft == "at" && dinucRight == "ac") {
        return INTRON_SEMICANONICAL;
    } else if (!forward && dinucLeft == "gt" && dinucRight == "ag") {
        return INTRON_SEMICANONICAL;
    } else
        return INTRON_NONINTRON;

}

// Naive implementation:
// Search for the best position left, and best position right independently
// then add an intron score.
void getBestSplicePositionA(int& bestRowL, int& bestColL, int& bestRowR, int& bestColR, // Used to store the best splicesite positions
                           const DynProg& dpLeft, const DynProg& dpRight,              // The dp objects for both ends of the gap
                           const Boundaries& offsetL, const Boundaries& offsetR,        // The offsets of the alligned strings
                           const string& ref,                                           // The reference string, needed to determine the introntype
                           const bool& fwStrand)
{
    int rowL        = 0;
    int rowR        = 0;
    int colBegin    = 0;
    int colEnd      = 0;
    int curBestColL = 0;
    int curBestColR = 0;

    int scoreBest   = MIN_INF;
    int scoreLeft   = MIN_INF;
    int scoreRight  = MIN_INF;
    int scoreIntron = MIN_INF;
    int score       = MIN_INF;

    string dinucLeft, dinucRight = "";

//    cout << "rowL bound: " << dpLeft.L2 <<endl;

    // try to put the splicesite between every possible nucleotide in the query
    for (rowL = 1, rowR = dpRight.L2 - 1; rowL < dpLeft.L2; rowL++, rowR--) {

        // Select the best position on this row in the left matrix
        scoreLeft = MIN_INF;
        colBegin = max(0, (int) rowL - dpLeft.bandLeft);
        colEnd = min(dpLeft.L1 + 1, (int) rowL + dpLeft.bandRight + 1);
        curBestColL = colBegin;
//        cout << "collBegin: " << colBegin << endl;
//        cout << "collEnd: " << colEnd << endl;
        
        for (int col = colBegin; col < colEnd; col++) {
            score = dpLeft.M[rowL][col];
//            cout << "score [" <<rowL <<", "<<col<<"]: " << score <<endl;
            // when equal score, choose the position as closely as possible to the main diagonal
            if (score > scoreLeft || (score == scoreLeft && abs(rowL - col) < abs(rowL - curBestColL))) {
                curBestColL = col;
                scoreLeft = score;
            }
        }
//                cout << "curBestColL:" << curBestColL <<endl;
        //    cout <<endl;

        // Select the best position on this row in the right matrix
        scoreRight = MIN_INF;
        colBegin = max(0, (int) rowR - dpRight.bandLeft);
        colEnd = min(dpRight.L1 + 1, (int) rowR + dpRight.bandRight + 1);
        curBestColR = colBegin;

        //    cout << "rowR: " << rowR << endl;
        //    cout << "colBegin: " << colBegin << "  colEnd: " << colEnd << endl;

        //    cout << "R: ";
        for (int col = colBegin; col < colEnd; col++) {
            //      cout << dpRight.M[rowR][col]<< " ";
            score = dpRight.M[rowR][col];
            // when equal score, choose the position that is closest to the main diagonal

            if (score > scoreRight || (score == scoreRight && abs(rowR - col) < abs(rowR - curBestColR))) {
                curBestColR = col;
                scoreRight = score;
            }
        }


        dinucLeft = ref.substr(offsetL.refB + curBestColL, 2);
        dinucRight = ref.substr(offsetR.refE - curBestColR - 2, 2);
        scoreIntron = intronScore(dinucLeft, dinucRight, fwStrand);
        //    debug(cerr << "colL: " << curBestColL << " dinuc: " << dinucLeft <<endl);
        //    debug(cerr << "colR: " << curBestColR << " dinuc: " << dinucRight <<endl);
        //    debug(cerr << "intronScore: " << scoreIntron << endl);
        //    cout << "curBestColL: " << curBestColL << " dinuc: " << dinucLeft <<endl;
        //    cout << rowR << " curBestColR: " << curBestColR << " dinuc: " << dinucRight <<endl;
        //    cout << dinucLeft << "-" << dinucRight << ": " << scoreIntron << endl;

        // Update the best positions if necessary
//        cerr << "scoreLeft:   " << scoreLeft <<endl;
//        cerr << "scoreRight:  " << scoreRight <<endl;
//        cerr << "scoreIntron: " << scoreIntron <<endl;
//        cerr << "scoreBest:   " << scoreBest <<endl;
        if (scoreLeft + scoreRight + scoreIntron > scoreBest) {
            bestRowL = rowL;
            bestColL = curBestColL;
            bestRowR = rowR;
            bestColR = curBestColR;
            scoreBest = scoreLeft + scoreRight + scoreIntron;
        }
        
//    cerr << "after rowL"<< rowL<<":" <<endl;
//    cerr << "bestRowL: "<< bestRowL << " bestColL: " << bestColL <<endl;
//    cerr << "bestRowR: "<< bestRowR << " bestColR: " << bestColR <<endl<<endl;
    }
//    cerr << "In:" <<endl;
//    cerr << "bestRowL: "<< bestRowL << " bestColL: " << bestColL <<endl;
//    cerr << "bestRowR: "<< bestRowR << " bestColR: " << bestColR <<endl;
}

// TODO: Less heuristic approach
void getBestSplicePositionB(int& bestRowL, int& bestColL, int& bestRowR, int& bestColR, // Used to store the best splicesite positions
                           const DynProg& dpLeft, const DynProg& dpRight,              // The dp objects for both ends of the gap
                           const Boundaries& offsetL, const Boundaries& offsetR,        // The offsets of the alligned strings
                           const string& ref,                                           // The reference string, needed to determine the introntype
                           const bool& fwStrand)
{
    int rowL, rowR, colBegin, colEnd;
    int curBestColL, curBestColR;

    int scoreBest, scoreLeft, scoreRight, scoreIntron, score = MIN_INF;

    string dinucLeft, dinucRight = "";

    // try to put the splicesite between every possible nucleotide in the query
    for (rowL = 1, rowR = dpRight.L2 - 1; rowL < dpLeft.L2; rowL++, rowR--) {

        // Select the best position on this row in the left matrix
        scoreLeft = MIN_INF;
        colBegin = max(0, (int) rowL - dpLeft.bandLeft);
        colEnd = min(dpLeft.L1 + 1, (int) rowL + dpLeft.bandRight + 1);
        curBestColL = colBegin;
        //    cout << "L: ";
        for (int col = colBegin; col < colEnd; col++) {
            //      cout << dpLeft.M[rowL][col]<< " ";
            score = dpLeft.M[rowL][col];
            // when equal score, choose the position as closely as possible to the main diagonal
            if (score > scoreLeft || (score == scoreLeft && abs(rowL - col) < abs(rowL - curBestColL))) {
                curBestColL = col;
                scoreLeft = score;
            }
        }
        //    cout <<endl;

        // Select the best position on this row in the right matrix
        scoreRight = MIN_INF;
        colBegin = max(0, (int) rowR - dpRight.bandLeft);
        colEnd = min(dpRight.L1 + 1, (int) rowR + dpRight.bandRight + 1);
        curBestColR = colBegin;

        //    cout << "rowR: " << rowR << endl;
        //    cout << "colBegin: " << colBegin << "  colEnd: " << colEnd << endl;

        //    cout << "R: ";
        for (int col = colBegin; col < colEnd; col++) {
            //      cout << dpRight.M[rowR][col]<< " ";
            score = dpRight.M[rowR][col];
            // when equal score, choose the position that is closest to the main diagonal

            if (score > scoreRight || (score == scoreRight && abs(rowR - col) < abs(rowR - curBestColR))) {
                curBestColR = col;
                scoreRight = score;
            }
        }
        //    cout << endl<<endl;

        dinucLeft = ref.substr(offsetL.refB + curBestColL, 2);
        dinucRight = ref.substr(offsetR.refE - curBestColR - 2, 2);
        scoreIntron = intronScore(dinucLeft, dinucRight, fwStrand);
        //    debug(cout << "colL: " << curBestColL << " dinuc: " << dinucLeft <<endl);
        //    debug(cout << "colR: " << curBestColR << " dinuc: " << dinucRight <<endl);
        //    debug(cout << "intronScore: " << scoreIntron << endl);
        //    cout << "curBestColL: " << curBestColL << " dinuc: " << dinucLeft <<endl;
        //    cout << rowR << " curBestColR: " << curBestColR << " dinuc: " << dinucRight <<endl;
        //    cout << dinucLeft << "-" << dinucRight << ": " << scoreIntron << endl;

        // Update the best positions if necessary
        if (scoreLeft + scoreRight + scoreIntron > scoreBest) {
            //        cout << "new bestRight," << score<<endl;
            //        cout << "col was: "<< bestColR << " now: " << curBestColR<<endl;
            bestRowL = rowL;
            bestColL = curBestColL;
            bestRowR = rowR;
            bestColR = curBestColR;
            scoreBest = scoreLeft + scoreRight + scoreIntron;
        }
    }
}

void getChainingA(Chains &chains, vector<Match_t> &matches, const long querylength, const SpliceOptions_t & spliceOptions) {
    int i = 0;
    //TODO: make sure the seeds in one cluster are from one chromosome
    int overlap; // length of the overlap between two consecutive seeds,
                 // will be negative when there is a gap between the two seeds
    const int max_overlap = 17; // the maximum value two seeds can overlap and peeled back
    while (i < matches.size()) {
        Chain chain;
        //build new cluster
        chain.push_back(matches[i]);
        long coverage = matches[i].len;
        int j = i + 1;
        while (j < matches.size() ) {  
            // Make sure the seeds are ordered on ref position as well
            if (matches[j].ref - chain.back().ref <= spliceOptions.maxSpliceLength 
            && chain.back().ref + chain.back().len - matches[j].ref < 0) {
                // overlap in query:
                overlap = chain.back().query + chain.back().len - matches[j].query;

                // The seeds do not overlap, so they can be added without any adjustment
                if (overlap <= 0) {
                    chain.push_back(matches[j]);
                    coverage += matches[j].len;
                }
                // Do not allow an overlap of greater then max_overlap
                else if (overlap < max_overlap) {
                    // peel back the two seeds a bit, so the overlap disappears if possible
                    if(chain.back().len - overlap > 2 && matches[j].len - overlap > 2) {
                        chain.back().len -= overlap;
                        chain.push_back(matches[j]);
                        chain.back().len -= overlap;
                        chain.back().ref += overlap;
                        chain.back().query += overlap;
                        coverage += matches[j].len;
                    }
                }
            }
            j++;
        }
        if ((coverage * 100) / querylength >= spliceOptions.minCoverage) {
            chains.push_back(chain);
            i = j + 1;
        } else {
            i++;
        }
    }    
}

void
init() {

}
