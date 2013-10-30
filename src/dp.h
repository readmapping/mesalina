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
 */

#ifndef DP_H
#define	DP_H

#include <iostream>
#include <string>
#include <algorithm>
#include <stdlib.h>
#include <stdio.h>

//#define DEBUG 0

using namespace std;

struct Boundaries{
    Boundaries(): refB(0), queryB(0), refE(0), queryE(0) {}
    Boundaries(long refBegin, long refEnd, long queryBegin, long queryEnd): refB(refBegin), refE(refEnd), queryB(queryBegin), queryE(queryEnd){}
    long refB;
    long refE;
    long queryB;
    long queryE;
    
    void setValues(long refBegin, long refEnd, long queryBegin, long queryEnd);
};

struct Dp_scores{
    Dp_scores(): match(0), mismatch(-1), extendGap(-1), openGap(0) {}
    Dp_scores(int m, int mm, int oG, int eG): match(m),mismatch(mm), openGap(oG), extendGap(eG) {}
    int match;
    int mismatch;
    int openGap;
    int extendGap;
    int scoreMatrix[127][127];
    void updateScoreMatrixDna();
};

struct Dp_type{
    Dp_type(): freeRefB(false),freeRefE(false),freeQueryB(false),freeQueryE(false),local(false) {}
    Dp_type(bool refB, bool refE, bool qB, bool qE, bool loc):freeRefB(refB),freeRefE(refE),freeQueryB(qB),freeQueryE(qE),local(loc) {}
    bool freeRefB;
    bool freeRefE;
    bool freeQueryB;
    bool freeQueryE;
    bool local;
};

struct Dp_output{
    Dp_output(): traceRef(""),traceQuery(""),editDist(0),dpScore(0), cigarChars(0), cigarLengths(0) {}
    string traceRef;//remove???
    string traceQuery;
    int editDist;
    int dpScore;
    vector<char> cigarChars;//static, fixed length field?
    vector<int> cigarLengths;
    void clear(){cigarChars.clear(); cigarLengths.clear(); editDist = dpScore = 0; traceRef = traceQuery = ""; };
};

enum outputType{
    DPSCORE,
    SCORE,
    ALIGNMENT,
    ERRORSTRING,
    ALL
};

struct DynProg {

//These fields should maybe be positioned somewhere else for multi-threading!!!
int ** M;
int ** UP;
int ** LEFT;
int DP_L2;
int DP_L1;
int L1;
int L2;
int bandSize;
int bandLeft;
int bandRight;
bool banded;
Dp_scores& scores;

// Constructor builds sparse suffix array.
DynProg(int dimension, bool affine, Dp_scores& scoreFunction);
~DynProg();

void initDPMatrix(int dimensionL2, int dimensionL1, bool affine);
void resizeDPMatrix(int dimensionL2, int dimensionL1, bool affine);
void deleteDPMatrix(bool affine);

int updateMatrix(const Dp_type& type);
int dpFillOptStatic(const string& ref,const string& query, bool forward, 
        const Boundaries& offset, const Dp_type& type);
int dpFillStatic(const string& ref,const string& query, bool forward,
        const Boundaries& offset, const Dp_type& type);
int dpFillStatic(const string& ref,const string& query, bool forward,
        const Boundaries& offset, const Dp_type& type, bool left);
int findTraceBackPosStatic(bool forward, int* const i, int* const j,const Dp_type& type);
int dpTraceBackStatic(int& i, int& j, const Dp_type& type, Dp_output& output,
        stringstream & ss, const Boundaries& offset, bool forward, const string& ref, const string& query);
int dpTraceBackStatic(int& i, int& j, const Dp_type& type, Dp_output& output,
        stringstream & ss, const Boundaries& offset, bool forward, const string& ref, const string& query, bool left);
int dpTraceBackOptStatic(int& i, int& j, const Dp_type& type, Dp_output& output,
        stringstream & ss, const Boundaries& offset, bool forward, const string& ref, const string& query);

void  print_matrices(const string& ref, const string& query, const Boundaries& offset, bool gapMatrices );
void  print_matrices(const string& ref, const string& query, const Boundaries& offset, bool gapMatrices, bool left );
void  print_seq( const string& ref,const string& query, Boundaries& offset );


//int dp( const string&, const string&, boundaries&, const dp_type&, 
//        const outputType&, dp_output&, bool print = false);

//int dpBand( const string&, const string&, boundaries&, const dp_type&, 
//        const outputType&, dp_output&, int bandSize, bool print = false);

int dpBandStatic( const string&, const string&, Boundaries&,
        const Dp_type&, const outputType&, Dp_output&, int bandSize, bool print = false);

int dpBandFull( const string&, const string&, Boundaries&,
        const Dp_type&, const outputType&, Dp_output&, int bandLeft, int bandRight, bool print = false);


};

#endif	/* DP_H */

