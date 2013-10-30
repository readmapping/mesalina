/* 
 * File:   types.h
 * Author: Dieter De Smedt
 *
 * Created on December 1, 2012, 9:52 PM
 */

#pragma once

#include "sparseSA.h"


// Type used for splicing, this will contain different seeds which may be
// extended to a full spliced mapping
typedef std::vector< Match_t >   Chain;
typedef std::vector< Chain > Chains;
typedef std::vector< Chain * > Chains_ptr;

enum IntronType {INTRON_GTAG, INTRON_GCAG, INTRON_ATAC, INTRON_NON};

typedef std::vector< IntronType > IntronTypes;
#define MIN_INF -200000

