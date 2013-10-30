/* 
 * File:   mate_t.h
 * Author: dieter
 *
 * Created on February 28, 2013, 11:05 AM
 */

#pragma once

#include <bitset>
#include <string>

using namespace std;

struct Mate_t {
    Mate_t(): pnextGlob(0), concordant(false), flag(0), tLength(0), rnext("*"),
    pnext(0) {}
    Mate_t(const Mate_t & o): pnextGlob(o.pnextGlob), concordant(o.concordant),
    flag(o.flag.to_ulong()), tLength(o.tLength), rnext(o.rnext), pnext(o.pnext) {}
    bitset<11> flag;
    string rnext;
    long pnext;
    long pnextGlob;
    int tLength;
    bool concordant;
};

