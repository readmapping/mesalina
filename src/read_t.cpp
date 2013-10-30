/* 
 * File:   read_t.cpp
 * Author: dieter
 * 
 * Created on February 28, 2013, 11:08 AM
 */

#include "read_t.h"

Read_t::Read_t() : qname(""), sequence(""), qual("*"), rcSequence(""), rQual(""), alignments(0), pairedAlignmentCount(0) {
}

Read_t::Read_t(string name, string &seq, string &qw, bool nucleotidesOnly) : qname(name),
sequence(seq), qual(qw), alignments(0), pairedAlignmentCount(0) {
    rcSequence = sequence; //copy
    Utils::reverse_complement(rcSequence, nucleotidesOnly);
    rQual = qual;
    reverse(rQual.begin(), rQual.end());
}

Read_t::~Read_t() {
    for (int i = 0; i < alignments.size(); i++) {
        alignments[i].~shared_ptr();
    }
}

void Read_t::init(bool nucleotidesOnly) {
    rcSequence = sequence; //copy
    Utils::reverse_complement(rcSequence, nucleotidesOnly);
    rQual = qual;
    reverse(rQual.begin(), rQual.end());
}

void Read_t::postprocess(const Dp_scores & scores, const SparseSA& sa) {
    int maxScore = scores.mismatch * sequence.length();
    assert(maxScore < 0); //works only for score<0 for now (to do: add sign switch to allow positive min scores)
    int secBestScore = scores.mismatch * sequence.length();
    for (int j = 0; j < alignments.size(); j++) {
        alignments[j]->setFieldsFromCigar(scores);
        if (alignments[j]->alignmentScore > maxScore) {
            secBestScore = maxScore;
            maxScore = alignments[j]->alignmentScore;
        } else if (alignments[j]->alignmentScore < maxScore && alignments[j]->alignmentScore > secBestScore)
            secBestScore = alignments[j]->alignmentScore;
    }
    int mapq = (secBestScore == scores.mismatch * sequence.length()) ? 255 :
            250 * (maxScore - secBestScore) / (maxScore - scores.mismatch * sequence.length());
    assert(mapq <= 255 && mapq >= 0);
    for (int j = 0; j < alignments.size(); j++) {
        alignments[j]->setLocalPos(sa);
        if (alignments[j]->alignmentScore == maxScore)
            alignments[j]->mapq = mapq;
        else
            alignments[j]->flag.set(8, true);
    }
}

string Read_t::emptyAlingment(bool paired, bool mateFailed, bool upstream) {
    stringstream ss;
    int flag = 4;
    if (paired) flag += 1;
    if (mateFailed) flag += 8;
    if (paired && upstream) flag += 64;
    if (paired && !upstream) flag += 128;
    ss << qname << "\t" << flag << "\t*\t0\t0\t*\t*\t0\t0\t" << sequence << "\t" << qual << endl;
    return ss.str();
}

int Read_t::alignmentCount() {//check the correct use of this
    return alignments.size();
}

void Read_t::removeLastAlignment() {
//    delete alignments[alignments.size() - 1];
    alignments[alignments.size() - 1].~shared_ptr();
    alignments.pop_back();
}

string Read_t::printUnpairedAlignment(int i) {
    stringstream ss;
    Alignment_t::Ptr a = alignments[i];
    ss << qname << "\t" << a->flag.to_ulong() << "\t"
            << a->rname << "\t" << a->pos << "\t" <<
            a->mapq << "\t" << a->cigar << "\t" <<
            "*" << "\t" << 0 << "\t" <<
            0 << "\t" << (a->flag.test(4) ? rcSequence : sequence) <<
            "\t" << (a->flag.test(4) ? rQual : qual) << "\tAS:i:" <<
            a->alignmentScore << "\tNM:i:" << a->editDist << "\tX0:Z:" <<
            a->NMtag << endl;
    return ss.str();
}

string Read_t::printUnpairedAlignmentSplice(int i) {
    stringstream ss;
    Alignment_t::Ptr a = alignments[i];
//    qname = qname.substr(0, qname.find(" "));
    ss << qname << "\t" << a->flag.to_ulong() << "\t"
            << a->rname << "\t" << a->pos << "\t" <<
            a->mapq << "\t" << a->cigar << "\t" <<
            "*" << "\t" << 0 << "\t" <<
            0 << "\t" << (a->flag.test(4) ? rcSequence : sequence) <<
            "\t" << (a->flag.test(4) ? rQual : qual) << "\tAS:i:" <<
            a->alignmentScore << "\tNM:i:" << a->editDist << "\tX0:Z:" <<
            a->NMtag << endl;
    return ss.str();
}

string Read_t::printPairedAlignments(int i) {
    stringstream ss;
    Alignment_t::Ptr a = alignments[i];
    if (a->paired()) {
        for (int j = 0; j < a->pairedCount(); j++) {
            //flag has to be changed
            ss << qname << "\t" << (a->flag.to_ulong() | a->mateInfo[j].flag.to_ulong()) << "\t"
                    << a->rname << "\t" << a->pos << "\t" <<
                    a->mapq << "\t" << a->cigar << "\t" <<
                    a->mateInfo[j].rnext << "\t" << a->mateInfo[j].pnext << "\t" <<
                    a->mateInfo[j].tLength << "\t" << (a->flag.test(4) ? rcSequence : sequence) <<
                    "\t" << (a->flag.test(4) ? rQual : qual) << "\tAS:i:" <<
                    a->alignmentScore << "\tNM:i:" << a->editDist << "\tX0:Z:" <<
                    a->NMtag << endl;
        }
        return ss.str();
    } else {
        return printUnpairedAlignment(i);
    }
}


void Read_t::addAlignment(Alignment_t::Ptr alignment) {
    this->alignments.push_back(alignment);
}

void Read_t::showAlignment(const size_t i, const size_t wraplength) const {
//    this->alignments.at(i)->print();
}