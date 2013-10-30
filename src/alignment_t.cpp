
#include "alignment_t.h"
#include <sstream>

Alignment_t::Alignment_t() : globPos(0), pos(0), cigar("*"), flag(0), rname("*"), mapq(0), editDist(0),
alignmentScore(0), cigarChars(0), cigarLengths(0), NMtag("*"), refLength(0), mateInfo(0), fwStrand(true) {
}

Alignment_t::Alignment_t(const Alignment_t & o) : globPos(o.globPos), pos(o.pos), cigar(o.cigar), flag(o.flag.to_ulong()),
rname(o.rname), mapq(o.mapq), editDist(o.editDist), alignmentScore(o.alignmentScore),
cigarChars(o.cigarChars), cigarLengths(o.cigarLengths), NMtag(o.NMtag),
refLength(o.refLength), mateInfo(o.mateInfo), fwStrand(o.fwStrand) {
}



bool Alignment_t::paired() {
    return !mateInfo.empty();
}

bool Alignment_t::concordant() const {
    return !mateInfo.empty() && mateInfo[0].concordant;
}

void Alignment_t::setLocalPos(const SparseSA& sa) {
    if (rname.empty() || rname == "*") {
        long descIndex;
        sa.from_set(globPos, descIndex, pos);
        rname = sa.descr[descIndex];
    }
}

int Alignment_t::pairedCount() {
    return mateInfo.size();
}

void Alignment_t::createCigar(bool newVersion) {
    stringstream ss;
    assert(cigarChars.size() == cigarLengths.size());
    if (newVersion)
        for (int i = 0; i < cigarChars.size(); i++)
            ss << cigarLengths[i] << cigarChars[i];
    else {
        int i = 0;
        while (i < cigarChars.size()) {
            if (cigarChars[i] == '=' || cigarChars[i] == 'X') {
                int tempLength = cigarLengths[i];
                while (i < cigarChars.size() - 1 && (cigarChars[i + 1] == '=' || cigarChars[i + 1] == 'X')) {
                    tempLength += cigarLengths[i + 1];
                    i++;
                }
                ss << tempLength << 'M';
            } else {
                ss << cigarLengths[i] << cigarChars[i];
            }
            i++;
        }
    }
    cigar = ss.str();
}

void Alignment_t::setFieldsFromCigar(const Dp_scores & scores) {//TODO: change to use special vectors
    stringstream sNM;
    stringstream sCig;
    assert(cigarChars.size() == cigarLengths.size());
    assert(cigarChars.size() > 0);
    int i = 0;
    while (i < cigarChars.size()) {
        if (cigarChars[i] == '=' || cigarChars[i] == 'X') {
            int tempLength = cigarLengths[i];
            sNM << cigarLengths[i] << cigarChars[i];
            alignmentScore += (cigarChars[i] == 'X' ? scores.mismatch * cigarLengths[i] : scores.match * cigarLengths[i]);
            while (i < cigarChars.size() - 1 && (cigarChars[i + 1] == '=' || cigarChars[i + 1] == 'X')) {
                i++;
                tempLength += cigarLengths[i];
                sNM << cigarLengths[i] << cigarChars[i];
                alignmentScore += (cigarChars[i] == 'X' ? scores.mismatch * cigarLengths[i] : scores.match * cigarLengths[i]);
            }
            sCig << tempLength << 'M';
        } else if(cigarChars[i] == 'N' ) {
            sNM << cigarLengths[i] << cigarChars[i];
            sCig << cigarLengths[i] << cigarChars[i];
        }
        
        else {
            sCig << cigarLengths[i] << cigarChars[i];
            sNM << cigarLengths[i] << cigarChars[i];
            if (cigarChars[i] == 'D' || cigarChars[i] == 'I') {
                alignmentScore += scores.openGap + scores.extendGap * cigarLengths[i];
            }
        }
        i++;
    }
    cigar = sCig.str();
    NMtag = sNM.str();
}

void Alignment_t::addMate(const Alignment_t::Ptr o, bool concordant, bool upstream) {
    Mate_t mate;
    mate.flag.set(0, true); //positions for mate to be important: 0,1,5,6,7
    mate.flag.set(1, true);
    mate.flag.set(5, o->flag.test(4));
    mate.flag.set(6, upstream);
    mate.flag.set(7, !upstream);
    mate.concordant = concordant;
    mate.pnext = o->pos;
    mate.pnextGlob = o->globPos;
    if (rname == o->rname) {
        mate.rnext = "=";
        mate.tLength = max(pos + (long) refLength - 1L, o->pos + (long) o->refLength - 1L) - min(pos, o->pos) + 1;
        long firstPos = flag.test(4) ? pos + (long) refLength - 1L : pos;
        long secondPos = o->flag.test(4) ? o->pos + (long) o->refLength - 1L : o->pos;
        if (firstPos > secondPos) mate.tLength *= -1;
    } else {
        mate.rnext = o->rname;
        mate.tLength = 0;
    }
    mateInfo.push_back(mate);
}

void Alignment_t::addMatch(const int length) {
    this->addSymbol('=', length);
}

void Alignment_t::addMisMatch(const int length) {
    this->addSymbol('X', length);
}

void Alignment_t::addDeletion(const int length) {
    this->addSymbol('D', length);
}

void Alignment_t::addInsertion(const int length) {
    this->addSymbol('I', length);
}

void Alignment_t::addSkippedRegion(const int length) {
    this->addSymbol('N', length);
}

void Alignment_t::addSymbol(const char symbol, const int length) {
    assert(length >= 0);
    assert(symbol == '=' || symbol == 'X' || symbol == 'D' ||
           symbol == 'I' || symbol == 'N' ||  symbol == 'M');

    this->cigarChars.push_back(symbol);
    this->cigarLengths.push_back(length);
}

void Alignment_t::deleteLastSymbol() {
    this->cigarChars.erase(this->cigarChars.end() - 1);
    this->cigarLengths.erase(this->cigarLengths.end() - 1);
}

void Alignment_t::addAlignmentForward(Alignment_t & alignment) {
    this->cigarChars.insert(this->cigarChars.end(),alignment.cigarChars.begin(),alignment.cigarChars.end());
    this->cigarLengths.insert(this->cigarLengths.end(),alignment.cigarLengths.begin(),alignment.cigarLengths.end());
}

void Alignment_t::addAlignmentReversed(Alignment_t & alignment) {
    this->cigarChars.insert(this->cigarChars.end(),alignment.cigarChars.rbegin(),alignment.cigarChars.rend());
    this->cigarLengths.insert(this->cigarLengths.end(),alignment.cigarLengths.rbegin(),alignment.cigarLengths.rend());        
}


void Alignment_t::addStringForward(string & s) {

    string::const_iterator it = s.cbegin();
    string::const_iterator end = s.cend();

    char c;
    size_t n = 1;
    while (it < end) {
        n = 1;
        c = *it;
        // Check for the same consecutive characters
        while (c == *(++it)) n++;     
        this->addSymbol(c, n);
    }
        
}

// Add the string of alignment characters to the alignment in reverse order
void Alignment_t::addStringReversed(string & s) {
  
    string::const_reverse_iterator rit = s.crbegin();
    string::const_reverse_iterator rend = s.crend();

    
    char c;
    size_t n = 1;
    while (rit < rend) {
        n = 1;
        c = *rit;
        // Check for the same consecutive characters
        while ( c == *(++rit)) n++;
        this->addSymbol(c, n);
    }
        
}


// Check if the cigarChars are correct, e.g. combine the same adjacent characters
void Alignment_t::postProcess() {
    CigarChars   newCigarChars   = CigarChars();
    CigarLengths newCigarLengths = CigarLengths();
    
    size_t n = 0;   
    char c;
    
    CigarChars::const_iterator   cit = this->cigarChars.cbegin();
    CigarLengths::const_iterator lit = this->cigarLengths.cbegin();
     
    while(lit < this->cigarLengths.cend()) {
        n = *lit; 
        c = *cit;
        
        while ( ++cit < this->cigarChars.cend() && c == *(cit) ) 
            n += *(++lit);
        
        lit++; //set the length pointer to the same position as the character pointer
        
        newCigarChars.push_back(c);
        newCigarLengths.push_back(n);
    }

    
    // Only copy when something has changed
    // Copy can be avoided by using (smart)pointers as fields
    if (newCigarChars.size() != this->cigarChars.size()) {
        this->cigarChars   = newCigarChars;
        this->cigarLengths = newCigarLengths;
    }
}


void Alignment_t::print(const string &ref, const string &fwQuery,const string &rcQuery, const size_t wraplength) {
    stringstream refstream;
    stringstream querystream;
    stringstream alignstream;
        
    const string query = this->fwStrand ? fwQuery : rcQuery;
    
    long refPos = this->globPos - 1; // one to zero based
    long qryPos = 0;
    
    CigarChars::const_iterator   cit = cigarChars.cbegin();
    CigarLengths::const_iterator lit = cigarLengths.cbegin();
    
    for(cit, lit; cit != cigarChars.cend(); cit++, lit++) {
        switch ( *cit ){
            case '=':
                refstream << ref.substr(refPos, *lit);
                querystream << query.substr(qryPos, *lit);
                alignstream << string(*lit, '|');
                refPos += *lit;
                qryPos += *lit;
                break;
            case 'X':
                refstream << ref.substr(refPos, *lit);
                querystream << query.substr(qryPos, *lit);
                alignstream << string(*lit, 'X');
                refPos += *lit;
                qryPos += *lit;
                break;
            case 'I':
                refstream << string(*lit, ' ');
                querystream << query.substr(qryPos, *lit);
                alignstream << string(*lit, 'I');
                qryPos += *lit;
                break;
            case 'D':
                refstream   << ref.substr(refPos, *lit);
                querystream << string(*lit, ' ');
                alignstream << string(*lit, 'D');
                refPos += *lit;
                break;
                break;
            case 'N':
                refstream << ref.substr(refPos, 3) << string(3, '.'); 
                refstream << ref.substr(refPos + *lit - 3,3);
                querystream << "  " << *lit << string(7 - to_string(*lit).size(), ' ') ;            
                alignstream << ">>>...>>>";
                refPos += *lit;
                break;               
        }
    }
    
    string refout  = refstream.str();
    string qryout  = querystream.str();
    string alignout = alignstream.str();
    size_t len;
    cout << endl;
    for(size_t i = 0; i < refout.size(); i+=wraplength) {
        len = (i + wraplength < refout.size() ) ? wraplength : i+2*wraplength - refout.size();
        cout << refout.substr(i, len) <<endl;
        cout << alignout.substr(i, len) <<endl;
        cout << qryout.substr(i, len) <<endl;
        cout << endl;
    }
    
    
}
