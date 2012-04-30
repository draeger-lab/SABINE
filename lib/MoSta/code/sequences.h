#ifndef SEQUENCES_H
#define SEQUENCES_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <iostream>
#include <string>
#include <cctype>
#include <map>

#include "exceptions.h"
#include "stringutils.h"

class CSequence {
 public:
  CSequence(std::string& asid,const std::string& asseq) throw(EWrongFormat);
  double get_gc();
  std::string get_sid() const;
  std::string get_sequence(bool breverse=false) const;
  long int get_length() const;
  static char basetrans(char cc,bool brev) throw(EWrongFormat);
  static std::map<char,char> dtrans;
  static std::map<char,char> dcompl;
 private:
  std::string sid;
  std::string sseq;
  std::string srseq;
  double gc;  
  friend std::ostream& operator<< (std::ostream& strm, const CSequence& obj);
};

class CSequences {
 public:
  CSequences(char* sf) throw(EFileNotFound,EWrongFormat);
  ~CSequences();
  bool eof();
  CSequence next();
 private:
  std::string sid;
  std::ifstream fin;
  bool beof;
};

//class for non-redundant sequences
class CNRSequences : public std::vector<CSequence> {
 public:
  void add(CSequence oseq);
  std::vector<int> vw; //vector for weights
};

//class for redundant sequences
class CMSequences : public std::vector<CSequence> {
 public:
  void add(CSequence oseq);
  std::vector<int> vw; //vector for weights
};



#endif

