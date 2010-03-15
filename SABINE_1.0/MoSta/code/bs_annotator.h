#ifndef BS_ANNOTATOR_H
#define BS_ANNOTATOR_H

#include <iostream>
#include <vector>
#include <string>
#include <map>

#include "pfm.h"
#include "sequences.h"

class CBSAnnotator {
 public:
  CBSAnnotator(CPfm& opfm);
  int annotate(const CSequence& oseq);
  static std::string get_header();
 private:
  std::map<char,int> dcompl;
  int t;
  int nlen;
  //scoring matrix
  std::vector<int> mpssm;
  //save hits
  std::string spfmid;
  std::vector<std::string> vsid;
  std::vector<int> vstrand;
  std::vector<long int> vpos_start;
  std::vector<long int> vpos_end;
  std::vector<long int> vseq_start;
  std::vector<long int> vseq_end;
  void _annotate(const CSequence& oseq,bool breverse);
  void addhit(std::string& sid,long int ipos,long int n,bool breverse);
  friend std::ostream& operator<< (std::ostream& strm, const CBSAnnotator& obj);
};


#endif

