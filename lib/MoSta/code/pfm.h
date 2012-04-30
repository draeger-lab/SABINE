#ifndef PFM_H
#define PFM_H

#include <vector>
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <map>
#include <cmath>

#include "pfm_helper.h"
#include "convolution.h"
#include "stringutils.h"
#include "sequences.h"

//Exceptions
class EThresholdNotSet{};

class CPfm {
 public:
  double nseqs;
  double ic_over;
  int n_over;
  std::string sid;
  int nlen;
  std::vector<int> mpssm;
  std::vector<double> mpcm;
  std::vector<long double> mpwm;
  std::vector<long double> mpfm;
  std::vector<double> vbg;
  bool bregularize;
  // given binding site in TRANSFAC file
  CNRSequences vbs;
  CMSequences vmbs;
  //constructs PSSM accepting only sequence sseq
  CPfm(CSequence oseq,double agc);
  //constructs and regularizes PSSM from PFM in a file
  CPfm(std::string sfn,double gc,bool bregularize,bool bsequences=false);
  CPfm(std::ifstream &f,double gc,bool bregularize,bool bsequences=false);
  //constructs PSSM by merging two PSSMs
  CPfm(CPfm& opfmA,CPfm& opfmB,int ishift,bool bcompl,bool bfakeA=false,bool bfakeB=false);
  //returns threshold
  int get_t() throw (EThresholdNotSet);
  int get_mint();
  int get_maxt();
  double get_alpha(int m=-1) throw (EThresholdNotSet);
  double get_beta() throw (EThresholdNotSet);
  void adjust_t(std::string smethod,double tp=-23880) throw(EThresholdNotSet);
  //returns information content
  double get_ic(int istart=0,int nr=-1);
  //returns gc content
  double get_gc();
  void reinit(double igc, bool bregularize);
  void get_cwords(std::vector<std::string> &vw);
  std::ostream& get_cwords(std::ostream& strm);
  CConvolution score_distr(bool brnd);
 private:
  //did we already compute convolutions
  bool bconv_rnd;
  bool bconv_sgn;
  //score distribution on random sequence
  CConvolution oconv_rnd;
  //score distribution on binding site model
  CConvolution oconv_sgn;
  void next_score(int i, std::vector<int> vj,int s,std::vector< std::vector<int> > &vw);
  double ic;
  double alpha;
  double beta;
  bool bt; //threshold set or not?
  int t;
  double gc;
  void read_pcm(std::ifstream &f,bool bsequences=false);
  void pcm2pssm(bool bregularize);
  friend std::ostream& operator<< (std::ostream& strm, const CPfm& obj);
};


#endif

