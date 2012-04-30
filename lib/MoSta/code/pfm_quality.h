#ifndef PFMQUALITY_H
#define PFMQUALITY_H

#include <cmath>

#include "pfm.h"
#include "simstat.h"
#include "sequences.h"

class CPfmQuality
{
 public:
  CPfmQuality(CPfm &opfm,double gc,std::string stmethod,double tp,bool bnr=false);
  int get_best_t();
  double get_best_q();
 private:
  std::string sid;
  std::vector<double> valpha;
  std::vector<double> vbeta;
  std::vector<double> vquality;
  std::vector<double> vquality_nr;
  std::vector<double> vsens;
  std::vector<double> vspec;
  std::vector<double> vprec;
  int tmin;
  int bestt;
  friend std::ostream& operator<< (std::ostream& strm, const CPfmQuality& obj);
};

#endif
