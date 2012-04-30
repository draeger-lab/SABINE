#ifndef COUNTSTAT_H
#define COUNTSTAT_H

#include <vector>
#include <cmath>
#include <ctime>

#include "pfm.h"
#include "countpars.h"
#include "mytimer.h"

class CCountStat
{
 public:
  long double theta1;
  long double alpha;
  long double xi;
  long double xip;
  long double xipnull;
  long double lambda1;
  long double lambda2;
  long double r;
  std::vector<long double> theta;
  std::vector<long double> p;
  CCountStat(CPfm &opfm,double seqgc);
  void calc_lambda();
  void calc_theta();
  long double calc_rate(long n);
  long double calc_p(long n,int x=0);
  double trun; //in s
  CCountPars opars;
 private:
  CTimer otimer;
  long double u;
  long double v;
  long double w;
  long double ww;
  long double precision;
  bool blambda;
};

#endif
