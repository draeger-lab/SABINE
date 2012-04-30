#ifndef SIMSTAT_H
#define SIMSTAT_H

#include <vector>
#include <cmath>
#include <ctime>

#include "pfm.h"
#include "countparshetero.h"
#include "mytimer.h"

class CSimStat
{
 public:
  CSimStat(CPfm &opfmA,CPfm &opfmB,double seqgc);
  double trun; //in s
  double smax;
  double ssum;
  double gABnull,gABpnull;
  int imax;
  bool bimaxp; //is on palindrome
 private:
  CTimer otimer;
  CCountParsHetero opars;
  long double precision;
};

#endif
