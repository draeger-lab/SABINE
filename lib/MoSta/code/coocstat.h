#ifndef COOCSTAT_H
#define COOCSTAT_H

#include <vector>
#include <cmath>
#include <ctime>

#include "pfm.h"
#include "countparshetero.h"
#include "mytimer.h"
#include "countstat.h"

class CCoocStat
{
 public:
  CCoocStat(CPfm &opfmA,CPfm &opfmB,double aalphaA=-1,double aalphaB=-1,double seqgc=-1);
  long double calc_winprob(long int nwinsize);
  double trun; //in s
  long double rA;
  long double rB;
  long double rAB;
  long double alphaA;
  long double alphaB;
 private:
  CTimer otimer;
  long double precision;
};

#endif
