//Class returns gamma's for the statistics for two different PFMs

#ifndef COUNTPARSHETERO_H
#define COUNTPARSHETERO_H

#include <vector>
#include <iostream>

#include "pfm.h"
#include "convolution.h"
#include "countpars.h"

class CCountParsHetero : public CCountPars
{
 public:
  CCountParsHetero(CPfm &opfmA,CPfm &opfmB,double seqgc=-1);
  double gAB,gBA,gApBp,gBpAp,gABp,gBpA,gApB,gBAp;
  double gABnull,gABpnull;
  bool bimaxp;
  int imax;
  double gmax;
  std::vector<double> gammaAB;
  std::vector<double> gammaABp;
 private:
  int tA;
  int tB;
  std::vector<int> mA;
  std::vector<int> mB;
  void get_simpars();
  void enlarge_right(std::vector<int> &m,int n);
  void enlarge_left(std::vector<int> &m,int n);
};

#endif
