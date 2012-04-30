//Class returns gamma's for the statistics

#ifndef COUNTPARS_H
#define COUNTPARS_H

#include <vector>
#include <iostream>

#include "pfm.h"
#include "convolution.h"

class CCountPars
{
 public:
  CCountPars(double seqgc);
  //calculate parameters
  void get_pars(CPfm &opfm);
  std::vector<double> gamma;
  std::vector<double> gammap;
 private:
  std::vector<double> vbg;
  void print_Qsub(std::vector<double> &Qsub,int nrangeA,int nrangeB);
 protected:
  void get_gamma(std::vector<int> &mA,int tA,std::vector<int> &mB,int tB,std::vector<double> &mygamma);
};

#endif
