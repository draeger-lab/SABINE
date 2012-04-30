//Class processes convolution and returns p-values

#ifndef CONVOLUTION_H
#define CONVOLUTION_H

#include <vector>
#include<iostream>
#include <cmath>
#include <string>
#include <sstream>

class CConvolution
{
 public:
  CConvolution();
  void convolute(std::vector<int> &m,std::vector<double> &vbg,int ipos);
  double pvalue(int t);
  double pvalue_region(int t,int n);
  int threshold(double p);
  int balanced(CConvolution oconv,int n);
  void probinregion(int n);
  std::string distr2r();
 private:
  std::vector<double> vdistr;
  std::vector<double> vdistr_cum;
  int nlo;
  void update_cum();
};

#endif
