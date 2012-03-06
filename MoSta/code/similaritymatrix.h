#ifndef SSIMILARITYMATRIX_H
#define SSIMILARITYMATRIX_H

#include <iostream>
#include <vector>
#include <string>

#ifdef OPENMP
#include "omp.h"
#endif

#include "pfm.h"
#include "simstat.h"
#include "stringutils.h"
#include "sge.h"
#include "pfmloader.h"

class CSimilarityMatrix : public ISGEClient
{
 public:
  CSimilarityMatrix(std::string &asflist,double &gc,std::string &stmethod,double &t,bool bdiagonal=false,bool bregularize=true);
  void fill_matrix(int akpartial=-2,bool bfilltwice=false);
  void sge_merge(const std::string &sfn,int i) throw (EFileNotFound,EWrongFormat);
 protected:
  std::vector<double> m;
  std::vector<int> mpos;
  CPfmLoader vopfm;
  std::vector<bool> bused;
  std::string stmethod;
  double t;
  int n;
  int nm;
  double igc;
  std::string sflist;
  std::vector<double> valpha;
  bool bfilldiagonal;
 private:
  std::string scmd;
  int kpartial;
  friend std::ostream& operator<< (std::ostream& strm, const CSimilarityMatrix& obj);
};

#endif
