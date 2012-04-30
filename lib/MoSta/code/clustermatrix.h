#ifndef SCLUSTERMATRIX_H
#define SCLUSTERMATRIX_H

#include <iostream>
#include <vector>
#include <string>

#ifdef OPENMP
#include "omp.h"
#endif

#include "pfm.h"
#include "simstat.h"
#include "similaritymatrix.h"
#include "stringutils.h"

class CClusterMatrix : public CSimilarityMatrix
{
 public:
  CClusterMatrix(std::string &asflist,double &gc,std::string &stmethod,double &t,bool absge, bool abloocv=false, bool abregularize=true);
  void do_cluster(double p=.95);
  void print_pfms(std::ostream& osout);
 private:
  bool bsge;
  bool bloocv;
};

#endif
