#ifndef SGE_H
#define SGE_H

#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <cstdlib>
#include <vector>

#include "stringutils.h"
#include "mytimer.h"
#include "exceptions.h"

class ESGEError : public EBase
{
 public:
  ESGEError(const std::string &ss) : EBase(ss) {};
  ESGEError() {};
};

class ISGEClient
{
 public:
  virtual void sge_merge(const std::string &sfn,int i) throw (EFileNotFound,EWrongFormat) = 0;
};

class CSGEMaster
{
 public:
  CSGEMaster(ISGEClient* aoclient,const char* csid,int n);
  void submit(std::string scmd);
  void finish() throw(EWrongFormat,ESGEError);
 private:
  ISGEClient *oclient;
  std::string sid;
  std::vector<long> vids;
  std::string sdir;
  int inr;
};

#endif
