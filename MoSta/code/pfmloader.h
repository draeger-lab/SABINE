#ifndef PFMLOADER_H
#define PFMLOADER_H

#include <vector>
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>

#include "stringutils.h"
#include "pfm.h"

class CPfmLoader : public std::vector<CPfm> {
 public:
  CPfmLoader(std::string sfn,double gc,bool bregularize,bool bsequences=false);
};

#endif
