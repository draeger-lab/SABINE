#ifndef HELPER_H
#define HELPER_H

#include <sstream>
#include <string>
#include <vector>

template <class T>
T strcast(std::string s);

//transforms matrix to R matrix
// parameter:
//   vector<T>: matrix
//   ncol: number of columns
template <class T>
std::string m2r(std::vector<T> &v,int ncol)
{
  //output string
  std::stringstream s("");
  s << "matrix(c(";
  //iterate over rows of matrix
  for (int i=0;i<(int)v.size()/ncol;i++)
    {
      //iterate over columns
      for (int j=0;j<ncol;j++)
	{
	  s << v.at(i*ncol+j) << ",";
	}
    }
  //add NULL for last comma
  s << "NULL),ncol=" << ncol << ",nrow=" << (int)v.size()/ncol << ",byrow=T)";
  return(s.str());
};

#endif
