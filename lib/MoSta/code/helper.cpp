#include "helper.h"

using namespace std;

template <class T>
T strcast(string s)
{
  T x;
  istringstream istr(s);
  istr >> x;
  return(x);
}

