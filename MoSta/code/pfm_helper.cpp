
#include "pfm_helper.h"

using namespace std;

//sorts vector and returns indices
// parameter:
//   v: vector
//   vi: vector of indices (is returned)
void pfm_sortindices(vector<int> &v,vector<int> &vi)
{
  //fill index vector
  vi.reserve(v.size());
  for (unsigned i = 0; i < v.size(); ++i)
    vi.push_back(i);
  //sort index vector
  sort(vi.begin(), vi.end(), index_cmp<vector<int>&>(v));
}

