#ifndef PFM_HELPER_H
#define PFM_HELPER_H

#include<vector>
#include<iostream>
#include<string>


//print matrix, assume that length is divisable by 4.
template <class T> 
void print_mat(std::vector<T> &m)
{
  for (int i=0;i<(m.size()/4);i++)
    {
      std::cout << i;
      for (int j=0;j<4;j++)
	std::cout << "\t" << m[4*i+j];
      std::cout << std::endl;
    }
};

//prints vector
template <class T>
void print_vec(std::vector<T> &v)
{
  std::cout << "(" << v[0];
  for (int i=1;i<v.size();i++)
    std::cout << "," << v[i];
  std::cout << ")\n";
};

//returns reverse complementary matrix
// parameter:
//   m: matrix
//   mc: reverse comlpementary matrix
template <class T>
void pfm_reversecompl(std::vector<T> &m,std::vector<T> &mc)
{
  int n = m.size()/4;
  mc.reserve(m.size());
  for (int i=0;i<n;i++)
    for (int j=0;j<4;j++)
      mc.push_back(m[(n-i-1)*4+(3-j)]);
};


void pfm_sortindices(std::vector<int> &v,std::vector<int> &vi);

//returns score for each part of each size of v
// parameter:
//   v: vector as input
//   m: output matrix
template <class T>
void pfm_maxscores(std::vector<T> &v,std::vector<std::vector<T> > &m)
{
  int n = v.size();
  
  //get maximum score
  T imax = v[0];
  for (int i=1;i<n;i++)
    {
      if (imax<v[i])
	imax = v[i];
    }

  //iterate over all sizes
  for (int isize=1;isize<n;isize++)
    {
      m[isize].reserve(n-isize+1);
      for (int ipos=0;ipos<n-isize+1;ipos++)
	{
	  m[isize].push_back(imax);
	  for (int i=ipos;i<ipos+isize;i++)
	    m[isize][ipos] -= v[i];
	}
    }
  m[n].push_back(0);
};

//returns minimum and maximum per row
// parameter:
//   m: matrix
//   vmin: first return vector for minimum
//   vmax: second return vector for maximum
template <class T>
void pfm_minmax(std::vector<T> &m,std::vector<T> &vmin,std::vector<T> &vmax)
{
  int n = m.size()/4;
  vmin.reserve(n);
  vmax.reserve(n);
  for (int i=0;i<n;i++)
    {
      T imin = m[i*4];
      T imax = m[i*4];
      for(int j=1;j<4;j++)
	{
	  if (m[i*4+j]<imin)
	    imin = m[i*4+j];
	  if (m[i*4+j]>imax)
	    imax = m[i*4+j];
	}
      vmin.push_back(imin);
      vmax.push_back(imax);
    }
};


//for the index sorting
template <class T> struct index_cmp {
  index_cmp(const T arr) : arr(arr) {}
  bool operator()(const size_t a, const size_t b) const
  { return arr[a] < arr[b]; }
  const T arr;
};

#endif
