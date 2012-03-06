#include "convolution.h"

using namespace std;

CConvolution::CConvolution() : nlo(0),vdistr(1,1) 
{
  update_cum();
}

//convolutes existing distribution with a column of m
// parameter:
//   m: matrix with 4 entries per row
//   vbg: background distribution
//   ipos: position of row
void CConvolution::convolute(vector<int> &m,vector<double> &vbg,int ipos)
{
  //get max and min value
  int nmin=m[ipos*4];
  int nmax=nmin;
  for (int j=1;j<4;j++)
    {
      nmin = min(m[ipos*4+j],nmin);
      nmax = max(m[ipos*4+j],nmax);
    }
  //get upper and lower bound
  int nnewlo = nlo+nmin;
  int nnewup = vdistr.size() - 1 + nlo + nmax;
  //create new vector
  vector<double> vdistrnew(nnewup-nnewlo+1,0);
  //iterate over old scores
  for (int iold=0;iold<vdistr.size();iold++)
    {
      int sold = nlo+iold;
      //iterate over the letters
      for (int j=0;j<4;j++)
	{
	  int inew = sold+m[ipos*4+j] - nnewlo;
	  vdistrnew[inew] += vdistr[iold]*vbg[j];
	}
    }
  nlo = nnewlo;
  vdistr = vdistrnew;
  //update cumulative distribution
  update_cum();
}

//updates cumulative distribution
void CConvolution::update_cum()
{
  int n= vdistr.size();
  vdistr_cum.clear();
  vdistr_cum.resize(n);
  vdistr_cum[n-1] = vdistr[n-1];
  for (int i=n-2;i>=0;i--)
    vdistr_cum[i] = vdistr[i]+vdistr_cum[i+1];
}

//returns pvalue for a certain threshold
// parameter:
//   t: threshold
double CConvolution::pvalue(int t)
{
  if (t-nlo<0)
    {
      return(1.0);
    } else if (t-nlo>=vdistr_cum.size())
    {
      return(0.0);
    } else
    {
      return(vdistr_cum[t-nlo]);
    }
}

//returns threshold for selected pvalue
// parameter:
//   p: alpha
int CConvolution::threshold(double p)
{
  for (int i=vdistr_cum.size()-1;i>=0;i--)
    if (vdistr_cum[i]>=p)
      return(i+nlo);
  return(nlo);
}

//returns balanced threshold between two convolutions
// parameter:
//   oconv: CConvolution object to compare with (has to be larger!!)
//   n: size of region
int CConvolution::balanced(CConvolution oconv,int n)
{
  int i=vdistr_cum.size()-1;
  while ((1-exp(-vdistr_cum[i]*(double)n))<1-oconv.pvalue(i+nlo))
    i--;
  return(i+nlo);
}

//transforms probability distribution such that prob is at least one in n
// parameter:
//   n: size of region
void CConvolution::probinregion(int n)
{
  for (int i=0;i<vdistr_cum.size();i++)
    vdistr_cum[i] = 1-exp(-vdistr_cum[i]*(double)n);
}

//returns pvalue for a certain threshold for a region
// parameter:
//   t: threshold
double CConvolution::pvalue_region(int t,int n)
{
  return(1-exp(-pvalue(t)*(double)n));
}

//returns R code to describe the convolution
string CConvolution::distr2r()
{
  stringstream s("");
  s << "matrix(c(";
  for (int i=0;i<vdistr.size();i++)
    {
      s << i+nlo << "," << vdistr.at(i) << "," << vdistr_cum.at(i) << ",\n";
    }
  s << "NULL),byrow=T,ncol=3,nrow=" << vdistr.size() << ")\n";
  return(s.str());
}
