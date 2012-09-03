#include <cstdlib>
#include "costat.h"

using namespace std;

int main(int argc, char* argv[])
{
  const string cswrongargs("Insufficient parameters\nCall costat <gc> <transfac-file> <threshold-method> <threshold-parameter> [<window-size>] [<bregularize>] [<file with empirical alphas>]\n");
  //check command line
  if (argc<5)
    {
      cout << cswrongargs;
      return(1);
    }

  bool bregularize=1;
  if (argc>6)
    {
      istringstream istr6(argv[6]);
      istr6 >> bregularize;
    }

  //window size
  long nseq = -1;
  if (argc>5)
    {
      istringstream istr5(argv[5]);
      istr5 >> nseq;
    }

  //get gc content
  double gc;
  istringstream istrgc(argv[1]);
  istrgc >> gc;
  if ((gc<=0) || (gc>=1))
    {
      cerr << "Error in input: given gc content " << gc << " is not between 0 and 1.\n";
      exit(1);
    }

  //get threshold method
  string stmethod;
  istringstream istr3(argv[3]);
  istr3 >> stmethod;
  //get threshold parameter
  double tp;
  istringstream istr4(argv[4]);
  istr4 >> tp;

  //matrix
  CPfmLoader vopfm(string(argv[2]),gc,bregularize,(stmethod=="nrwords"));
  
  //adjust matrices
  vector<double> valpha;
  valpha.reserve(vopfm.size());
  for (int i=0;i<vopfm.size();i++)
    {
      //create PFM object
      vopfm[i].adjust_t(stmethod,tp);
      //get alpha
      valpha.push_back(vopfm[i].get_alpha());
    }

  //do we also have empirical alphas?
  if (argc>7)
    {
      //lets use them
      string sfalpha(argv[7]);
      //read file
      ifstream falpha;
      falpha.open(sfalpha.c_str());
      if (!falpha)
	{
	  cerr << "Error in costat.cpp::main: Cannot open file " 
	       << sfalpha 
	       << " for input.\n";
	  exit(1);
	}
      string sl;
      int i = 0;
      getline(falpha,sl);
      while ((sl!="") && (i<valpha.size()))
	{
	  double a;
	  istringstream ia(sl);
	  ia >> a;
	  valpha[i] = a;
	  getline(falpha,sl);
	  i++;
	}
      if (i<valpha.size())
	{
	  cerr << "Error in costat.cpp::main: Not enough values in " 
	       << sfalpha 
	       << " for input.\n";
	  exit(1);
	}
      falpha.close();
    }

  //pair-wise similarity probabilities
  if (nseq>0)
    cout << "matrixA\tmatrixB\tp\n";
  else
    cout << "matrixA\tmatrixB\trA\trB\trAB\talphaA\talphaB\n";
  for (int i=0;i<vopfm.size();i++)
    {
      for (int j=i+1;j<vopfm.size();j++)
	{
	  //create object for the cooccurences
	  CCoocStat ocstat(vopfm[i],vopfm[j],valpha[i],valpha[j],gc);
	  if (nseq>0)
	    {
	      long double p = ocstat.calc_winprob(nseq);
	      //output
	      cout << vopfm[i].sid << "\t" << vopfm[j].sid << "\t" << p << endl;
	    } 
	  else
	    cout << vopfm[i].sid << "\t" << vopfm[j].sid << "\t" << ocstat.rA << "\t" << ocstat.rB << "\t" << ocstat.rAB << "\t" << ocstat.alphaA << "\t" << ocstat.alphaB << endl;
	}
    }

  return(0);
}
