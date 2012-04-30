
#include "pfmic.h"

using namespace std;

int main(int argc, char* argv[])
{
  const string cswrongargs("Insufficient parameters\nCall pfmic <gc> <[list:]transfac-file> <threshold-method> <threshold-parameter> [<bregularize>]\n");
  //check command line
  if (argc<5)
    {
      cout << cswrongargs;
      return(1);
    }

  //get gc content
  double gc;
  istringstream istrgc(argv[1]);
  istrgc >> gc;

  //get threshold method
  string stmethod;
  istringstream istr3(argv[3]);
  istr3 >> stmethod;
  //get threshold parameter
  double tp;
  istringstream istr4(argv[4]);
  istr4 >> tp;

  bool bregularize=1;
  if (argc>5)
    {
      istringstream istr5(argv[5]);
      istr5 >> bregularize;
    }

  //matrix
  CPfmLoader vopfm(string(argv[2]),gc,bregularize,(stmethod=="nrwords"));

  //prepare output
  cout << "matrix\tic\talpha\tbeta\tlen\tt\n";

  //iterate over matrices
  for (int i=0;i<vopfm.size();i++)
    {
      //create PFM object
      vopfm[i].adjust_t(stmethod,tp);
      /*
      //print pssm
      cerr << "AC " << vopfm[i].sid << endl;
      for (int j=0;j<vopfm[i].nlen;j++)
	{
	  cerr << j << "\t";
	  for (int k=0;k<4;k++)
	    {
	      cerr << "\t" << vopfm[i].mpssm[j*4+k];
	    }
	  cerr << endl;
	}
      */
      cout << vopfm[i].sid << "\t" ;
      cout << vopfm[i].get_ic() << "\t";
      cout << vopfm[i].get_alpha(500) << "\t";
      cout << vopfm[i].get_beta() << "\t";
      cout << vopfm[i].nlen << "\t";
      cout << vopfm[i].get_t() << endl;
    }
  return(0);
}
