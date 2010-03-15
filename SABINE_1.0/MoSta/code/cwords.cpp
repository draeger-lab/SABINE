
#include "pfmic.h"

using namespace std;

int main(int argc, char* argv[])
{
  const string cswrongargs("Insufficient parameters\nCall cwords <gc> <[list:]transfac-file> <threshold-method> <threshold-parameter> [<bregularize>]\n");
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
  CPfmLoader vopfm(string(argv[2]),gc,bregularize);

  //prepare outpout
  cout << "matrix\tseq\n";

  //iterate over matrices
  for (int i=0;i<vopfm.size();i++)
    {
      //create PFM object
      vopfm[i].adjust_t(stmethod,tp);
      //print words
      vopfm[i].get_cwords(cout);
    }
  return(0);
}
