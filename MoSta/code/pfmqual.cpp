
#include "pfmqual.h"

using namespace std;

int main(int argc, char* argv[])
{
  const string cswrongargs("Insufficient parameters\nCall pfmqual <gc> <[list:]transfac-file> <threshold-method> <threshold-parameter> [<bregularize>] [<bunweighted>]\n");
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

  bool bnr=false;
  if (argc>6)
    {
      istringstream istr6(argv[6]);
      istr6 >> bnr;
    }

  //matrix
  CPfmLoader vopfm(string(argv[2]),gc,bregularize,true);

  //prepare output
  cout << "matrix\tquality\tquality.nr\tsens\tspec\tprec\tthreshold\talpha\tbeta\n";

  //iterate over matrices
  for (int i=0;i<vopfm.size();i++)
    {
      //compute quality
      CPfmQuality oq(vopfm[i],gc,stmethod,tp,bnr);
      cout << oq;
    }
  return(0);
}
