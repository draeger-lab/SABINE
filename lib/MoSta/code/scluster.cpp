
#include "scluster.h"
#include "pfm_helper.h"
#include "clustermatrix.h"

using namespace std;

int main(int argc, char* argv[])
{
  //check command line
  if (argc<5)
    {
      cout << "Insufficient parameters\nCall scluster <gc> <[list:]transfac-file> <threshold-method> <threshold-parameter> [<use-sge>] [<p=.95>] [<LOOCV=0>] [<regularize>]\n";
      return(1);
    }

  //read parameters
  //get gc content
  double gc;
  istringstream istr1(argv[1]);
  istr1 >> gc;
  if ((gc<=0) || (gc>=1))
    {
      cerr << "Error in input: given gc content " << gc << " is not between 0 and 1.\n";
      exit(1);
    }

  //get threshold method
  string stmethod;
  istringstream istr3(argv[3]);
  istr3 >> stmethod;
  //and corresponding parameter
  double t;
  istringstream istr4(argv[4]);
  istr4 >> t;

  //use SGE engine?
  bool bsge = false;
  if (argc>5)
    {
      istringstream istr5(argv[5]);
      istr5 >> bsge;
    }

  //define p?
  double p = .95;
  if (argc>6)
    {
      istringstream istr6(argv[6]);
      istr6 >> p;
    }

  //perform LeaveOneOutCrossValidation
  bool bloocv = 0;
  if (argc>7)
    {
      istringstream istr7(argv[7]);
      istr7 >> bloocv;
    }

  //regularize matrix?
  bool bregularize = true;
  if (argc>8)
    {
      istringstream istr8(argv[8]);
      istr8 >> bregularize;
    }

   //get matrices
  string sflist(argv[2]);

  //prepare matrix to hold similarities
  CClusterMatrix ocm(sflist,gc,stmethod,t,bsge,bloocv,bregularize);
  ocm.do_cluster(p);

  //create file for matrix output
  string sfout(sflist);
  if (sflist.substr(0,5)=="list:")
    sfout = sflist.substr(5);
  stringstream sfmatout;
  sfmatout << sfout << ".cluster";
  ofstream fmatout(sfmatout.str().c_str());  
  ocm.print_pfms(fmatout);
  fmatout.close();
  
  return(0);
}
