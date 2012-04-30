
#include "sstat.h"
#include "similaritymatrix.h"
#include "pfm.h"
#include "pfm_helper.h"

using namespace std;

int main(int argc, char* argv[])
{
  const string cswrongargs("Insufficient parameters\nCall sstat <gc> <[list:]transfac-file> <threshold-method> <threshold-parameter> [<partial-execution>] [<return diagonal>] [<regularize>]\n");
  //check command line
  if (argc<5)
    {
      cout << cswrongargs;
      return(0);
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


  //partial execution
  int kpartial=-2;
  if (argc>5)
    {
      istringstream istr5(argv[5]);
      istr5 >> kpartial;
    }

  //return diagonal
  int bdiagonal=false;
  if (argc>6)
    {
      istringstream istr6(argv[6]);
      istr6 >> bdiagonal;
    }

  //regularization?
  int bregularize=true;
  if (argc>7)
    {
      istringstream istr7(argv[7]);
      istr7 >> bregularize;
    }
  
  //get similarity matrix
  string sflist(argv[2]);
  CSimilarityMatrix osm(sflist,gc,stmethod,t,bdiagonal,bregularize);
  osm.fill_matrix(kpartial);

  cout << osm;

  return(1);
}
