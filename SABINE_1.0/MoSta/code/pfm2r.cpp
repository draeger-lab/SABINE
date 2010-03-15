
#include "pfm2r.h"

using namespace std;

int main(int argc, char* argv[])
{
  const string cswrongargs("Insufficient parameters\nCall pfm2r <gc> <[list:]transfac-file> <threshold-parameter> [<bregularize>]\n");
  //check command line
  if (argc<4)
    {
      cout << cswrongargs;
      return(1);
    }

  //get gc content
  double gc;
  istringstream istrgc(argv[1]);
  istrgc >> gc;

  //get threshold parameter
  double tp;
  istringstream istr3(argv[3]);
  istr3 >> tp;

  bool bregularize=1;
  if (argc>4)
    {
      istringstream istr4(argv[4]);
      istr4 >> bregularize;
    }

  //matrix
  CPfmLoader vopfm(string(argv[2]),gc,bregularize,true);

  //output list
  cout << "pfm2r.o <- list()\n";
  //iterate over matrices
  for (int i=0;i<vopfm.size();i++)
    {
      cout << "pfm2r.o[[" << i+1 << "]] <- list()\n";
      //define PCM
      cout << "pfm2r.o[[" << i+1 << "]]$mpcm <- " 
	   << m2r(vopfm[i].mpcm,4) << endl;
      //define PFM
      cout << "pfm2r.o[[" << i+1 << "]]$mpfm <- " 
	   << m2r(vopfm[i].mpfm,4) << endl;
      //define PwM
      cout << "pfm2r.o[[" << i+1 << "]]$mpwm <- " 
	   << m2r(vopfm[i].mpwm,4) << endl;
      //define PSSM
      cout << "pfm2r.o[[" << i+1 << "]]$mpssm <- " 
	   << m2r(vopfm[i].mpssm,4) << endl;
      //get BS
      cout << "pfm2r.o[[" << i+1 << "]]$vbs <- vector(mode='character')\n";
      //iterate over binding sites
      for (int j=0;j<vopfm[i].vbs.size();j++)
	{
	  cout << "pfm2r.o[[" << i+1 << "]]$vbs[" << j+1 << "] <- '" 
	       << vopfm[i].vbs.at(j) << "'\n";
	}
      //get score distributions
      cout << "pfm2r.o[[" << i+1 << "]]$mbdistr <- "
	   << vopfm[i].score_distr(true).distr2r() << endl;
      cout << "pfm2r.o[[" << i+1 << "]]$msdistr <- "
	   << vopfm[i].score_distr(false).distr2r() << endl;
      //get balanced threshold
      vopfm[i].adjust_t("balanced",tp);
      cout << "pfm2r.o[[" << i+1 << "]]$tbalanced <- " 
	   << vopfm[i].get_t() << endl;
      //get typeI threshold
      vopfm[i].adjust_t("typeI",tp);
      cout << "pfm2r.o[[" << i+1 << "]]$ttypeI <- " 
	   << vopfm[i].get_t() << endl;
      //get typeII threshold
      vopfm[i].adjust_t("typeII",tp);
      cout << "pfm2r.o[[" << i+1 << "]]$ttypeII <- " 
	   << vopfm[i].get_t() << endl;
      //get typeIext threshold
      vopfm[i].adjust_t("typeIext",tp);
      cout << "pfm2r.o[[" << i+1 << "]]$ttypeIext <- " 
	   << vopfm[i].get_t() << endl;
    }
  return(0);
}
