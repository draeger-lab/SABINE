
#include "bsanno.h"

using namespace std;

int main(int argc, char* argv[])
{
  const string cswrongargs("Insufficient parameters\nCall bsanno <sequence file> <[list:]transfac-file> <threshold-method> <threshold-parameter> [<bregularize>] [<statistic>] [<gc-content for global option>]\n");
  //check command line
  if (argc<5)
    {
      cout << cswrongargs;
      return(1);
    }

  //read parameters

  //regularization
  bool bregularize=1;
  if (argc>5)
    {
      istringstream istr(argv[5]);
      istr >> bregularize;
    }

  //report pvalues?
  bool bstat=false;
  if (argc>6)
    {
      istringstream istr6(argv[6]);
      istr6 >> bstat;
    }

  //global gc content
  double gc=-1;
  if (argc>7)
    {
      istringstream istr7(argv[7]);
      istr7 >> gc;

      if ((gc<=0) || (gc>=1))
	{
	  cerr << "Error in input: given gc content " << gc << " is not between 0 and 1.\n";
	  exit(1);
	}

    }

  //get threshold method
  string stmethod;
  istringstream istr3(argv[3]);
  istr3 >> stmethod;
  //get threshold parameter
  double tp;
  istringstream istr4(argv[4]);
  istr4 >> tp;

  //PFM, set gc content arbitrarily since it'll be reinit!
  CPfmLoader vopfm(string(argv[2]),.5,bregularize,(stmethod=="nrwords"));

  //header of output
  if (bstat)
    {
      cout << "matrix\tgene\tgc\thits\tpvalue\n";
    }
  else
    {
      cout << CBSAnnotator::get_header();
    }

  //iterate over sequences
  CSequences oseqs(argv[1]);
  while (!oseqs.eof())
    {
      CSequence oseq = oseqs.next();

      //iterate over matrices
      for (int i=0;i<vopfm.size();i++)
	{
	  double igc=gc;
	  if (igc<0)
	    igc = oseq.get_gc();

	  //create PFM object
	  vopfm[i].reinit(igc,bregularize);
	  vopfm[i].adjust_t(stmethod,tp);

	  //create annotator
	  CBSAnnotator obs(vopfm[i]);
	  //annotate
	  int x = obs.annotate(oseq);

	  if (bstat)
	    {
	      //print statistic for fuond binding sites
	      CCountStat ocstat(vopfm[i],igc);
	      long double pval = ocstat.calc_p(oseq.get_length(),x);
	      //output
	      cout << vopfm[i].sid << "\t";
	      cout << oseq.get_sid() << "\t";
	      cout << igc << "\t" << x << "\t" << pval << endl;
	    } else
	    {
	      //print binding site without statistic
	      cout << obs;
	    }
	}
    }
  return(0);
}
