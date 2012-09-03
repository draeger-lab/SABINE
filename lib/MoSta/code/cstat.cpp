#include <cstdlib>
#include "cstat.h"

using namespace std;

int main(int argc, char* argv[])
{
  const string cswrongargs("Insufficient parameters\nCall cstat <gc> <[list:]transfac-file> <threshold-method> <threshold-parameter> [bregularize] [parameter|lambda|theta|rate|cpd] [sequence-length]\n");
  //check command line
  if (argc<5)
    {
      cout << cswrongargs;
      return(1);
    }

  bool bregularize=1;
  if (argc>5)
    {
      istringstream istr5(argv[5]);
      istr5 >> bregularize;
    }

  bool bverbose = (argc>6);
  bool blongoutput = 0;
  long nseq = -1;
  string soutput("parameter");

  if (argc>6)
    {
      istringstream istr6(argv[6]);
      istr6 >> soutput;
      blongoutput = (soutput != "parameter");
      if ((soutput == "rate") || (soutput=="cpd"))
	{
	  if (argc<8)
	    {
	      cout << cswrongargs << "\noutput parameter was set to rate or cpd, thus, we need a sequence length in the next parameter!\n\n";
	      return(1);
	    }
	  istringstream istr7(argv[7]);
	  istr7 >> nseq;
	}
    }

  //read parameters

  //get gc content
  double gc;
  istringstream istrgc(argv[1]);
  istrgc >> gc;
  if ((gc<=0) || (gc>=1))
    {
      cerr << "Error in input: given gc content " << gc << " is not between 0 and 1.\n";
      exit(1);
    }

  //is there a file with more gc contents
  vector<double> vgc;
  if (argc>8)
    {
      string sgc(argv[8]);

      //is it a list with many gc contents?
      if (sgc.substr(0,5)=="file:")
	{
	  string sfgclist(sgc.substr(5));
	  ifstream fgc;
	  fgc.open(sfgclist.c_str());
	  if (!fgc)
	    {
	      cerr << "Error in cstat.cpp::main(): Cannot open file " 
		   << sfgclist 
		   << " for input.\n";
	      exit(1);
	    }
	  string sk;
	  while (getline(fgc,sk))
	    {
	      if (sk!="")
		{
		  double kgc;
		  istringstream istrgc(sk);
		  istrgc >> kgc;
		  if ((kgc<=0) || (kgc>=1))
		    {
		      cerr << "Error in input: given gc content " << kgc << " is not between 0 and 1.\n";
		      exit(1);
		    }
		  vgc.push_back(kgc);
		}
	    }
	} else
	{
	  double gcseq;
	  istringstream istrgcd(sgc);
	  istrgcd >> gcseq;
	  if ((gcseq<=0) || (gcseq>=1))
	    {
	      cerr << "Error in input: given gc content " << gcseq << " is not between 0 and 1.\n";
	      exit(1);
	    }
	  vgc.push_back(gcseq);
	}
    } else 
    {
      vgc.push_back(gc);
    }

  //get threshold method
  string stmethod;
  istringstream istr3(argv[3]);
  istr3 >> stmethod;
  //get threshold parameter
  double tp;
  istringstream istr4(argv[4]);
  istr4 >> tp;

  //load PFMs
  CPfmLoader vopfm(string(argv[2]),gc,bregularize,(stmethod=="nrwords"));
  
  //prepare output
  if (!blongoutput)
    {
      //print header
      cout << "matrix\talpha\ttheta1\txi\txip\txipnull\tt\tlength\tseqgc";
      if (bverbose)
	cout << "\ttime";
      cout << endl;
    } else 
    {
      cout << "matrix\tparameter\tindex\tvalue\n";
    }
  //iterate over matrices
  for (int i=0;i<vopfm.size();i++)
    {
      //create PFM object
      vopfm[i].adjust_t(stmethod,tp);

      for (int j=0;j<vgc.size();j++)
	{
	  //create object for the count statistic
	  CCountStat ocstat(vopfm[i],vgc[j]);

	  //output
	  if (!blongoutput)
	    {
	      cout << vopfm[i].sid << "\t";
	      cout << ocstat.alpha << "\t";
	      cout << ocstat.theta1 << "\t";
	      cout << ocstat.xi << "\t";
	      cout << ocstat.xip << "\t";
	      cout << ocstat.xipnull << "\t";
	      cout << vopfm[i].get_t() << "\t";
	      cout << vopfm[i].nlen << "\t";
	      cout << vgc[j];
	      if (bverbose)
		cout << "\t" << ocstat.trun;
	      cout << endl;
	    } else 
	    {
	      cout << vopfm[i].sid << "\talpha\t0\t" << ocstat.alpha << endl;
	      cout << vopfm[i].sid << "\ttheta1\t0\t" << ocstat.theta1 << endl;
	      cout << vopfm[i].sid << "\txi\t0\t" << ocstat.xi << endl;
	      cout << vopfm[i].sid << "\txip\t0\t" << ocstat.xip << endl;
	      cout << vopfm[i].sid << "\txipnull\t0\t" << ocstat.xipnull << endl;
	      cout << vopfm[i].sid << "\tt\t0\t" << vopfm[i].get_t() << endl;
	      cout << vopfm[i].sid << "\tlength\t0\t" << vopfm[i].nlen << endl;
	      cout << vopfm[i].sid << "\tseqgc\t0\t" << vgc[j] << endl;
	      
	      //lambda calculation and output
	      ocstat.calc_lambda();
	      cout << vopfm[i].sid << "\tlambda1\t0\t" << ocstat.lambda1 << endl;
	      cout << vopfm[i].sid << "\tlambda2\t0\t" << ocstat.lambda2 << endl;
	      
	      if (soutput != "lambda")
		{
		  if (soutput == "rate")
		    {
		      //calculate the rate
		      ocstat.calc_rate(nseq);
		      cout << vopfm[i].sid << "\tr\t0\t" << ocstat.r << endl;
		    } else 
		    {
		      //calculate theta and output
		      ocstat.calc_theta();
		      for (int j =0;j<ocstat.theta.size();j++)
			cout << vopfm[i].sid << "\ttheta\t" << j << "\t" << ocstat.theta[j] << endl;
		      
		      if (soutput == "cpd")
			{
			  ocstat.calc_p(nseq);
			  for (int j =0;j<ocstat.p.size();j++)
			    cout << vopfm[i].sid << "\tp\t" << j << "\t" << ocstat.p[j] << endl;
			}
		    }
		}

	      cout << vopfm[i].sid << "\ttime\t0\t" << ocstat.trun << endl;
	    }
	}
    }
  return(0);
}
