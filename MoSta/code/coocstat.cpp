#include "coocstat.h"

using namespace std;

//type for random variable set
struct RV{
  int position;
  bool bp;
  bool bA;
};

CCoocStat::CCoocStat(CPfm &opfmA,CPfm &opfmB,double aalphaA,double aalphaB,double seqgc) : otimer(),precision(.000000001), rA(-1.0), rB(-1.0), rAB(-1.0), alphaA((long double)aalphaA),alphaB((long double)aalphaB)
{
  //get timer
  otimer.start();
  
  //constant definition and variable declarations
  long double epsilon = .0000000001;

  //get length and alpha of matrices
  int nA=opfmA.nlen;
  int nB=opfmB.nlen;
  if (alphaA<0)
    alphaA=opfmA.get_alpha();
  if (alphaB<0)
    alphaB=opfmB.get_alpha();
  //get parameter
  CCountStat ostatA(opfmA,seqgc);
  CCountStat ostatB(opfmB,seqgc);
  CCountParsHetero oparsAB(opfmA,opfmB,seqgc);

  //length of gamma
  int nAB = oparsAB.gammaAB.size();

  //get gamma(A and B) while A has the fixed position 
  vector<long double> gammaAaB;
  vector<long double> gammapAaB;
  gammaAaB.reserve(nAB);
  gammapAaB.reserve(nAB);
  for (int k=0;k<nAB;k++)
    {
      if (nA>=nB)
	{
	  gammaAaB.push_back(min((long double)oparsAB.gammaAB.at(k),1-epsilon));
	  gammapAaB.push_back(min((long double)oparsAB.gammaABp.at(k),1-epsilon));
	}
      else
	{
	  gammaAaB.push_back(min((long double)oparsAB.gammaAB.at(nAB-k-1),1-epsilon));
	  gammapAaB.push_back(min((long double)oparsAB.gammaABp.at(nAB-k-1),1-epsilon));
	}
    }

  //get start position for gamma for A and B
  int offsetA=nB-1;

  //create set of random variables
  vector<RV> vRVs;
  vRVs.reserve(2*nA+2*nB);
  int nmAB = max(nA,nB);
  //add RVs
  for (int k=0;k<nmAB;k++)
    {
      //As
      //normal strand
      struct RV XA;
      XA.position=k;
      XA.bp=false;
      XA.bA=true;
      vRVs.push_back(XA);
      //alternative strand
      struct RV YA;
      YA.position=k;
      YA.bp=true;
      YA.bA=true;
      vRVs.push_back(YA);
      //add Bs
      //normal strand
      struct RV XB;
      XB.position=k;
      XB.bp=false;
      XB.bA=false;
      vRVs.push_back(XB);
      //alternative strand
      struct RV YB;
      YB.position=k;
      YB.bp=true;
      YB.bA=false;
      vRVs.push_back(YB);
    }

  
  //get probability
  long double Z = 1-alphaA;
  //Z- without actual position
  long double Zm=Z;
  long double ZA = Z;
  long double ZmA = Z;
  long double ZB = 1-alphaB;
  long double ZmB = ZB;

  int nr=vRVs.size();
  for (int j=1;j<nr;j++)
    {
      long double s=1;
      long double sA=1;
      long double sB=1;
      RV X = vRVs[j];
      for (int i=0;i<j;i++)
	{
	  //we have P(Y|X)
	  RV Y = vRVs[i];
	  /*
	  if (X.bA)
	    cerr << "A";
	  else
	    cerr << "B";
	  if (X.bp)
	    cerr << "'";
	  cerr << "(" << X.position << ")\t";
	  if (Y.bA)
	    cerr << "A";
	  else
	    cerr << "B";
	  if (Y.bp)
	    cerr << "'";
	  cerr << "(" << Y.position << ")\t";
	  */
	  //get correct gamma
	  long double g=1;
	  int ipos = X.position-Y.position;
	  if (X.bA==Y.bA)
	    {
	      if (X.bA)
		{
		  if (ipos<nA)
		    {
		      if (X.bp==Y.bp)
			g = ostatA.opars.gamma.at(ipos);
		      else
			g = ostatA.opars.gammap.at(ipos);
		    }
		  else
		    g = alphaA;
		}
	      else
		{
		  if (ipos<nB)
		    {
		      if (X.bp==Y.bp)
			g = ostatB.opars.gamma.at(ipos);
		      else
			g = ostatB.opars.gammap.at(ipos);
		    }
		  else
		    g = alphaB;
		}
	    }
	  else
	    {
	      if (X.bA)
		{
		  if (offsetA-ipos>=0)
		    {
		      if (X.bp==Y.bp)
			{
			  if (!X.bp)
			    g = gammaAaB.at(nAB-(offsetA-ipos)-1)/alphaA;
			  else
			    g = gammaAaB.at(offsetA-ipos)/alphaA;
			}
		      else
			{
			  if (!X.bp)
			    g = gammapAaB.at(nAB-(offsetA-ipos)-1)/alphaA;
			  else
			    g = gammapAaB.at(offsetA-ipos)/alphaA;
			}
		    }
		  else
		    g = alphaB;
		}
	      else
		{
		  if (offsetA+ipos<nAB)
		    {
		      if (X.bp==Y.bp)
			{
			  if (!X.bp)
			    g = gammaAaB.at(offsetA+ipos)/alphaB;
			  else
			    g = gammaAaB.at(nAB-(offsetA+ipos)-1)/alphaB;
			}
		      else
			{
			  if (!X.bp)
			    g = gammapAaB.at(offsetA+ipos)/alphaB;
			  else
			    g = gammapAaB.at(nAB-(offsetA+ipos)-1)/alphaB;
			}
		    } 
		  else
		    g = alphaA;
		}
	    }
	  //	  cerr << g << endl;
	  //add to product
	  s *= (1-g);
	  //and for marginal
	  if (X.bA & Y.bA)
	    sA *= (1-g);
	  if ((!X.bA) & (!Y.bA))
	    sB *= (1-g);
	}
      if (X.bA)
	s *= alphaA;
      else
	s *= alphaB;
      //compute Z
      Z -= s;
      if (j<nr-4)
	Zm -= s;
      //compute marginal
      if (X.bA)
	{
	  ZA -= (sA*alphaA);
	  if (j<nr-4)
	    ZmA -= (sA*alphaA);
	}
      if (!X.bA)
	{
	  ZB -= (sB*alphaB);
	  if (j<nr-4)
	    ZmB -= (sB*alphaB);
	}
    }

  rAB = 1-Z/Zm;
  rA = 1-ZA/ZmA;
  rB = 1-ZB/ZmB;

  /*
  cerr << "ZA=" << ZA << "\tZA-=" << ZmA << endl;
  cerr << "ZB=" << ZB << "\tZB-=" << ZmB << endl;
  cerr << "ZAB=" << Z << "\tZAB-=" << Zm << endl;
  cerr << "rA=" << rA << endl;
  cerr << "rB=" << rB << endl;
  cerr << "rAB=" << rAB << endl;
  */

  //stop time measuring
  trun = otimer.get_time();

}


//computes probability for a windows size
// parameter:
//   nwinsize: size of the window
long double CCoocStat::calc_winprob(long int nwinsize)
{
  return(1-exp(-rA*nwinsize)-exp(-rB*nwinsize)+exp(-rAB*nwinsize));
}


