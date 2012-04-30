#include "countparshetero.h"

using namespace std;


//ATTENTION: we use bg distribution from A!
CCountParsHetero::CCountParsHetero(CPfm &opfmA,CPfm &opfmB,double seqgc) : CCountPars(seqgc),tA(opfmA.get_t()),tB(opfmB.get_t()),mA(opfmA.mpssm),mB(opfmB.mpssm),bimaxp(false),imax(0),gmax(0.0)
{
  //  cerr << "\n\n" << opfmA.sid << "\t" << opfmB.sid << "\n";

  get_simpars();
}

// enlarges given matrix to the left
// parameter:
//   n: final size!
void CCountParsHetero::enlarge_left(vector<int> &m,int n)
{
  if (m.size()<n)
    {
      int nloc = m.size();
      vector<int> mt=m;
      m.clear();
      m.reserve(n);
      for (int i=0;i<(n-nloc);i++)
	m.push_back(0);
      for (int i=0;i<nloc;i++)
	m.push_back(mt[i]);
    }
}
// enlarges given matrix to the right
// parameter:
//   n: final size!
void CCountParsHetero::enlarge_right(vector<int> &m,int n)
{
  if (m.size()<n)
    {
      int nloc = m.size();
      m.reserve(n);
      for (int i=0;i<(n-nloc);i++)
	{
	  m.push_back(0);
	}
    }
}

void print_m(vector<int> &m)
{
  cerr << "i\tA\tC\tG\tT\n";
  for (int i=0;i<m.size()/4;i++)
    {
      cerr << i << "\t";
      for (int j=0;j<4;j++)
	cerr << m[i*4+j] << "\t" ;
      cerr << endl;
    }
  cerr << "\n";
}

//computes gamma values for similarity measure
void CCountParsHetero::get_simpars()
{
  //get reverse complementary matrix
  vector<int> mAc;
  pfm_reversecompl(mA,mAc);
  vector<int> mBc;
  pfm_reversecompl(mB,mBc);

  int nlen = mA.size()+mB.size()-4;
  int nA=mA.size()/4;
  int nB=mB.size()/4;

  //make matrices equally large, first A
  vector<int> mAl(mA);
  //now for B
  vector<int> mBr(mB);
  vector<int> mBcr(mBc);

  //enlarge matrices, first A
  enlarge_left(mAl,nlen);
  //now for B
  enlarge_right(mBr,nlen);
  enlarge_right(mBcr,nlen);

  //get necessary gammas
  get_gamma(mAl,tA,mBr,tB,gammaAB);
  get_gamma(mAl,tA,mBcr,tB,gammaABp);

  //initialize gXXs
  gAB =0;
  gBA=0;
  gApBp=0;
  gBpAp =0;
  gABp=0;
  gBpA=0;
  gApB=0;
  gBAp=0;
  //compute gXX
  //  cerr << "i\tgammaAB\tgammaABp" << endl;
  for (int i=0;i<nlen/4;i++)
    {
      //      cerr << i << "\t" << gammaAB[i] << "\t" << gammaABp[i] << endl;
      gAB += gammaAB[i];
      gABp += gammaABp[i];
      if (gammaAB[i]>gmax)
	{
	  gmax = gammaAB[i];
	  imax = i-(nB-1);
	  bimaxp=false;
	}
      if (gammaABp[i]>gmax)
	{
	  gmax = gammaABp[i];
	  imax = i-(nB-1);
	  bimaxp=true;
	}
    }

  //use symmetrie
  gApBp=gAB;
  gBpAp=gAB;
  gBA=gAB;
  gApB=gABp;
  gBpA=gABp;
  gBAp=gABp;

  //get probs for shift zero
  gABnull = gammaAB[nB-1];
  gABpnull = gammaABp[nB-1];

  /*  
  //define remaining matrices (for testing purposes)
  vector<int> mAr(mA);
  vector<int> mAcl(mAc);
  vector<int> mAcr(mAc);
  vector<int> mBl(mB);
  vector<int> mBcl(mBc);

  enlarge_right(mAcr,nlen);
  enlarge_left(mBl,nlen);
  enlarge_left(mBcl,nlen);
  enlarge_right(mAr,nlen);
  enlarge_left(mAcl,nlen);
  
  //get remaining gammas
  vector<double> gammaApBp;
  vector<double> gammaBA;
  vector<double> gammaApB;
  vector<double> gammaBpA;
  vector<double> gammaBpAp;
  vector<double> gammaBAp;
  get_gamma(mBcl,tB,mAcr,tA,gammaBpAp);
  get_gamma(mBl,tB,mAcr,tA,gammaBAp);
  get_gamma(mAcl,tA,mBcr,tB,gammaApBp);
  get_gamma(mBl,tB,mAr,tA,gammaBA);
  get_gamma(mAcl,tA,mBr,tB,gammaApB);
  get_gamma(mBcl,tB,mAr,tA,gammaBpA);

  //initialize gXXs
  gAB =0;
  gBA=0;
  gApBp=0;
  gBpAp =0;
  gABp=0;
  gBpA=0;
  gApB=0;
  gBAp=0;

  //compute gXX
  cerr << "i\tgammaAB\tgammaApBp\tgammaBA\tgammaBpAp\tgammaABp\tgammaApB\tgammaBpA\tgammaBAp" << endl;
  for (int i=0;i<nlen/4;i++)
    {
      cerr << i << "\t" << gammaAB[i] << "\t" << gammaApBp[i] << "\t" << gammaBA[i] << "\t" << gammaBpAp[i] << "\t" << gammaABp[i] << "\t" << gammaApB[i] << "\t" << gammaBpA[i] << "\t" << gammaBAp[i] << endl;
      gAB += gammaAB[i];
      gBA += gammaBA[i];
      gApBp += gammaApBp[i];
      gBpAp += gammaBpAp[i];
      gABp += gammaABp[i];
      gBAp += gammaBAp[i];
      gApB += gammaApB[i];
      gBpA += gammaBpA[i];
    }
  cerr << "\t" << gAB << "\t" << gApBp << "\t" << gBA << "\t" << gBpAp << "\t" << gABp << "\t" << gApB << "\t" << gBpA << "\t" << gBAp << endl;
  if ((gAB!=gBA) || (gBA!=gApBp) || (gApBp!=gBpAp))
    cerr << "***** STOP for (i)\n";
  if ((gABp!=gBpA) || (gBpA!=gApB) || (gApB!=gBAp))
    cerr << "***** STOP for (ii):\t"  << gABp << "==" << gBpA  << "==" << gApB << "==" << gBAp << endl;

  //get probs for shift zero
  gABnull = gammaAB[nB-1];
  gABpnull = gammaABp[nB-1];
  gApBnull = gammaBAp[nA-1];
  gApBpnull = gammaBpAp[nA-1];

  cerr << "gAB0=" << gABnull << "\tgABp0=" << gABpnull << "\tgApB0=" << gApBnull << "\tgApBpnull=" << gApBpnull << endl;

  */

}
