#include "simstat.h"

using namespace std;

//imax is always B shifted against A, and bimaxp is true
//if B has to be the reverse complementary pfm.
CSimStat::CSimStat(CPfm &opfmA,CPfm &opfmB,double seqgc) : otimer(),opars(opfmA,opfmB,seqgc),precision(.000000001),ssum(0),smax(0),imax(0),bimaxp(true)
{
  //get timer
  otimer.start();
  
  //get length of matrices
  int nA=opfmA.nlen;
  int nB=opfmB.nlen;

  double alpha = opfmA.get_alpha();
  double beta = opfmB.get_alpha();

  //get max score
  smax=log(opars.gmax/(alpha*beta));
  imax=opars.imax;
  bimaxp=opars.bimaxp;

  //get covariance
  double cov = 2 * (opars.gAB + opars.gABp);
  cov -= 4 * (nA+nB) * alpha * beta;
  cov += 4*alpha*beta;

  //make joint for complete overlap public
  gABnull = opars.gABnull;
  gABpnull = opars.gABpnull;

  //we don't need that term sind gAB and so on only contain g(0) term once.
  //  cov -= alpha * (opars.gABnull+opars.gApBnull+opars.gABpnull+opars.gApBpnull-4*beta);

  ssum = cov;

  //stop time measuring
  trun = otimer.get_time();
}
