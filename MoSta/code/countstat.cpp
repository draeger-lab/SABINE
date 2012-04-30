#include "countstat.h"

using namespace std;

CCountStat::CCountStat(CPfm &opfm,double seqgc) : otimer(),opars(seqgc),precision(.000000001),r(-1),blambda(false)
{
  //get timer
  otimer.start();
  
  //get parameter
  opars.get_pars(opfm);

  //constant definition and variable declarations
  long double epsilon = .0000000001;
  int n = opfm.nlen;
  //read parameters from CCountPars
  long double gammapnull = min((long double)opars.gammap[0],1-epsilon);

  theta1 = 1-gammapnull;
  for (int k=1;k<n;k++)
    {
      theta1 *= pow(1-(long double)opars.gamma[k],2);
      theta1 *= pow(1-(long double)opars.gammap[k],2);
    }
  
  //compute xi's
  xi=0;
  xip=0;
  xipnull=gammapnull/(1-gammapnull);

  for (int k=1;k<n;k++)
    {
      //substitute terms
      long double factor=1;
      for (int kp=1;kp<n-k;kp++)
	{
	  factor *= (1-opars.gamma[kp]) / (1-opars.gamma[kp+k]);
	  factor *= (1-opars.gammap[kp]) / (1-opars.gammap[kp+k]);
	}
      //extend terms
      for (int kp=n-k;kp<n;kp++)
	{
	  factor *= (1-opars.gamma[kp]) * (1-opars.gammap[kp]);
	}
      //compute xiks
      long double xik = opars.gamma[k] / (1-opars.gamma[k]);
      xik *= (1-gammapnull) / (1-opars.gammap[k]);
      long double xipk = opars.gammap[k] / (1-opars.gammap[k]);
      xi += xik*factor;
      xip += xipk*factor;
    }
  
  alpha = opfm.get_alpha();
  //stop time measuring
  trun = otimer.get_time();
}

//computes lambdas
void CCountStat::calc_lambda()
{
  otimer.start();
  w = pow(xi+xip,2)+4*xi*xipnull;
  ww = pow(w,(long double).5);
  //eigenvalues
  lambda1 = (xi+xip)/2 + .5 * ww;
  lambda2 = (xi+xip)/2 - .5 * ww;
  if (!((lambda1==0) && (lambda2==0)))
    {
      //constants
      u = (w + (xi+xip+2*xipnull)*ww)/2/w;
      v = (w - (xi+xip+2*xipnull)*ww)/2/w;
    } else 
    {
      u =1; //such that rate calculation is not na
      v = 0;
    }
  blambda=true;
  trun = otimer.get_time();
}

//get theta
void CCountStat::calc_theta()
{
  otimer.start();
  //lambda already computed?
  if (!blambda)
    calc_lambda();
  //thetas
  theta.push_back(0);
  theta.push_back(theta1);
  int i = 1;
  long double s = theta1;
  while (theta[i]>precision)
    {
      long double tt = (u*pow(lambda1,i) + v*pow(lambda2,i))*theta1;
      theta.push_back(tt);
      i++;
      s += theta[i];
    }
  //normalize theta due to numerical instabilities
  for (i=0;i<theta.size();i++)
    theta[i] /= s;
  trun = otimer.get_time();
}

//get rate
// parameter:
//   n: length of sequence
long double CCountStat::calc_rate(long n)
{
  otimer.start();
  if (!blambda)
    calc_lambda();
  r = 2*alpha*n/(u/pow(1-lambda1,2) + v/pow(1-lambda2,2));
  trun = otimer.get_time();
  return(r);
}


//get p-values
// parameter:
//   n: length of sequence
//   x: number of hits
long double CCountStat::calc_p(long n,int x)
{
  otimer.start();
  //get rate
  if (r<0)
    calc_rate(n);
  //get theta
  if (theta.size()==0)
    calc_theta();
  //get some memory for the vectors
  vector<long double> prob;
  prob.reserve(100);
  p.reserve(100);
  //initialize vectors
  p.push_back((long double)1);
  prob.push_back(exp(-r));
  int i=0;
  long double ssum = prob[0];
  //compute CPD
  while ((p[i]>precision) || (i<x))
    {
      long double tt = 0;
      for (int j=0;j<=i;j++)
	{
	  long double ttt = (i+1-j<theta.size()) ? theta[i+1-j] : 0;
	  tt += (i+1-j)*ttt*prob[j];
	}
      prob.push_back(tt*r/(i+1));
      //p-value
      p.push_back(1-ssum);
      ssum += prob.back();
      i++;
    }
  trun = otimer.get_time();
  return(p[x]);
}
